#!/usr/bin/env python3
"""Attach differential-expression evidence from an Excel workbook to dark candidates.

The script reads all sheets by default. Each sheet is treated as a contrast unless
--contrast-column is supplied and present. It auto-detects ID, log2FC, and adjusted
P-value columns, then writes one row per candidate with compact DE support fields.
"""

import argparse
import csv
import math
import os
import re
import sys
from collections import Counter, defaultdict

import pandas as pd

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    limit = sys.maxsize
    while True:
        limit //= 10
        try:
            csv.field_size_limit(limit)
            break
        except OverflowError:
            continue


def parse_args():
    p = argparse.ArgumentParser(description="Add DE evidence from Excel workbook to dark candidates.")
    p.add_argument("--candidates", required=True, help="Candidate TSV")
    p.add_argument("--excel", required=True, help="DE Excel workbook")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--prefix", default="equina_dark_candidates")
    p.add_argument("--candidate-id-column", default="protein_id")
    p.add_argument("--candidate-extra-id-columns", default="transcript_id,gene_id,source_fasta_id", help="Comma-separated optional candidate ID columns for matching")
    p.add_argument("--de-id-column", default=None)
    p.add_argument("--logfc-column", default=None)
    p.add_argument("--padj-column", default=None)
    p.add_argument("--contrast-column", default=None)
    p.add_argument("--sheets", default="ALL", help="Comma-separated sheet list or ALL")
    p.add_argument("--padj-threshold", type=float, default=0.05)
    p.add_argument("--abs-logfc-threshold", type=float, default=1.0)
    p.add_argument("--expressed-padj-threshold", type=float, default=1.0, help="Rows present in DE table are considered tested; this is retained for summary compatibility")
    return p.parse_args()


def clean(value):
    value = str(value or "").strip().strip('"').strip("'")
    return "" if value in {"", ".", "NA", "NaN", "nan", "None", "none"} else value


def to_float(value):
    value = clean(value)
    if not value:
        return None
    try:
        return float(value.replace(",", ""))
    except Exception:
        return None


def detect_col(columns, explicit, candidates, label, required=True):
    cols = [str(c).strip() for c in columns]
    if explicit:
        for col in cols:
            if col == explicit:
                return col
        raise SystemExit(f"ERROR: requested {label} column {explicit!r} not found. Available: {', '.join(cols)}")
    lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    # regex fallback
    for col in cols:
        col_l = col.lower()
        for cand in candidates:
            if re.fullmatch(cand.lower(), col_l):
                return col
    if required:
        raise SystemExit(f"ERROR: could not detect {label}. Available: {', '.join(cols)}")
    return None


def read_candidates(path, candidate_id_column, extra_cols):
    with open(path, "r", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []
        rows = list(reader)
    if candidate_id_column not in headers:
        raise SystemExit(f"ERROR: candidate ID column {candidate_id_column!r} not found")
    id_cols = [candidate_id_column] + [c for c in extra_cols if c in headers]
    candidate_ids = defaultdict(set)
    for i, row in enumerate(rows):
        for col in id_cols:
            value = clean(row.get(col))
            if value:
                candidate_ids[value].add(i)
    return headers, rows, candidate_ids, id_cols


def selected_sheets(workbook, sheets_arg):
    if sheets_arg.upper() == "ALL":
        return workbook.sheet_names
    wanted = [s.strip() for s in sheets_arg.split(",") if s.strip()]
    missing = [s for s in wanted if s not in workbook.sheet_names]
    if missing:
        raise SystemExit(f"ERROR: requested sheet(s) not found: {', '.join(missing)}")
    return wanted


def classify_de(logfc, padj, padj_threshold, abs_logfc_threshold):
    if padj is None or logfc is None:
        return "tested_no_numeric_de_values"
    if padj <= padj_threshold and abs(logfc) >= abs_logfc_threshold:
        return "de_up" if logfc > 0 else "de_down"
    if padj <= padj_threshold:
        return "significant_padj_below_logfc_threshold"
    return "not_de"


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    extra_cols = [x.strip() for x in args.candidate_extra_id_columns.split(",") if x.strip()]
    cand_headers, cand_rows, candidate_ids, candidate_id_cols = read_candidates(args.candidates, args.candidate_id_column, extra_cols)

    xls = pd.ExcelFile(args.excel)
    sheets = selected_sheets(xls, args.sheets)

    candidate_hits = defaultdict(list)
    parse_summary = []
    all_de_rows = []

    for sheet in sheets:
        df = pd.read_excel(args.excel, sheet_name=sheet)
        df.columns = [str(c).strip() for c in df.columns]
        if df.empty:
            parse_summary.append({"sheet": sheet, "rows": 0, "matched_rows": 0, "status": "empty_sheet"})
            continue

        id_col = detect_col(
            df.columns,
            args.de_id_column,
            ["protein_id", "transcript_id", "gene_id", "id", "ID", "target_id", "Geneid", "gene", "transcript", "Unnamed: 0"],
            "DE ID",
        )
        logfc_col = detect_col(
            df.columns,
            args.logfc_column,
            ["log2FoldChange", "log2fc", "logFC", "lfc", "LFC", "log2_fold_change"],
            "logFC",
        )
        padj_col = detect_col(
            df.columns,
            args.padj_column,
            ["padj", "FDR", "fdr", "qvalue", "q.value", "adj.P.Val", "p_fdr", "p.adjust"],
            "adjusted P-value",
        )
        contrast_col = detect_col(df.columns, args.contrast_column, ["contrast", "comparison", "experiment"], "contrast", required=False)

        matched = 0
        for _, row in df.iterrows():
            de_id = clean(row.get(id_col))
            if not de_id or de_id not in candidate_ids:
                continue
            matched += 1
            logfc = to_float(row.get(logfc_col))
            padj = to_float(row.get(padj_col))
            contrast = clean(row.get(contrast_col)) if contrast_col else sheet
            status = classify_de(logfc, padj, args.padj_threshold, args.abs_logfc_threshold)
            record = {
                "de_id": de_id,
                "contrast": contrast,
                "sheet": sheet,
                "logfc": logfc,
                "padj": padj,
                "de_status": status,
                "logfc_column": logfc_col,
                "padj_column": padj_col,
                "id_column": id_col,
            }
            all_de_rows.append(record)
            for idx in candidate_ids[de_id]:
                candidate_hits[idx].append(record)

        parse_summary.append({
            "sheet": sheet,
            "rows": len(df),
            "matched_rows": matched,
            "id_column": id_col,
            "logfc_column": logfc_col,
            "padj_column": padj_col,
            "contrast_column": contrast_col or "sheet_name",
            "status": "parsed",
        })

    out_rows = []
    status_counts = Counter()
    contrast_counts = Counter()
    de_candidate_count = 0
    tested_candidate_count = 0

    extra_out = [
        "expression_de_tested_status",
        "expression_de_any_significant",
        "expression_de_any_up",
        "expression_de_any_down",
        "expression_de_n_records",
        "expression_de_n_significant_records",
        "expression_de_significant_contrasts",
        "expression_de_top_hit",
    ]
    out_headers = cand_headers + [x for x in extra_out if x not in cand_headers]

    for idx, row in enumerate(cand_rows):
        hits = candidate_hits.get(idx, [])
        out = dict(row)
        if not hits:
            out.update({
                "expression_de_tested_status": "not_found_in_de_workbook",
                "expression_de_any_significant": "no",
                "expression_de_any_up": "no",
                "expression_de_any_down": "no",
                "expression_de_n_records": 0,
                "expression_de_n_significant_records": 0,
                "expression_de_significant_contrasts": "NA",
                "expression_de_top_hit": "NA",
            })
            status_counts["not_found_in_de_workbook"] += 1
            out_rows.append(out)
            continue

        tested_candidate_count += 1
        sig = [h for h in hits if h["de_status"] in {"de_up", "de_down", "significant_padj_below_logfc_threshold"}]
        sig_directional = [h for h in hits if h["de_status"] in {"de_up", "de_down"}]
        if sig_directional:
            de_candidate_count += 1
        for h in sig_directional:
            contrast_counts[h["contrast"]] += 1

        top_candidates = [h for h in hits if h["padj"] is not None]
        top = sorted(top_candidates, key=lambda h: (h["padj"], -abs(h["logfc"] or 0)))[0] if top_candidates else hits[0]
        sig_contrasts = ";".join(f"{h['contrast']}:{h['de_status']}:logFC={h['logfc']}:padj={h['padj']}" for h in sig_directional) if sig_directional else "NA"
        status = "differentially_expressed" if sig_directional else "tested_not_differentially_expressed"
        status_counts[status] += 1
        out.update({
            "expression_de_tested_status": status,
            "expression_de_any_significant": "yes" if sig else "no",
            "expression_de_any_up": "yes" if any(h["de_status"] == "de_up" for h in hits) else "no",
            "expression_de_any_down": "yes" if any(h["de_status"] == "de_down" for h in hits) else "no",
            "expression_de_n_records": len(hits),
            "expression_de_n_significant_records": len(sig_directional),
            "expression_de_significant_contrasts": sig_contrasts,
            "expression_de_top_hit": f"{top['contrast']}:logFC={top['logfc']}:padj={top['padj']}:status={top['de_status']}",
        })
        out_rows.append(out)

    out_tsv = os.path.join(args.outdir, f"{args.prefix}.expression_de_context.tsv")
    long_tsv = os.path.join(args.outdir, f"{args.prefix}.expression_de_matches.long.tsv")
    summary_txt = os.path.join(args.outdir, f"{args.prefix}.expression_de_context.summary.txt")

    with open(out_tsv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=out_headers, extrasaction="ignore")
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)

    long_headers = ["de_id", "contrast", "sheet", "logfc", "padj", "de_status", "id_column", "logfc_column", "padj_column"]
    with open(long_tsv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=long_headers)
        writer.writeheader()
        for row in all_de_rows:
            writer.writerow(row)

    with open(summary_txt, "w", encoding="utf-8") as summary:
        summary.write("Dark candidate differential-expression context summary\n")
        summary.write("======================================================\n")
        summary.write(f"Candidate TSV: {args.candidates}\n")
        summary.write(f"DE workbook: {args.excel}\n")
        summary.write(f"Sheets analysed: {','.join(sheets)}\n")
        summary.write(f"Candidate rows read: {len(cand_rows)}\n")
        summary.write(f"Candidate ID columns used: {','.join(candidate_id_cols)}\n")
        summary.write(f"padj threshold: {args.padj_threshold}\n")
        summary.write(f"abs logFC threshold: {args.abs_logfc_threshold}\n")
        summary.write("\nSheet parse summary:\n")
        for item in parse_summary:
            summary.write("- " + ", ".join(f"{k}={v}" for k, v in item.items()) + "\n")
        summary.write("\nCandidate DE status counts:\n")
        for key, value in status_counts.most_common():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nDirectional DE candidate counts by contrast:\n")
        for key, value in contrast_counts.most_common():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nOutputs:\n")
        summary.write(f"- Candidate DE context table: {out_tsv}\n")
        summary.write(f"- Long matched DE table: {long_tsv}\n")
        summary.write(f"- Summary: {summary_txt}\n")

    print("Done.")
    print(f"Summary: {summary_txt}")
    print(f"Candidate DE context table: {out_tsv}")
    print(f"Long matched DE table: {long_tsv}")


if __name__ == "__main__":
    main()
