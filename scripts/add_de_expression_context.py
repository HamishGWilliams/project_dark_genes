#!/usr/bin/env python3
"""Integrate RNA-seq differential-expression evidence with dark candidates.

The DE input may be either:
  1. Wide multi-contrast CSV/TSV, with one row per gene/transcript and many
     logFC / padj columns, or
  2. Long format, with columns for ID, contrast, logFC, and padj.

Outputs:
  <prefix>.de_context.tsv
  <prefix>.de_significant_long.tsv
  <prefix>.de_context.summary.txt

The script is intentionally conservative: a candidate is called differentially
expressed only when padj <= --padj-threshold and abs(logFC) >= --lfc-threshold.
"""

import argparse
import csv
import math
import os
import re
import sys
from collections import Counter, defaultdict

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
    p = argparse.ArgumentParser(description="Add DE/RNA-seq context to dark candidates.")
    p.add_argument("--candidates", required=True)
    p.add_argument("--de", required=True, help="DE CSV/TSV table")
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", default="equina_dark_candidates")
    p.add_argument("--candidate-id-columns", default="protein_id,transcript_id,gene_id,source_fasta_id")
    p.add_argument("--de-id-column", default=None)
    p.add_argument("--de-format", choices=["auto", "wide", "long"], default="auto")
    p.add_argument("--contrast-column", default=None)
    p.add_argument("--lfc-column", default=None)
    p.add_argument("--padj-column", default=None)
    p.add_argument("--padj-threshold", type=float, default=0.05)
    p.add_argument("--lfc-threshold", type=float, default=1.0)
    p.add_argument("--delimiter", choices=["auto", "comma", "tab"], default="auto")
    return p.parse_args()


def clean(x):
    x = str(x or "").strip().strip('"').strip("'")
    return "" if x in {"", ".", "NA", "N/A", "none", "None", "nan", "NaN"} else x


def as_float(x):
    x = clean(x)
    if not x:
        return None
    try:
        value = float(x.replace(",", ""))
        if math.isnan(value):
            return None
        return value
    except Exception:
        return None


def detect_delimiter(path, forced="auto"):
    if forced == "comma":
        return ","
    if forced == "tab":
        return "\t"
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                return "," if line.count(",") > line.count("\t") else "\t"
    return ","


def read_table(path, delimiter):
    with open(path, "r", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = list(reader)
        headers = reader.fieldnames or []
    return headers, rows


def normalise_id(value):
    value = clean(value)
    if not value:
        return ""
    # Keep the exact ID first elsewhere; this normalisation helps common FASTA/GFF suffix differences.
    value = value.split()[0]
    value = re.sub(r"^(gene:|transcript:|protein:|cds:)", "", value, flags=re.IGNORECASE)
    return value


def detect_de_id_column(headers, explicit=None):
    if explicit:
        if explicit not in headers:
            raise SystemExit(f"ERROR: requested DE ID column {explicit!r} not found")
        return explicit
    lower = {h.lower(): h for h in headers}
    for cand in [
        "protein_id", "transcript_id", "gene_id", "id", "ID", "Name", "target_id",
        "feature_id", "row", "rownames", "seqid", "locus", "gene", "transcript",
    ]:
        if cand.lower() in lower:
            return lower[cand.lower()]
    # R/CSV exports sometimes write row names as an empty first column or X.
    for cand in ["", "X", "x", "...1"]:
        if cand in headers:
            return cand
    raise SystemExit(f"ERROR: could not detect DE ID column. Available columns: {', '.join(headers)}")


def classify_lfc_col(name):
    n = name.lower()
    return any(token in n for token in ["log2foldchange", "log2fc", "logfc", "lfc"])


def classify_padj_col(name):
    n = name.lower()
    return any(token in n for token in ["padj", "fdr", "qvalue", "q_value", "adj.p", "adj_p", "p_adj"])


def strip_stat_tokens(name):
    n = name
    replacements = [
        "log2FoldChange", "log2foldchange", "log2FC", "log2fc", "logFC", "logfc", "lfc", "LFC",
        "padj", "FDR", "fdr", "qvalue", "q_value", "adj.P.Val", "adj_p", "p_adj", "qval",
    ]
    for r in replacements:
        n = n.replace(r, "")
    n = re.sub(r"(^[._:-]+|[._:-]+$)", "", n)
    n = re.sub(r"[._:-]+", "_", n)
    return n or name


def detect_wide_pairs(headers):
    lfc_cols = [h for h in headers if classify_lfc_col(h)]
    padj_cols = [h for h in headers if classify_padj_col(h)]
    pairs = []

    used_padj = set()
    for lfc in lfc_cols:
        lbase = strip_stat_tokens(lfc).lower()
        best = None
        best_score = -1
        for padj in padj_cols:
            if padj in used_padj:
                continue
            pbase = strip_stat_tokens(padj).lower()
            score = 0
            if lbase == pbase:
                score = 100
            elif lbase and pbase and (lbase in pbase or pbase in lbase):
                score = 50
            else:
                # Token overlap fallback.
                ltoks = set(lbase.split("_"))
                ptoks = set(pbase.split("_"))
                score = len(ltoks & ptoks)
            if score > best_score:
                best = padj
                best_score = score
        if best is not None:
            used_padj.add(best)
            contrast = strip_stat_tokens(lfc)
            pairs.append((contrast, lfc, best))

    # If there is one generic LFC and one generic padj, include that too.
    if not pairs and len(lfc_cols) == 1 and len(padj_cols) == 1:
        pairs.append(("contrast", lfc_cols[0], padj_cols[0]))

    return pairs


def detect_long_columns(headers, args):
    lower = {h.lower(): h for h in headers}
    def pick(explicit, options, label):
        if explicit:
            if explicit not in headers:
                raise SystemExit(f"ERROR: requested {label} column {explicit!r} not found")
            return explicit
        for opt in options:
            if opt.lower() in lower:
                return lower[opt.lower()]
        raise SystemExit(f"ERROR: could not detect {label} column. Available columns: {', '.join(headers)}")
    contrast = pick(args.contrast_column, ["contrast", "comparison", "condition", "experiment"], "contrast")
    lfc = pick(args.lfc_column, ["log2FoldChange", "log2FC", "logFC", "lfc"], "logFC")
    padj = pick(args.padj_column, ["padj", "FDR", "qvalue", "q_value", "adj_p"], "padj/FDR")
    return contrast, lfc, padj


def build_de_index(headers, rows, args):
    de_id_col = detect_de_id_column(headers, args.de_id_column)

    wide_pairs = detect_wide_pairs(headers)
    use_long = args.de_format == "long" or (args.de_format == "auto" and not wide_pairs and args.contrast_column)

    index = defaultdict(list)
    parse_counts = Counter()
    detected = {"de_id_column": de_id_col}

    if use_long:
        contrast_col, lfc_col, padj_col = detect_long_columns(headers, args)
        detected.update({"de_format_used": "long", "contrast_column": contrast_col, "lfc_column": lfc_col, "padj_column": padj_col})
        for row in rows:
            rid = normalise_id(row.get(de_id_col))
            if not rid:
                parse_counts["missing_id"] += 1
                continue
            lfc = as_float(row.get(lfc_col))
            padj = as_float(row.get(padj_col))
            if lfc is None or padj is None:
                parse_counts["missing_lfc_or_padj"] += 1
                continue
            index[rid].append({"contrast": clean(row.get(contrast_col)) or "contrast", "logfc": lfc, "padj": padj})
            parse_counts["parsed_de_record"] += 1
    else:
        if args.lfc_column and args.padj_column:
            wide_pairs = [("contrast", args.lfc_column, args.padj_column)]
        if not wide_pairs:
            raise SystemExit("ERROR: could not detect any logFC/padj column pairs in wide DE table")
        detected.update({"de_format_used": "wide", "wide_pairs_detected": len(wide_pairs)})
        detected["wide_pairs"] = ";".join(f"{c}:{l}:{p}" for c, l, p in wide_pairs[:50])
        for row in rows:
            rid = normalise_id(row.get(de_id_col))
            if not rid:
                parse_counts["missing_id"] += 1
                continue
            for contrast, lfc_col, padj_col in wide_pairs:
                lfc = as_float(row.get(lfc_col))
                padj = as_float(row.get(padj_col))
                if lfc is None or padj is None:
                    parse_counts["missing_lfc_or_padj"] += 1
                    continue
                index[rid].append({"contrast": contrast, "logfc": lfc, "padj": padj})
                parse_counts["parsed_de_record"] += 1

    return index, parse_counts, detected


def find_candidate_de_records(row, candidate_cols, de_index):
    tried = []
    for col in candidate_cols:
        if col not in row:
            continue
        raw = clean(row.get(col))
        if not raw:
            continue
        ids = [normalise_id(raw)]
        # Try common transcript/protein suffix variants.
        if ".p" in raw:
            ids.append(normalise_id(raw.split(".p")[0]))
        if ".t" in raw:
            ids.append(normalise_id(raw.split(".t")[0]))
        for rid in ids:
            if not rid or rid in tried:
                continue
            tried.append(rid)
            if rid in de_index:
                return rid, col, de_index[rid]
    return "", "", []


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    cand_headers, cand_rows = read_table(args.candidates, "\t")
    de_delim = detect_delimiter(args.de, args.delimiter)
    de_headers, de_rows = read_table(args.de, de_delim)

    candidate_cols = [x.strip() for x in args.candidate_id_columns.split(",") if x.strip()]
    available_candidate_cols = [c for c in candidate_cols if c in cand_headers]
    if not available_candidate_cols:
        raise SystemExit(f"ERROR: none of candidate ID columns are present: {candidate_cols}")

    de_index, de_parse_counts, detected = build_de_index(de_headers, de_rows, args)
    detected["de_delimiter"] = "comma" if de_delim == "," else "tab"

    extra = [
        "de_match_status", "de_matched_id", "de_matched_candidate_column",
        "de_tested_contrast_count", "de_any_significant", "de_significant_contrast_count",
        "de_significant_contrasts", "de_min_padj", "de_max_abs_logfc",
    ]
    out_headers = cand_headers + [x for x in extra if x not in cand_headers]
    out_rows = []
    sig_long_rows = []
    status_counts = Counter()
    contrast_sig_counts = Counter()

    for row in cand_rows:
        out = dict(row)
        matched_id, matched_col, records = find_candidate_de_records(row, available_candidate_cols, de_index)
        if not records:
            status = "no_de_record_matched"
            out.update({
                "de_match_status": status,
                "de_matched_id": "NA",
                "de_matched_candidate_column": "NA",
                "de_tested_contrast_count": 0,
                "de_any_significant": "no",
                "de_significant_contrast_count": 0,
                "de_significant_contrasts": "NA",
                "de_min_padj": "NA",
                "de_max_abs_logfc": "NA",
            })
            status_counts[status] += 1
            out_rows.append(out)
            continue

        significant = []
        min_padj = None
        max_abs_lfc = 0.0
        for rec in records:
            lfc = rec["logfc"]
            padj = rec["padj"]
            min_padj = padj if min_padj is None else min(min_padj, padj)
            max_abs_lfc = max(max_abs_lfc, abs(lfc))
            if padj <= args.padj_threshold and abs(lfc) >= args.lfc_threshold:
                direction = "up" if lfc > 0 else "down" if lfc < 0 else "zero"
                sig = {"contrast": rec["contrast"], "logfc": lfc, "padj": padj, "direction": direction}
                significant.append(sig)
                contrast_sig_counts[rec["contrast"]] += 1

        status = "de_significant" if significant else "de_record_matched_not_significant"
        status_counts[status] += 1
        sig_text = ";".join(f"{s['contrast']}|{s['direction']}|logFC={s['logfc']:.4g}|padj={s['padj']:.4g}" for s in significant) if significant else "NA"
        out.update({
            "de_match_status": status,
            "de_matched_id": matched_id,
            "de_matched_candidate_column": matched_col,
            "de_tested_contrast_count": len(records),
            "de_any_significant": "yes" if significant else "no",
            "de_significant_contrast_count": len(significant),
            "de_significant_contrasts": sig_text,
            "de_min_padj": f"{min_padj:.6g}" if min_padj is not None else "NA",
            "de_max_abs_logfc": f"{max_abs_lfc:.6g}",
        })
        out_rows.append(out)

        candidate_id = clean(row.get("protein_id")) or clean(row.get("transcript_id")) or matched_id
        for sig in significant:
            sig_long_rows.append({
                "candidate_id": candidate_id,
                "matched_de_id": matched_id,
                "matched_candidate_column": matched_col,
                "contrast": sig["contrast"],
                "direction": sig["direction"],
                "logfc": sig["logfc"],
                "padj": sig["padj"],
            })

    out_tsv = os.path.join(args.outdir, f"{args.prefix}.de_context.tsv")
    sig_tsv = os.path.join(args.outdir, f"{args.prefix}.de_significant_long.tsv")
    summary_txt = os.path.join(args.outdir, f"{args.prefix}.de_context.summary.txt")

    with open(out_tsv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=out_headers, extrasaction="ignore")
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)

    sig_headers = ["candidate_id", "matched_de_id", "matched_candidate_column", "contrast", "direction", "logfc", "padj"]
    with open(sig_tsv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=sig_headers)
        writer.writeheader()
        for row in sig_long_rows:
            writer.writerow(row)

    with open(summary_txt, "w", encoding="utf-8") as summary:
        summary.write("Dark candidate differential-expression context summary\n")
        summary.write("======================================================\n")
        summary.write(f"Candidate TSV: {args.candidates}\n")
        summary.write(f"DE table: {args.de}\n")
        summary.write(f"Output directory: {args.outdir}\n")
        summary.write(f"Candidate rows read: {len(cand_rows)}\n")
        summary.write(f"DE rows read: {len(de_rows)}\n")
        summary.write(f"DE IDs indexed: {len(de_index)}\n")
        summary.write(f"padj threshold: {args.padj_threshold}\n")
        summary.write(f"abs(logFC) threshold: {args.lfc_threshold}\n")
        summary.write("\nDetected DE settings:\n")
        for key, value in detected.items():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nDE parse counts:\n")
        for key, value in de_parse_counts.most_common():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nCandidate DE status counts:\n")
        for key, value in status_counts.most_common():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nSignificant dark-candidate counts by contrast:\n")
        for key, value in contrast_sig_counts.most_common():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nOutputs:\n")
        summary.write(f"- Candidate DE-context table: {out_tsv}\n")
        summary.write(f"- Significant DE long table: {sig_tsv}\n")
        summary.write(f"- Summary: {summary_txt}\n")

    print("Done.")
    print(f"Summary: {summary_txt}")
    print(f"DE-context table: {out_tsv}")
    print(f"Significant long table: {sig_tsv}")


if __name__ == "__main__":
    main()
