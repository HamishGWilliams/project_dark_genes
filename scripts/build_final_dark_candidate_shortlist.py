#!/usr/bin/env python3
"""Build final integrated and shortlist tables for Project Dark Genes.

Inputs:
  --de-context : candidate-level table after annotation/genome/repeat/DMR/DE integration
  --sig-long   : long table of significant DE dark-candidate records

Outputs:
  <prefix>.final_integrated.tsv
  <prefix>.final_shortlist.tsv
  <prefix>.final_summary.txt

The shortlist score is deliberately transparent and conservative. It rewards:
  - BUSCO-backed high priority tier
  - non-assembly-redundant duplication interpretation
  - DE significance in focal stress contrasts
  - diesel-added significance specifically
  - strong absolute logFC / adjusted p-value evidence
  - no repeat/TE overlap as a modest confidence bonus

Repeat/TE overlap is retained as biological context, not used as an automatic exclusion.
"""

import argparse
import csv
import math
import os
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
    p = argparse.ArgumentParser(description="Build final integrated dark-candidate shortlist.")
    p.add_argument("--de-context", required=True, help="Candidate-level DE context TSV")
    p.add_argument("--sig-long", required=True, help="Significant DE long TSV")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--prefix", default="equina_dark_candidates")
    p.add_argument("--candidate-id-column", default="protein_id")
    p.add_argument(
        "--focal-contrasts",
        default="diesel_added_wald,combined_wald,full_model_LRT,interactive_only_wald,interactive_only_LRT",
        help="Comma-separated contrasts considered focal stress-response evidence.",
    )
    p.add_argument("--top-n", type=int, default=250, help="Number of rows to write to the shortlist table.")
    return p.parse_args()


def clean(value):
    value = str(value or "").strip()
    return "" if value in {"", ".", "NA", "None", "none", "nan", "NaN"} else value


def as_float(value):
    value = clean(value)
    if not value:
        return None
    try:
        x = float(value.replace(",", ""))
        if math.isnan(x):
            return None
        return x
    except Exception:
        return None


def read_tsv(path):
    with open(path, "r", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        headers = reader.fieldnames or []
    return headers, rows


def write_tsv(path, rows, headers):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=headers, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def detect_col(headers, candidates):
    lower = {h.lower(): h for h in headers}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def candidate_id(row, preferred):
    for col in [preferred, "protein_id", "candidate_id", "transcript_id", "gene_id", "source_fasta_id"]:
        if col and col in row and clean(row.get(col)):
            return clean(row.get(col))
    return ""


def bool_is_yes(value):
    return clean(value).lower() in {"yes", "true", "1", "de_significant"}


def load_sig_records(path, focal_contrasts):
    headers, rows = read_tsv(path)
    by_candidate = defaultdict(list)
    for row in rows:
        cid = clean(row.get("candidate_id"))
        if not cid:
            continue
        contrast = clean(row.get("contrast"))
        logfc = as_float(row.get("logfc"))
        padj = as_float(row.get("padj"))
        direction = clean(row.get("direction")) or ("up" if logfc and logfc > 0 else "down" if logfc and logfc < 0 else "NA")
        rec = {
            "contrast": contrast,
            "direction": direction,
            "logfc": logfc,
            "padj": padj,
            "is_focal": contrast in focal_contrasts,
            "is_diesel_added": contrast == "diesel_added_wald",
        }
        by_candidate[cid].append(rec)
    return by_candidate


def score_candidate(row, records, cols, focal_contrasts):
    score = 0
    reasons = []

    priority = clean(row.get(cols.get("priority") or ""))
    duplication = clean(row.get(cols.get("duplication") or ""))
    repeat_status = clean(row.get(cols.get("repeat") or ""))
    de_status = clean(row.get(cols.get("de_status") or ""))

    if priority == "high_priority":
        score += 25
        reasons.append("high_priority")
    elif priority == "medium_priority":
        score += 10
        reasons.append("medium_priority")
    elif priority:
        reasons.append(priority)

    if duplication == "single_copy_or_no_near_identical_dark_duplicate":
        score += 20
        reasons.append("single_copy_or_no_near_identical_dark_duplicate")
    elif duplication == "possible_biological_gene_family_expansion":
        score += 12
        reasons.append("possible_biological_gene_family_expansion")
    elif duplication == "possible_assembly_redundancy_or_haplotig_duplication":
        score -= 20
        reasons.append("assembly_redundancy_or_haplotig_caution")

    if de_status == "de_significant" or records:
        score += 20
        reasons.append("DE_significant_any_contrast")

    focal_records = [r for r in records if r["is_focal"]]
    diesel_records = [r for r in records if r["is_diesel_added"]]

    if focal_records:
        score += 20
        reasons.append("focal_stress_DE")
        score += min(10, 2 * len(set(r["contrast"] for r in focal_records)))

    if diesel_records:
        score += 15
        reasons.append("diesel_added_DE")

    if records:
        best_abs_lfc = max(abs(r["logfc"] or 0.0) for r in records)
        best_padj = min((r["padj"] for r in records if r["padj"] is not None), default=None)
        if best_abs_lfc >= 2:
            score += 6
            reasons.append("abs_logFC_ge_2")
        elif best_abs_lfc >= 1.5:
            score += 3
            reasons.append("abs_logFC_ge_1.5")
        if best_padj is not None and best_padj <= 0.001:
            score += 6
            reasons.append("padj_le_0.001")
        elif best_padj is not None and best_padj <= 0.01:
            score += 3
            reasons.append("padj_le_0.01")

    if repeat_status == "no_repeat_overlap":
        score += 5
        reasons.append("no_repeat_TE_overlap")
    elif repeat_status == "repeat_overlap":
        reasons.append("repeat_TE_overlap_context")

    # DMR overlap was negative in both experiments. Capture if columns are present.
    dmr_status_values = [clean(row.get(col)) for col in cols.get("dmr_cols", []) if clean(row.get(col))]
    if dmr_status_values and all(v == "no_dmr_overlap" for v in dmr_status_values):
        score += 3
        reasons.append("no_DMR_overlap")

    return score, reasons


def summarise_records(records):
    if not records:
        return {
            "final_focal_stress_responsive": "no",
            "final_diesel_added_responsive": "no",
            "final_n_significant_contrasts": 0,
            "final_n_focal_significant_contrasts": 0,
            "final_best_abs_logfc": "NA",
            "final_best_padj": "NA",
            "final_significant_contrast_summary": "NA",
        }

    contrasts = sorted(set(r["contrast"] for r in records))
    focal_contrasts = sorted(set(r["contrast"] for r in records if r["is_focal"]))
    best_abs_lfc = max(abs(r["logfc"] or 0.0) for r in records)
    best_padj = min((r["padj"] for r in records if r["padj"] is not None), default=None)
    sig_summary = ";".join(
        f"{r['contrast']}|{r['direction']}|logFC={r['logfc']:.4g}|padj={r['padj']:.4g}"
        for r in sorted(records, key=lambda x: (x["padj"] if x["padj"] is not None else 1.0, -abs(x["logfc"] or 0.0)))
    )
    return {
        "final_focal_stress_responsive": "yes" if focal_contrasts else "no",
        "final_diesel_added_responsive": "yes" if any(r["is_diesel_added"] for r in records) else "no",
        "final_n_significant_contrasts": len(contrasts),
        "final_n_focal_significant_contrasts": len(focal_contrasts),
        "final_best_abs_logfc": f"{best_abs_lfc:.6g}",
        "final_best_padj": f"{best_padj:.6g}" if best_padj is not None else "NA",
        "final_significant_contrast_summary": sig_summary,
    }


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    focal_contrasts = [x.strip() for x in args.focal_contrasts.split(",") if x.strip()]
    headers, rows = read_tsv(args.de_context)
    sig_records = load_sig_records(args.sig_long, set(focal_contrasts))

    cols = {
        "priority": detect_col(headers, ["priority_tier"]),
        "duplication": detect_col(headers, ["duplication_interpretation", "duplication_status"]),
        "repeat": detect_col(headers, ["repeat_overlap_status"]),
        "de_status": detect_col(headers, ["de_match_status"]),
        "dmr_cols": [h for h in headers if h.endswith("dmr_overlap_status") or h == "dmr_overlap_status"],
    }

    out_rows = []
    status_counts = Counter()
    priority_counts = Counter()
    shortlist_reason_counts = Counter()

    for row in rows:
        cid = candidate_id(row, args.candidate_id_column)
        records = sig_records.get(cid, [])
        # If candidate_id differs from protein_id, try transcript/gene IDs as a fallback.
        if not records:
            for alt_col in ["transcript_id", "gene_id", "source_fasta_id"]:
                alt = clean(row.get(alt_col))
                if alt and alt in sig_records:
                    records = sig_records[alt]
                    break

        score, reasons = score_candidate(row, records, cols, set(focal_contrasts))
        rec_summary = summarise_records(records)

        final_category = "retain_context_only"
        if score >= 75 and rec_summary["final_focal_stress_responsive"] == "yes":
            final_category = "top_stress_responsive_dark_candidate"
        elif score >= 55 and rec_summary["final_focal_stress_responsive"] == "yes":
            final_category = "strong_stress_responsive_dark_candidate"
        elif rec_summary["final_focal_stress_responsive"] == "yes":
            final_category = "stress_responsive_dark_candidate"
        elif records:
            final_category = "DE_dark_candidate_non_focal_contrast"

        out = dict(row)
        out.update(rec_summary)
        out.update({
            "final_candidate_id": cid,
            "final_shortlist_score": score,
            "final_shortlist_category": final_category,
            "final_shortlist_reasons": ";".join(reasons) if reasons else "NA",
        })
        out_rows.append(out)
        status_counts[final_category] += 1
        priority_counts[clean(row.get(cols.get("priority") or "")) or "NA"] += 1
        for reason in reasons:
            shortlist_reason_counts[reason] += 1

    out_rows.sort(key=lambda r: (
        -int(float(r.get("final_shortlist_score", 0))),
        as_float(r.get("final_best_padj")) if as_float(r.get("final_best_padj")) is not None else 1.0,
        -as_float(r.get("final_best_abs_logfc")) if as_float(r.get("final_best_abs_logfc")) is not None else 0.0,
        r.get("final_candidate_id", ""),
    ))

    extra_headers = [
        "final_candidate_id",
        "final_shortlist_score",
        "final_shortlist_category",
        "final_shortlist_reasons",
        "final_focal_stress_responsive",
        "final_diesel_added_responsive",
        "final_n_significant_contrasts",
        "final_n_focal_significant_contrasts",
        "final_best_abs_logfc",
        "final_best_padj",
        "final_significant_contrast_summary",
    ]
    out_headers = extra_headers + [h for h in headers if h not in extra_headers]

    integrated_path = os.path.join(args.outdir, f"{args.prefix}.final_integrated.tsv")
    shortlist_path = os.path.join(args.outdir, f"{args.prefix}.final_shortlist.tsv")
    summary_path = os.path.join(args.outdir, f"{args.prefix}.final_summary.txt")

    write_tsv(integrated_path, out_rows, out_headers)
    write_tsv(shortlist_path, out_rows[: args.top_n], out_headers)

    with open(summary_path, "w", encoding="utf-8") as summary:
        summary.write("Final dark-candidate synthesis summary\n")
        summary.write("======================================\n")
        summary.write(f"DE context table: {args.de_context}\n")
        summary.write(f"Significant DE long table: {args.sig_long}\n")
        summary.write(f"Output directory: {args.outdir}\n")
        summary.write(f"Candidate rows read: {len(rows)}\n")
        summary.write(f"Significant-DE candidate IDs indexed: {len(sig_records)}\n")
        summary.write(f"Shortlist rows written: {min(args.top_n, len(out_rows))}\n")
        summary.write(f"Focal contrasts: {','.join(focal_contrasts)}\n")
        summary.write("\nFinal shortlist category counts:\n")
        for key, value in status_counts.most_common():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nPriority tier counts in input:\n")
        for key, value in priority_counts.most_common():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nMost common shortlist evidence reasons:\n")
        for key, value in shortlist_reason_counts.most_common(30):
            summary.write(f"- {key}: {value}\n")
        summary.write("\nOutputs:\n")
        summary.write(f"- Final integrated table: {integrated_path}\n")
        summary.write(f"- Final shortlist table: {shortlist_path}\n")
        summary.write(f"- Summary: {summary_path}\n")

    print("Done.")
    print(f"Summary: {summary_path}")
    print(f"Final integrated table: {integrated_path}")
    print(f"Final shortlist table: {shortlist_path}")


if __name__ == "__main__":
    main()
