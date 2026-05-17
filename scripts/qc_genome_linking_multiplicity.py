#!/usr/bin/env python3
"""
QC genome-linked dark-candidate mappings and collapse inflated multi-match output.

The genome lookup step may produce multiple GFF3 matches per dark-candidate ID when
protein/transcript aliases are shared, truncated, or too permissive. This script
summarises the multiplicity and writes a deterministic one-row-per-candidate
primary mapping table for downstream figure generation and prioritisation.

It does not delete the full multi-match table. It creates:
  - equina_dark_candidates.genome_linked.primary.tsv
  - equina_dark_candidates.genome_linking_multiplicity.tsv
  - equina_dark_candidates.genome_linking_multiplicity.summary.txt
"""

import argparse
import csv
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
    p = argparse.ArgumentParser(
        description="Summarise genome-linking multiplicity and write a primary one-row-per-candidate mapping table."
    )
    p.add_argument("--linked", required=True, help="Genome-linked candidate TSV from build_genome_lookup.py")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--candidate-id-column", default="protein_id", help="Candidate ID column in the linked TSV")
    p.add_argument("--prefix", default="equina_dark_candidates", help="Output filename prefix")
    return p.parse_args()


def clean(value):
    value = str(value or "").strip()
    return value if value and value not in {".", "NA", "None", "none", "nan"} else ""


def as_int(value):
    try:
        return int(value)
    except Exception:
        return None


def rank_row(row, candidate_id_column):
    """Lower rank is better. Designed for deterministic primary match selection."""
    candidate_id = clean(row.get(candidate_id_column))
    matched_alias = clean(row.get("genome_matched_by_alias"))
    matched_col = clean(row.get("genome_matched_by_column"))
    gff3_protein_id = clean(row.get("gff3_protein_id"))
    transcript_id = clean(row.get("transcript_id"))
    gene_id = clean(row.get("gene_id"))

    if candidate_id and gff3_protein_id and candidate_id == gff3_protein_id:
        primary_rank = 0
        primary_basis = "candidate_id_equals_gff3_protein_id"
    elif candidate_id and transcript_id and candidate_id == transcript_id:
        primary_rank = 1
        primary_basis = "candidate_id_equals_transcript_id"
    elif candidate_id and gene_id and candidate_id == gene_id:
        primary_rank = 2
        primary_basis = "candidate_id_equals_gene_id"
    elif candidate_id and matched_alias and candidate_id == matched_alias:
        primary_rank = 3
        primary_basis = "candidate_id_equals_matched_alias"
    elif matched_col and matched_col == candidate_id_column:
        primary_rank = 4
        primary_basis = "matched_by_candidate_id_column"
    else:
        primary_rank = 10
        primary_basis = "fallback_deterministic_choice"

    # Tie-breakers prefer rows with usable CDS/transcript coordinates, then shorter transcript span,
    # then stable lexical IDs. This avoids non-deterministic output when many aliases tie.
    cds_total_bp = as_int(row.get("cds_total_bp"))
    transcript_span_bp = as_int(row.get("transcript_span_bp"))
    has_cds = 0 if cds_total_bp is not None and cds_total_bp > 0 else 1
    span = transcript_span_bp if transcript_span_bp is not None else 10**18

    tie_key = (
        primary_rank,
        has_cds,
        span,
        clean(row.get("contig")),
        as_int(row.get("transcript_start")) if as_int(row.get("transcript_start")) is not None else 10**18,
        transcript_id,
        gff3_protein_id,
    )

    return tie_key, primary_rank, primary_basis


def write_tsv(path, rows, headers):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=headers, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if not os.path.exists(args.linked):
        raise SystemExit(f"ERROR: linked TSV not found: {args.linked}")

    candidate_stats = defaultdict(lambda: {
        "rows": 0,
        "genes": set(),
        "transcripts": set(),
        "gff3_proteins": set(),
        "contigs": set(),
        "raw_statuses": Counter(),
        "matched_columns": Counter(),
        "matched_aliases": Counter(),
    })

    best_rows = {}
    best_keys = {}
    best_ranks = {}
    best_bases = {}
    headers = None
    total_rows = 0

    with open(args.linked, "r", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames or []

        if args.candidate_id_column not in headers:
            raise SystemExit(
                f"ERROR: candidate column {args.candidate_id_column!r} not present. "
                f"Available columns: {', '.join(headers)}"
            )

        for row in reader:
            total_rows += 1
            candidate_id = clean(row.get(args.candidate_id_column))
            if not candidate_id:
                candidate_id = f"__missing_candidate_id_row_{total_rows}"

            stats = candidate_stats[candidate_id]
            stats["rows"] += 1

            for field, target in [
                ("gene_id", "genes"),
                ("transcript_id", "transcripts"),
                ("gff3_protein_id", "gff3_proteins"),
                ("contig", "contigs"),
            ]:
                value = clean(row.get(field))
                if value:
                    stats[target].add(value)

            stats["raw_statuses"][clean(row.get("genome_match_status")) or "NA"] += 1
            stats["matched_columns"][clean(row.get("genome_matched_by_column")) or "NA"] += 1
            stats["matched_aliases"][clean(row.get("genome_matched_by_alias")) or "NA"] += 1

            key, rank, basis = rank_row(row, args.candidate_id_column)
            if candidate_id not in best_keys or key < best_keys[candidate_id]:
                best_keys[candidate_id] = key
                best_rows[candidate_id] = dict(row)
                best_ranks[candidate_id] = rank
                best_bases[candidate_id] = basis

    if headers is None:
        raise SystemExit(f"ERROR: no header found in {args.linked}")

    primary_headers = list(headers)
    for extra in [
        "primary_mapping_status",
        "primary_mapping_rank",
        "primary_mapping_basis",
        "candidate_total_genome_matches",
        "candidate_n_matched_genes",
        "candidate_n_matched_transcripts",
        "candidate_n_matched_gff3_proteins",
        "candidate_n_matched_contigs",
    ]:
        if extra not in primary_headers:
            primary_headers.append(extra)

    primary_rows = []
    multiplicity_rows = []
    primary_status_counts = Counter()
    rank_counts = Counter()

    for candidate_id in sorted(candidate_stats):
        stats = candidate_stats[candidate_id]
        best = dict(best_rows[candidate_id])
        rank = best_ranks[candidate_id]
        basis = best_bases[candidate_id]

        if stats["rows"] == 1:
            primary_status = "matched_unique"
        elif rank <= 3:
            primary_status = "matched_primary_from_multiple_strict"
        elif rank == 4:
            primary_status = "matched_primary_from_multiple_candidate_column"
        else:
            primary_status = "ambiguous_primary_from_multiple"

        primary_status_counts[primary_status] += 1
        rank_counts[str(rank)] += 1

        best.update({
            "primary_mapping_status": primary_status,
            "primary_mapping_rank": rank,
            "primary_mapping_basis": basis,
            "candidate_total_genome_matches": stats["rows"],
            "candidate_n_matched_genes": len(stats["genes"]),
            "candidate_n_matched_transcripts": len(stats["transcripts"]),
            "candidate_n_matched_gff3_proteins": len(stats["gff3_proteins"]),
            "candidate_n_matched_contigs": len(stats["contigs"]),
        })
        primary_rows.append(best)

        top_aliases = ";".join(f"{k}:{v}" for k, v in stats["matched_aliases"].most_common(10))
        top_columns = ";".join(f"{k}:{v}" for k, v in stats["matched_columns"].most_common(10))
        raw_statuses = ";".join(f"{k}:{v}" for k, v in stats["raw_statuses"].most_common())

        multiplicity_rows.append({
            "candidate_id": candidate_id,
            "n_linked_rows": stats["rows"],
            "n_matched_genes": len(stats["genes"]),
            "n_matched_transcripts": len(stats["transcripts"]),
            "n_matched_gff3_proteins": len(stats["gff3_proteins"]),
            "n_matched_contigs": len(stats["contigs"]),
            "primary_mapping_status": primary_status,
            "primary_mapping_rank": rank,
            "primary_mapping_basis": basis,
            "raw_genome_match_statuses": raw_statuses,
            "top_matched_by_columns": top_columns,
            "top_matched_by_aliases": top_aliases,
        })

    primary_path = os.path.join(args.outdir, f"{args.prefix}.genome_linked.primary.tsv")
    multiplicity_path = os.path.join(args.outdir, f"{args.prefix}.genome_linking_multiplicity.tsv")
    summary_path = os.path.join(args.outdir, f"{args.prefix}.genome_linking_multiplicity.summary.txt")

    multiplicity_headers = [
        "candidate_id",
        "n_linked_rows",
        "n_matched_genes",
        "n_matched_transcripts",
        "n_matched_gff3_proteins",
        "n_matched_contigs",
        "primary_mapping_status",
        "primary_mapping_rank",
        "primary_mapping_basis",
        "raw_genome_match_statuses",
        "top_matched_by_columns",
        "top_matched_by_aliases",
    ]

    write_tsv(primary_path, primary_rows, primary_headers)
    write_tsv(multiplicity_path, multiplicity_rows, multiplicity_headers)

    with open(summary_path, "w", encoding="utf-8") as summary:
        summary.write("Genome-linking multiplicity QC\n")
        summary.write("==============================\n")
        summary.write(f"Input linked TSV: {args.linked}\n")
        summary.write(f"Candidate ID column: {args.candidate_id_column}\n")
        summary.write(f"Total linked rows read: {total_rows}\n")
        summary.write(f"Unique candidate IDs: {len(candidate_stats)}\n")
        summary.write(f"Primary rows written: {len(primary_rows)}\n")
        summary.write("\nPrimary mapping status counts:\n")
        for key, value in primary_status_counts.most_common():
            summary.write(f"- {key}: {value}\n")
        summary.write("\nPrimary rank counts:\n")
        for key, value in sorted(rank_counts.items(), key=lambda x: int(x[0])):
            summary.write(f"- rank_{key}: {value}\n")
        summary.write("\nMultiplicity distribution:\n")
        row_count_counter = Counter(stats["rows"] for stats in candidate_stats.values())
        for n_rows, n_candidates in sorted(row_count_counter.items())[:25]:
            summary.write(f"- {n_rows} linked rows: {n_candidates} candidates\n")
        if len(row_count_counter) > 25:
            summary.write("- additional multiplicity bins omitted from text summary; see TSV.\n")
        summary.write("\nOutputs:\n")
        summary.write(f"- Primary one-row-per-candidate table: {primary_path}\n")
        summary.write(f"- Multiplicity table: {multiplicity_path}\n")
        summary.write(f"- Summary: {summary_path}\n")

    print("Done.")
    print(f"Primary table: {primary_path}")
    print(f"Multiplicity table: {multiplicity_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
