#!/bin/bash
#SBATCH --job-name=diamond_rep_full
#SBATCH --mem=128G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=1-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

module load diamond/2.1.10

THREADS="${SLURM_NTASKS:-8}"
EVALUE_THRESHOLD="1e-5"

mkdir -p 02_annotation/db
mkdir -p 02_annotation/diamond/raw
mkdir -p 02_annotation/diamond/filtered
mkdir -p logs/outputs
mkdir -p logs/errors

# Full representative no-stop proteome
QUERY="02_annotation/input/equina_representative_longest_per_gene.no_stop.fa"

SWISSPROT_FASTA="02_annotation/reference_dbs/swissprot_all.fasta"
TREMBL_FASTA="02_annotation/reference_dbs/trembl_cnidaria_selected.fasta"

if [[ ! -f "$QUERY" ]]; then
    echo "ERROR: query FASTA not found: $QUERY" >&2
    exit 1
fi

run_diamond_search () {
    local label="$1"
    local ref_fasta="$2"
    local db_prefix="$3"
    local raw_out_tsv="$4"
    local filtered_out_tsv="$5"
    local filtered_top_hits="$6"

    if [[ ! -f "$ref_fasta" ]]; then
        echo "ERROR: reference FASTA not found for ${label}: $ref_fasta" >&2
        exit 1
    fi

    echo "========================================"
    echo "Running DIAMOND for: ${label}"
    echo "Query FASTA:      $QUERY"
    echo "Reference FASTA:  $ref_fasta"
    echo "Database prefix:  $db_prefix"
    echo "Raw output TSV:   $raw_out_tsv"
    echo "Filtered TSV:     $filtered_out_tsv"
    echo "Filtered top hit: $filtered_top_hits"
    echo "E-value cutoff:   $EVALUE_THRESHOLD"
    echo "Threads:          $THREADS"
    echo "========================================"

    diamond makedb \
      --in "$ref_fasta" \
      --db "$db_prefix"

    diamond blastp \
      --query "$QUERY" \
      --db "$db_prefix" \
      --out "$raw_out_tsv" \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle \
      --evalue "$EVALUE_THRESHOLD" \
      --threads "$THREADS"

    # Post-processing threshold filter.
    # Column 11 is evalue in this outfmt:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle
    awk -F '\t' -v e="$EVALUE_THRESHOLD" 'BEGIN {OFS="\t"} $11 <= e {print}' \
      "$raw_out_tsv" > "$filtered_out_tsv"

    # Keep the first retained hit per query after e-value filtering.
    # DIAMOND output is already ordered by best hit, so this gives the filtered top hit.
    awk -F '\t' '!seen[$1]++' "$filtered_out_tsv" > "$filtered_top_hits"

    echo
    echo "Line counts for ${label}:"
    wc -l "$raw_out_tsv" "$filtered_out_tsv" "$filtered_top_hits"

    echo
    echo "Check max e-value in filtered output for ${label}:"
    awk -F '\t' 'BEGIN {max=0} {if ($11 > max) max=$11} END {print max}' "$filtered_out_tsv"

    echo
    echo "First few filtered top hits for ${label}:"
    head "$filtered_top_hits"
    echo
}

run_diamond_search \
    "Swiss-Prot" \
    "$SWISSPROT_FASTA" \
    "02_annotation/db/swissprot_all" \
    "02_annotation/diamond/raw/equina_vs_swissprot_all_representative.diamond.tsv" \
    "02_annotation/diamond/filtered/equina_vs_swissprot_all_representative.evalue_1e-5.diamond.tsv" \
    "02_annotation/diamond/filtered/equina_vs_swissprot_all_representative.evalue_1e-5.top_hits.tsv"

run_diamond_search \
    "Cnidarian TrEMBL" \
    "$TREMBL_FASTA" \
    "02_annotation/db/trembl_cnidaria_selected" \
    "02_annotation/diamond/raw/equina_vs_trembl_cnidaria_selected_representative.diamond.tsv" \
    "02_annotation/diamond/filtered/equina_vs_trembl_cnidaria_selected_representative.evalue_1e-5.diamond.tsv" \
    "02_annotation/diamond/filtered/equina_vs_trembl_cnidaria_selected_representative.evalue_1e-5.top_hits.tsv"

echo "All full representative DIAMOND searches completed successfully."