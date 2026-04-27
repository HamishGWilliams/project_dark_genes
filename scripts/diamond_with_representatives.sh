#!/bin/bash
#SBATCH --job-name=diamond_representative
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

mkdir -p 02_annotation/db
mkdir -p 02_annotation/diamond
mkdir -p logs/outputs
mkdir -p logs/errors

# representative proteome
QUERY="02_annotation/input/equina_representative_longest_per_gene.no_stop.fa"

SWISSPROT_FASTA="02_annotation/reference_dbs/swissprot_all.fasta"
TREMBL_FASTA="02_annotation/reference_dbs/trembl_cnidaria_selected.fasta"

run_diamond_search () {
    local label="$1"
    local ref_fasta="$2"
    local db_prefix="$3"
    local out_tsv="$4"
    local top_hits="$5"

    if [[ ! -f "$ref_fasta" ]]; then
        echo "ERROR: reference FASTA not found for ${label}: $ref_fasta" >&2
        exit 1
    fi

    echo "========================================"
    echo "Running DIAMOND for: ${label}"
    echo "Reference FASTA: $ref_fasta"
    echo "Database prefix: $db_prefix"
    echo "Output TSV:      $out_tsv"
    echo "Threads:         $THREADS"
    echo "========================================"

    diamond makedb \
      --in "$ref_fasta" \
      --db "$db_prefix"

    diamond blastp \
      --query "$QUERY" \
      --db "$db_prefix" \
      --out "$out_tsv" \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
      --threads "$THREADS"

    awk -F '\t' '!seen[$1]++' "$out_tsv" > "$top_hits"

    echo
    echo "Line counts for ${label}:"
    wc -l "$out_tsv" "$top_hits"

    echo
    echo "First few top hits for ${label}:"
    head "$top_hits"
    echo
}

run_diamond_search \
    "Swiss-Prot" \
    "$SWISSPROT_FASTA" \
    "02_annotation/db/swissprot_all" \
    "02_annotation/diamond/equina_vs_swissprot_all_representative.diamond.tsv" \
    "02_annotation/diamond/equina_vs_swissprot_all_representative.top_hits.tsv"

run_diamond_search \
    "TrEMBL" \
    "$TREMBL_FASTA" \
    "02_annotation/db/trembl_cnidaria_selected" \
    "02_annotation/diamond/equina_vs_trembl_cnidaria_selected_representative.diamond.tsv" \
    "02_annotation/diamond/equina_vs_trembl_cnidaria_selected_representative.top_hits.tsv"

echo "All DIAMOND searches completed successfully."