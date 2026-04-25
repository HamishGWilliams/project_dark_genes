#!/bin/bash
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

: "${DB_NAME:?Need DB_NAME}"
: "${REF_FASTA:?Need REF_FASTA}"
: "${QUERY:?Need QUERY}"

mkdir -p 02_annotation/db 02_annotation/diamond

if [[ ! -f "$QUERY" ]]; then
    echo "ERROR: Query FASTA not found: $QUERY" >&2
    exit 1
fi
if [[ ! -f "$REF_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found: $REF_FASTA" >&2
    exit 1
fi

DIAMOND_DB="02_annotation/db/${DB_NAME}"
DIAMOND_OUT="02_annotation/diamond/equina_vs_${DB_NAME}.diamond.tsv"
DIAMOND_TOP="02_annotation/diamond/equina_vs_${DB_NAME}.diamond.top_hits.tsv"

echo "Running DIAMOND for $DB_NAME"
diamond makedb --in "$REF_FASTA" --db "$DIAMOND_DB"

diamond blastp \
  --query "$QUERY" \
  --db "$DIAMOND_DB" \
  --out "$DIAMOND_OUT" \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
  --threads "$THREADS"

awk -F '\t' '!seen[$1]++' "$DIAMOND_OUT" > "$DIAMOND_TOP"

wc -l "$DIAMOND_OUT" "$DIAMOND_TOP"
head "$DIAMOND_TOP" || true
