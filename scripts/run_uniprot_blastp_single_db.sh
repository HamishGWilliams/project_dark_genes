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

module load blast+/2.14.1 || true

THREADS="${SLURM_NTASKS:-8}"

: "${DB_NAME:?Need DB_NAME}"
: "${REF_FASTA:?Need REF_FASTA}"
: "${QUERY_NOHIT:?Need QUERY_NOHIT}"

mkdir -p 02_annotation/db 02_annotation/blastp

if [[ ! -f "$QUERY_NOHIT" ]]; then
    echo "ERROR: No-hit query FASTA not found: $QUERY_NOHIT" >&2
    exit 1
fi
if [[ ! -f "$REF_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found: $REF_FASTA" >&2
    exit 1
fi
if ! command -v makeblastdb >/dev/null 2>&1 || ! command -v blastp >/dev/null 2>&1; then
    echo "ERROR: BLAST+ tools not available in PATH" >&2
    exit 1
fi

BLAST_DB="02_annotation/db/${DB_NAME}_blast"
BLAST_OUT="02_annotation/blastp/equina_no_swissprot_diamond_hit_vs_${DB_NAME}.blastp.tsv"
BLAST_TOP="02_annotation/blastp/equina_no_swissprot_diamond_hit_vs_${DB_NAME}.blastp.top_hits.tsv"

echo "Running BLASTp fallback for $DB_NAME"
makeblastdb -in "$REF_FASTA" -dbtype prot -out "$BLAST_DB"

blastp \
  -query "$QUERY_NOHIT" \
  -db "$BLAST_DB" \
  -out "$BLAST_OUT" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -num_threads "$THREADS"

awk -F '\t' '!seen[$1]++' "$BLAST_OUT" > "$BLAST_TOP"

wc -l "$BLAST_OUT" "$BLAST_TOP"
head "$BLAST_TOP" || true
