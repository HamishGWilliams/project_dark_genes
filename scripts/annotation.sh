#!/bin/bash

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

mkdir -p logs/outputs logs/errors
mkdir -p 02_annotation/db
mkdir -p 02_annotation/diamond
mkdir -p 02_annotation/blastp

QUERY="00_raw/equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa"

# Toggle these if you only want one search type
RUN_DIAMOND="true"
RUN_BLASTP="false"

# Database definitions
DB_NAMES=(
  "swissprot_all"
  "trembl_cnidaria_selected"
  "metazoa_selected"
)

DB_FASTAS=(
  "02_annotation/reference_dbs/swissprot_all.fasta"
  "02_annotation/reference_dbs/trembl_cnidaria_selected.fasta"
  "02_annotation/reference_dbs/metazoa_selected.fasta"
)

# Write a reusable worker script once
WORKER_SCRIPT="scripts/run_uniprot_search_single_db.sh"

cat > "$WORKER_SCRIPT" <<'EOF'
#!/usr/bin/env bash
#SBATCH --mem=256G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=7-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

module load diamond/2.1.10
module load blast+/2.14.1 || true

THREADS="${SLURM_NTASKS:-8}"

: "${DB_NAME:?Need DB_NAME}"
: "${REF_FASTA:?Need REF_FASTA}"
: "${QUERY:?Need QUERY}"
: "${RUN_DIAMOND:?Need RUN_DIAMOND}"
: "${RUN_BLASTP:?Need RUN_BLASTP}"

if [[ ! -f "$QUERY" ]]; then
    echo "ERROR: Query FASTA not found: $QUERY" >&2
    exit 1
fi

if [[ ! -f "$REF_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found: $REF_FASTA" >&2
    exit 1
fi

mkdir -p 02_annotation/db
mkdir -p 02_annotation/diamond
mkdir -p 02_annotation/blastp

echo "========================================"
echo "DB_NAME      : $DB_NAME"
echo "REF_FASTA    : $REF_FASTA"
echo "QUERY        : $QUERY"
echo "THREADS      : $THREADS"
echo "RUN_DIAMOND  : $RUN_DIAMOND"
echo "RUN_BLASTP   : $RUN_BLASTP"
echo "========================================"

if [[ "$RUN_DIAMOND" == "true" ]]; then
    DIAMOND_DB="02_annotation/db/${DB_NAME}"
    DIAMOND_OUT="02_annotation/diamond/equina_vs_${DB_NAME}.diamond.tsv"
    DIAMOND_TOP="02_annotation/diamond/equina_vs_${DB_NAME}.diamond.top_hits.tsv"

    echo "Building DIAMOND DB for $DB_NAME"
    diamond makedb \
      --in "$REF_FASTA" \
      --db "$DIAMOND_DB"

    echo "Running DIAMOND blastp for $DB_NAME"
    diamond blastp \
      --query "$QUERY" \
      --db "$DIAMOND_DB" \
      --out "$DIAMOND_OUT" \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
      --threads "$THREADS"

    awk -F '\t' '!seen[$1]++' "$DIAMOND_OUT" > "$DIAMOND_TOP"

    echo "DIAMOND outputs for $DB_NAME"
    wc -l "$DIAMOND_OUT" "$DIAMOND_TOP"
    head "$DIAMOND_TOP" || true
fi

if [[ "$RUN_BLASTP" == "true" ]]; then
    if ! command -v makeblastdb >/dev/null 2>&1 || ! command -v blastp >/dev/null 2>&1; then
        echo "ERROR: BLAST+ tools not available in PATH" >&2
        exit 1
    fi

    BLAST_DB="02_annotation/db/${DB_NAME}_blast"
    BLAST_OUT="02_annotation/blastp/equina_vs_${DB_NAME}.blastp.tsv"
    BLAST_TOP="02_annotation/blastp/equina_vs_${DB_NAME}.blastp.top_hits.tsv"

    echo "Building BLAST DB for $DB_NAME"
    makeblastdb \
      -in "$REF_FASTA" \
      -dbtype prot \
      -out "$BLAST_DB"

    echo "Running BLASTp for $DB_NAME"
    blastp \
      -query "$QUERY" \
      -db "$BLAST_DB" \
      -out "$BLAST_OUT" \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
      -num_threads "$THREADS"

    awk -F '\t' '!seen[$1]++' "$BLAST_OUT" > "$BLAST_TOP"

    echo "BLASTp outputs for $DB_NAME"
    wc -l "$BLAST_OUT" "$BLAST_TOP"
    head "$BLAST_TOP" || true
fi

echo "Completed searches for $DB_NAME"
EOF

chmod +x "$WORKER_SCRIPT"

# Submit one SLURM job per database
for i in "${!DB_NAMES[@]}"; do
    DB_NAME="${DB_NAMES[$i]}"
    REF_FASTA="${DB_FASTAS[$i]}"

    if [[ ! -f "$REF_FASTA" ]]; then
        echo "WARNING: Skipping $DB_NAME because FASTA not found: $REF_FASTA" >&2
        continue
    fi

    echo "Submitting job for $DB_NAME"
    sbatch \
      --job-name="search_${DB_NAME}" \
      --export=ALL,DB_NAME="$DB_NAME",REF_FASTA="$REF_FASTA",QUERY="$QUERY",RUN_DIAMOND="$RUN_DIAMOND",RUN_BLASTP="$RUN_BLASTP" \
      "$WORKER_SCRIPT"
done