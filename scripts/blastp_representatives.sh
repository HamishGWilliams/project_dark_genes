#!/bin/bash
#SBATCH --job-name=blastp_representative
#SBATCH --mem=300G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=3-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

module load blast-plus/2.16.0

THREADS="${SLURM_NTASKS:-8}"
TOP_N=10

mkdir -p 02_annotation/db
mkdir -p 02_annotation/blastp
mkdir -p logs/outputs
mkdir -p logs/errors

# representative proteome
QUERY="02_annotation/input/equina_representative_longest_per_gene.no_stop.fa"

SWISSPROT_FASTA="02_annotation/reference_dbs/swissprot_all.fasta"
TREMBL_FASTA="02_annotation/reference_dbs/trembl_cnidaria_selected.fasta"

run_blastp_search () {
    local label="$1"
    local ref_fasta="$2"
    local blast_db="$3"
    local out_tsv="$4"
    local top_hits="$5"

    if [[ ! -f "$ref_fasta" ]]; then
        echo "ERROR: reference FASTA not found for ${label}: $ref_fasta" >&2
        exit 1
    fi

    echo "========================================"
    echo "Running BLASTp for: ${label}"
    echo "Reference FASTA: $ref_fasta"
    echo "BLAST DB:        $blast_db"
    echo "Raw top 10 TSV:  $out_tsv"
    echo "Best hit TSV:    $top_hits"
    echo "Threads:         $THREADS"
    echo "Top N:           $TOP_N"
    echo "========================================"

    makeblastdb \
      -in "$ref_fasta" \
      -dbtype prot \
      -out "$blast_db"

    blastp \
      -query "$QUERY" \
      -db "$blast_db" \
      -out "$out_tsv" \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
      -num_threads "$THREADS" \
      -max_target_seqs "$TOP_N" \
      -max_hsps 1

    awk -F '\t' '!seen[$1]++' "$out_tsv" > "$top_hits"

    echo
    echo "Line counts for ${label}:"
    wc -l "$out_tsv" "$top_hits"

    echo
    echo "First few best hits for ${label}:"
    head "$top_hits"
    echo
}

run_blastp_search \
    "Swiss-Prot" \
    "$SWISSPROT_FASTA" \
    "02_annotation/db/swissprot_all_blast" \
    "02_annotation/blastp/equina_vs_swissprot_all.blastp.top10_raw_representative.tsv" \
    "02_annotation/blastp/equina_vs_swissprot_all.best_hit_representative.tsv"

run_blastp_search \
    "Cnidarian TrEMBL" \
    "$TREMBL_FASTA" \
    "02_annotation/db/trembl_cnidaria_selected_blast" \
    "02_annotation/blastp/equina_vs_trembl_cnidaria_selected.blastp.top10_raw_representative.tsv" \
    "02_annotation/blastp/equina_vs_trembl_cnidaria_selected.best_hit_representative.tsv"

echo "All BLASTp searches completed successfully."