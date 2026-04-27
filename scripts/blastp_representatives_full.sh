#!/bin/bash
#SBATCH --job-name=blastp_rep_full
#SBATCH --mem=300G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=3-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

module load blast-plus/2.16.0

THREADS="${SLURM_NTASKS:-16}"
TOP_N=10
EVALUE_THRESHOLD="1e-5"

mkdir -p 02_annotation/db
mkdir -p 02_annotation/blastp/raw
mkdir -p 02_annotation/blastp/filtered
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

make_blast_db_if_needed () {
    local ref_fasta="$1"
    local blast_db="$2"
    local label="$3"

    if [[ -s "${blast_db}.pin" && -s "${blast_db}.psq" && -s "${blast_db}.phr" ]]; then
        echo "BLAST database already exists for ${label}: ${blast_db}"
        echo "Skipping makeblastdb."
    else
        echo "Building BLAST database for ${label}: ${blast_db}"

        makeblastdb \
          -in "$ref_fasta" \
          -dbtype prot \
          -out "$blast_db"
    fi
}

run_blastp_search () {
    local label="$1"
    local ref_fasta="$2"
    local blast_db="$3"
    local raw_out_tsv="$4"
    local filtered_out_tsv="$5"
    local filtered_best_hit="$6"

    if [[ ! -f "$ref_fasta" ]]; then
        echo "ERROR: reference FASTA not found for ${label}: $ref_fasta" >&2
        exit 1
    fi

    echo "========================================"
    echo "Running BLASTp for: ${label}"
    echo "Query FASTA:       $QUERY"
    echo "Reference FASTA:   $ref_fasta"
    echo "BLAST DB:          $blast_db"
    echo "Raw top ${TOP_N} TSV:    $raw_out_tsv"
    echo "Filtered TSV:      $filtered_out_tsv"
    echo "Filtered best hit: $filtered_best_hit"
    echo "E-value cutoff:    $EVALUE_THRESHOLD"
    echo "Threads:           $THREADS"
    echo "Top N:             $TOP_N"
    echo "========================================"

    make_blast_db_if_needed "$ref_fasta" "$blast_db" "$label"

    # Speed-related choices:
    # - Use explicit e-value filtering during BLASTp itself.
    # - Keep only top N targets per query rather than unlimited hits.
    # - Keep only one HSP per target.
    # - Reuse existing BLAST databases instead of rebuilding them every run.
    blastp \
      -query "$QUERY" \
      -db "$blast_db" \
      -out "$raw_out_tsv" \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
      -evalue "$EVALUE_THRESHOLD" \
      -num_threads "$THREADS" \
      -max_target_seqs "$TOP_N" \
      -max_hsps 1

    # Post-processing threshold filter.
    # Column 11 is evalue in this outfmt:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
    awk -F '\t' -v e="$EVALUE_THRESHOLD" 'BEGIN {OFS="\t"} $11 <= e {print}' \
      "$raw_out_tsv" > "$filtered_out_tsv"

    # Keep the first retained hit per query after e-value filtering.
    # BLASTp output is ordered by best hit per query under this setup.
    awk -F '\t' '!seen[$1]++' "$filtered_out_tsv" > "$filtered_best_hit"

    echo
    echo "Line counts for ${label}:"
    wc -l "$raw_out_tsv" "$filtered_out_tsv" "$filtered_best_hit"

    echo
    echo "Check max e-value in filtered output for ${label}:"
    awk -F '\t' 'BEGIN {max=0} {if ($11 > max) max=$11} END {print max}' "$filtered_out_tsv"

    echo
    echo "First few filtered best hits for ${label}:"
    head "$filtered_best_hit"
    echo
}

run_blastp_search \
    "Swiss-Prot" \
    "$SWISSPROT_FASTA" \
    "02_annotation/db/swissprot_all_blast" \
    "02_annotation/blastp/raw/equina_vs_swissprot_all.blastp.top10_raw_representative.tsv" \
    "02_annotation/blastp/filtered/equina_vs_swissprot_all.evalue_1e-5.blastp.tsv" \
    "02_annotation/blastp/filtered/equina_vs_swissprot_all.evalue_1e-5.best_hit_representative.tsv"

run_blastp_search \
    "Cnidarian TrEMBL" \
    "$TREMBL_FASTA" \
    "02_annotation/db/trembl_cnidaria_selected_blast" \
    "02_annotation/blastp/raw/equina_vs_trembl_cnidaria_selected.blastp.top10_raw_representative.tsv" \
    "02_annotation/blastp/filtered/equina_vs_trembl_cnidaria_selected.evalue_1e-5.blastp.tsv" \
    "02_annotation/blastp/filtered/equina_vs_trembl_cnidaria_selected.evalue_1e-5.best_hit_representative.tsv"

echo "All full representative BLASTp searches completed successfully."