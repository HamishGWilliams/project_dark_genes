#!/bin/bash
#SBATCH --job-name=interproscan
#SBATCH --mem=128G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=7-00:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

module load python/3.9.12
module load openjdk/11.0.20.1_1

export JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))
export PATH="$JAVA_HOME/bin:$PATH"

# Edit this if your InterProScan installation is elsewhere.
IPS="${IPS:-$HOME/software/interproscan/interproscan-5.77-108.0/interproscan.sh}"

THREADS="${SLURM_NTASKS:-16}"

QUERY="${QUERY:-00_raw/equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa}"
RUN_NAME="${RUN_NAME:-equina_full}"

OUTDIR="02_annotation/interproscan/raw"
SUMMARY_DIR="02_annotation/interproscan/summary"
TMPDIR_IPS="/uoa/scratch/users/r02hw22/project_dark_genes/interproscan_tmp_${SLURM_JOB_ID:-manual}"

OUTBASE="${OUTDIR}/${RUN_NAME}.interproscan"

mkdir -p "$OUTDIR" "$SUMMARY_DIR" "$TMPDIR_IPS" logs/outputs logs/errors

echo "========================================"
echo "Project dir : $PROJECT_DIR"
echo "InterProScan: $IPS"
echo "Query FASTA : $QUERY"
echo "Run name    : $RUN_NAME"
echo "Threads     : $THREADS"
echo "JAVA_HOME   : $JAVA_HOME"
echo "Temp dir    : $TMPDIR_IPS"
echo "Output base : $OUTBASE"
echo "========================================"

if [[ ! -x "$IPS" ]]; then
    echo "ERROR: InterProScan script not found or not executable: $IPS" >&2
    exit 1
fi

if [[ ! -f "$QUERY" ]]; then
    echo "ERROR: Query protein FASTA not found: $QUERY" >&2
    exit 1
fi

java -version
python3 --version

"$IPS" \
  -i "$QUERY" \
  -t p \
  -f TSV,GFF3,XML \
  -goterms \
  -pa \
  -cpu "$THREADS" \
  -T "$TMPDIR_IPS" \
  -b "$OUTBASE" \
  -dp

echo "InterProScan finished."

TSV="${OUTBASE}.tsv"

if [[ -f "$TSV" ]]; then
    echo "Summarising TSV: $TSV"

    awk -F '\t' '{print $1}' "$TSV" | sort -u > "${SUMMARY_DIR}/${RUN_NAME}.interproscan_matched_protein_ids.txt"

    awk -F '\t' '{print $4}' "$TSV" | sort | uniq -c | sort -nr \
      > "${SUMMARY_DIR}/${RUN_NAME}.member_database_counts.txt"

    awk -F '\t' 'BEGIN{OFS="\t"} $12 != "-" {print $1,$4,$5,$6,$7,$8,$12,$13}' "$TSV" | sort -u \
      > "${SUMMARY_DIR}/${RUN_NAME}.interpro_domain_summary.tsv"

    awk -F '\t' 'BEGIN{OFS="\t"} $14 != "-" {
        n=split($14, go, "|");
        for (i=1; i<=n; i++) print $1, go[i]
    }' "$TSV" | sort -u \
      > "${SUMMARY_DIR}/${RUN_NAME}.go_terms.tsv"

    awk -F '\t' 'BEGIN{OFS="\t"} $15 != "-" {
        n=split($15, pw, "|");
        for (i=1; i<=n; i++) print $1, pw[i]
    }' "$TSV" | sort -u \
      > "${SUMMARY_DIR}/${RUN_NAME}.pathways.tsv"

    echo "Raw TSV rows:"
    wc -l "$TSV"

    echo "Proteins with at least one InterProScan match:"
    wc -l "${SUMMARY_DIR}/${RUN_NAME}.interproscan_matched_protein_ids.txt"

    echo "Top member databases:"
    head "${SUMMARY_DIR}/${RUN_NAME}.member_database_counts.txt"

    echo "Summary files:"
    ls -lh "$SUMMARY_DIR"/"${RUN_NAME}".*
else
    echo "ERROR: Expected TSV not found: $TSV" >&2
    exit 1
fi

rm -rf "$TMPDIR_IPS"

echo "Done."