#!/bin/bash
#SBATCH --job-name=signalp
#SBATCH --mem=32G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=12:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

set -euo pipefail

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

# Initialise module command if needed
if [[ -f /etc/profile.d/lmod.sh ]]; then
    source /etc/profile.d/lmod.sh
elif [[ -f /usr/share/lmod/lmod/init/bash ]]; then
    source /usr/share/lmod/lmod/init/bash
fi

module load signalp/5.0b

QUERY="${QUERY:-02_annotation/input/equina_representative_longest_per_gene.no_stop.fa}"
RUN_NAME="${RUN_NAME:-equina_representative_full}"

OUTDIR="02_annotation/signalp/raw/${RUN_NAME}"
SUMMARY_DIR="02_annotation/signalp/summary"

mkdir -p "$OUTDIR" "$SUMMARY_DIR" logs/outputs logs/errors

echo "========================================"
echo "Project dir : $PROJECT_DIR"
echo "Query FASTA : $QUERY"
echo "Run name    : $RUN_NAME"
echo "Output dir  : $OUTDIR"
echo "Summary dir : $SUMMARY_DIR"
echo "Loaded modules:"
module list || true
echo "SignalP path:"
which signalp || true
echo "========================================"

if [[ ! -f "$QUERY" ]]; then
    echo "ERROR: Query FASTA not found: $QUERY" >&2
    exit 1
fi

if ! command -v signalp >/dev/null 2>&1; then
    echo "ERROR: signalp command not found after loading signalp/5.0b" >&2
    exit 1
fi

signalp -h 2>&1 | head -n 30 || true

signalp \
  -fasta "$QUERY" \
  -org euk \
  -format short \
  -prefix "${OUTDIR}/${RUN_NAME}"

echo "SignalP command finished."

echo "Raw output files:"
find "$OUTDIR" -maxdepth 1 -type f -print -exec ls -lh {} \;

SUMMARY_FILE="${OUTDIR}/${RUN_NAME}_summary.signalp5"
PRED_FILE="${OUTDIR}/${RUN_NAME}_prediction_results.txt"

if [[ -s "$SUMMARY_FILE" ]]; then
    echo "Detected SignalP 5 summary file: $SUMMARY_FILE"
    cp "$SUMMARY_FILE" "${SUMMARY_DIR}/${RUN_NAME}.signalp5_summary.tsv"
elif [[ -s "$PRED_FILE" ]]; then
    echo "Detected SignalP 5 prediction file: $PRED_FILE"
    cp "$PRED_FILE" "${SUMMARY_DIR}/${RUN_NAME}.signalp5_summary.tsv"
else
    echo "ERROR: Could not find expected SignalP 5 output." >&2
    find "$OUTDIR" -maxdepth 1 -type f -print
    exit 1
fi

awk '
BEGIN { OFS="\t" }
$0 ~ /^#/ { next }
tolower($0) ~ /signal peptide| sp\(|\tsp\t|^.*\tsp/ { print }
' "${SUMMARY_DIR}/${RUN_NAME}.signalp5_summary.tsv" \
  > "${SUMMARY_DIR}/${RUN_NAME}.signalp5_positive_predictions.tsv" || true

echo "Summary files:"
ls -lh "$SUMMARY_DIR"/"${RUN_NAME}".*

echo "Done."