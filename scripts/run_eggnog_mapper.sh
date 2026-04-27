#!/bin/bash
#SBATCH --job-name=eggnog_mapper
#SBATCH --mem=128G
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

EGGNOG_ENV="/uoa/scratch/users/r02hw22/project_dark_genes/conda_envs/eggnog-mapper-2.1.13"
EGGNOG_DATA_DIR="${EGGNOG_DATA_DIR:-/uoa/scratch/users/r02hw22/project_dark_genes/eggnog_data}"

export EGGNOG_DATA_DIR
export PATH="$EGGNOG_ENV/bin:$PATH"

THREADS="${SLURM_NTASKS:-16}"

QUERY="${QUERY:-02_annotation/input/equina_representative_longest_per_gene.no_stop.fa}"
RUN_NAME="${RUN_NAME:-equina_representative_full}"

OUTDIR="02_annotation/eggnog/raw"
SUMMARY_DIR="02_annotation/eggnog/summary"
TMPDIR_EGGNOG="/uoa/scratch/users/r02hw22/project_dark_genes/eggnog_tmp_${SLURM_JOB_ID:-manual}"

mkdir -p "$OUTDIR" "$SUMMARY_DIR" "$TMPDIR_EGGNOG" logs/outputs logs/errors

echo "========================================"
echo "Project dir      : $PROJECT_DIR"
echo "Query FASTA      : $QUERY"
echo "Run name         : $RUN_NAME"
echo "Threads          : $THREADS"
echo "eggNOG env       : $EGGNOG_ENV"
echo "eggNOG data dir  : $EGGNOG_DATA_DIR"
echo "Output directory : $OUTDIR"
echo "Temp directory   : $TMPDIR_EGGNOG"
echo "========================================"

if [[ ! -x "$EGGNOG_ENV/bin/emapper.py" ]]; then
    echo "ERROR: emapper.py not found or not executable: $EGGNOG_ENV/bin/emapper.py" >&2
    exit 1
fi

if [[ ! -f "$QUERY" ]]; then
    echo "ERROR: Query FASTA not found: $QUERY" >&2
    exit 1
fi

if [[ ! -s "$EGGNOG_DATA_DIR/eggnog.db" ]]; then
    echo "ERROR: eggnog.db missing or empty: $EGGNOG_DATA_DIR/eggnog.db" >&2
    exit 1
fi

if [[ ! -s "$EGGNOG_DATA_DIR/eggnog_proteins.dmnd" ]]; then
    echo "ERROR: eggnog_proteins.dmnd missing or empty: $EGGNOG_DATA_DIR/eggnog_proteins.dmnd" >&2
    exit 1
fi

which emapper.py
emapper.py --version

emapper.py \
  -i "$QUERY" \
  --itype proteins \
  -m diamond \
  --data_dir "$EGGNOG_DATA_DIR" \
  --cpu "$THREADS" \
  --output "$RUN_NAME" \
  --output_dir "$OUTDIR" \
  --temp_dir "$TMPDIR_EGGNOG" \
  --override

ANNOT="${OUTDIR}/${RUN_NAME}.emapper.annotations"

if [[ ! -s "$ANNOT" ]]; then
    echo "ERROR: eggNOG-mapper did not produce a non-empty annotations file: $ANNOT" >&2
    exit 1
fi

echo "eggNOG-mapper completed successfully."

echo "Raw output files:"
ls -lh "$OUTDIR"/"${RUN_NAME}".emapper.*

echo "Annotation rows excluding comments:"
grep -v '^#' "$ANNOT" | wc -l

# Lightweight downstream summaries
grep -v '^#' "$ANNOT" > "${SUMMARY_DIR}/${RUN_NAME}.annotations.no_comments.tsv"

grep -v '^#' "$ANNOT" | cut -f1 | sort -u \
  > "${SUMMARY_DIR}/${RUN_NAME}.annotated_query_ids.txt"

# COG category counts; standard eggNOG-mapper v2 annotations usually store COG_category in field 7.
awk -F'\t' '
BEGIN { OFS="\t" }
!/^#/ {
    cog=$7
    if (cog != "" && cog != "-" && cog != "NA") {
        split(cog, chars, "")
        for (i in chars) print chars[i]
    }
}' "$ANNOT" | sort | uniq -c | sort -nr \
  > "${SUMMARY_DIR}/${RUN_NAME}.cog_category_counts.txt" || true

# GO terms; standard eggNOG-mapper v2 annotations usually store GOs in field 10.
awk -F'\t' '
BEGIN { OFS="\t" }
!/^#/ {
    go=$10
    if (go != "" && go != "-" && go != "NA") {
        n=split(go, arr, ",")
        for (i=1; i<=n; i++) print $1, arr[i]
    }
}' "$ANNOT" | sort -u \
  > "${SUMMARY_DIR}/${RUN_NAME}.go_terms.tsv" || true

echo "Summary files:"
ls -lh "$SUMMARY_DIR"/"${RUN_NAME}".* || true

rm -rf "$TMPDIR_EGGNOG"

echo "Done."