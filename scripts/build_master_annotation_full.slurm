#!/bin/bash
#SBATCH --job-name=master_annot_full
#SBATCH --mem=64G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=06:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

mkdir -p logs/outputs
mkdir -p logs/errors
mkdir -p 02_annotation/master/full

echo "========================================"
echo "Building full representative master annotation table"
echo "Project dir: $PROJECT_DIR"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "========================================"

# Load Python if needed.
# Use whichever Python module works on Maxwell.
module load python/3.11.9 || module load python/3.9.12 || true

PYTHON_BIN="$(command -v python3 || command -v python)"

echo "Using Python:"
"$PYTHON_BIN" --version
echo "$PYTHON_BIN"

echo
echo "Checking expected input files:"
INPUTS=(
  "02_annotation/input/equina_representative_longest_per_gene.no_stop.fa"
  "01_qc/isoforms/representative_longest_protein_per_gene.tsv"
  "02_annotation/diamond/filtered/equina_vs_swissprot_all_representative.evalue_1e-5.top_hits.tsv"
  "02_annotation/diamond/filtered/equina_vs_trembl_cnidaria_selected_representative.evalue_1e-5.top_hits.tsv"
  "02_annotation/blastp/filtered/equina_vs_swissprot_all.evalue_1e-5.best_hit_representative.tsv"
  "02_annotation/blastp/filtered/equina_vs_trembl_cnidaria_selected.evalue_1e-5.best_hit_representative.tsv"
  "02_annotation/interproscan/raw/equina_representative_longest_per_gene.interproscan.tsv"
  "02_annotation/eggnog/raw/equina_representative_full.emapper.annotations"
  "02_annotation/signalp/summary/equina_representative_full.signalp5_summary.tsv"
)

for file in "${INPUTS[@]}"; do
    if [[ -s "$file" ]]; then
        echo "FOUND:   $file"
    else
        echo "MISSING OR EMPTY: $file"
    fi
done

echo
echo "Running full master annotation compiler..."

"$PYTHON_BIN" scripts/build_master_annotation_full.py \
  --fasta 02_annotation/input/equina_representative_longest_per_gene.no_stop.fa \
  --rep-lookup 01_qc/isoforms/representative_longest_protein_per_gene.tsv \
  --diamond-swissprot 02_annotation/diamond/filtered/equina_vs_swissprot_all_representative.evalue_1e-5.top_hits.tsv \
  --diamond-trembl-cnidaria 02_annotation/diamond/filtered/equina_vs_trembl_cnidaria_selected_representative.evalue_1e-5.top_hits.tsv \
  --blastp-swissprot 02_annotation/blastp/filtered/equina_vs_swissprot_all.evalue_1e-5.best_hit_representative.tsv \
  --blastp-trembl-cnidaria 02_annotation/blastp/filtered/equina_vs_trembl_cnidaria_selected.evalue_1e-5.best_hit_representative.tsv \
  --interproscan 02_annotation/interproscan/raw/equina_representative_longest_per_gene.interproscan.tsv \
  --eggnog 02_annotation/eggnog/raw/equina_representative_full.emapper.annotations \
  --signalp 02_annotation/signalp/summary/equina_representative_full.signalp5_summary.tsv \
  --out 02_annotation/master/full/equina_representative_full.master_annotation.tsv

echo
echo "Output files:"
ls -lh 02_annotation/master/full/equina_representative_full.master_annotation.tsv
ls -lh 02_annotation/master/full/equina_representative_full.master_annotation.summary.txt

echo
echo "Line counts:"
wc -l 02_annotation/master/full/equina_representative_full.master_annotation.tsv
cat 02_annotation/master/full/equina_representative_full.master_annotation.summary.txt

echo
echo "Finished: $(date)"
echo "Full representative master annotation table completed successfully."