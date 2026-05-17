#!/bin/bash
#SBATCH --job-name=busco_dup_validation
#SBATCH --mem=32G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

set -euo pipefail

cd /uoa/home/r02hw22/sharedscratch/project_dark_genes

PROJECT_DIR="${PROJECT_DIR:-/uoa/scratch/users/r02hw22/project_dark_genes}"
LINEAGE="${BUSCO_LINEAGE:-metazoa_odb10}"
BUSCO_MODE="proteins"
THREADS="${SLURM_CPUS_ON_NODE:-16}"

# Prefer the representative no-stop protein set used throughout the annotation workflow.
PROTEOME="${BUSCO_PROTEOME:-${PROJECT_DIR}/02_annotation/input/equina_representative_longest_per_gene.no_stop.fa}"
if [[ ! -s "${PROTEOME}" ]]; then
    PROTEOME="${PROJECT_DIR}/02_annotation/eggnog/input/equina_representative_longest_per_gene.no_stop.fa"
fi

PRIMARY_LINKED="${PRIMARY_LINKED:-${PROJECT_DIR}/06_genome_lookup/equina_dark_candidates.genome_linked.primary.tsv}"
LOOKUP_TSV="${LOOKUP_TSV:-${PROJECT_DIR}/06_genome_lookup/equina_gff3_gene_transcript_protein_contig_lookup.tsv}"
DARK_FASTA="${DARK_FASTA:-${PROJECT_DIR}/03_dark_candidates/equina_dark_candidates.all.fa}"
DUP_OUTDIR="${DUP_OUTDIR:-${PROJECT_DIR}/07_duplication_context}"
BUSCO_OUTDIR="${BUSCO_OUTDIR:-${PROJECT_DIR}/05_busco}"
BUSCO_RUN_NAME="${BUSCO_RUN_NAME:-equina_representative_longest_per_gene_${LINEAGE}}"
BUSCO_STABLE_DIR="${BUSCO_OUTDIR}/${BUSCO_RUN_NAME}"
BUSCO_FULL_TABLE="${BUSCO_FULL_TABLE:-}"

mkdir -p \
    "${PROJECT_DIR}/logs/outputs" \
    "${PROJECT_DIR}/logs/errors" \
    "${BUSCO_OUTDIR}" \
    "${BUSCO_STABLE_DIR}" \
    "${DUP_OUTDIR}"

echo "Project directory: ${PROJECT_DIR}"
echo "Proteome: ${PROTEOME}"
echo "BUSCO lineage: ${LINEAGE}"
echo "BUSCO output directory: ${BUSCO_OUTDIR}"
echo "BUSCO stable directory: ${BUSCO_STABLE_DIR}"
echo "Primary linked candidates: ${PRIMARY_LINKED}"
echo "Duplication output directory: ${DUP_OUTDIR}"

if [[ ! -s "${PROTEOME}" ]]; then
    echo "ERROR: BUSCO proteome FASTA not found: ${PROTEOME}" >&2
    exit 1
fi

if [[ ! -s "${PRIMARY_LINKED}" ]]; then
    echo "ERROR: primary genome-linked candidate table not found: ${PRIMARY_LINKED}" >&2
    exit 1
fi

if [[ ! -s "${LOOKUP_TSV}" ]]; then
    echo "ERROR: genome lookup table not found: ${LOOKUP_TSV}" >&2
    exit 1
fi

if [[ ! -s "${DARK_FASTA}" ]]; then
    echo "ERROR: dark-candidate FASTA not found: ${DARK_FASTA}" >&2
    exit 1
fi

# If the user supplied BUSCO_FULL_TABLE, use it. Otherwise try to reuse an existing BUSCO full table.
if [[ -z "${BUSCO_FULL_TABLE}" ]]; then
    BUSCO_FULL_TABLE="$(find "${BUSCO_OUTDIR}" "${PROJECT_DIR}" -type f \( -name "full_table.tsv" -o -name "full_table*.tsv" -o -name "full_table*.txt" \) 2>/dev/null | head -n 1 || true)"
fi

if [[ -z "${BUSCO_FULL_TABLE}" || ! -s "${BUSCO_FULL_TABLE}" ]]; then
    echo "No existing BUSCO full_table found. Running BUSCO..."

    if ! command -v busco >/dev/null 2>&1; then
        echo "ERROR: busco command not found. Load BUSCO first, for example:" >&2
        echo "  module load busco" >&2
        echo "or activate the conda environment containing BUSCO, then rerun this script." >&2
        exit 1
    fi

    # BUSCO creates run-specific directories inside --out_path.
    busco \
        --in "${PROTEOME}" \
        --out "${BUSCO_RUN_NAME}" \
        --out_path "${BUSCO_OUTDIR}" \
        --mode "${BUSCO_MODE}" \
        --lineage_dataset "${LINEAGE}" \
        --cpu "${THREADS}" \
        --force

    BUSCO_FULL_TABLE="$(find "${BUSCO_OUTDIR}/${BUSCO_RUN_NAME}" -type f \( -name "full_table.tsv" -o -name "full_table*.tsv" -o -name "full_table*.txt" \) | head -n 1 || true)"
fi

if [[ -z "${BUSCO_FULL_TABLE}" || ! -s "${BUSCO_FULL_TABLE}" ]]; then
    echo "ERROR: BUSCO full table still not found after search/run." >&2
    exit 1
fi

echo "BUSCO full table: ${BUSCO_FULL_TABLE}"

# Copy concise BUSCO files to a stable path for downstream reproducibility.
cp -f "${BUSCO_FULL_TABLE}" "${BUSCO_STABLE_DIR}/full_table.tsv"

BUSCO_SHORT_SUMMARY="$(find "$(dirname "${BUSCO_FULL_TABLE}")" "${BUSCO_OUTDIR}/${BUSCO_RUN_NAME}" -maxdepth 3 -type f -name "short_summary*.txt" 2>/dev/null | head -n 1 || true)"
if [[ -n "${BUSCO_SHORT_SUMMARY}" && -s "${BUSCO_SHORT_SUMMARY}" ]]; then
    cp -f "${BUSCO_SHORT_SUMMARY}" "${BUSCO_STABLE_DIR}/short_summary.txt"
fi

BUSCO_FULL_TABLE_STABLE="${BUSCO_STABLE_DIR}/full_table.tsv"

# Reuse existing cluster file if present. The downstream script can regenerate exact clusters if absent.
CLUSTER_FILE="${DUP_OUTDIR}/equina_dark_candidate_clusters.exact.tsv"

echo "Rerunning duplication/prioritisation with BUSCO full table..."
python3 scripts/add_busco_duplication_context.py \
    --linked "${PRIMARY_LINKED}" \
    --lookup "${LOOKUP_TSV}" \
    --dark-fasta "${DARK_FASTA}" \
    --busco-full-table "${BUSCO_FULL_TABLE_STABLE}" \
    --cluster-file "${CLUSTER_FILE}" \
    --cluster-method exact_python_fallback \
    --outdir "${DUP_OUTDIR}"

echo "Rerunning figures with BUSCO-backed prioritisation..."
GENOME_LINKED_TSV="${PRIMARY_LINKED}" \
PRIORITISED_TSV="${DUP_OUTDIR}/equina_dark_candidates.prioritised.tsv" \
BUSCO_FULL_TABLE="${BUSCO_FULL_TABLE_STABLE}" \
bash scripts/make_figures.sh

echo
echo "BUSCO-backed validation complete. Check these files:"
echo "  ${BUSCO_FULL_TABLE_STABLE}"
echo "  ${BUSCO_STABLE_DIR}/short_summary.txt"
echo "  ${DUP_OUTDIR}/equina_duplication_context.summary.txt"
echo "  ${PROJECT_DIR}/08_figures/figure_tables/figure_14_priority_tier_counts.tsv"
