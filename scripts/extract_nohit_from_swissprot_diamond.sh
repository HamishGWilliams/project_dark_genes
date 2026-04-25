#!/bin/bash
#SBATCH --mem=32G
#SBATCH --partition=uoa-compute
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=r02hw22@abdn.ac.uk
#SBATCH --time=02:00:00
#SBATCH --output=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/outputs/%x_%j.out
#SBATCH --error=/uoa/home/r02hw22/sharedscratch/project_dark_genes/logs/errors/%x_%j.err

PROJECT_DIR="/uoa/home/r02hw22/sharedscratch/project_dark_genes"
cd "$PROJECT_DIR"

: "${QUERY:?Need QUERY}"
: "${SWISSPROT_TOP:?Need SWISSPROT_TOP}"
: "${NOHIT_FASTA:?Need NOHIT_FASTA}"
: "${NOHIT_LIST:?Need NOHIT_LIST}"

python3 - <<'PY'
from pathlib import Path

query = Path("00_raw/equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa")
top_hits = Path("02_annotation/diamond/equina_vs_swissprot_all.diamond.top_hits.tsv")
nohit_fasta = Path("02_annotation/intermediate/equina_no_swissprot_diamond_hit.fa")
nohit_list = Path("02_annotation/intermediate/equina_no_swissprot_diamond_hit.ids.txt")

hit_ids = set()
with top_hits.open() as fh:
    for line in fh:
        if line.strip():
            hit_ids.add(line.split('\t', 1)[0])

with query.open() as fh_in, nohit_fasta.open('w') as fh_out, nohit_list.open('w') as fh_ids:
    keep = False
    current_id = None
    for line in fh_in:
        if line.startswith('>'):
            current_id = line[1:].strip().split()[0]
            keep = current_id not in hit_ids
            if keep:
                fh_out.write(line)
                fh_ids.write(current_id + '\n')
        else:
            if keep:
                fh_out.write(line)

print(f"Swiss-Prot DIAMOND hits: {len(hit_ids)}")
print(f"No-hit FASTA written to: {nohit_fasta}")
print(f"No-hit ID list written to: {nohit_list}")
PY

wc -l "$NOHIT_LIST"
grep -c '^>' "$NOHIT_FASTA"
head "$NOHIT_LIST" || true
