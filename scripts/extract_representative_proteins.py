#!/usr/bin/env python3

from pathlib import Path

fasta = Path("02_annotation/interproscan/input/equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.no_stop.fa")
ids_file = Path("01_qc/isoforms/representative_longest_protein_ids.txt")
out = Path("02_annotation/interproscan/input/equina_representative_longest_per_gene.no_stop.fa")

keep_ids = set(line.strip() for line in ids_file.open() if line.strip())

written = 0
keep = False

with fasta.open() as fin, out.open("w") as fout:
    for line in fin:
        if line.startswith(">"):
            seq_id = line[1:].strip().split()[0]
            keep = seq_id in keep_ids
            if keep:
                written += 1
                fout.write(line)
        else:
            if keep:
                fout.write(line)

print(f"Representative IDs requested: {len(keep_ids)}")
print(f"Representative FASTA entries written: {written}")
print(f"Wrote: {out}")

if written != len(keep_ids):
    raise SystemExit("ERROR: Not all representative IDs were found in the FASTA.")
