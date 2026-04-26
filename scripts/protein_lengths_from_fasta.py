#!/usr/bin/env python3

from pathlib import Path

inp = Path("02_annotation/interproscan/input/equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.no_stop.fa")
out = Path("01_qc/isoforms/interproscan_cleaned_protein_lengths.tsv")

seq_id = None
chunks = []

with inp.open() as fin, out.open("w") as fout:
    fout.write("mrna_id\tprotein_length\n")

    def flush():
        if seq_id is None:
            return
        seq = "".join(chunks).replace(" ", "").replace("\t", "")
        fout.write(f"{seq_id}\t{len(seq)}\n")

    for line in fin:
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            seq_id = line[1:].strip().split()[0]
            chunks = []
        else:
            chunks.append(line)

    flush()

print(f"Wrote: {out}")
