#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict

lookup_file = Path("01_qc/isoforms/gff3_mrna_gene_lookup.tsv")
length_file = Path("01_qc/isoforms/interproscan_cleaned_protein_lengths.tsv")

out_reps = Path("01_qc/isoforms/representative_longest_protein_per_gene.tsv")
out_ids = Path("01_qc/isoforms/representative_longest_protein_ids.txt")
out_excluded = Path("01_qc/isoforms/excluded_nonrepresentative_isoforms.tsv")

mrna_to_gene = {}
with lookup_file.open() as fh:
    next(fh)
    for line in fh:
        gene_id, mrna_id, seqid, start, end, strand, feature_type = line.rstrip("\n").split("\t")
        mrna_to_gene[mrna_id] = gene_id

lengths = {}
with length_file.open() as fh:
    next(fh)
    for line in fh:
        mrna_id, length = line.rstrip("\n").split("\t")
        lengths[mrna_id] = int(length)

gene_to_mrnas = defaultdict(list)

# Include all proteins. If an ID is not in the GFF3 lookup, treat it as its own gene-like unit.
for mrna_id, length in lengths.items():
    gene_id = mrna_to_gene.get(mrna_id, mrna_id)
    gene_to_mrnas[gene_id].append((mrna_id, length))

representatives = {}
excluded = []

for gene_id, records in gene_to_mrnas.items():
    # Sort by longest length, then stable ID order.
    records_sorted = sorted(records, key=lambda x: (-x[1], x[0]))
    rep_id, rep_len = records_sorted[0]
    representatives[gene_id] = (rep_id, rep_len, len(records_sorted))

    for mrna_id, length in records_sorted[1:]:
        excluded.append((gene_id, mrna_id, length, rep_id, rep_len))

with out_reps.open("w") as out:
    out.write("gene_id\trepresentative_mrna_id\trepresentative_protein_length\tisoform_count\n")
    for gene_id, (rep_id, rep_len, iso_count) in sorted(representatives.items()):
        out.write(f"{gene_id}\t{rep_id}\t{rep_len}\t{iso_count}\n")

with out_ids.open("w") as out:
    for gene_id, (rep_id, rep_len, iso_count) in sorted(representatives.items()):
        out.write(rep_id + "\n")

with out_excluded.open("w") as out:
    out.write("gene_id\texcluded_mrna_id\texcluded_protein_length\trepresentative_mrna_id\trepresentative_protein_length\n")
    for row in sorted(excluded):
        out.write("\t".join(map(str, row)) + "\n")

print(f"Genes/loci represented: {len(representatives)}")
print(f"Representative proteins: {len(representatives)}")
print(f"Excluded non-representative isoforms: {len(excluded)}")
print(f"Wrote: {out_reps}")
print(f"Wrote: {out_ids}")
print(f"Wrote: {out_excluded}")
