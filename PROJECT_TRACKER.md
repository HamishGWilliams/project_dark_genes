# Project tracker

## Current project state

### Inputs currently available
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.nucl.fa`
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa`

### What these files appear to be
The current input files are most likely an existing **paired CDS/protein gene set**, not a raw genome assembly.

Evidence collected so far:
- Sequence count in nucleotide FASTA: **55,607**
- Sequence count in amino acid FASTA: **55,607**
- Matching headers observed in both files, e.g. `evm.utg4.1`
- Mean protein length: **452.05 aa**
- Minimum protein length: **32 aa**
- Maximum protein length: **10,458 aa**

Interpretation:
- `nucl.fa` is currently being treated as **coding nucleotide sequences (CDS/transcripts)**
- `aa.fa` is currently being treated as the corresponding **predicted proteome**
- The project should therefore start with **proteome/gene-set QC and functional annotation**, not repeat masking from a raw assembly

---

## Completed so far
- [x] Defined the conceptual dark-gene analysis workflow for cnidarian genomes
- [x] Adapted the workflow for a single-species HPC implementation
- [x] Identified that chromosome-level assembly is not required to begin the project
- [x] Confirmed that available files are likely an existing gene set rather than a raw genome assembly
- [x] Verified that nucleotide and amino acid FASTA files contain the same number of entries (**55,607**)
- [x] Verified that FASTA headers match between nucleotide and amino acid files
- [x] Calculated basic protein-length summary statistics
- [x] Reframed phase 1 to focus on **gene-set/proteome QC** instead of genome repeat-masking

---

## In progress
- [X] Confirm CDS/protein pairing more rigorously using shared ID and translation-length checks
- [X] Generate QC tables for CDS lengths and protein lengths
- [X] Flag short proteins for review rather than immediate removal
- [X] Run BUSCO in **protein mode** on the amino acid FASTA

---

## Immediate next tasks

### 1. Confirm sequence pairing and translation plausibility
- [X] Check shared IDs across nucleotide and protein FASTA files
- [X] Check whether CDS lengths are divisible by 3
- [X] Check whether protein lengths are compatible with CDS lengths

### 2. Generate QC tables
- [X] `01_qc/protein_lengths.tsv`
- [X] `01_qc/cds_lengths.tsv`
- [X] `01_qc/proteins_lt50aa.tsv`
- [X] `01_qc/proteins_lt100aa.tsv`
- [X] `01_qc/proteins_lt150aa.tsv`

### 3. Run BUSCO on proteins
- [X] Run BUSCO with `-m proteins`
- [ ] Record lineage dataset used
- [ ] Save `short_summary` output
- [ ] Save the full BUSCO result directory for downstream interpretation

---

## Future annotation tasks

### Functional annotation
- [ ] Run DIAMOND or BLASTp against UniProt / reference proteomes
- [ ] Run InterProScan
- [ ] Run eggNOG-mapper
- [ ] Run SignalP for unresolved candidates
- [ ] Build a master annotation table with one row per gene/protein

### Dark-gene classification
- [ ] Classify proteins as annotated, sequence-dark, function-dark, or high-confidence dark candidates
- [ ] Separate genuinely dark proteins from low-confidence short ORFs and fragmented models

### Comparative and validation steps
- [ ] Download comparison proteomes from cnidarians and suitable outgroups
- [ ] Run OrthoFinder
- [ ] Evaluate lineage restriction versus broader conservation
- [ ] Run homology-detection-failure control where appropriate
- [ ] Prioritise top dark-gene candidates for structure/function follow-up

---

## Data still needed for stronger final conclusions
These are **not required to continue right now**, but they will strengthen later interpretation.

- [ ] Genome assembly FASTA (scaffolds/contigs)
- [ ] GFF3 annotation file
- [ ] RNA-seq read evidence or alignments
- [ ] Information on how the EVM-derived gene set was generated and filtered

These additional files will allow:
- repeat-overlap assessment
- contamination and genome-context checks
- genomic neighbourhood analysis
- stronger filtering of artefactual dark-gene calls

---

## BUSCO outputs needed to proceed
To interpret the current gene set and decide how strict the downstream dark-gene filtering should be, collect the following from the BUSCO run:

### Minimum required
- [ ] BUSCO version used
- [ ] BUSCO mode used (`proteins`)
- [ ] Lineage dataset used (for example `metazoa_odb10`, `metazoa_odb12`, or another lineage)
- [ ] Short summary values:
  - [ ] Complete (C)
  - [ ] Single-copy (S)
  - [ ] Duplicated (D)
  - [ ] Fragmented (F)
  - [ ] Missing (M)
  - [ ] Total BUSCO groups searched (`n`)

### Very useful
- [ ] Path to the BUSCO output directory
- [ ] `short_summary*.txt` or `short_summary*.json`
- [ ] `full_table.tsv`
- [ ] Any warnings about unusually high duplication, fragmentation, or failed searches

### Why this matters
BUSCO will help determine whether the current proteome is:
- reasonably complete
- strongly duplicated
- fragmented
- potentially over-predicted

That directly affects how cautious we need to be when interpreting proteins with no annotation as candidate dark genes.

---

## Suggested working directories on the HPC

```text
project_dark_genes/
├── 00_raw/
├── 01_qc/
├── 02_annotation/
├── 03_dark_candidates/
├── 04_reports/
├── logs/
└── scripts/
```

---

## Notes
This tracker should be updated after every major step so that the distinction between:
- completed work
- active work
- blocked steps
- future validation

remains clear throughout the project.
