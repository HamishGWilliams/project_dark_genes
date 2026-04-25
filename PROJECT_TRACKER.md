# Project tracker

## Current project state

### Inputs currently available
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.nucl.fa`
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa`
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.trans.fa.gz`
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.proteins.gff3.gz`
- `equina_smartden.arrow4.noredun.fa.gz`
- RNA-seq data from a multi-stressor experiment including control treatments

### What these files appear to be
The current input files now form a coherent annotation set plus genome assembly:
- `nucl.fa` is being treated as CDS-like nucleotide sequences
- `aa.fa` is being treated as the corresponding predicted proteome
- `trans.fa.gz` is transcript FASTA matching the EVM-style IDs
- `proteins.gff3.gz` is the annotation GFF3 linking genes/transcripts/proteins to genomic coordinates
- `arrow4.noredun.fa.gz` is the genome assembly FASTA
- RNA-seq data can provide expression evidence, splice support, and treatment responsiveness

Evidence collected so far:
- Sequence count in nucleotide FASTA: **55,607**
- Sequence count in amino acid FASTA: **55,607**
- Matching headers observed in nucleotide, amino acid, and transcript FASTA files, e.g. `evm.utg4.1`
- GFF3 mRNA IDs match the FASTA IDs, e.g. `ID=evm.utg4.1`
- GFF3 seqids match genome contig names, e.g. `utg4` in the GFF3 and `>utg4`-style headers in the genome FASTA
- Mean protein length: **452.05 aa**
- Minimum protein length: **32 aa**
- Maximum protein length: **10,458 aa**

Interpretation:
- The project can now proceed as a **genome-aware and expression-aware dark-gene workflow** rather than a proteome-only survey
- Candidate dark genes can later be mapped back to genomic coordinates and evaluated in structural context
- RNA-seq can be used to distinguish expressed dark candidates from unsupported predictions and to prioritise stress-responsive candidates

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
- [x] Confirmed transcript FASTA IDs match the protein/CDS IDs
- [x] Confirmed GFF3 transcript IDs match the FASTA IDs
- [x] Confirmed GFF3 seqids correspond to genome assembly contig names
- [x] Located genome assembly FASTA
- [x] Located GFF3 annotation file
- [x] Located transcript FASTA
- [x] Established that RNA-seq data are available for downstream validation and prioritisation

---

## In progress
- [x] Confirm CDS/protein pairing more rigorously using shared ID and translation-length checks
- [x] Generate QC tables for CDS lengths and protein lengths
- [x] Flag short proteins for review rather than immediate removal
- [x] Run BUSCO in **protein mode** on the amino acid FASTA

---

## BUSCO checkpoint

### Run details
- [x] BUSCO version used: **v5.3.2**
- [x] BUSCO mode used: **proteins**
- [x] Lineage dataset used: **metazoa_odb10** (2024-01-08)
- [x] BUSCO output directory identified: `/uoa/scratch/users/r02hw22/project_dark_genes/equina_busco_proteins`

### Short summary values
- [x] Complete (C): **95.6%**
- [x] Single-copy (S): **52.9%**
- [x] Duplicated (D): **42.7%**
- [x] Fragmented (F): **2.6%**
- [x] Missing (M): **1.8%**
- [x] Total BUSCO groups searched (`n`): **954**
- [x] Complete BUSCOs: **912**
- [x] Single-copy BUSCOs: **505**
- [x] Duplicated BUSCOs: **407**
- [x] Fragmented BUSCOs: **25**
- [x] Missing BUSCOs: **17**

### Interpretation
- The proteome appears **highly complete**
- The very high duplicated BUSCO fraction suggests caution about gene-set redundancy, isoform retention, haplotig effects, or overprediction
- Downstream dark-gene calling should therefore be conservative, especially for short proteins and singletons

---

## Immediate next tasks

### 1. Archive BUSCO outputs for interpretation
- [x] Record lineage dataset used
- [ ] Save `short_summary` output into the project repository as a lightweight text summary
- [ ] Save or inspect `full_table.tsv` for duplicated/fragmented BUSCO patterns

### 2. Start functional annotation
- [ ] Run DIAMOND or BLASTp against UniProt / reference proteomes
- [ ] Run InterProScan
- [ ] Run eggNOG-mapper
- [ ] Run SignalP for unresolved candidates
- [ ] Build a master annotation table with one row per gene/protein

### 3. Start genome-aware validation
- [ ] Parse the GFF3 into a gene-to-transcript-to-protein lookup table
- [ ] Map candidate dark genes back to genome coordinates
- [ ] Summarise exon count, CDS span, transcript span, and scaffold location for each candidate
- [ ] Use genome context later to help distinguish plausible genes from suspicious models

### 4. Integrate RNA-seq evidence
- [ ] Organise RNA-seq metadata table with sample IDs, treatment, control/stressor labels, replicate IDs, and file paths
- [ ] Run raw read QC and adapter/quality trimming if needed
- [ ] Align reads to the genome or transcriptome
- [ ] Quantify expression at transcript/gene level
- [ ] Identify which candidate dark genes are detectably expressed in any condition
- [ ] Test for differential expression between controls and stressor treatments
- [ ] Prioritise dark candidates with reproducible expression and/or stress responsiveness

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

### Expression-aware validation
- [ ] Require evidence thresholds for expression-supported dark genes
- [ ] Distinguish constitutive expression from condition-specific expression
- [ ] Evaluate whether stress-responsive dark genes are enriched among high-confidence candidates

---

## Data still useful for stronger final conclusions
These are not blocking now, but would still help later.

- [ ] Information on how the EVM-derived gene set was generated and filtered
- [ ] If available, splice-aware read alignments or existing transcript assemblies from the same RNA-seq experiment

These additional data will allow:
- stronger filtering of artefactual dark-gene calls
- more confident interpretation of duplicated BUSCOs and short ORFs
- direct comparison between original annotation evidence and current expression evidence

---

## Suggested working directories on the HPC

```text
project_dark_genes/
├── 00_raw/
├── 01_qc/
├── 02_annotation/
├── 03_dark_candidates/
├── 04_reports/
├── 05_rnaseq/
│   ├── 00_metadata/
│   ├── 01_qc/
│   ├── 02_trimmed/
│   ├── 03_alignment/
│   ├── 04_quant/
│   └── 05_de/
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
