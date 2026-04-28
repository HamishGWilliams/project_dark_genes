# Project tracker

Last updated: 2026-04-28

## Current project state

The project has now completed the first full representative master annotation compilation for the *Actinia equina* predicted proteome. The current annotation outputs are considered robust and consistent with the methods defined so far: representative protein selection, no-stop FASTA preparation, threshold-filtered DIAMOND/BLASTp homology evidence, InterProScan, eggNOG-mapper, SignalP, and hierarchical master-table classification.

The project is now ready to move from annotation compilation into candidate dark-gene extraction, genome-context annotation, and RNA-seq-supported prioritisation.

---

## Inputs currently available

- `equina_smart.rnam-trna.merged.ggf.curated.remredun.nucl.fa`
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa`
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.trans.fa.gz`
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.proteins.gff3.gz`
- `equina_smartden.arrow4.noredun.fa.gz`
- RNA-seq data from a multi-stressor experiment including control treatments
- Experimental design notes file in the repository `notes/` directory describing treatment conditions and total sample size

### Interpretation of current input files

The current input files form a coherent annotation set plus genome assembly:

- `nucl.fa` is being treated as CDS-like nucleotide sequences
- `aa.fa` is being treated as the corresponding predicted proteome
- `trans.fa.gz` is transcript FASTA matching the EVM-style IDs
- `proteins.gff3.gz` is the annotation GFF3 linking genes/transcripts/proteins to genomic coordinates
- `arrow4.noredun.fa.gz` is the genome assembly FASTA
- RNA-seq data can provide expression evidence, splice support, and treatment responsiveness

### Evidence collected so far

- Sequence count in nucleotide FASTA: **55,607**
- Sequence count in amino acid FASTA: **55,607**
- Matching headers observed in nucleotide, amino acid, and transcript FASTA files, e.g. `evm.utg4.1`
- GFF3 mRNA IDs match the FASTA IDs, e.g. `ID=evm.utg4.1`
- GFF3 seqids match genome contig names, e.g. `utg4` in the GFF3 and `>utg4`-style headers in the genome FASTA
- Mean protein length: **452.05 aa**
- Minimum protein length: **32 aa**
- Maximum protein length: **10,458 aa**
- RNA-seq libraries are **paired-end**
- RNA-seq libraries are **unstranded**
- Experimental design information is available in a text file under `notes/`

### Current interpretation

- The project is now a **genome-aware and expression-aware dark-gene workflow**, not a proteome-only survey.
- Candidate dark genes can be mapped back to genomic coordinates and evaluated in structural context.
- RNA-seq can be used to distinguish expressed dark candidates from unsupported predictions and to prioritise stress-responsive candidates.
- Because BUSCO duplication is high, downstream dark-gene calls remain conservative and are based on the representative protein set rather than the redundant original proteome.

---

## Completed so far

### Core setup and QC

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
- [x] Confirmed RNA-seq library layout is paired-end and unstranded
- [x] Located an experiment notes file describing treatment conditions and sample size

### Representative proteome preparation

- [x] Generated a representative protein set by selecting the longest representative protein per gene/isoform group
- [x] Removed terminal stop characters from representative protein sequences for compatibility with InterProScan and downstream tools
- [x] Used the representative no-stop proteome for the annotation workflow

### Annotation software setup

- [x] Installed and tested InterProScan
- [x] Installed and tested eggNOG-mapper after manually resolving database download issues
- [x] Confirmed SignalP module availability and successfully tested `signalp/5.0b`

### Homology annotation revision

- [x] Revised DIAMOND/BLASTp homology handling to use `e-value <= 1e-5`
- [x] Updated DIAMOND output fields to retain `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen`
- [x] Updated BLASTp output fields to retain `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen`
- [x] Generated raw DIAMOND outputs locally/HPC-only
- [x] Generated filtered DIAMOND outputs locally/HPC-only
- [x] Generated filtered DIAMOND top-hit files after applying `e-value <= 1e-5`
- [x] Generated raw BLASTp outputs locally/HPC-only
- [x] Generated filtered BLASTp outputs locally/HPC-only
- [x] Generated filtered BLASTp best-hit files after applying `e-value <= 1e-5`
- [x] Confirmed filtered homology files are suitable for master annotation compilation
- [x] Treated Swiss-Prot and TrEMBL evidence differently in interpretation, with TrEMBL as lower-confidence support

### Test100 annotation checkpoint

- [x] Built a 100-protein test master annotation table
- [x] Built and ran a 100-protein master-table QC script
- [x] Updated and QC-tested the 100-protein workflow after introducing the `e-value <= 1e-5` threshold
- [x] Confirmed the 100-protein test workflow passed after adding expanded homology fields and filtered source files

### Full annotation checkpoint

- [x] Ran full representative DIAMOND Swiss-Prot annotation using the revised threshold strategy
- [x] Ran full representative DIAMOND Cnidaria-TrEMBL annotation using the revised threshold strategy
- [x] Ran full representative BLASTp Swiss-Prot annotation using the revised threshold strategy
- [x] Ran full representative BLASTp Cnidaria-TrEMBL annotation using the revised threshold strategy
- [x] Confirmed full representative InterProScan outputs are available for master annotation compilation
- [x] Confirmed full representative eggNOG-mapper outputs are available for master annotation compilation
- [x] Confirmed full representative SignalP outputs are available for master annotation compilation
- [x] Built the first full representative master annotation table
- [x] Ran and checked full master annotation QC
- [x] Confirmed that the current full master annotation results are robust and true to the methods outlined so far

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

- The proteome appears **highly complete**.
- The very high duplicated BUSCO fraction suggests caution about gene-set redundancy, isoform retention, haplotig effects, or overprediction.
- Downstream dark-gene calling should therefore be conservative, especially for short proteins and singletons.
- This concern is partly addressed by using a representative longest-protein-per-gene/isoform-group set for the annotation workflow.

---

## Annotation threshold revision checkpoint

Detailed working note: `notes/annotation_threshold_revision_plan.md`

### Rationale

The original DIAMOND/BLASTp workflow generated top-hit files before applying a formal hit-quality threshold. This was acceptable for pipeline testing but not for final dark-gene classification, because weak sequence-similarity hits could incorrectly rescue proteins from the dark/unresolved category.

The revised homology threshold is aligned with the coral dark-gene approach from Stephens et al. and uses:

```text
e-value <= 1e-5
```

This threshold is now applied to DIAMOND and BLASTp hits before selecting top/best hits for the master annotation table.

### Revised sequence-similarity handling

- [x] Use `e-value <= 1e-5` as the primary sequence-similarity inclusion threshold.
- [x] Retain percent identity, alignment length, bitscore, query length, subject length, query coverage, and subject coverage as QC/reporting fields where possible.
- [x] Do not use identity or coverage as hard filters for the main annotation classes unless a separate high-confidence subcategory is defined later.
- [x] Where subject descriptions are available, especially for TrEMBL, flag ambiguous descriptions such as:
  - `uncharacterized protein`
  - `uncharacterised protein`
  - `hypothetical protein`
  - `predicted protein`
  - `expressed protein`
  - `unnamed protein product`
  - `unknown protein`
- [x] Treat TrEMBL evidence as lower-confidence than Swiss-Prot.
- [x] Avoid treating ambiguous TrEMBL hits as equivalent to reviewed Swiss-Prot annotation.

### Revised homology-output format

DIAMOND/BLASTp runs now retain:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
```

Where possible, subject titles/descriptions should also be retained or mapped back later for ambiguity review, especially for Cnidaria-TrEMBL hits.

### Revised annotation hierarchy

The master annotation table assigns one primary class per representative protein in this order:

1. `annotated_swissprot_supported`
2. `sequence_supported_trembl_cnidaria`
3. `domain_supported_interpro`
4. `orthology_supported_eggnog`
5. `function_dark_but_signalp_secretory_candidate`
6. `function_dark_no_current_annotation`

Only threshold-passing DIAMOND/BLASTp evidence is allowed to assign categories 1 or 2. TrEMBL support is interpreted conservatively and should not be treated as equivalent to reviewed Swiss-Prot evidence.

---

## Current annotation deliverables

Generated locally/HPC-only unless explicitly added as lightweight summaries:

```text
02_annotation/diamond/raw/
02_annotation/diamond/filtered/
02_annotation/blastp/raw/
02_annotation/blastp/filtered/
02_annotation/interproscan/raw/
02_annotation/eggnog/raw/
02_annotation/signalp/raw/
02_annotation/signalp/summary/
02_annotation/master/test100/equina_representative_test100.master_annotation.tsv
02_annotation/master/test100/equina_representative_test100.master_annotation.summary.txt
02_annotation/master/test100/qc/
02_annotation/master/full/equina_representative_full.master_annotation.tsv
02_annotation/master/full/equina_representative_full.master_annotation.summary.txt
02_annotation/master/full/qc/
```

Large generated outputs should remain local/HPC-only and should not be committed to GitHub. Scripts, notes, lightweight summaries, and tracker updates should be committed.

---

## In progress / active follow-up

- [ ] Preserve the current full annotation checkpoint in the repository using lightweight summaries only
- [ ] Update QC scripts to explicitly report TrEMBL ambiguous/uncharacterised hit counts
- [ ] Confirm that ambiguous TrEMBL hits do not independently drive `sequence_supported_trembl_cnidaria` classification
- [ ] Record final full master annotation summary counts and proportions in a lightweight report
- [ ] Decide whether to rerun or remap homology outputs to retain `stitle`/subject descriptions for all relevant TrEMBL hits

---

## Immediate next tasks

### 1. Preserve and document the current full annotation checkpoint

- [ ] Record full master annotation class counts and proportions in `PROJECT_TRACKER.md` or a lightweight report
- [ ] Record source-coverage summary counts from the full QC report
- [ ] Record TrEMBL ambiguity summary counts once the updated QC script has been run
- [ ] Confirm `.gitignore` excludes large generated annotation outputs
- [ ] Commit only scripts, notes, tracker updates, and lightweight summaries

### 2. Extract candidate dark genes

- [ ] Extract proteins classified as `function_dark_no_current_annotation`
- [ ] Extract proteins classified as `function_dark_but_signalp_secretory_candidate`
- [ ] Create dark-candidate TSV files
- [ ] Create dark-candidate FASTA files
- [ ] Retain non-classifying or ambiguous evidence for manual review
- [ ] Summarise candidate counts and proportions

### 3. Start genome-aware validation

- [ ] Parse the GFF3 into a gene-to-transcript-to-protein lookup table
- [ ] Map candidate dark genes back to genome coordinates
- [ ] Summarise exon count, CDS span, transcript span, and scaffold location for each candidate
- [ ] Use genome context to help distinguish plausible genes from suspicious models
- [ ] Run or integrate RepeatModeler/RepeatMasker outputs to identify repeats
- [ ] Identify TE-associated candidates and flag them conservatively

### 4. Integrate RNA-seq evidence

- [ ] Convert the `notes/` experiment file into a tabular RNA-seq sample sheet with sample IDs, treatment/control labels, replicate IDs, and file paths
- [ ] Perform RNA-seq read QC
- [ ] Align paired-end unstranded reads to the genome or transcriptome
- [ ] Quantify expression at transcript/gene level
- [ ] Identify which candidate dark genes are detectably expressed in any condition
- [ ] Test for differential expression between controls and stressor treatments
- [ ] Prioritise dark candidates with reproducible expression and/or stress responsiveness

### 5. Comparative and validation steps

- [ ] Download comparison proteomes from cnidarians and suitable outgroups
- [ ] Run OrthoFinder
- [ ] Evaluate lineage restriction versus broader conservation
- [ ] Run homology-detection-failure controls where appropriate
- [ ] Prioritise top dark-gene candidates for structure/function follow-up

---

## Dark-gene classification tasks

- [ ] Classify proteins as annotated, sequence-supported, function-dark, or high-confidence dark candidates based on the full master table
- [ ] Separate genuinely dark proteins from low-confidence short ORFs and fragmented models
- [ ] Extract `function_dark_no_current_annotation` candidates
- [ ] Extract `function_dark_but_signalp_secretory_candidate` candidates
- [ ] Create dark-candidate TSV and FASTA files
- [ ] Summarise dark-candidate counts and proportions
- [ ] Record best subthreshold or ambiguous sequence hits, if any, for later manual review

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
│   ├── diamond/
│   │   ├── raw/
│   │   └── filtered/
│   ├── blastp/
│   │   ├── raw/
│   │   └── filtered/
│   ├── interproscan/
│   ├── eggnog/
│   ├── signalp/
│   └── master/
│       ├── test100/
│       └── full/
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
