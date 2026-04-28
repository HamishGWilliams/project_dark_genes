# project_dark_genes

Exploring and understanding how to identify novel **“dark” genes** in a non-model cnidarian dataset using *Actinia equina* sequence resources, genome annotation files, and RNA-seq data on an HPC workflow.

## What I hope to develop

Although this project is focused on *Actinia equina* for the time being, once I have validated this pipeline and shown it works for this species, I would like to transform this repository into a space to walk people through various bioinformatic methods and techniques towards dark-gene novelty. The secondary goal is to provide a clear scaffold for others, especially early-career researchers, to develop their own workflows while understanding what each method, threshold, and filtering step contributes.

## Project aim

The goal of this project is to identify and prioritise **candidate dark genes** in *Actinia equina* by combining:

- an existing predicted gene/protein set
- genome assembly and GFF3 annotation
- functional annotation tools
- expression evidence from RNA-seq
- downstream comparative and validation steps

In this project, a **dark gene** is treated as a gene/protein that remains unresolved after layered annotation and quality control, rather than simply a sequence with “no BLAST hit”.

---

## Current project state

### Core input files currently available

- `equina_smart.rnam-trna.merged.ggf.curated.remredun.nucl.fa`
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa`
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.trans.fa.gz`
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.proteins.gff3.gz`
- `equina_smartden.arrow4.noredun.fa.gz`

### RNA-seq resources currently available

- paired-end RNA-seq data from a **multi-stressor experiment**
- includes **control treatments**
- libraries are **unstranded**
- experimental design notes are available in the `notes/` directory

### What has been established so far

- the nucleotide and amino acid FASTA files each contain **55,607** entries
- IDs match across the CDS/protein/transcript files, e.g. `evm.utg4.1`
- GFF3 `mRNA` IDs match the FASTA IDs
- GFF3 seqids match genome contig names
- mean protein length is **452.05 aa**
- protein length range is **32 aa** to **10,458 aa**
- the available data support a **genome-aware and expression-aware dark-gene workflow**
- a representative no-stop protein set has been generated to reduce redundancy from isoforms
- the 100-protein test group successfully passed the revised annotation compilation workflow
- the first full representative master annotation table has now been compiled and checked
- the current full master annotation results are considered robust and consistent with the methods defined so far

---

## BUSCO checkpoint

BUSCO was run on the predicted proteome in **protein mode**.

### BUSCO run details

- **BUSCO version:** v5.3.2
- **mode:** proteins
- **lineage:** `metazoa_odb10`
- **output directory:** `/uoa/scratch/users/r02hw22/project_dark_genes/equina_busco_proteins`

### BUSCO summary

- **C:** 95.6%
- **S:** 52.9%
- **D:** 42.7%
- **F:** 2.6%
- **M:** 1.8%
- **n:** 954

Detailed counts:

- **912** complete BUSCOs
- **505** complete and single-copy
- **407** complete and duplicated
- **25** fragmented
- **17** missing

### Interpretation

The proteome appears **highly complete**, which supports downstream annotation and candidate dark-gene discovery.

However, the **high duplicated BUSCO fraction (42.7%)** suggests caution about:

- annotation redundancy
- retained isoforms
- haplotig-related duplication
- overprediction

This means downstream dark-gene classification should be **conservative**, particularly for:

- short proteins
- singletons
- low-support loci
- proteins lacking expression support

---

## Annotation checkpoint

The first full representative master annotation compilation has been completed.

The current annotation strategy integrates:

- representative no-stop protein FASTA as the row universe
- representative gene/isoform lookup information
- DIAMOND Swiss-Prot homology evidence
- DIAMOND Cnidaria-TrEMBL homology evidence
- BLASTp Swiss-Prot homology evidence
- BLASTp Cnidaria-TrEMBL homology evidence
- InterProScan domain/family/GO/pathway evidence
- eggNOG-mapper orthology and functional annotation evidence
- SignalP secretory-signal prediction

### Homology filtering currently used

DIAMOND and BLASTp evidence is filtered using:

```text
E-value <= 1e-5
```

The homology output format has been expanded to retain useful interpretive fields:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
```

Where available, subject-title/description information is also retained for TrEMBL ambiguity checks.

### TrEMBL interpretation

TrEMBL evidence is treated conservatively and is **not considered equivalent to Swiss-Prot**.

Cnidaria-TrEMBL hits are used as lower-confidence sequence-support evidence only when they are threshold-passing and informative. Ambiguous descriptions such as hypothetical, uncharacterized/uncharacterised, predicted, unnamed, or unknown proteins are flagged for review and should not independently rescue a protein from dark-gene status.

---

## What has been completed so far

- defined the conceptual dark-gene workflow for cnidarian genomes
- adapted the workflow for a single-species HPC implementation
- confirmed that chromosome-level assembly is **not required** to begin
- confirmed that the available FASTA files represent an existing **gene/protein set**
- verified ID consistency across CDS, protein, transcript, and GFF3 resources
- confirmed linkage between GFF3 coordinates and genome contigs
- generated basic QC summaries for the proteome
- flagged short proteins for later review
- completed a BUSCO checkpoint on the protein set
- generated a representative no-stop proteome for downstream annotation
- identified RNA-seq resources for downstream validation and prioritisation
- defined a **hierarchical UniProt strategy**: Swiss-Prot first, Cnidaria-TrEMBL second
- installed and tested InterProScan, eggNOG-mapper, and SignalP
- built and QC-tested a 100-protein master annotation table
- successfully updated and QC-tested the 100-protein DIAMOND/BLASTp homology workflow with `e-value <= 1e-5`, expanded outfmt 6 fields including `qlen` and `slen`, post-processing e-value filtering, and coverage-field propagation into the master table
- regenerated full representative DIAMOND and BLASTp outputs using the revised threshold-filtered homology strategy
- compiled and checked the first full representative master annotation table
- confirmed that the current full master annotation results are robust and consistent with the agreed annotation hierarchy
- recorded the revised annotation-threshold plan in `notes/annotation_threshold_revision_plan.md`

---

## Immediate next steps

### 1. Preserve the current full annotation checkpoint

- record the full master annotation summary counts in `PROJECT_TRACKER.md`
- retain lightweight QC summaries and interpretation notes in the repository
- keep large raw/filtered annotation output files excluded from Git
- document any manual interpretation decisions made during this checkpoint

### 2. Extract candidate dark genes

- extract proteins classified as `function_dark_no_current_annotation`
- extract proteins classified as `function_dark_but_signalp_secretory_candidate`
- create dark-candidate TSV files
- create dark-candidate FASTA files
- retain non-classifying or ambiguous evidence for manual review

### 3. Build genome-aware lookup tables

- parse the GFF3 into a gene → transcript → protein → contig lookup table
- link candidate dark genes back to genomic coordinates
- summarise exon count, CDS span, transcript span, and scaffold location
- evaluate repeat/TE context using RepeatModeler and RepeatMasker outputs when available

### 4. Integrate RNA-seq evidence

- convert the experiment notes into a structured sample sheet
- perform RNA-seq QC
- align paired-end unstranded reads
- quantify expression
- identify expressed candidate dark genes
- test for differential expression between controls and stressor treatments

### 5. Prioritise high-confidence dark candidates

For each candidate gene/protein, integrate:

- length/QC information
- absence of threshold-passing Swiss-Prot similarity
- absence of threshold-passing informative Cnidaria-TrEMBL similarity
- absence of InterProScan domain/family support
- absence of eggNOG orthology-based functional support
- SignalP status
- genomic context
- repeat/TE context
- expression support
- differential expression status

---

## Current workflow overview

```mermaid
flowchart TD
    A[Project start] --> B[Inventory available files]

    B --> C[Existing annotation resources identified]
    C --> C1[CDS FASTA<br/>nucl.fa]
    C --> C2[Protein FASTA<br/>aa.fa]
    C --> C3[Transcript FASTA<br/>trans.fa.gz]
    C --> C4[GFF3 annotation<br/>proteins.gff3.gz]
    C --> C5[Genome assembly<br/>arrow4.noredun.fa.gz]
    C --> C6[RNA-seq data<br/>paired-end and unstranded]
    C --> C7[Experiment notes<br/>treatments and sample size]

    C1 --> D[Confirm paired CDS and protein set]
    C2 --> D
    C3 --> E[Confirm transcript IDs match protein and CDS IDs]
    C4 --> F[Confirm GFF3 mRNA IDs match FASTA IDs]
    C5 --> G[Confirm GFF3 contigs match genome contig names]

    D --> H[Generate proteome QC tables]
    H --> H1[CDS lengths]
    H --> H2[Protein lengths]
    H --> H3[Flag short proteins for review]

    H --> I[BUSCO proteins mode]
    I --> I1[C 95.6 percent]
    I --> I2[S 52.9 percent]
    I --> I3[D 42.7 percent]
    I --> I4[F 2.6 percent]
    I --> I5[M 1.8 percent]

    I --> J{Proteome suitable for annotation?}
    J -->|Yes, but high duplication requires caution| K[Generate representative proteome]
    K --> K1[Select representative protein per gene or isoform group]
    K1 --> K2[Remove terminal stop characters]
    K2 --> K3[Representative no-stop proteome]

    K3 --> L1[DIAMOND Swiss-Prot search]
    K3 --> L2[DIAMOND Cnidaria-TrEMBL search]
    K3 --> L3[BLASTp Swiss-Prot search]
    K3 --> L4[BLASTp Cnidaria-TrEMBL search]

    L1 --> M1[Filter hits<br/>e-value <= 1e-5<br/>retain qlen and slen]
    L2 --> M2[Filter hits<br/>e-value <= 1e-5<br/>retain qlen, slen and descriptions]
    L3 --> M3[Filter hits<br/>e-value <= 1e-5<br/>retain qlen and slen]
    L4 --> M4[Filter hits<br/>e-value <= 1e-5<br/>retain qlen, slen and descriptions]

    M1 --> N1[Filtered Swiss-Prot top hits]
    M2 --> N2[Filtered Cnidaria-TrEMBL top hits]
    M3 --> N3[Filtered BLASTp Swiss-Prot top hits]
    M4 --> N4[Filtered BLASTp Cnidaria-TrEMBL top hits]

    N2 --> N5[Flag ambiguous TrEMBL descriptions]
    N4 --> N6[Flag ambiguous TrEMBL descriptions]

    K3 --> O1[InterProScan]
    K3 --> O2[eggNOG-mapper]
    K3 --> O3[SignalP 5.0b]

    O1 --> P1[Domain, family, GO and pathway evidence]
    O2 --> P2[Orthology, COG, GO, KEGG and Pfam evidence]
    O3 --> P3[Signal peptide evidence]

    N1 --> Q[Master annotation compiler]
    N2 --> Q
    N3 --> Q
    N4 --> Q
    N5 --> Q
    N6 --> Q
    P1 --> Q
    P2 --> Q
    P3 --> Q

    Q --> Q1[Test100 master table]
    Q1 --> Q2[Test100 QC checks]
    Q2 --> Q3{QC passed?}
    Q3 -->|Yes| Q4[Full representative master table]
    Q3 -->|No| Q5[Revise scripts or source paths]
    Q5 --> Q
    Q4 --> Q6[Full master-table QC and interpretation]

    Q6 --> X{Hierarchical classification}
    X -->|1 Swiss-Prot threshold-passing hit| X1[annotated_swissprot_supported]
    X -->|2 informative Cnidaria-TrEMBL threshold-passing hit| X2[sequence_supported_trembl_cnidaria]
    X -->|3 InterProScan support only| X3[domain_supported_interpro]
    X -->|4 eggNOG support only| X4[orthology_supported_eggnog]
    X -->|5 SignalP positive only| X5[function_dark_but_signalp_secretory_candidate]
    X -->|6 No current annotation support| X6[function_dark_no_current_annotation]

    X5 --> Y[Dark-candidate outputs]
    X6 --> Y
    Y --> Y1[Candidate TSV files]
    Y --> Y2[Candidate FASTA files]
    Y --> Y3[Genome-context filters]
    Y --> Y4[Expression-support filters]
    Y4 --> Z[Prioritised high-confidence dark genes]

    C4 --> R[Parse GFF3 lookup table]
    C5 --> R
    R --> R1[Map proteins and transcripts to contigs and coordinates]
    R1 --> R2[Summarise exon count, CDS span and scaffold context]
    R2 --> Y3

    C5 --> R3[RepeatModeler and RepeatMasker]
    R3 --> R4[Flag TE and repeat overlap]
    R4 --> Y3

    C6 --> S[Organise RNA-seq sample sheet]
    C7 --> S
    S --> T[RNA-seq QC and metadata validation]
    T --> U[Align paired-end unstranded reads]
    U --> V[Quantify expression]
    V --> W[Differential expression<br/>control vs stressors]
    W --> Y4

    Z --> ZA[Later comparative analysis<br/>OrthoFinder and lineage restriction]
    ZA --> ZB[Final shortlist]
```

---

## Dark-gene classification strategy

The working classification logic is:

### Annotated
Any representative gene/protein with convincing support from:
- threshold-passing reviewed Swiss-Prot similarity
- conserved domain/family annotation
- orthology-based functional annotation
- corroborated sequence evidence across multiple sources

### Sequence-supported but lower-confidence
A protein with:
- no threshold-passing Swiss-Prot hit
- but a threshold-passing Cnidaria-TrEMBL hit with an informative, non-ambiguous subject description

These should be described conservatively and not treated as equivalent to reviewed Swiss-Prot functional assignments.

### Not fully dark
A protein with:
- no convincing sequence hit
- but additional support such as InterPro domains or eggNOG orthology

### Function-dark
No convincing:
- threshold-passing Swiss-Prot hit
- threshold-passing informative Cnidaria-TrEMBL hit
- conserved domain hit
- orthology-based annotation

### SignalP-positive dark candidate
A function-dark candidate that is SignalP-positive. This remains functionally unresolved, but may be biologically interesting as a possible secreted or signal-peptide-bearing protein.

### High-confidence dark candidate
A function-dark candidate that also has:
- structurally plausible annotation
- acceptable QC status
- supportive genomic context
- ideally expression support from RNA-seq

---

## Role of RNA-seq in this project

RNA-seq is used as a **validation and prioritisation layer**, not as a replacement for annotation.

It will help determine whether candidate dark genes are:

- actually expressed
- reproducibly detected
- structurally supported by splice-aware mapping
- constitutively expressed or condition-specific
- responsive to stress treatments

Dark genes with reproducible expression and/or treatment responsiveness are stronger candidates than unannotated proteins with no expression support.

---

## Suggested project structure

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

## Tracker

A task-by-task record of progress is maintained in:

- [PROJECT_TRACKER.md](PROJECT_TRACKER.md)

The current threshold-revision plan is maintained in:

- [notes/annotation_threshold_revision_plan.md](notes/annotation_threshold_revision_plan.md)

---

## Notes

This repository is intended to track:

- workflow development
- project documentation
- lightweight summary outputs
- scripts
- annotation logic
- interpretation notes

Large raw bioinformatics files are intentionally excluded from version control via `.gitignore`.
