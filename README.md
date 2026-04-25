# project_dark_genes

Exploring and understanding how to identify novel **“dark” genes** in a non-model cnidarian dataset using *Actinia equina* sequence resources, genome annotation files, and RNA-seq data on an HPC workflow.

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
- the available data now support a **genome-aware and expression-aware dark-gene workflow**

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
- identified RNA-seq resources for downstream validation and prioritisation
- defined a **hierarchical UniProt strategy**: Swiss-Prot first, curated TrEMBL second

---

## Immediate next steps

### 1. Archive BUSCO outputs
- save `short_summary`
- inspect/copy `full_table.tsv`

### 2. Start functional annotation
- run **DIAMOND or BLASTp** against **UniProtKB/Swiss-Prot** first
- run a second-pass **DIAMOND or BLASTp** against a **curated TrEMBL set**, preferably drawn from **UniProt Reference Proteomes** and filtered to biologically relevant taxa
- run **InterProScan**
- run **eggNOG-mapper**
- run **SignalP** later for unresolved candidates

### 3. Build genome-aware lookup tables
- parse the GFF3 into a gene → transcript → protein → contig lookup table
- link candidate dark genes back to genomic coordinates
- summarise exon count, CDS span, transcript span, and scaffold location

### 4. Integrate RNA-seq evidence
- convert the experiment notes into a structured sample sheet
- align paired-end unstranded reads
- quantify expression
- identify expressed candidate dark genes
- test for differential expression between controls and stressor treatments

### 5. Build the master annotation table
For each gene/protein, integrate:

- length/QC information
- Swiss-Prot similarity results
- curated TrEMBL similarity results
- InterPro domain results
- eggNOG annotation
- genomic context
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
    C --> C6[RNA-seq data<br/>paired-end, unstranded]
    C --> C7[Experiment notes<br/>treatments + sample size]

    %% QC and validation already completed
    C1 --> D[Confirm paired CDS/protein set]
    C2 --> D
    C3 --> E[Confirm transcript IDs match protein/CDS IDs]
    C4 --> F[Confirm GFF3 mRNA IDs match FASTA IDs]
    C5 --> G[Confirm GFF3 contigs match genome contig names]

    D --> H[Generate QC tables]
    H --> H1[CDS lengths]
    H --> H2[Protein lengths]
    H --> H3[Flag short proteins]

    H --> I[BUSCO proteins mode]
    I --> I1[C:95.6%]
    I --> I2[S:52.9%]
    I --> I3[D:42.7%]
    I --> I4[F:2.6%]
    I --> I5[M:1.8%]

    I --> J{Proteome suitable for annotation?}
    J -->|Yes, but caution due duplication| K[Proceed with conservative annotation workflow]

    %% Hierarchical functional annotation
    K --> L1[DIAMOND / BLASTp vs UniProtKB/Swiss-Prot]
    L1 --> L2{Swiss-Prot hit found?}
    L2 -->|Yes| L3[Assign highest-confidence sequence annotation]
    L2 -->|No| L4[Search curated UniProtKB/TrEMBL<br/>Reference Proteomes / relevant taxa]
    L4 --> L5{Curated TrEMBL hit found?}
    L5 -->|Yes| L6[Record homology support<br/>use conservative functional wording]
    L5 -->|No| L7[Retain as unresolved for deeper screening]

    K --> M[InterProScan]
    K --> N[eggNOG-mapper]
    K --> O[SignalP for unresolved proteins]

    L3 --> P[Master annotation table]
    L6 --> P
    L7 --> P
    M --> P
    N --> P
    O --> P

    %% Genome-aware validation
    C4 --> Q[Parse GFF3 lookup table]
    C5 --> Q
    Q --> Q1[Map proteins/transcripts to contigs and coordinates]
    Q1 --> P

    %% RNA-seq integration
    C6 --> R[Organise RNA-seq sample sheet]
    C7 --> R
    R --> S[RNA-seq QC / metadata validation]
    S --> T[Align paired-end unstranded reads]
    T --> U[Quantify expression]
    U --> V[Differential expression<br/>control vs stressors]
    V --> P

    %% Candidate classification
    P --> W{Classification}
    W -->|Swiss-Prot and/or strong corroborated annotation| X[Annotated / rescued proteins]
    W -->|TrEMBL-only homology support| Y[Sequence-supported but lower-confidence proteins]
    W -->|No sequence hit, but domain/orthology support| Z1[Not fully dark]
    W -->|No sequence hit + no domain + no function| Z2[Function-dark candidates]

    %% Prioritisation
    Z2 --> AA[Filter using genome context]
    AA --> AB[Filter using expression support]
    AB --> AC[Prioritise high-confidence dark genes]
    AC --> AD[Later comparative analysis<br/>OrthoFinder / lineage restriction]
    AD --> AE[Final shortlist]
```

---

## Dark-gene classification strategy

The working classification logic is:

### Annotated
Any gene/protein with convincing support from:
- reviewed Swiss-Prot similarity
- conserved domain/family annotation
- orthology-based functional annotation
- corroborated sequence evidence across multiple sources

### Sequence-supported but lower-confidence
A protein with:
- no Swiss-Prot hit
- but a plausible hit to a **curated TrEMBL** entry

These should be described conservatively and not treated as equivalent to reviewed Swiss-Prot functional assignments.

### Not fully dark
A protein with:
- no convincing sequence hit
- but some additional support such as InterPro domains or eggNOG orthology

### Function-dark
No convincing:
- Swiss-Prot hit
- curated TrEMBL hit
- conserved domain hit
- orthology-based annotation

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
