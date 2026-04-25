# project_dark_genes

Exploring and understanding how to identify novel **"dark" genes** in a non-model cnidarian dataset, using *Actinia equina* sequence resources and an HPC workflow.

## Current status

At the moment, the project is starting **from an existing paired CDS/protein gene set**, rather than from a raw genome assembly.

### Current files in hand
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.nucl.fa`
- `equina_smart.rnam-trna.merged.ggf.curated.remredun.aa.fa`

### What has been established so far
- Both files contain **55,607** entries
- Headers match between files (example: `evm.utg4.1`)
- Mean protein length is **452.05 aa**
- Minimum protein length is **32 aa**
- Maximum protein length is **10,458 aa**

Interpretation so far:
- `nucl.fa` is currently being treated as **CDS/transcript-like nucleotide sequences**
- `aa.fa` is currently being treated as the corresponding **predicted proteome**
- The immediate workflow therefore begins with **proteome QC and functional annotation**

A task-by-task record is maintained in [PROJECT_TRACKER.md](PROJECT_TRACKER.md).

---

## Workflow overview

```mermaid
flowchart TD
    A[Project start] --> B[Assess available input files]

    B --> C{What is available?}

    C -->|Current route| D[Existing CDS FASTA + protein FASTA]
    C -->|Future/expanded route| E[Genome assembly FASTA + optional GFF3 + RNA-seq]

    %% Current route
    D --> D1[Confirm FASTA counts match]
    D1 --> D2[Inspect headers and infer file type]
    D2 --> D3[Check CDS/protein pairing<br/>shared IDs and translation logic]
    D3 --> D4[Generate QC tables<br/>CDS lengths and protein lengths]
    D4 --> D5[Flag very short proteins<br/>do not remove yet]
    D5 --> D6[Run BUSCO in protein mode]
    D6 --> D7{Proteome quality acceptable?}

    D7 -->|Yes| D8[Run similarity search<br/>DIAMOND or BLASTp vs UniProt]
    D7 -->|Borderline / poor| D9[Increase caution<br/>treat dark calls as provisional]
    D9 --> D8

    D8 --> D10[Run InterProScan]
    D10 --> D11[Run eggNOG-mapper]
    D11 --> D12[Run SignalP / TM prediction for unresolved proteins]
    D12 --> D13[Build master annotation table]
    D13 --> D14{Classification}

    D14 -->|Sequence hit and/or domain/function| D15[Annotated or rescued proteins]
    D14 -->|No sequence hit only| D16[Sequence-dark candidates]
    D14 -->|No sequence hit + no domain + no function| D17[Function-dark candidates]

    D17 --> D18[Download comparator cnidarian proteomes]
    D18 --> D19[Run OrthoFinder]
    D19 --> D20[Assess lineage restriction]
    D20 --> D21[Control for homology-detection failure]
    D21 --> D22[Prioritise high-confidence dark genes]
    D22 --> D23[Optional structure modelling<br/>ColabFold / AlphaFold-style follow-up]

    %% Future route from raw assembly
    E --> E1[Assembly QC and contamination checks]
    E1 --> E2[Repeat discovery and masking]
    E2 --> E3{Evidence available?}
    E3 -->|RNA-seq + proteins| E4[BRAKER3]
    E3 -->|Proteins only| E5[BRAKER2]
    E3 -->|Genome only| E6[Ab initio fallback]
    E4 --> E7[Predicted genes and proteins]
    E5 --> E7
    E6 --> E7
    E7 --> D6
```

---

## What has been completed so far
- Defined the conceptual dark-gene analysis workflow for cnidarian genomes
- Adapted the workflow for an HPC-based single-species analysis
- Confirmed that chromosome-level assembly is not required to begin the project
- Determined that the currently available files are most likely an **existing EVM-like gene set** rather than a raw genome assembly
- Verified matching sequence counts between nucleotide and protein FASTA files
- Verified matching example headers between nucleotide and protein FASTA files
- Calculated first-pass protein length summary statistics
- Reframed phase 1 toward **gene-set QC**, rather than assembly repeat masking

---

## Immediate next steps
1. Confirm CDS/protein pairing rigorously using shared IDs and translation-length checks
2. Generate QC tables for CDS lengths and protein lengths
3. Flag short proteins for review
4. Run **BUSCO in protein mode**
5. Use BUSCO to decide how strict downstream filtering must be before dark-gene classification

---

## Information needed from BUSCO next
To continue to the annotation stage with confidence, record the following from the BUSCO run:

### Essential
- BUSCO version
- BUSCO mode (`proteins`)
- lineage dataset used
- short summary values:
  - Complete (C)
  - Single-copy (S)
  - Duplicated (D)
  - Fragmented (F)
  - Missing (M)
  - Total BUSCO groups searched (`n`)

### Very useful
- BUSCO output directory path
- `short_summary` text or JSON file
- `full_table.tsv`
- any warnings about duplication, fragmentation, or failed searches

BUSCO is the first major checkpoint because it tells us whether this proteome looks broadly complete, over-duplicated, fragmented, or suspiciously inflated. That directly changes how conservatively we should interpret proteins with no annotation.
