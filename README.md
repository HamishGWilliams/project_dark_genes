# project_dark_genes
Exploring and understanding how to identify novel "dark" genes in genomes in non-model species. using my Actinia equina RNA-seq data from my PhD and Dr Craig Wilding's (et al.)  genome scaffolds, nucleotide and amino acid fasta data files.

```mermaid
flowchart TD
    A[Genome FASTA] --> B[Assembly QC<br/>stats + BUSCO genome]
    B --> C[Repeat discovery<br/>RepeatModeler2]
    C --> D[Repeat masking<br/>RepeatMasker]
    D --> E{Evidence available?}

    E -->|RNA-seq + related proteins| F[BRAKER3]
    E -->|Related proteins only| G[BRAKER2]
    E -->|Genome only| H[Ab initio fallback<br/>high false-dark risk]

    F --> I[Predicted genes + proteins]
    G --> I
    H --> I

    I --> J[Proteome QC<br/>BUSCO proteins + OMArk]
    J --> K[Filter weak/repeat/contaminant models]

    K --> L[BLASTp / DIAMOND vs UniProt]
    K --> M[InterProScan]
    K --> N[eggNOG-mapper]
    K --> O[SignalP / TM prediction]

    L --> P[Master annotation table]
    M --> P
    N --> P
    O --> P

    P --> Q{Classification}
    Q -->|annotated| R[Known / rescued genes]
    Q -->|no sequence hit only| S[Sequence-dark]
    Q -->|no sequence + no domain + no function| T[Function-dark candidates]

    T --> U[OrthoFinder with cnidarian proteomes]
    U --> V[GenEra / HDF control]
    V --> W[High-confidence dark genes]
```
    W --> X[Structure modelling<br/>ColabFold / AlphaFold DB]
    X --> Y[Priority shortlist]
