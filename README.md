README
================

G4 in stem cells and differentiation
====================================

This repository collect code used for the analysis of stem cells by G4-ChIP-seq, by RNA-seq and by ATAC-seq. <br /> <br /> A total of 62 samples have been processed and they are deposited on [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161531).

The goal of the study was to obtain new insights about G4 in a stem cell system. Human embryonic stem cells have been used here and they have been differentiated towards two different states: CNCC cells (Cranial Neural Crest Cells) and NSC (Neuro Stem cells). <br /> <br />

The cells have been profiled by G4-ChIP-seq, RNA-seq and ATAC-seq and we processed and analysed them.

Sequencing processing and data analysis
---------------------------------------

Data analysis:

-   [**Processing of ATAC-seq**](./processing_ATAC_seq.md)

-   [**Processing of RNA-seq**](./processing_RNA_seq.md)

    -   [**Differential gene expression analysis**](./differential_gene_expression_analysis.md)

-   [**Processing of external histone mark ChIP-seq**](./processing_external_chip.md)

-   [**Processing of G4-ChIP-seq**](./processing_G4_chip.md)

    -   [**Differential biding and G4 integration with other data (RNA-seq, ATAC, external, etc..)**](./additional_G4_processing.md)
    -   [**Transcriptional variability**](./transcriptional_variability.md)

List of files used for this study:

-   bg4 peaks (all individual peaks, consensus peaks and cell-specific peaks)
-   atac-seq peaks
-   tpm matrix (all libraries)
-   tpm matrix (median by rep, cells)
-   histone mark peaks
-   promoter coordinates
-   annotation of proromoters
-   G4 signals
