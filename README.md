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

-   [**Image processing protocol**](.//Users/simeon01/Documents/git_repositories/G4_in_stem_cells_diff/microscopy_image_quant/ICY_microscope_image_analysis_file.protocol)


List of files used for this study:

-   bg4 peaks (all individual peaks, consensus peaks and cell-specific peaks)
-   atac-seq peaks
-   tpm matrix (all libraries)
-   tpm matrix (median by rep, cells)
-   histone mark peaks
-   promoter coordinates
-   annotation of proromoters
-   G4 signals

#### Software, tools and environment used

<table style="width:100%;">
<colgroup>
<col width="50%" />
<col width="45%" />
<col width="4%" />
</colgroup>
<thead>
<tr class="header">
<th>Step</th>
<th>Software name and version</th>
<th></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>fastqc</td>
<td>FastQC v0.11.7</td>
<td></td>
</tr>
<tr class="even">
<td>adapter trimmin</td>
<td>cutadapt version 1.16</td>
<td></td>
</tr>
<tr class="odd">
<td>alignment</td>
<td>bwa mem 0.7.17-r1188</td>
<td></td>
</tr>
<tr class="even">
<td>duplicate marking</td>
<td>picard-2.20.3</td>
<td></td>
</tr>
<tr class="odd">
<td>bam file indexing, sorting and handeling</td>
<td>samtools Version: 1.8 (using htslib 1.8)</td>
<td></td>
</tr>
<tr class="even">
<td>bigWig track generation</td>
<td>bamCoverage 3.3.0 (deepTools 3.3.0)</td>
<td></td>
</tr>
<tr class="odd">
<td>peak calling</td>
<td>macs2 2.1.2</td>
<td></td>
</tr>
<tr class="even">
<td>bed files processing and manipulation</td>
<td>bedtools v2.27.1</td>
<td></td>
</tr>
<tr class="odd">
<td></td>
<td></td>
<td></td>
</tr>
<tr class="even">
<td>System info:</td>
<td></td>
<td></td>
</tr>
<tr class="odd">
<td>Linux kernel version</td>
<td>3.10.0-1127.18.2.el7.x86_64</td>
<td></td>
</tr>
<tr class="even">
<td>cluster management and job scheduling system</td>
<td>slurm 20.02.4</td>
<td></td>
</tr>
</tbody>
</table>

R session info:

``` bash
print(sessionInfo())
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ComplexHeatmap_2.0.0 dplyr_1.0.2          edgeR_3.26.8         limma_3.40.6        

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.5          pillar_1.4.6        compiler_3.6.1      RColorBrewer_1.1-2  tools_3.6.1         digest_0.6.25       evaluate_0.14      
 [8] lifecycle_0.2.0     tibble_3.0.3        lattice_0.20-41     clue_0.3-57         pkgconfig_2.0.3     png_0.1-7           rlang_0.4.7        
[15] rstudioapi_0.11     yaml_2.2.1          parallel_3.6.1      xfun_0.17           knitr_1.29          cluster_2.1.0       generics_0.0.2     
[22] GlobalOptions_0.1.2 vctrs_0.3.4         locfit_1.5-9.4      tidyselect_1.1.0    glue_1.4.2          R6_2.4.1            GetoptLong_1.0.2   
[29] rmarkdown_2.3       purrr_0.3.4         magrittr_1.5        ellipsis_0.3.1      htmltools_0.5.0     splines_3.6.1       shape_1.4.5        
[36] circlize_0.4.10     colorspace_1.4-1    crayon_1.3.4        rjson_0.2.20   
```

