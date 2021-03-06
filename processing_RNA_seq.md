Katie RNA-seq
================
Angela Simeone

hESC (ESC), CNCC and NSC total RNA have been profiled by RNA-seq.

In the following is the number of biological and technical replicates used in this study.

| sample\_num | sample | oigin | cell\_type | origin | rep               | cell\_type | pool (tech. rep) |
|-------------|--------|-------|------------|--------|-------------------|------------|------------------|
| 5           | sample | H9    | hESC       | <NA>   | H9\_hESC\_NA\_p2  | H9\_hESC   | 2                |
| 6           | sample | H9    | hESC       | <NA>   | H9\_hESC\_NA\_p2  | H9\_hESC   | 2                |
| 7           | sample | H9    | hESC       | <NA>   | H9\_hESC\_NA\_p3  | H9\_hESC   | 3                |
| 8           | sample | H9    | hESC       | <NA>   | H9\_hESC\_NA\_p3  | H9\_hESC   | 3                |
| 35          | sample | H9    | hESC       | <NA>   | H9\_hESC\_NA\_p3  | H9\_hESC   | 3                |
| 36          | sample | H9    | hESC       | <NA>   | H9\_hESC\_NA\_p3  | H9\_hESC   | 3                |
| 19          | sample | H9    | hESC       | <NA>   | H9\_hESC\_NA\_p10 | H9\_hESC   | 10               |
| 20          | sample | H9    | hESC       | <NA>   | H9\_hESC\_NA\_p11 | H9\_hESC   | 11               |
| 21          | sample | H9    | hESC       | <NA>   | H9\_hESC\_NA\_p12 | H9\_hESC   | 12               |
| 25          | sample | H9    | CNCC       | A2     | H9\_CNCC\_A2\_p16 | H9\_CNCC   | 16               |
| 26          | sample | H9    | CNCC       | B4     | H9\_CNCC\_B4\_p17 | H9\_CNCC   | 17               |
| 28          | sample | H9    | CNCC       | B3     | H9\_CNCC\_B3\_p19 | H9\_CNCC   | 19               |
| 10          | sample | H9    | NSC        | <NA>   | H9\_NSC\_NA\_p4.2 | H9\_NSC    | 4.2              |
| 11          | sample | H9    | NSC        | <NA>   | H9\_NSC\_NA\_p4.3 | H9\_NSC    | 4.3              |
| 12          | sample | H9    | NSC        | <NA>   | H9\_NSC\_NA\_p4.4 | H9\_NSC    | 4.4              |
| 9           | sample | H9    | NSC        | <NA>   | H9\_NSC\_NA\_p4.1 | H9\_NSC    | 4.1              |

Basic processing of the sequencing files included:

-   fastQC quality control
-   adaptor trimming
-   aligment to transcriptome obtained from hg38

Prior processing, hg38 have been downloaded together with annotations and the reference transcritome has been built from it.

Processing
----------

Download the genome

``` bash
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
#Release 28 (GRCh38.p12)
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz
```

Prepare reference genome hg38

``` bash
gtf='rna_seq/reference_genomes/hg38/gencode.v28.annotation.gtf'
fasta_genome='stem_cells_rna/reference_genomes/hg38/hg38.fa'
out_dir='stem_cells_rna/reference_genomes/hg38'
cmd_rsem_ref="rsem-prepare-reference --bowtie2 --bowtie2-path /Users/simeon01/.linuxbrew/bin/ --gtf $gtf --bowtie $fasta_genome stem_cells_rna/reference_genomes/hg38/rsem_bowtie2_hg38"
sbatch --mem 32G -o $out_dir/out_rsem_ref.out -e $out_dir/out_rsem_ref.err --wrap="${cmd_rsem_ref}"
```

Prepare reference genome hg19

``` bash
gtf='stem_cells_rna/reference_genomes/hg19_gen/gencode.v19.annotation.gtf_withproteinids'
fasta_genome='stem_cells_rna/reference_genomes/hg19_gen/hsa.hg19.fa'
out_dir='stem_cells_rna/reference_genomes/hg19_gen'
cmd_rsem_ref="/Users/simeon01/bin/rsem_1.3.1/RSEM-1.3.1/rsem-prepare-reference --bowtie2 --bowtie2-path /Users/simeon01/.linuxbrew/bin/ --gtf $gtf --bowtie $fasta_genome stem_cells_rna/reference_genomes/hg19_gen/rsem_bowtie2_hg19"
sbatch --mem 32G -o $out_dir/out_rsem_ref.out -e $out_dir/out_rsem_ref.err --wrap="${cmd_rsem_ref}"
```

fastQC check

``` bash
# === fastqc quality check QC ===
cd /stem_cells_rna/
mkdir fastqc
for file in *.fq.gz
do
  sbatch -o %j.out -e %j.err --mem 4G --wrap "fastqc --noextract --nogroup -q -o fastqc/ $file"
done
```

Adapter trimming

``` bash
## === cut illumina adapters ===
cd stem_cells_rna/Data/rnaseq_diff_cells
mkdir trimmed
for f in *r_1.fq.gz; 
do    
sbatch -o %j.$f.tmp.log --mem 8000 -J cutadapt --wrap="cutadapt -q 20 -m 20 -a AGATCGGAAGAGC -o trimmed/${f%%.fq.gz}.trimmed.fq.gz $f" 
done
```

Align to hg38 (transcriptome)

``` bash
## === alignment ===
cd stem_cells_rna/Data/rnaseq_diff_cells/trimmed/
out_dir='stem_cells_rna/Data/rnaseq_diff_cells/trimmed/rsem_out'
for f in *gz
do
  sbatch --mem 60G -o out_${f}.out -J RSEM_map -e out_${f}.err --wrap "rsem-calculate-expression --phred33-quals -p 40 --bowtie2 --append-names --sort-bam-by-coordinate $f stem_cells_rna/reference_genomes/hg38/rsem_bowtie2_hg38 $out_dir/${f%%.fastq.gz*}_rsem"
done
```

The alignment produced 2 type of files: `genes.results` and `isoform.results`. For the differential analysis we used the pseudocounts from the first type of result file. For other analysis we used the TPM quantification store in `genes.results`.
