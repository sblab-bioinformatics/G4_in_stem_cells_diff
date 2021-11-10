Processing G4-ChIP-seq
================
Angela Simeone

Three distinct cell types have been profiled here: hESC, CNCC and NSC. Often we refer to hESC as simply ESC. <br /> <br />

For each cell, 3 biological replicates have been performed; each biological replicate have been validated with 3 technical replicates. For each biological replicate an input control library have been prepared and sequenced.

General basic processing of the data include:

-   quality controls
-   adaptor trimming
-   alignment
-   stats generation and genome-wide track (bigWig)
-   peak calling
-   identification of consensus regions per cell type

For peak calling, each pull-down library have been used together with the paired input library. For each biological replicate there is one input library 3 technical replicates.

<table>
<colgroup>
<col width="5%" />
<col width="14%" />
<col width="17%" />
<col width="7%" />
<col width="27%" />
<col width="26%" />
</colgroup>
<thead>
<tr class="header">
<th>sample</th>
<th>cell</th>
<th>cell-extended</th>
<th>derived_from</th>
<th>peaks</th>
<th>run1</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Sample 1</td>
<td>CNCC_Rep1_3xBG4_rep1</td>
<td>cranial neural crest cells (CNCC)</td>
<td></td>
<td>CNCC_Rep1_3xBG4_rep1_hg38_broad_peaks.broadPeak</td>
<td>CNCC_Rep1_3xBG4_ChIP1_SLX-18496.i707_i505.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 2</td>
<td>CNCC_Rep1_3xBG4_rep2</td>
<td>cranial neural crest cells (CNCC)</td>
<td></td>
<td>CNCC_Rep1_3xBG4_rep2_hg38_broad_peaks.broadPeak</td>
<td>CNCC_Rep1_3xBG4_ChIP2_SLX-18496.i708_i505.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 3</td>
<td>CNCC_Rep1_3xBG4_rep4</td>
<td>cranial neural crest cells (CNCC)</td>
<td></td>
<td>CNCC_Rep1_3xBG4_rep4_hg38_broad_peaks.broadPeak</td>
<td>CNCC_Rep1_3xBG4_ChIP4_SLX-18496.i709_i505.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 4</td>
<td>CNCC_Rep2_3xBG4_rep1</td>
<td>cranial neural crest cells (CNCC)</td>
<td></td>
<td>CNCC_Rep2_3xBG4_rep1_hg38_broad_peaks.broadPeak</td>
<td>CNCC_Rep2_3xBG4_ChIP1_SLX-18355.i703_i506.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 5</td>
<td>CNCC_Rep2_3xBG4_rep3</td>
<td>cranial neural crest cells (CNCC)</td>
<td></td>
<td>CNCC_Rep2_3xBG4_rep3_hg38_broad_peaks.broadPeak</td>
<td>CNCC_Rep2_3xBG4_ChIP3_SLX-18355.i704_i506.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 6</td>
<td>CNCC_Rep2_3xBG4_rep4</td>
<td>cranial neural crest cells (CNCC)</td>
<td></td>
<td>CNCC_Rep2_3xBG4_rep4_hg38_broad_peaks.broadPeak</td>
<td>CNCC_Rep2_3xBG4_ChIP4_SLX-18355.i705_i506.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 7</td>
<td>CNCC_Rep3_3xBG4_rep2</td>
<td>cranial neural crest cells (CNCC)</td>
<td></td>
<td>CNCC_Rep3_3xBG4_rep2_hg38_broad_peaks.broadPeak</td>
<td>CNCC_Rep3_3xBG4_ChIP2_SLX-18495.i708_i508.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 8</td>
<td>CNCC_Rep3_3xBG4_rep3</td>
<td>cranial neural crest cells (CNCC)</td>
<td></td>
<td>CNCC_Rep3_3xBG4_rep3_hg38_broad_peaks.broadPeak</td>
<td>CNCC_Rep3_3xBG4_ChIP3_SLX-18495.i709_i508.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 9</td>
<td>CNCC_Rep3_3xBG4_rep4</td>
<td>cranial neural crest cells (CNCC)</td>
<td></td>
<td>CNCC_Rep3_3xBG4_rep4_hg38_broad_peaks.broadPeak</td>
<td>CNCC_Rep3_3xBG4_ChIP4_SLX-18495.i711_i508.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 10</td>
<td>ESC_Rep1_3xBG4_rep1</td>
<td>human embryonic stem cells (ESC)</td>
<td></td>
<td>ESC_Rep1_3xBG4_rep1_hg38_broad_peaks.broadPeak</td>
<td>ESC_Rep1_3xBG4_ChIP1_SLX-18496.i707_i503.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 11</td>
<td>ESC_Rep1_3xBG4_rep2</td>
<td>human embryonic stem cells (ESC)</td>
<td></td>
<td>ESC_Rep1_3xBG4_rep2_hg38_broad_peaks.broadPeak</td>
<td>ESC_Rep1_3xBG4_ChIP2_SLX-18496.i708_i503.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 12</td>
<td>ESC_Rep1_3xBG4_rep4</td>
<td>human embryonic stem cells (ESC)</td>
<td></td>
<td>ESC_Rep1_3xBG4_rep4_hg38_broad_peaks.broadPeak</td>
<td>ESC_Rep1_3xBG4_ChIP4_SLX-18496.i709_i503.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 13</td>
<td>ESC_Rep2_3xBG4_rep2</td>
<td>human embryonic stem cells (ESC)</td>
<td></td>
<td>ESC_Rep2_3xBG4_rep2_hg38_broad_peaks.broadPeak</td>
<td>ESC_Rep2_3xBG4_ChIP2_SLX-18355.i703_i517.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 14</td>
<td>ESC_Rep2_3xBG4_rep3</td>
<td>human embryonic stem cells (ESC)</td>
<td></td>
<td>ESC_Rep2_3xBG4_rep3_hg38_broad_peaks.broadPeak</td>
<td>ESC_Rep2_3xBG4_ChIP3_SLX-18355.i704_i517.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 15</td>
<td>ESC_Rep2_3xBG4_rep4</td>
<td>human embryonic stem cells (ESC)</td>
<td></td>
<td>ESC_Rep2_3xBG4_rep4_hg38_broad_peaks.broadPeak</td>
<td>ESC_Rep2_3xBG4_ChIP4_SLX-18355.i705_i517.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 16</td>
<td>ESC_Rep3_3xBG4_rep1</td>
<td>human embryonic stem cells (ESC)</td>
<td></td>
<td>ESC_Rep3_3xBG4_rep1_hg38_broad_peaks.broadPeak</td>
<td>ESC_Rep3_3xBG4_ChIP1_SLX-18495.i708_i504.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 17</td>
<td>ESC_Rep3_3xBG4_rep2</td>
<td>human embryonic stem cells (ESC)</td>
<td></td>
<td>ESC_Rep3_3xBG4_rep2_hg38_broad_peaks.broadPeak</td>
<td>ESC_Rep3_3xBG4_ChIP2_SLX-18495.i709_i504.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 18</td>
<td>ESC_Rep3_3xBG4_rep4</td>
<td>human embryonic stem cells (ESC)</td>
<td></td>
<td>ESC_Rep3_3xBG4_rep4_hg38_broad_peaks.broadPeak</td>
<td>ESC_Rep3_3xBG4_ChIP4_SLX-18495.i711_i504.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 19</td>
<td>NSC_Rep1_3xBG4_rep2</td>
<td>neural stem cells (NSC)</td>
<td></td>
<td>NSC_Rep1_3xBG4_rep2_hg38_broad_peaks.broadPeak</td>
<td>NSC_Rep1_3xBG4_ChIP2_SLX-18496.i707_i504.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 20</td>
<td>NSC_Rep1_3xBG4_rep3</td>
<td>neural stem cells (NSC)</td>
<td></td>
<td>NSC_Rep1_3xBG4_rep3_hg38_broad_peaks.broadPeak</td>
<td>NSC_Rep1_3xBG4_ChIP3_SLX-18496.i708_i504.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 21</td>
<td>NSC_Rep1_3xBG4_rep4</td>
<td>neural stem cells (NSC)</td>
<td></td>
<td>NSC_Rep1_3xBG4_rep4_hg38_broad_peaks.broadPeak</td>
<td>NSC_Rep1_3xBG4_ChIP4_SLX-18496.i709_i504.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 22</td>
<td>NSC_Rep2_3xBG4_rep1</td>
<td>neural stem cells (NSC)</td>
<td></td>
<td>NSC_Rep2_3xBG4_rep1_hg38_broad_peaks.broadPeak</td>
<td>NSC_Rep2_3xBG4_ChIP1_SLX-18355.i703_i503.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 23</td>
<td>NSC_Rep2_3xBG4_rep2</td>
<td>neural stem cells (NSC)</td>
<td></td>
<td>NSC_Rep2_3xBG4_rep2_hg38_broad_peaks.broadPeak</td>
<td>NSC_Rep2_3xBG4_ChIP2_SLX-18355.i704_i503.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 24</td>
<td>NSC_Rep2_3xBG4_rep4</td>
<td>neural stem cells (NSC)</td>
<td></td>
<td>NSC_Rep2_3xBG4_rep4_hg38_broad_peaks.broadPeak</td>
<td>NSC_Rep2_3xBG4_ChIP4_SLX-18355.i705_i503.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 25</td>
<td>NSC_Rep3_3xBG4_rep2</td>
<td>neural stem cells (NSC)</td>
<td></td>
<td>NSC_Rep3_3xBG4_rep2_hg38_broad_peaks.broadPeak</td>
<td>NSC_Rep3_3xBG4_ChIP2_SLX-18495.i708_i505.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 26</td>
<td>NSC_Rep3_3xBG4_rep3</td>
<td>neural stem cells (NSC)</td>
<td></td>
<td>NSC_Rep3_3xBG4_rep3_hg38_broad_peaks.broadPeak</td>
<td>NSC_Rep3_3xBG4_ChIP3_SLX-18495.i709_i505.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 27</td>
<td>NSC_Rep3_3xBG4_rep4</td>
<td>neural stem cells (NSC)</td>
<td></td>
<td>NSC_Rep3_3xBG4_rep4_hg38_broad_peaks.broadPeak</td>
<td>NSC_Rep3_3xBG4_ChIP4_SLX-18495.i711_i505.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 28</td>
<td>NSC_Rep3b_and_3_3xBG4_rep3</td>
<td>neural stem cells (NSC)</td>
<td></td>
<td>NSC_Rep3b_and_3_3xBG4_rep3_hg38_broad_peaks.broadPeak</td>
<td>NSC_Rep3b_3xBG4_ChIP2_SLX-18355.i708_i505.r_1.fq.g</td>
</tr>
<tr class="odd">
<td>Sample 54</td>
<td>CNCC_Rep1_3xBG4_input</td>
<td>cranial neural crest cells (CNCC)</td>
<td></td>
<td></td>
<td>CNCC_Rep1_3xBG4_input_SLX-18496.i712_i505.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 55</td>
<td>CNCC_Rep2_3xBG4_input</td>
<td>cranial neural crest cells (CNCC)</td>
<td></td>
<td></td>
<td>CNCC_Rep2_3xBG4_input_SLX-18355.i706_i506.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 56</td>
<td>CNCC_Rep3_3xBG4_input</td>
<td>cranial neural crest cells (CNCC)</td>
<td></td>
<td></td>
<td>CNCC_Rep3_3xBG4_input_SLX-18495.i712_i508.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 57</td>
<td>ESC_Rep1_3xBG4_input</td>
<td>human embryonic stem cells (ESC)</td>
<td></td>
<td></td>
<td>ESC_Rep1_3xBG4_input_SLX-18496.i712_i503.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 58</td>
<td>ESC_Rep2_3xBG4_input</td>
<td>human embryonic stem cells (ESC)</td>
<td></td>
<td></td>
<td>ESC_Rep2_3xBG4_input_SLX-18355.i706_i517.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 59</td>
<td>ESC_Rep3_3xBG4_input</td>
<td>human embryonic stem cells (ESC)</td>
<td></td>
<td></td>
<td>ESC_Rep3_3xBG4_input_SLX-18495.i712_i504.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 60</td>
<td>NSC_Rep1_3xBG4_input</td>
<td>neural stem cells (NSC)</td>
<td></td>
<td></td>
<td>NSC_Rep1_3xBG4_input_SLX-18496.i712_i504.r_1.fq.gz</td>
</tr>
<tr class="even">
<td>Sample 61</td>
<td>NSC_Rep2_3xBG4_input</td>
<td>neural stem cells (NSC)</td>
<td></td>
<td></td>
<td>NSC_Rep2_3xBG4_input_SLX-18355.i706_i503.r_1.fq.gz</td>
</tr>
<tr class="odd">
<td>Sample 62</td>
<td>NSC_Rep3_3xBG4_input</td>
<td>neural stem cells (NSC)</td>
<td></td>
<td></td>
<td>NSC_Rep3_3xBG4_input_SLX-18495.i712_i505.r_1.fq.gz</td>
</tr>
</tbody>
</table>

Processing
----------

fastQC check

``` bash
# === fastqc quality check QC ===
cd /stem_cells/
mkdir fastqc
for file in *.fq.gz
do
  sbatch -o %j.out -e %j.err --mem 4G --wrap "fastqc --noextract --nogroup -q -o fastqc/ $file"
done
```

Cut TruSeq Illumina adapters

``` bash
## === cut illumina adapters ===

out_dir=/stem_cells
cd $out_dir
for f in *.fq.gz 
do
  sbatch --mem 8G --wrap "cutadapt -m 10 -q 20 -O 3 -a CTGTCTCTTATACACATCT -o ${f%%.fastq.gz}.trimmed.fastq.gz $f"
done
```

Align to hg38

``` bash
## === alignment ===
cd /stem_cells
ls
path=/stem_cells
g='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa'
w='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.sorted.bed'
mkdir aligned
for f in *trimmed.fastq.gz
do
    sbatch --mem 16000 --wrap "bwa mem -M $g $f | samtools view -S -u -F2304 -q 10 -L $w - | samtools sort -@ 8 - > $path/aligned/${f%%.trimmed.fastq.gz}.hg38.sort.bam"
done
```

Remove duplicates

``` bash
## === remove duplicates ===
cd /stem_cells/aligned
for file in *.hg38.sort.bam
do
  #/Users/marsic01/applications/picard-2.8.2.jar MarkDuplicates ==> this was the one used before
  # now picard in home
  #picard MarkDuplicates --version
  #2.18.12-SNAPSHOT
  sbatch --mem 16G  --wrap "java -Xmx7g -jar ~/picard-2.20.3.jar MarkDuplicates INPUT=$file OUTPUT=${file%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${file%%.bam}.markduplicates_metrics.txt"
  
  echo "==========="
  echo "          "
done
```

Generate stats

``` bash
## === generate stats ===
for file in *.hg38.sort.bam
  do
  
  bam_hg38=$file
  #bam_f_nodup=${bam%%.hg38.sort.bam}.hg38.sort.markduplicates.bam
  bam_hg38_nodup=${file%%.bam}.markduplicates.bam
  
  ls -lth $bam_hg38
  echo "====="
  ls -lth $bam_hg38_nodup
  
  sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_hg38 >  ${file%%.bam}.stat2"
  sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $bam_hg38_nodup >  ${file%%.bam}.stat5"
done
```

Generating tracks

``` bash
# ===  === 
cd 
for bam in *.markduplicates.bam
    do
  sbatch -o %j.out.log --mem 16G --wrap "samtools index $bam"
done

# === generate tracks === 
for bam in *.markduplicates.bam
    do
        tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`
        scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }')
        echo $scal_factor_hg38
      echo sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
        sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
        echo "=============="
done
```

Call peaks

``` bash
# === call peaks === 
cd /stem_cells/aligned
out_dir_narrow='macs2_individual_rep'
out_dir_broad='macs2_broad_individual_rep'
mkdir $out_dir_narrow
mkdir $out_dir_broad

#CNCC case
rep=(1 3 4)
cases=(CNCC_Rep2_3xBG4)
#StemCell3_NT_BG4 # this one has rep2-3-4 and not 1-2-3
for i in ${rep[@]}
do
  #echo $i
  for j in ${cases[@]}
    do
    echo "rep:$i"
    echo $j
    t=`ls ${j}_ChIP${i}*.hg38.sort.markduplicates.bam`
    c=`ls ${j}_input*.hg38.sort.markduplicates.bam`
    ls -lht $t
    ls -lht $c
    echo "=========================="
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $t -c $c -n $out_dir_narrow/${j}_rep${i}" # human is default
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all --broad -t $t -c $c -n $out_dir_broad/${j}_rep${i}_hg38_broad" # human is default
  done
done

#ESC case
rep=(2 3 4)
cases=(ESC_Rep2_3xBG4)
for i in ${rep[@]}
do
  #echo $i
  for j in ${cases[@]}
    do
    echo "rep:$i"
    echo $j
    t=`ls ${j}_ChIP${i}*.hg38.sort.markduplicates.bam`
    c=`ls ${j}_input*.hg38.sort.markduplicates.bam`
    ls -lht $t
    ls -lht $c
    echo "=========================="
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $t -c $c -n $out_dir_narrow/${j}_rep${i}" # human is default
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all --broad -t $t -c $c -n $out_dir_broad/${j}_rep${i}_hg38_broad" # human is default
  done
done

#NSC case
rep=(1 2 4)
cases=(NSC_Rep2_3xBG4)
for i in ${rep[@]}
do
  #echo $i
  for j in ${cases[@]}
    do
    echo "rep:$i"
    echo $j
    t=`ls ${j}_ChIP${i}*.hg38.sort.markduplicates.bam`
    c=`ls ${j}_input*.hg38.sort.markduplicates.bam`
    ls -lht $t
    ls -lht $c
    echo "=========================="
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $t -c $c -n $out_dir_narrow/${j}_rep${i}" # human is default
    sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all --broad -t $t -c $c -n $out_dir_broad/${j}_rep${i}_hg38_broad" # human is default
  done
done
```

In R, generate overlap across technical replicates and generate consensus for each biological replicate. Addditionally, check overlap of individual libraries (peaks) with OQs.

``` r
compare_tech_rep_overlap_and_OQs <- function(conditions,list_FOI,OQs_map,search_thresholds) {
  # list_FOI is the list of files of interest
  summary_OQs <- c()
  counter <- 1
  temp_label <- c()
  for (j in 1:length(search_thresholds)) { 
    
    for (i in 1:length(conditions)){
      
      curr_condition <- conditions[i]
      
      if (search_thresholds[j] != 'default') {
        
        matches_temp1 <- unique(grep(curr_condition,list_FOI, value=TRUE))
        matches_temp2 <- unique(grep(search_thresholds[j],matches_temp1, value = TRUE))
      
        }  else   {
        
        matches_temp1 <- unique(grep(curr_condition,list_FOI, value=TRUE))
        matches_temp2 <- unique(matches_temp1[!grepl(paste(search_thresholds[(1+1):length(search_thresholds)], collapse ="|"),matches_temp1)])
        }
      
      out_file_name_multi2 <- paste0(curr_condition,'_',search_thresholds[j],'_multi2.bed')
      out_file_name_multi3 <- paste0(curr_condition,'_',search_thresholds[j],'_multi3.bed')
        
      system(paste("multiIntersectBed -i", paste(matches_temp2, collapse=' ') ,"| awk '{if($4>=2) print $0}'| sortBed -i - | mergeBed -i - >",out_file_name_multi2))
      system(paste("multiIntersectBed -i", paste(matches_temp2, collapse=' ') ,"| awk '{if($4>=3) print $0}'| sortBed -i - | mergeBed -i - >",out_file_name_multi3))
      multi2_overlap_OQs <- as.numeric(system(paste("intersectBed -a ",out_file_name_multi2, "-b ",OQs_map," -wa | sortBed -i - | mergeBed -i - | wc -l"),intern = TRUE))
      multi3_overlap_OQs <- as.numeric(system(paste("intersectBed -a ",out_file_name_multi3, "-b ",OQs_map," -wa | sortBed -i - | mergeBed -i - | wc -l"),intern = TRUE)) 
      len_multi2 <- as.numeric(system(paste("wc -l",out_file_name_multi2,"| awk '{print $1}'"),intern = TRUE))
      len_multi3 <- as.numeric(system(paste("wc -l",out_file_name_multi3,"| awk '{print $1}'"),intern = TRUE))
      summ <- rbind(len_multi2,len_multi3,multi2_overlap_OQs,multi3_overlap_OQs,multi2_overlap_OQs/len_multi2,multi3_overlap_OQs/len_multi3)
      print(j)
      print(i)
      summary_OQs <- cbind(summary_OQs,summ)
      temp_label <- c(temp_label,paste0(curr_condition,"_",search_thresholds[j]))
      counter <- counter+1
      
      }
  }
  rownames(summary_OQs) <- c('len_multi2','len_multi3','multi2_overlap_OQs','multi3_overlap_OQs','perc_multi2_overlap_OQs','perc_multi3_overlap_OQs')
  colnames(summary_OQs) <- temp_label
  return(summary_OQs)
} 


## ************************************************
# broad peaks
path_peaks <- '/stem_cells/aligned/macs2_broad_individual_rep'
setwd(path_peaks)
# get filenames of narrowPeaks - all -  even various conditions tested.
list_files_broad_peaks <- list.files(path = path_peaks, pattern= '_peaks.broadPeak') 
# for human print peak number and save it into a structure that we can plot/export
N_peaks_hg38 <- vector(mode="numeric", length=length(list_files_broad_peaks))
for (i in 1:length(list_files_broad_peaks)){
  N_peaks_hg38[i] <- as.numeric(system(paste("wc -l",list_files_broad_peaks[i],"|  awk '{print $1}' "),intern = TRUE))
}
names(N_peaks_hg38) <- list_files_broad_peaks

# write output to file
write.table(N_peaks_hg38,'SLX-18496_broadPeak_summary_number_peaks.csv',col.names=NA,row.names = T,quote = F, sep = ",")

# for each condition check multiintersect (2ou3 and 3out3) and then generate dthe intersection with oqs map
conditions_hg38 <- c("ESC_Rep1_3xBG4_rep", "NSC_Rep1_3xBG4_rep","CNCC_Rep1_3xBG4_rep")
pval_hg38 <- c('default')
OQs_hg38 <- '/Users/simeon01/Documents/Katie/OQ_hits.lifted_hg19_to_hg38.bed'
sumary_OQs_hg38 <- compare_tech_rep_overlap_and_OQs(conditions_hg38,list_files_broad_peaks,OQs_hg38,pval_hg38)
write.table(t(sumary_OQs_hg38),'SLX-18496_summary_confirmed_broadpeaks_overlap_OQS_hg38.csv',col.names=NA,row.names = T,quote = F, sep = ",")
```

One library from NSC had a very low sequencing detph therefore it has been resequenced since they had to low depth. In order to incorporate them in the analysis they will be trated as explained in the following:

1.  merge the aligned bam (prior dup-removal)
2.  remove duplicates
3.  call broad peaks
4.  regenerate stats and tracks starting from merged files

The steps above are follower whenev one library have been sequenced more than one time (re-sequenced)

``` bash
bs list runs
#++----------------------------------------------------------------+----------+
#||                         ExperimentName                         |  Status  |
#++----------------------------------------------------------------+----------+
#|           SLX18495_StemCell_Rep3_3xBG4                          | Complete |
#|          SLX-18351_PDTX                                         | Complete |



NSC_bam_seq1=/stem_cells/SLX-18495/aligned/NSC_Rep3_3xBG4_ChIP2.all.hg38.merged.bam 
NSC_bam_seq2=/stem_cells/aligned/SLX-18355/aligned/NSC_Rep3b_3xBG4_ChIP2_SLX-18355.i708_i505.hg38.sort.bam
NSC_bam_merged=/stem_cells/aligned/SLX-18355/aligned/NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.bam
NSC_bam_merged_nodup=/stem_cells/aligned/SLX-18355/aligned/NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.markduplicates.bam

## == mergeBam == 
cd /stem_cells/aligned/SLX-18355/aligned
sbatch --mem 16G --wrap="samtools merge NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.bam $NSC_bam_seq1 $NSC_bam_seq2 && samtools index NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.bam"

## == remove dup == 
sbatch --mem 16G  --wrap "java -Xmx7g -jar ~/picard-2.20.3.jar MarkDuplicates INPUT=NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.bam OUTPUT=NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=NSC_Rep3b_and_3_3xBG4_ChIP2.hg38.sort.markduplicates_metrics.txt"

## == generate stats ==
sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $NSC_bam_merged >  ${NSC_bam_merged%%.bam}.stat2"
sbatch -o %j.out.stat --mem 8G --wrap "samtools view -c -F 260 $NSC_bam_merged_nodup >  ${NSC_bam_merged_nodup%%.bam}.stat5"

## == call peaks only for this rep ==
# note: input has been sequenced before! /stem_cells/aligned/NSC_Rep3_3xBG4_input.all.hg38.merged.markduplicates.bam
cd /stem_cells/aligned/SLX-18355/aligned
out_dir_narrow='macs2_individual_rep'
out_dir_broad='macs2_broad_individual_rep'
NSC_rep3_input=/stem_cells/aligned/NSC_Rep3_3xBG4_input.all.hg38.merged.markduplicates.bam
sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all -t $NSC_bam_merged_nodup -c $NSC_rep3_input -n $out_dir_narrow/NSC_Rep3b_and_3_3xBG4_rep3" # human is default
sbatch -o %j.tmp.log --mem 8000 --wrap "macs2 callpeak --keep-dup all --broad -t $NSC_bam_merged_nodup -c $NSC_rep3_input -n $out_dir_broad/NSC_Rep3b_and_3_3xBG4_rep3_hg38_broad" # human is default
```
