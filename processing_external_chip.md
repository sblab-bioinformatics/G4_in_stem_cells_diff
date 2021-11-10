Processing external ChIP-seq (Histone marks)
================
Angela Simeone

We used data for H3K4me3 and H3K27me3 for each of the 3 cells. We donwloaded the data and reprocessed them for hg38.

Annotations ESC

    H3K4me3
    ESC_H3K4me3 ==  Sample GSM602296   == 
    SRR067375
    ESC_H3K27me3 ==  Sample GSM602293 == 
    SRR067372

     #### input
    ESC_input == Sample GSM602292 == 
    SRR067371   

Annotation CNCC

    H3K27me3_human1_rep1 ==  Sample GSM1817179  == 
    SRR2096422  
    CNCC_H3K27me3_human1_rep1

    H3K4me3_human1_rep1  == Sample GSM1817174 == 
    SRR2096417
    CNCC_H3K4me3_human1_rep1

    H3K4me3_human2_rep1 ==  Sample GSM1817175 == 
    SRR2096418  
    CNCC_H3K4me3_human2_rep1

    H3K4me3_human3_rep1  == Sample GSM1817176 == 
    SRR2096419  
    CNCC_H3K4me3_human3_rep1

     ### input
    Input_human_rep1 == Sample GSM1817222 == 
    SRR2096456  
    CNCC_Input_human_rep1

Annoations NPC (neuro stem cells)

    H3K27me3 in Neural Progenitor Cells H1 derived renlab.H3K27me3.NPC.02.01 ==  Sample GSM818033 == 
    SRR353319   
    NSC_H3K27me3.NPC.02.01

    H3K27me3 in Neural Progenitor Cells; renlab.H3K27me3.NPC.01.01  ==  Sample GSM818032 == 
    SRR353318
    NSC_H3K27me3.NPC.01.01


    H3K4me3 in Neural Progenitor Cells; renlab.H3K4me3.NPC.01.01 ==  Sample GSM767350 == 
    SRR316168   
    NSC_H3K4me3.NPC.01.01

    H3K4me3 in Neural Progenitor Cells; renlab.H3K4me3.NPC.02.01 ==  Sample GSM767351 == 
    SRR316169
    NSC_H3K4me3.NPC.02.01


     ### inputs
     Reference Epigenome: ChIP-Seq Input from Neural Progenitor Cells; renlab.Input.NPC.01.01 == Sample GSM767355 == 
     SRR316173
     NSC_Input.NPC.01.01
     
     Reference Epigenome: ChIP-Seq Input from Neural Progenitor Cells; renlab.Input.NPC.02.01 == Sample GSM767356 == 
     SRR316174
     NSC_Input.NPC.02.01
     

Download from GEO
-----------------

``` bash

mkdir /stem_cells/extern_stem_cells_data
cd /stem_cells/extern_stem_cells_data

cd /stem_cells/extern_stem_cells_data/rada_iglesia_ESC
# download file: prefetch will download and save SRA file related to SRR accession in 
sbatch --mem 4G --wrap "fastq-dump SRR067375"
sbatch --mem 4G --wrap "fastq-dump SRR067372"
cd /stem_cells/extern_stem_cells_data/rada_iglesia_ESC
SRR_id=(SRR067371)
for curr_id in ${SRR_id[@]}
  do
  sbatch --mem 4G --wrap "fastq-dump --split-files $curr_id && gzip -1 ${curr_id}_1.fastq && gzip -1 ${curr_id}_2.fastq"
done

mkdir /stem_cells/extern_stem_cells_data/prescott_CNCC
cd /stem_cells/extern_stem_cells_data/prescott_CNCC
#H3K27me3_human1_rep1 ==  Sample GSM1817179  == 
sbatch --mem 4G --wrap "fastq-dump SRR2096422"
#H3K4me3_human1_rep1  == Sample GSM1817174 == 
sbatch --mem 4G --wrap "fastq-dump SRR2096417"
#H3K4me3_human2_rep1 ==  Sample GSM1817175 == 
sbatch --mem 4G --wrap "fastq-dump SRR2096418"
#H3K4me3_human3_rep1  == Sample GSM1817176 == 
sbatch --mem 4G --wrap "fastq-dump SRR2096419"
cd /stem_cells/extern_stem_cells_data/prescott_CNCC
SRR_id=(SRR2096456)
for curr_id in ${SRR_id[@]}
  do
  sbatch --mem 4G --wrap "fastq-dump --split-files $curr_id && gzip -1 ${curr_id}_1.fastq && gzip -1 ${curr_id}_2.fastq"
done


mkdir /stem_cells/extern_stem_cells_data/xie_NSC
cd /stem_cells/extern_stem_cells_data/xie_NSC
#H3K27me3 in Neural Progenitor Cells H1 derived renlab.H3K27me3.NPC.02.01 ==  Sample GSM818033 == 
sbatch --mem 4G --wrap "fastq-dump SRR353319"   
#H3K27me3 in Neural Progenitor Cells; renlab.H3K27me3.NPC.01.01  ==  Sample GSM818032 == 
sbatch --mem 4G --wrap "fastq-dump SRR353318"   

#H3K4me3 in Neural Progenitor Cells; renlab.H3K4me3.NPC.01.01 ==  Sample GSM767350 == 
sbatch --mem 4G --wrap "fastq-dump SRR316168"   
#H3K4me3 in Neural Progenitor Cells; renlab.H3K4me3.NPC.02.01 ==  Sample GSM767351 == 
sbatch --mem 4G --wrap "fastq-dump SRR316169"   
# download inputs
cd /stem_cells/extern_stem_cells_data/xie_NSC
SRR_id=(SRR316173 SRR316174)
for curr_id in ${SRR_id[@]}
  do
  sbatch --mem 4G --wrap "fastq-dump --split-files $curr_id && gzip -1 ${curr_id}_1.fastq && gzip -1 ${curr_id}_2.fastq"
done
```

Data processing
---------------

``` bash
cd /stem_cells/extern_stem_cells_data

path_to_process=(prescott_CNCC rada_iglesia_ESC xie_NSC)

for dir in ${path_to_process[@]}
do
  cd /stem_cells/extern_stem_cells_data/$dir
  mkdir trimmed
  pwd
  for fq in *.fastq.gz
  do
    bname=`basename $fq`
    sbatch --mem 16000 --wrap "cutadapt -q 20 -O 3 -a n -o trimmed/${bname%%.fastq.gz}.trimmed.fq.gz $fq"
  done
  echo "======"
done


cd /stem_cells/extern_stem_cells_data

path_to_process=(prescott_CNCC rada_iglesia_ESC xie_NSC)

g='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa'
w='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.sorted.bed'

for dir in ${path_to_process[@]}
do
  cd /stem_cells/extern_stem_cells_data/$dir/trimmed
  mkdir aligned
  pwd
  
  for f in *gz
  do
      sbatch -o %j.$f.tmp.out -e %j.$f.tmp.err --mem 16000 --wrap "bwa mem -M $g $f | samtools view -S -u -F2304 -q 10 -L $w - | samtools sort -@ 8 - > aligned/${f%%.trimmed.fq.gz}.hg38.sort.bam"
  done
  echo " ====== "
done


cd /stem_cells/extern_stem_cells_data
path_to_process=(prescott_CNCC rada_iglesia_ESC)

for dir in ${path_to_process[@]}
do
  cd /stem_cells/extern_stem_cells_data/$dir/trimmed/aligned
  for bam in *hg38.sort.bam
  do
    sbatch -o %j.$bam.tmp.out -e %j.$bam.tmp.err --mem 32000  --wrap "java -Xmx7g -jar /home/simeon01/applications/picard-2.20.3.jar MarkDuplicates INPUT=$bam OUTPUT=${bam%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${bam%%.bam}.markduplicates_metrics.txt"
  done
done

for file in *.sort.bam
  do
  bam_hg38=$file
  #bam_f_nodup=${bam%%.hg38.sort.bam}.hg38.sort.markduplicates.bam
  bam_hg38_nodup=${file%%.bam}.markduplicates.bam
  
  ls -lth $bam_hg38
  echo "====="
  ls -lth $bam_hg38_nodup
  
  sbatch  --mem 4G --wrap "samtools view -c -F 260 $bam_hg38 >  ${file%%.bam}.stat2"
  sbatch  --mem 4G --wrap "samtools view -c -F 260 $bam_hg38_nodup >  ${file%%.bam}.stat5"
done

for bam in *.markduplicates.bam
    do
  sbatch -o %j.out.log --mem 16G --wrap "samtools index $bam"
done

for bam in *.markduplicates.bam
    do
        tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`
        scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }')
        echo $scal_factor_hg38
      echo sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
        sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
        echo "=============="
done


## PEAK CALLING ESC
cd /stem_cells/extern_stem_cells_data/rada_iglesia_ESC/trimmed/aligned
mkdir macs2_peaks
sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR067375.hg38.sort.markduplicates.bam -c SRR067371_1.hg38.sort.markduplicates.bam -n macs2_peaks/H3K4me3_SRR067375_vs_SRR067371"
sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR067372.hg38.sort.markduplicates.bam -c SRR067371_1.hg38.sort.markduplicates.bam -n macs2_peaks/H3K27me3_SRR067372_vs_SRR067371"
sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR067372.hg38.sort.markduplicates.bam -c SRR067371_1.hg38.sort.markduplicates.bam --broad -n macs2_peaks/H3K27me3_SRR067372_vs_SRR067371"

#H3K4me3
#ESC_H3K4me3 ==  Sample GSM602296   == 
#SRR067375
#ESC_H3K27me3 ==  Sample GSM602293 == 
#SRR067372
#### input
#ESC_input == Sample GSM602292 == 
#SRR067371  

## PEAK CALLING CNCC
cd /stem_cells/extern_stem_cells_data/prescott_CNCC/trimmed/aligned
mkdir macs2_peaks

sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR2096422.hg38.sort.markduplicates.bam -c SRR2096456_1.hg38.sort.markduplicates.bam -n macs2_peaks/H3K27me3_SRR2096422_vs_SRR2096456"

sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR2096417.hg38.sort.markduplicates.bam -c SRR2096456_1.hg38.sort.markduplicates.bam -n macs2_peaks/H3K4me3_human1_rep1_SRR2096417_vs_SRR2096456"

sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR2096418.hg38.sort.markduplicates.bam -c SRR2096456_1.hg38.sort.markduplicates.bam -n macs2_peaks/H3K4me3_human2_rep1_SRR2096418_vs_SRR2096456"

sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR2096419.hg38.sort.markduplicates.bam -c SRR2096456_1.hg38.sort.markduplicates.bam -n macs2_peaks/H3K4me3_human3_rep1_SRR2096419_vs_SRR2096456"

## call peaks with broad option for K27 only
sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR2096422.hg38.sort.markduplicates.bam -c SRR2096456_1.hg38.sort.markduplicates.bam --broad -n macs2_peaks/H3K27me3_SRR2096422_vs_SRR2096456"

#H3K27me3_human1_rep1 ==  Sample GSM1817179  == 
#SRR2096422 

#H3K4me3_human1_rep1  == Sample GSM1817174 == 
#SRR2096417

#H3K4me3_human2_rep1 ==  Sample GSM1817175 == 
#SRR2096418 

#H3K4me3_human3_rep1  == Sample GSM1817176 == 
#SRR2096419 

 ### input
#Input_human_rep1 == Sample GSM1817222 == 
#SRR2096456 


### CALL PEAKS FOR NSC

cd /stem_cells/extern_stem_cells_data/xie_NSC/trimmed/aligned
mkdir macs2_peaks

sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR353319.hg38.sort.markduplicates.bam -c SRR316174_1.hg38.sort.markduplicates.bam  -n macs2_peaks/renlab.H3K27me3.NPC.02.01_SRR353319_vs_SRR316174"

sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR353318.hg38.sort.markduplicates.bam -c SRR316173_1.hg38.sort.markduplicates.bam -n macs2_peaks/renlab.H3K27me3.NPC.01.01_SRR353319_vs_SRR316173"

sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR316168.hg38.sort.markduplicates.bam -c SRR316173_1.hg38.sort.markduplicates.bam -n macs2_peaks/renlab.H3K4me3.NPC.01.01_SRR316168_vs_SRR316173"

sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR316169.hg38.sort.markduplicates.bam -c SRR316174_1.hg38.sort.markduplicates.bam -n macs2_peaks/renlab.H3K4me3.NPC.02.01_SRR316169_vs_SRR316174"

## call peaks with broad option for K27 only
sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR353319.hg38.sort.markduplicates.bam -c SRR316174_1.hg38.sort.markduplicates.bam  --broad -n macs2_peaks/renlab.H3K27me3.NPC.02.01_SRR353319_vs_SRR316174"

sbatch --mem 16G --wrap="macs2 callpeak --keep-dup all -t SRR353318.hg38.sort.markduplicates.bam -c SRR316173_1.hg38.sort.markduplicates.bam --broad -n macs2_peaks/renlab.H3K27me3.NPC.01.01_SRR353319_vs_SRR316173"

#H3K27me3 in Neural Progenitor Cells H1 derived renlab.H3K27me3.NPC.02.01 ==  Sample GSM818033 == 
#SRR353319  
#H3K27me3 in Neural Progenitor Cells; renlab.H3K27me3.NPC.01.01  ==  Sample GSM818032 == 
#SRR353318  

#H3K4me3 in Neural Progenitor Cells; renlab.H3K4me3.NPC.01.01 ==  Sample GSM767350 == 
#SRR316168  
#H3K4me3 in Neural Progenitor Cells; renlab.H3K4me3.NPC.02.01 ==  Sample GSM767351 == 
#SRR316169  

### inputs
#Reference Epigenome: ChIP-Seq Input from Neural Progenitor Cells; renlab.Input.NPC.01.01 == Sample GSM767355 == 
#SRR316173
#Reference Epigenome: ChIP-Seq Input from Neural Progenitor Cells; renlab.Input.NPC.02.01 == Sample GSM767356 == 
#SRR316174
```

curate annotations for bivalent promoters
-----------------------------------------

For the annotation of the bivalent promoters we used only the 2 histone marks: H3K4me3 and H3K27me3.

``` bash
#locally /Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw
cd /Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw

multiIntersectBed -i H3K4me3_human*_rep*.narrowPeak | awk '{if ($4>=2) print $0 }' | sortBed -i - |mergeBed -i - > H3K4me3_human1.multi2.bed
ls 
intersectBed -a H3K27me3_SRR2096422_vs_SRR2096456_peaks.broadPeak -b H3K4me3_human1.multi2.bed | sort | uniq | wc -l
intersectBed -a H3K27me3_SRR2096422_vs_SRR2096456_peaks.narrowPeak -b H3K4me3_human1.multi2.bed | sort | uniq | wc -l
intersectBed -a H3K27me3_SRR2096422_vs_SRR2096456_peaks.narrowPeak -b H3K4me3_human1.multi2.bed | sort | uniq | wc -l

intersectBed -a H3K27me3_SRR2096422_vs_SRR2096456_peaks.narrowPeak -b H3K4me3_human1.multi2.bed | sort | uniq > CNCC_bival_regions.H3K27me3_SRR2096422.H3K4me3_human1.multi2.intersect.bed

intersectBed -a H3K27me3_SRR2096422_vs_SRR2096456_peaks.broadPeak -b H3K4me3_human1.multi2.bed | sort | uniq > CNCC_bival_regions.H3K27me3_SRR2096422.broad.H3K4me3_human1.multi2.intersect.bed

#K27 narrowPeak
intersectBed -a renlab.H3K27me3.NPC.01.01_SRR353319_vs_SRR316173_peaks.narrowPeak -b renlab.H3K27me3.NPC.02.01_SRR353319_vs_SRR316174_peaks.narrowPeak | sortBed -i - | mergeBed -i - > renlab.H3K27me3.NPC.intersect.bed

intersectBed -a renlab.H3K4me3.NPC.01.01_SRR316168_vs_SRR316173_peaks.narrowPeak -b renlab.H3K4me3.NPC.02.01_SRR316169_vs_SRR316174_peaks.narrowPeak | sortBed -i - | mergeBed -i - >renlab.H3K4me3.NPC.intersect.bed 
intersectBed -a renlab.H3K27me3.NPC.intersect.bed -b renlab.H3K4me3.NPC.intersect.bed > NSC_bival_regions.renlab.H3K27me3.H3K4me3.intersect.bed

#K27 broadPeak
intersectBed -a renlab.H3K27me3.NPC.01.01_SRR353319_vs_SRR316173_peaks.broadPeak -b renlab.H3K27me3.NPC.02.01_SRR353319_vs_SRR316174_peaks.broadPeak | sortBed -i - | mergeBed -i - > renlab.H3K27me3.broad.NPC.intersect.bed

intersectBed -a renlab.H3K4me3.NPC.01.01_SRR316168_vs_SRR316173_peaks.narrowPeak -b renlab.H3K4me3.NPC.02.01_SRR316169_vs_SRR316174_peaks.narrowPeak | sortBed -i - | mergeBed -i - >renlab.H3K4me3.NPC.intersect.bed 
intersectBed -a renlab.H3K27me3.broad.NPC.intersect.bed -b renlab.H3K4me3.NPC.intersect.bed > NSC_bival_regions.renlab.H3K27me3.broad.H3K4me3.intersect.bed
```

``` bash

output_biv_data=/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/bivalent_promoters

promoter1kb=/scratcha/sblab/simeon01/reference_genomes/hg38/gencode.v28.TSS.slop1000.bed

epig_annotations=/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells

cd $epig_annotations

# for ESC intersect the 2 marks
intersectBed -a Rada_Iglesias_lifted_hg38.H3K27me3.bed -b Rada_Iglesias_lifted_hg38.H3K4me3.bed > $output_biv_data/ESC_bival_regions.Rada_Iglesias_lifted_hg38.H3K27me3.H3K4me3.intersect.bed

# for CNCC intersect the 2 marks
intersectBed -a CNCC_annotations_hg38_no_.H3K27me3.bed -b CNCC_annotations_hg38_no_.H3K4me3.bed > $output_biv_data/CNCC_bival_regions.CNCC_annotations_hg38_no_.H3K27me3.H3K4me3.intersect.bed

# for NSC intersect the 2 marks
intersectBed -a encodeH1_histone_released_H3histMark_hg38_no_.H3K27me3.bed -b encodeH1_histone_released_H3histMark_hg38_no_.H3K4me3.bed > $output_biv_data/NSC_bival_regions.encodeH1_histone_released_H3histMark_hg38.H3K27me3.H3K4me3.intersect.bed


# i have also copied my new bivalent domains (based on peak calling I have performed on selected data)

cd $output_biv_data
cp /stem_cells/extern_stem_cells_data/final_reference_sets_stem_cells_published_data/*bival_regions*broad*bed .

# generate promoters regions (segments) that are bivalent in each pairwise comparison
for file in *bival_regions*bed
do
intersectBed -a $promoter1kb -b $file -wa | sort | uniq > bivalent_promoter.${file}
done

# extract promoters that are bivalent in each pairwise comparison
#ESC_CNCC
intersectBed -a ESC_bival_regions.Rada_Iglesias_lifted_hg38.H3K27me3.H3K4me3.intersect.bed -b CNCC_bival_regions.H3K27me3_SRR2096422.broad.H3K4me3_human1.multi2.intersect.bed | intersectBed -a $promoter1kb -b - -wa | sortBed -i  | uniq  > ESC_bival.CNCC_bival.bed

#ESC_NSC
intersectBed -a ESC_bival_regions.Rada_Iglesias_lifted_hg38.H3K27me3.H3K4me3.intersect.bed -b NSC_bival_regions.renlab.H3K27me3.broad.H3K4me3.intersect.bed | intersectBed -a $promoter1kb -b - -wa | sortBed -i  | uniq  > ESC_bival.NSC_bival.bed

#CNCC_NSC
intersectBed -a CNCC_bival_regions.H3K27me3_SRR2096422.broad.H3K4me3_human1.multi2.intersect.bed -b NSC_bival_regions.renlab.H3K27me3.broad.H3K4me3.intersect.bed | intersectBed -a $promoter1kb -b - -wa | sortBed -i  | uniq  > CNCC_bival.NSC_bival.bed

# promoter bivalent in all 3 cellsESC_epig_annot <- read.delim('ESC_annotated_promoters.Rada_Iglesias_lifted_hg38.bed',sep = "\t", header = T, stringsAsFactors = F)

sortBed -i CNCC_bival_regions.H3K27me3_SRR2096422.broad.H3K4me3_human1.multi2.intersect.bed > CNCC_bival_regions.H3K27me3_SRR2096422.broad.H3K4me3_human1.multi2.intersect.sorted.bed

sortBed -i ESC_bival_regions.Rada_Iglesias_lifted_hg38.H3K27me3.H3K4me3.intersect.bed > ESC_bival_regions.Rada_Iglesias_lifted_hg38.H3K27me3.H3K4me3.intersect.sorted.bed

sortBed -i NSC_bival_regions.renlab.H3K27me3.broad.H3K4me3.intersect.bed > NSC_bival_regions.renlab.H3K27me3.broad.H3K4me3.intersect.sorted.bed

multiIntersectBed -i ESC_bival_regions.Rada_Iglesias_lifted_hg38.H3K27me3.H3K4me3.intersect.sorted.bed CNCC_bival_regions.H3K27me3_SRR2096422.broad.H3K4me3_human1.multi2.intersect.sorted.bed NSC_bival_regions.renlab.H3K27me3.broad.H3K4me3.intersect.sorted.bed | awk '{if ($4==3) print $0}' | sort | uniq > ESC_CNCC_NSC.bival3.bed
```

process external pax6 data
--------------------------

``` bash

g='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa'
w='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.sorted.bed'
cd /scratcha/sblab/simeon01/Data/20200729_external_NSC_pax6_data


# transform qseq into fastq

for file in *.qseq.gz
do
sbatch --mem 12G --time 12:00:00 --wrap "zcat $file | /scratcha/sblab/simeon01/Data/20200729_external_NSC_pax6_data/qseq2fastq.pl > ${file%%.qseq.gz}.fastq && gzip ${file%%.qseq.gz}.fastq && rm ${file%%.qseq.gz}.fastq"
done

sbatch --mem 12G --time 12:00:00 --wrap "zcat HS003-SR-R00035_BD0E2FACXX.s_2_qseq.CHN021.PAX6.qseq.gz| /scratcha/sblab/simeon01/Data/20200729_external_NSC_pax6_data/qseq2fastq.pl > HS003-SR-R00035_BD0E2FACXX.s_2_qseq.CHN021.PAX6.fastq && gzip HS003-SR-R00035_BD0E2FACXX.s_2_qseq.CHN021.PAX6.fastq"



for file in *.qseq.gz
do
echo sbatch --mem 12G --time 12:00:00 --wrap "zcat $file |python qseq_to_fastq.py > ${file%%.qseq.gz}.converted.fastq"
done


mkdir aligned
pwd
g='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa'
w='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.sorted.bed'
cd /scratcha/sblab/simeon01/Data/20200729_external_NSC_pax6_data

for f in *fastq
do
    sbatch --time 12:00:00 --mem 12G --wrap "bwa mem -M $g $f | samtools view -S -u -F2304 -q 10 -L $w - | samtools sort -@ 8 - > aligned/${f%%.fastq}.hg38.sort.bam"
done


# markduplicates picard
for bam in *hg38.sort.bam
  do
    sbatch --time 6:00:00 --mem 12000  --wrap "java -Xmx7g -jar /home/simeon01/applications/picard-2.20.3.jar MarkDuplicates INPUT=$bam OUTPUT=${bam%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${bam%%.bam}.markduplicates_metrics.txt"
done


# call peaks
cd /scratcha/sblab/simeon01/Data/20200729_external_NSC_pax6_data/aligned
mkdir macs2_output
sbatch --time 04:00:00 --mem 16G --wrap="macs2 callpeak --keep-dup all -t H3K27ac_GA003-SR-R00089.s_1_qseq.filtered.CHN013.hg38.sort.markduplicates.bam -c Input_GA003-SR-R00089.s_3_qseq.filtered.CHN015.hg38.sort.markduplicates.bam -n macs2_output/H3K27ac"

sbatch --time 04:00:00 --mem 16G --wrap="macs2 callpeak --keep-dup all -t H3K27me3_GA003-SR-R00089.s_2_qseq.filtered.CHN014.hg38.sort.markduplicates.bam -c Input_GA003-SR-R00089.s_3_qseq.filtered.CHN015.hg38.sort.markduplicates.bam -n macs2_output/H3K27me3"

sbatch --time 04:00:00 --mem 16G --wrap="macs2 callpeak --keep-dup all -t PAX6_HS003-SR-R00035_BD0E2FACXX.s_2_qseq.CHN021.hg38.sort.markduplicates.bam -c Input_GA003-SR-R00089.s_3_qseq.filtered.CHN015.hg38.sort.markduplicates.bam -q 0.1 -n macs2_output/PAX6_FDR0.1 "

# compute stats (tot reads) and use them to generate bedgraph
for file in *.sort.bam
  do
  bam_hg38=$file
  #bam_f_nodup=${bam%%.hg38.sort.bam}.hg38.sort.markduplicates.bam
  bam_hg38_nodup=${file%%.bam}.markduplicates.bam
  
  ls -lth $bam_hg38
  echo "====="
  ls -lth $bam_hg38_nodup
  
  sbatch --time 00:15:00 --mem 4G --wrap "samtools view -c -F 260 $bam_hg38 >  ${file%%.bam}.stat2"
  sbatch --time 00:15:00 --mem 4G --wrap "samtools view -c -F 260 $bam_hg38_nodup >  ${file%%.bam}.stat5"
done

for bam in *.markduplicates.bam
    do
  sbatch --time 00:15:00 -o %j.out.log --mem 16G --wrap "samtools index $bam"
done

for bam in *.markduplicates.bam
    do
        tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`
        scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }')
        echo $scal_factor_hg38
      echo sbatch --time 03:00:00 --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
        sbatch --time 05:00:00 --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
        echo "=============="
done


# prepare annotations for running GAT enrichment analysis
cd /scratcha/sblab/simeon01/Data/20200729_external_NSC_pax6_data/aligned/macs2_output/

touch Pax6_RepreocessedStudy.hg38.NEC.anno 
for file in *narrowPeak
do
label=${file%%_peaks.narrowPeak}
echo $label
awk -v lab=$label '{print $1"\t"$2"\t"$3"\t"lab}' $file >> Pax6_RepreocessedStudy.hg38.NEC.anno 
done

cp Pax6_RepreocessedStudy.hg38.NEC.anno /scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells

cd /scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells
ls

cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2

# run enrichment
./script_for_count_feat_overlaps_pax6nec_reprocessed_pariwise_interaction.sh

#output folder is under /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2
#Fold_enrichment_scenario_additional_PAX6Nec_reprocessed
tar -zcvf Fold_enrichment_scenario_additional_PAX6Nec_reprocessed.tar.gz Fold_enrichment_scenario_additional_PAX6Nec_reprocessed/
```

add H3K9me3 Homo sapiens neuronal stem cell originated from H9
--------------------------------------------------------------

Data downloaded from ENCODE

``` bash
# go to annotations

cd /scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells
ls

#download hg38 H3K9me3 for NSC
wget https://www.encodeproject.org/files/ENCFF452NFM/@@download/ENCFF452NFM.bed.gz
gunzip ENCFF452NFM.bed.gz
#rename
mv ENCFF452NFM.bed NSC_H3K9me3_ENCFF452NFM.bed
# download hg38 H3K9me3 for ESC
wget https://www.encodeproject.org/files/ENCFF112ULZ/@@download/ENCFF112ULZ.bed.gz
gunzip ENCFF112ULZ.bed.gz
#rename
mv ENCFF112ULZ.bed ESC_H3K9me3_ENCFF112ULZ.bed

# create annotation file fpr ESC_NSC_encode
touch ESC_ENCFF112ULZ_NSC_ENCFF452NFM_encode_H3K9me3.anno
awk '{print $1"\t"$2"\t"$3"\tNSC_H3K9me3"}' NSC_H3K9me3_ENCFF452NFM.bed >> ESC_ENCFF112ULZ_NSC_ENCFF452NFM_encode_H3K9me3.anno 
awk '{print $1"\t"$2"\t"$3"\tESC_H3K9me3"}' ESC_H3K9me3_ENCFF112ULZ.bed >> ESC_ENCFF112ULZ_NSC_ENCFF452NFM_encode_H3K9me3.anno
```

Super enhancers CNCC
--------------------

``` bash
# locally
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/super_enhancers
ls *txt
wc -l *txt
#    1048 Wilderman_2018_CNCC_enriched_superenhancers.txt
#    4345 Wilderman_2018_CNCC_superenhancers.txt

# when you type the command, for ^M you type Ctrl-V Ctrl-M
sed -e "s/^M//" Wilderman_2018_CNCC_enriched_superenhancers.txt > Wilderman_2018_CNCC_enriched_superenhancers_mod.txt
sed -e "s/^M//" Wilderman_2018_CNCC_superenhancers.txt > Wilderman_2018_CNCC_superenhancers_mod.txt
cat Wilderman_2018_CNCC_enriched_superenhancers_mod.txt Wilderman_2018_CNCC_superenhancers_mod.txt > CNCC_superenhancers.anno


scp CNCC_superenhancers.anno simeon01@10.20.236.34:/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells
```

RUN gat on superenhancers and H3K9me3 (heterochromatine)
--------------------------------------------------------

``` bash

cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2

./script_for_count_feat_overlaps_H3K9me3_superEnhancers.sh
tar -zcvf Fold_enrichment_scenario_additional_H3K9me3_and_cncc_superenhan.tar.gz Fold_enrichment_scenario_additional_H3K9me3_and_cncc_superenhan/
```
