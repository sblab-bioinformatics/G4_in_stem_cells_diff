Additional processing G4 data
================
Angela Simeone

G4 structures have been profiled in a complex system: ESC cells and 2 cell types derived from ESCs that possess 2 different differentation steps (Neuro crest -CNCC - and neural stem cells).

Katie's has also profiled the transcriptome of those cells.

Epigenetic profiles and other transcription factor chip-seq data used here for comparative analysis have been previously reported by other lab.

Data
----

Consensus voting strategy: \* 2 out 3 technical rep \* 2 out 3 bio rep

!!! NOTE: I would like to test as biological confirmed peaks the UNION OF THE INTERSECTION (as Giovanni's done previously) and see how the final results change. THIS IS NOT BEEN TESTE HERE!!!

To evaluate if/not merge peaks in a certain proximity, explore interpeak distances.

``` bash

cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/confirmed_bio_peaks
for file in *.bio2.bed
  do
  complementBed -i $file -g /Users/simeon01/Documents/genomes/hg38/hg38.sorted.genome > ${file%%.bed}.interpeak_distance.bed
  done
```

Load resulting bed\_files in R and check

``` r
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/confirmed_bio_peaks')

## == actual confirmed peaks ==
ESC_table <- read.table('ESC_multi2.bio2.bed',stringsAsFactors = F)
NSC_table <- read.table('NSC_multi2.bio2.bed',stringsAsFactors = F)
CNCC_table <- read.table('CNCC_multi2.bio2.bed',stringsAsFactors = F)

## == inter-peak space ==
ESC_interpeaks <- read.table('ESC_multi2.bio2.interpeak_distance.bed',stringsAsFactors = F)
NSC_interpeaks <- read.table('NSC_multi2.bio2.interpeak_distance.bed',stringsAsFactors = F)
CNCC_interpeaks <- read.table('CNCC_multi2.bio2.interpeak_distance.bed',stringsAsFactors = F)

## == create a dataframe with everything for analyasis and visualization purposes ==

Cells_features <- data.frame(cell_type= c(rep('ESC',length(ESC_table$V1)),rep('NSC',length(NSC_table$V1)),rep('CNCC',length(CNCC_table$V1))),
                             peak_sizes= c(ESC_table$V3-ESC_table$V2,NSC_table$V3-NSC_table$V2,CNCC_table$V3-CNCC_table$V2))
Cell_interpeak_features <- data.frame(cell_type= c(rep('ESC',length(ESC_interpeaks$V1)),rep('NSC',length(NSC_interpeaks$V1)),rep('CNCC',length(CNCC_interpeaks$V1))),
                                      inter_peak_space = c(ESC_interpeaks$V3-ESC_interpeaks$V2,NSC_interpeaks$V3-NSC_interpeaks$V2,CNCC_interpeaks$V3-CNCC_interpeaks$V2))

ggplot(Cells_features,aes(x=cell_type,y=peak_sizes)) + geom_boxplot()

ggplot(Cell_interpeak_features,aes(x=cell_type,y=inter_peak_space)) + geom_boxplot()
```

Main scenario to test - generate peak files
-------------------------------------------

3 are the scenario I am testing on the ESC and differentiated cells.

1.  peaks confimerd in 2/3 of the technical replicates and then 2/3 of the biological replicate.
2.  peaks in 5/9 of all peaks
3.  2/3 extended of 1kb

``` bash
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad
mkdir scenario1_multi2
mkdir scenario2_6out9
mkdir scenario2_5out9
mkdir scenario3_multi2_extended_1kb


## == SCENARIO 1 ==
cells=(ESC NSC CNCC)
for c in ${cells[@]}
do
echo $c
multi2_files=`ls $c*multi2.bed`
multiIntersectBed -i $multi2_files | awk '{if ($4>=2) print $0}' | sortBed -i - | mergeBed -i - > ./scenario1_multi2/${c}.bio2out3.bed
done


## == SCENARIO 2 ==
for c in ${cells[@]}
do
echo $c
broad_files=`ls $c*peaks.broadPeak`
multiIntersectBed -i $broad_files | awk '{if ($4>=5) print $0}' | sortBed -i - | mergeBed -i - > ./scenario2_5out9/${c}.tech5out9.bed
done

## == SCENARIO 3 ==
cells=(ESC NSC CNCC)
for c in ${cells[@]}
do
echo $c
multi2_files=`ls $c*multi2.bed`
multiIntersectBed -i $multi2_files | awk '{if ($4>=2) print $0}' | sortBed -i - | mergeBed -i - -d 1000 | sortBed -i - | mergeBed -i - > ./scenario3_multi2_extended_1kb/${c}.bio2out3.merged_1000.bed
done


##== summary of file lengths ==
#    9456 ./scenario1_multi2/CNCC.bio2out3.bed
#   17950 ./scenario1_multi2/ESC.bio2out3.bed
#    4436 ./scenario1_multi2/NSC.bio2out3.bed
#    9423 ./scenario2_5out9/CNCC.tech5out9.bed
#    4465 ./scenario2_5out9/NSC.tech5out9.bed
#   17705 ./scenario2_5out9/ESC.tech5out9.bed
#    8581 ./scenario3_multi2_extended_1kb/CNCC.bio2out3_1000.bed
#   15509 ./scenario3_multi2_extended_1kb/ESC.bio2out3_1000.bed
#    4175 ./scenario3_multi2_extended_1kb/NSC.bio2out3_1000.bed
```

annotate Freire\_Pritchett\_2017\_PIR regions
---------------------------------------------

``` r
library(dplyr)
supp_table1 <- read.table('/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/Freire_Pritchett_2017_PIR_supp_table1/SupplementaryTable1.txt', sep = "\t",stringsAsFactors = F, header = T)
# this contains the clusters
supp_table2 <- read.table('/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/Freire_Pritchett_2017_PIR_supp_table1/SupplementaryTable2.txt', sep = "\t",stringsAsFactors = F, header = T)
supp_table3 <- read.table('/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/Freire_Pritchett_2017_PIR_supp_table1/SupplementaryTable3.txt', sep = "\t",stringsAsFactors = F, header = T)

unique(supp_table1$PIRStateESC)
unique(supp_table1$BaitStateESC)

###  gerenate an annotation for various cases (ESC- specific only)

# == generate PIR bait lists annotation list == 
# extract CRU as they define it in the paper (see figure3)
# collection of all PIR plus associated promoter
PIR_bait_anno <- supp_table1 %>% group_by(bait_ID) %>% select(bait_ID,bait_start,bait_end,BaitStateESC,Associated.Gene.Name,bait_chr) %>% distinct()
A <- supp_table1 %>% group_by(bait_ID) %>% filter(ESConly=="TRUE" | Both == "TRUE") 
A_unique_bait_ID <- unique(A$bait_ID)
CRUs <- matrix(nrow = length(A_unique_bait_ID),ncol = 5)
for ( i in 1:10){ #length(A_unique_bait_ID)){
 curr_id  <-  A_unique_bait_ID[i]
 subset_A <- A[which(A$bait_ID == curr_id),]
 CRUs[i,] <-c(curr_id,unique(subset_A$bait_chr),as.numeric(min(subset_A[,c(2,3,6,7)])),as.numeric(max(subset_A[,c(2,3,6,7)])),unique(subset_A$Associated.Gene.Name))
}

CRUs <- data.frame(matrix(nrow = length(A_unique_bait_ID),ncol = 5))
colnames(CRUs) <- c("cru_chr","cru_start","cru_end","bait_ID","Gene")
 for ( i in 1:length(A_unique_bait_ID)){
  curr_id  <-  A_unique_bait_ID[i]
  #print(curr_id)
  subset_A <- A[which(A$bait_ID == curr_id),]
  CRUs$cru_chr[i] <- unique(subset_A$bait_chr)
  CRUs$cru_start[i] <- min(subset_A[,c(2,3,6,7)])
  CRUs$cru_end[i] <- max(subset_A[,c(2,3,6,7)])
  CRUs$bait_ID[i] <- curr_id
  CRUs$Gene[i] <- unique(subset_A$Associated.Gene.Name)
  #print(CRUs[i,])
}


PIR_anno <- supp_table1 %>% group_by(PIR_ID) %>% select(PIR_chr,PIR_start,PIR_end,PIRStateESC,Associated.Gene.Name,ESConly,Both,PIRStateESC)
PIR_anno$combined_ESConly_Both <- gsub("\\s", "",paste(PIR_anno$ESConly,'-',PIR_anno$Both))
colnames(PIR_anno)
#[1] "PIR_ID"                "PIR_chr"               "PIR_start"             "PIR_end"               "PIRStateESC"           "Associated.Gene.Name"  "ESConly"             [8] "Both"                  "combined_ESConly_Both"
# print annotation of PIR to file
write.table(PIR_bait_anno[,c(6,2,3,4)],file='/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/Freire_Pritchett_2017_PIR_supp_table1/PIR_bait.anno.bed',col.names = F, row.names=F, quote =F, sep = "\t")

write.table(PIR_anno[,c(2,3,4,5)],file='/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/Freire_Pritchett_2017_PIR_supp_table1/PIR.anno.bed',col.names = F, row.names=F, quote =F, sep = "\t")

#unique(PIR_anno$combined_ESConly_Both)
#[1] "FALSE-TRUE"  "TRUE-FALSE"  "FALSE-FALSE"
PIR_anno.F_T <- PIR_anno[which(PIR_anno$combined_ESConly_Both=="FALSE-TRUE"),]
PIR_anno.T_F <- PIR_anno[which(PIR_anno$combined_ESConly_Both=="TRUE-FALSE"),]
PIR_anno.F_F <- PIR_anno[which(PIR_anno$combined_ESConly_Both=="FALSE-FALSE"),]
write.table(PIR_anno[,c(2,3,4,5)],file='/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/Freire_Pritchett_2017_PIR_supp_table1/PIR.anno.bed',col.names = F, row.names=F, quote =F, sep = "\t")

write.table(PIR_anno.F_T[,c(2,3,4,5)],file='/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/Freire_Pritchett_2017_PIR_supp_table1/PIR.FALSE_TRUE.anno.bed',col.names = F, row.names=F, quote =F, sep = "\t")

write.table(PIR_anno.T_F[,c(2,3,4,5)],file='/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/Freire_Pritchett_2017_PIR_supp_table1/PIR.TRUE_FALSE.anno.bed',col.names = F, row.names=F, quote =F, sep = "\t")

write.table(PIR_anno.F_F[,c(2,3,4,5)],file='/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/Freire_Pritchett_2017_PIR_supp_table1/PIR.FALSE_FALSE.anno.bed',col.names = F, row.names=F, quote =F, sep = "\t")

##== now export info about clusters
CRUs_clusters <- left_join(CRUs,supp_table3,by="Gene") %>% filter(!is.na(cluster_ESC)) %>% group_by("Gene") %>% select(cru_chr,cru_start,cru_end,cluster_ESC,bait_ID,Gene,MedianRPKM_ESC)

write.table(CRUs_clusters[,c(2,3,4,5)],file='/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/Freire_Pritchett_2017_PIR_supp_table1/PCRUs_clusters.anno.bed',col.names = F, row.names=F, quote =F, sep = "\t")
```

annotation from HiC data have been generated and are moved to cluster
---------------------------------------------------------------------

``` bash
cd /scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells
ls -l 
#-rw-r--r-- 1 simeon01 sblab  173471 Oct  9 16:44 PCRUs_clusters.anno.bed
#-rw-r--r-- 1 simeon01 sblab 2755144 Oct  9 16:44 PIR.anno.bed
#-rw-r--r-- 1 simeon01 sblab  402585 Oct  9 16:44 PIR_bait.anno.bed
#-rw-r--r-- 1 simeon01 sblab  841074 Oct  9 16:44 PIR.FALSE_FALSE.anno.bed
#-rw-r--r-- 1 simeon01 sblab 1085318 Oct  9 16:44 PIR.FALSE_TRUE.anno.bed
#-rw-r--r-- 1 simeon01 sblab  828752 Oct  9 16:44 PIR.TRUE_FALSE.anno.bed

chain_hg19_to_hg38=/scratcha/sblab/simeon01/reference_genomes/liftover_chain/hg19ToHg38.over.chain.gz
for file in *anno.bed
do
  echo $file
  echo ${file%%.bed}.hg38.bed
  /Users/simeon01/applications/liftOver $file $chain_hg19_to_hg38 ${file%%.bed}.hg38.bed  ${file%%.bed}unmapped
  #/Users/simeon01/applications/liftOver $file $chain_hg19_to_hg38 ${file%%.bed}.hg38.bed
done

#remove Chr[d]_ eg. chr7_KI270803v1_alt
for file in *anno.hg38.bed
do 
  echo $file
  grep -v "_" $file > ${file/hg38.bed/hg38_no_.bed}
done
```

re-define genomic annotations files
-----------------------------------

``` bash
tail -n +6  gencode.v28.annotation.sorted.gtf |grep "exon_number 1;"  | awk '{print $1"\t"$4"\t"$5}'| sortBed -i -  > gencode.v28.annotation.exon1.bed
tail -n +6  gencode.v28.annotation.gtf |grep "exon_number 1;"  |awk '{print $1"\t"$4"\t"$5"\t"$10}' | sed 's/\";//g' |sed 's/\"//g' > gencode.v28.annotation.exon1.ensID.bed
sortBed -i gencode.v28.genebody.bed | mergeBed -i -  |bedtools complement -i - -g hg38.sorted.genome > gencode.v28.intergenic.bed
bedtools subtract -a gencode.v28.intergenic.bed -b gencode_v28_exon_merged.bed > gencode_v28_intergenic.exon_subtract.bed
# generate set of exons without first exons
bedtools subtract -a gencode_v28_exon_merged.bed -b gencode.v28.annotation.exon1.merged.bed | sortBed -i - | mergeBed -i -  > gencode_v28_exon_minus_first_exons.bed
bedtools subtract -a gencode.v28.genebody.bed  -b gencode.v28.annotation.exon1.ensID.bed | head

exon=gencode_v28_exon_merged.bed
intron=gencode.v28.intron.bed
firstExon=gencode.v28.annotation.exon1.merged.bed
intergenic=gencode_v28_intergenic.exon_subtract.bed
gene_body=gencode.v28.genebody.bed
NoFirstExon=gencode_v28_exon_minus_first_exons.bed
TSS_single_coord=gencode.v28.TSS.bed
TSS_1kb_slop=gencode.v28.TSS.slop1000.bed


rm annotations_genomic_regions_gencode_v28.anno.bed
touch annotations_genomic_regions_gencode_v28.anno.bed
awk '{print $1"\t"$2"\t"$3"\tAllExons"}' $exon >> annotations_genomic_regions_gencode_v28.anno.bed
awk '{print $1"\t"$2"\t"$3"\tFirstExons"}' $firstExon >> annotations_genomic_regions_gencode_v28.anno.bed
awk '{print $1"\t"$2"\t"$3"\tNOFirstExons"}' $NoFirstExon >> annotations_genomic_regions_gencode_v28.anno.bed
awk '{print $1"\t"$2"\t"$3"\tIntrons"}' $intron >> annotations_genomic_regions_gencode_v28.anno.bed
awk '{print $1"\t"$2"\t"$3"\tGeneBody"}' $gene_body >> annotations_genomic_regions_gencode_v28.anno.bed
awk '{print $1"\t"$2"\t"$3"\tIntergenic"}' $intergenic >> annotations_genomic_regions_gencode_v28.anno.bed
awk '{print $1"\t"$2"\t"$3"\tTSS_single_coord"}' $TSS_single_coord >> annotations_genomic_regions_gencode_v28.anno.bed
awk '{print $1"\t"$2"\t"$3"\tTSS_1kb"}' $TSS_1kb_slop >> annotations_genomic_regions_gencode_v28.anno.bed

# after producin 3prime and 5 prime, include them in the genomic regions annotations
cat annotations_genomic_regions_gencode_v28.anno.bed > annotations_genomic_regions_gencode_v28.anno2.bed
three_prime=gencode.v28.3_UTR.bed
five_prime=gencode.v28.5_UTR.bed
awk '{print $1"\t"$2"\t"$3"\t3_prime"}' $three_prime >>annotations_genomic_regions_gencode_v28.anno2.bed
awk '{print $1"\t"$2"\t"$3"\t5_prime"}' $five_prime >>annotations_genomic_regions_gencode_v28.anno2.bed
```

``` bash
# define annotations to run fold change enrichments 

scenario=(scenario1_multi2 scenario2_5out9 scenario3_multi2_extended_1kb)

workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed

#genomic_annotation=genomic_regions_UCSC_hg38_igenomes.bed6.bed
genomic_annotation=annotations_genomic_regions_gencode_v28.anno.bed
repetitive_annotations=repeatMasker_hg38.elements_no.anno.bed

ESC_anno=Rada_Iglesias_lifted_hg38.bed
CNCC_anno=CNCC_annotations_hg38_no_.bed
NEC_anno=NEC_GSE24447_annotations_hg38_no_.bed
NSC_anno=encodeH1_histone_released_H3histMark_hg38_no_.bed

#annotations derived from HiC paper Freire_Pritchett
FP_CRUs_anno=PCRUs_clusters.anno.hg38_no_.bed
FP_PIR_anno=PIR.anno.hg38_no_.bed
FP_PIR_bait_anno=PIR_bait.anno.hg38_no_.bed
FP_PIR_FALSE_FALSE_anno=PIR.FALSE_FALSE.anno.hg38_no_.bed
FP_PIR_FALSE_TRUE_anno=PIR.FALSE_TRUE.anno.hg38_no_.bed
FP_PIR_TRUE_FALSE_anno=PIR.TRUE_FALSE.anno.hg38_no_.bed

annotation_list1=($genomic_annotation $repetitive_annotations $ESC_anno $NSC_anno $CNCC_anno $FP_CRUs_anno $FP_PIR_anno $FP_PIR_bait_anno $FP_PIR_FALSE_FALSE_anno $FP_PIR_FALSE_TRUE_anno $FP_PIR_TRUE_FALSE_anno)

annotation_list1=($genomic_annotation)

cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
# check overlap in terms of counts and fractio of peaks
for case in ${scenario[@]}
do
  cd $case
  # == check overlap (counts and fraction)
  mkdir intervene_output
  sbatch --mem 12G --wrap "intervene pairwise -i ./*bed --type genomic  --compute count --htype pie --filenames --output ./intervene_output"
  sbatch --mem 12G --wrap "intervene pairwise -i ./*bed --type genomic  --compute frac --htype pie --filenames --output ./intervene_output"
  cd ..
done
  

# for each scenario generate common and specific bed_files from pairwise intersection
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
scenario=(scenario1_multi2 scenario2_5out9 scenario3_multi2_extended_1kb)
for case in ${scenario[@]}
do
  cd ./$case
  esc=`ls ESC*.bed`
  nsc=`ls NSC*.bed`
  cncc=`ls CNCC*.bed`
  mkdir bed_intersections_common
  
  sbatch --mem 6G --wrap "intersectBed -a $esc -b $nsc -v | sort | uniq > ./bed_intersections_common/specific_ESC_vs_NSC.bed"
  sbatch --mem 6G --wrap "intersectBed -b $esc -a $nsc -v | sort | uniq > ./bed_intersections_common/specific_NSC_vs_ESC.bed"
  
  sbatch --mem 6G --wrap "intersectBed -a $esc -b $cncc -v | sort | uniq > ./bed_intersections_common/specific_ESC_vs_CNCC.bed"
  sbatch --mem 6G --wrap "intersectBed -b $esc -a $cncc -v | sort | uniq > ./bed_intersections_common/specific_CNCC_vs_ESC.bed"
  
  sbatch --mem 6G --wrap "intersectBed -a $nsc -b $cncc -v | sort | uniq > ./bed_intersections_common/specific_NSC_vs_CNCC.bed"
  sbatch --mem 6G --wrap "intersectBed -b $nsc -a $cncc -v | sort | uniq > ./bed_intersections_common/specific_CNCC_vs_NSC.bed"
  
  #common sets from both sides
  sbatch --mem 6G --wrap "intersectBed -a $esc -b $nsc -wa | sort | uniq > ./bed_intersections_common/common_ESC_vs_NSC.bed"
  sbatch --mem 6G --wrap "intersectBed -a $esc -b $nsc | sortBed -i - | mergeBed -i -> ./bed_intersections_common/common_ESC_vs_NSC.strict.merged.bed"
  sbatch --mem 6G --wrap "intersectBed -b $esc -a $nsc -wa | sort | uniq > ./bed_intersections_common/common_NSC_vs_ESC.bed"
  
  
  sbatch --mem 6G --wrap "intersectBed -a $esc -b $cncc -wa | sort | uniq > ./bed_intersections_common/common_ESC_vs_CNCC.bed"
  sbatch --mem 6G --wrap "intersectBed -a $esc -b $cncc| sortBed -i - | mergeBed -i -> ./bed_intersections_common/common_ESC_vs_CNCC.strict.merged.bed"
  sbatch --mem 6G --wrap "intersectBed -b $esc -a $cncc -wa | sort | uniq > ./bed_intersections_common/common_CNCC_vs_ESC.bed"
  
  sbatch --mem 6G --wrap "intersectBed -a $nsc -b $cncc -wa | sort | uniq > ./bed_intersections_common/common_NSC_vs_CNCC.bed"
  sbatch --mem 6G --wrap "intersectBed -a $nsc -b $cncc| sortBed -i - | mergeBed -i -> ./bed_intersections_common/common_NSC_vs_CNCC.strict.merged.bed"
  sbatch --mem 6G --wrap "intersectBed -b $nsc -a $cncc -wa | sort | uniq > ./bed_intersections_common/common_CNCC_vs_NSC.bed"
  
  cd ..
done

# create the mergePeaks files for all scenatio. This is used then for differential binding analysis
# for each scenario generate common and specific bed_files from pairwise intersection
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
scenario=(scenario1_multi2 scenario2_5out9 scenario3_multi2_extended_1kb)
for case in ${scenario[@]}
do
  cd ./$case
  esc=`ls ESC*.bed`
  nsc=`ls NSC*.bed`
  cncc=`ls CNCC*.bed`
  mkdir bed_merge_all
  sbatch --mem 6G --wrap "cat $esc $nsc $cncc | sortBed -i - | mergeBed -i - > ./bed_merge_all/${scenario}.all_peaks_merged.bed"
  cd ..
done

Bam_repA=/scratcha/sblab/simeon01/Data/190829_Katie_stemcells/SLX-18496/aligned #.hg38.sort.markduplicates.bam
Bam_repB=/scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4/aligned/merged # markduplicates.bam
Bam_repC=/scratcha/sblab/simeon01/Data/20190917_SLX-18355_Katie_StemCells_Rep2_3xBG4/SLX-18355/aligned #.hg38.sort.markduplicates.bam

#./peaks_coverages
#./promoter_coverages

#loop over scenario in order to run coverages -  for the moment I run it only for scenario1
bam_folders=($Bam_repA $Bam_repB $Bam_repC)
scenario=(scenario1_multi2)
path_analysis=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
for case in ${scenario[@]}
  do
  cd $path_analysis/$case
  bed_all_peaks=`ls $path_analysis/$case/bed_merge_all/*.all_peaks_merged.bed`
  out_dir=$path_analysis/$case/bed_merge_all/coverage_peaks
  mkdir $out_dir
  echo $bed_all_peaks
  for folder_bam in ${bam_folders[@]}
    do
    cd $folder_bam
    pwd
    for bam in `ls *.markduplicates.bam`
      do
      echo $bam
      pwd
      echo sbatch --mem 8G -o %j.log.out -e %j.log.err --wrap "bedtools coverage -a $bed_all_peaks -b $bam  -counts > $out_dir/${bam%%.bam}.union_bg4_peaks"
      sbatch --mem 8G -o %j.log.out -e %j.log.err --wrap "bedtools coverage -a $bed_all_peaks -b $bam  -counts > $out_dir/${bam%%.bam}.union_bg4_peaks"
      echo "   "
      done
  
    done
  
  done
    
  
# perform fold enrichmeent analysis on any bed file in each subfolder
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
path_anno='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells'
scenario=(scenario1_multi2 scenario2_5out9 scenario3_multi2_extended_1kb)
genomic_annotation=annotations_genomic_regions_gencode_v28.anno.bed
annotation_list1=($genomic_annotation)
#scenario=(scenario3_multi2_extended_1kb)
for case in ${scenario[@]}
do
  cd ./$case
  # run gat analysis on full bed files 
  for file in *bed
  do
    echo $file
    for anno in ${annotation_list1[@]}
    do
      echo $anno
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$path_anno/${anno} --segments=${file} -S ${file%%.bed}_${anno%%.bed}.dat -t 4"
      echo ${file%%.bed}_${anno%%.bed}.dat
      echo $cmd_gat
      sbatch --mem 6G --wrap "$cmd_gat"
      echo " "
      
    #sbatch --mem 6G --wrap "$cmd_gat"
    done
  done

  #run the enrichment for subsets (common and specific ones)
  cd ./bed_intersections_common
  pwd
  for file in *bed
  do
    echo $file
    for anno in ${annotation_list1[@]}
    do
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$path_anno/${anno} --segments=${file} -S ${file%%.bed}_${anno%%.bed}.dat -t 4"
      echo ${file%%.bed}.${anno%%.bed}.dat
      pwd
      echo $cmd_gat
      sbatch --mem 6G --wrap "$cmd_gat"
      echo " "
    #sbatch --mem 6G --wrap "$cmd_gat"
    done
  done
  echo " ** == ** =="
  cd ../..
done
  
```

RNA-seq data processed
----------------------

    TPM_hg38_all_bio_rep_meanTechReps.csv

compute BG4 coverage at promoters/ union of BG4 peaks
-----------------------------------------------------

geeneratee graphs on peaks sizes and N. peaks
---------------------------------------------

``` r
rm(list = ls())
library(ggplot2)
library(dplyr)
library(VennDiagram)
# setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2')
# esc <- read.table(file='ESC.bio2out3.bed',stringsAsFactors = F, header = F)
# nsc <- read.table('NSC.bio2out3.bed',stringsAsFactors = F)
# cncc <- read.table('CNCC.bio2out3.bed',stringsAsFactors = F)

setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario2_5out9/')
esc <- read.table(file='ESC.tech5out9.bed',stringsAsFactors = F, header = F)
nsc <- read.table('NSC.tech5out9.bed',stringsAsFactors = F)
cncc <- read.table('CNCC.tech5out9.bed',stringsAsFactors = F)

# setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario3_multi2_extended_1kb/')
# esc <- read.table(file='ESC.bio2out3.merged_1000.bed',stringsAsFactors = F, header = F)
# nsc <- read.table('NSC.bio2out3.merged_1000.bed',stringsAsFactors = F)
# cncc <- read.table('CNCC.bio2out3.merged_1000.bed',stringsAsFactors = F)


stem_cells_df <- data.frame(cell_type=c(rep('1_ESC',length(esc$V1)),rep('2_NSC',length(nsc$V1)),rep('3_CNCC',length(cncc$V1))),
                            peak_sizes=c(esc$V3-esc$V2,nsc$V3-nsc$V2,cncc$V3-cncc$V2))

stem_colors=c('coral1','lightslateblue','mediumseagreen')
pdf("./Peak_sizes.confirmedPeaks.pdf")
stem_cells_df %>% ggplot(aes(y=peak_sizes,x=cell_type,fill=cell_type)) + geom_boxplot( outlier.size =  -1) + ylim(0,600)
dev.off()

pdf("./Peak_numbers.confirmedPeaks.pdf")
stem_cells_df %>% ggplot(aes(cell_type,fill=cell_type)) + geom_bar(colour="black") 
dev.off()

boxplot(esc$V3-esc$V2,nsc$V3-nsc$V2,cncc$V3-cncc$V2, outline = F)

# generate venn diagrams with different overlaps
setwd('./intervene_output')
matrix_counts <- read.table('Intervene_pairwise_count_matrix.txt',stringsAsFactors = F)
set1_cncc <- matrix_counts[1,1]
set2_esc <- matrix_counts[2,2]
set3_nsc <- matrix_counts[3,3]
set1_set2_common <- matrix_counts[1,2]
set1_set3_common <- matrix_counts[1,3]
set2_set3_common <- matrix_counts[2,3]


pdf('../overlap_CNCC_ESC.pdf')
grid.newpage()
draw.pairwise.venn(set1_cncc,set2_esc,set1_set2_common,
                   category=c("CNCC","ESC"),
                   lty = rep("blank",2),
                   fill = c("light blue", "red"), alpha = rep(0.5, 2),
                 cat.pos = c(0,0), cat.dist = rep(0.025, 2))
dev.off()

pdf('../overlap_NSC_ESC.pdf')
grid.newpage()
draw.pairwise.venn(set3_nsc,set2_esc,set2_set3_common,
                   category=c("NSC","ESC"),
                   lty = rep("blank",2),
                   fill = c("light green", "red"), alpha = rep(0.5, 2),
                 cat.pos = c(0,0), cat.dist = rep(0.025, 2))
dev.off()

pdf('../overlap_NSC_CNCC.pdf')
grid.newpage()
draw.pairwise.venn(set3_nsc,set1_cncc,set1_set3_common,
                   category=c("NSC","CNCC"),
                   lty = rep("blank",2),
                   fill = c("light green", "light blue"), alpha = rep(0.5, 2),
                 cat.pos = c(0,0), cat.dist = rep(0.025, 2))
dev.off()
```

Differential binding analysis at peaks for the various stem cells
-----------------------------------------------------------------

``` r
## == fuctions  and packages ==
library(dendextend)
library(dplyr)
library(ggdendro)
library(ggplot2)
tss_file <- '/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS.bed'
source('/Users/simeon01/Documents/Karen/20190613_Karen_continue_pol2_G4_maps/G4_signal_DRB/function_to_export_tables.R')

norm_sign_to_tot_map_reads <- function(matrix,list_stat5_mapped_reads){
  Library_sizes <- c()
  matrix_norm <- matrix
  for (i in 1:dim(matrix)[2]){
    print(i)
    tmp_cond <- colnames(matrix)[i]
    tmp_cond <- gsub('.union_bg4_peaks','',tmp_cond)
    print(tmp_cond)
    tmp_stat <- grep(substr(tmp_cond,1,20),list_stat5_mapped_reads,value = T)
    tot_reads <- as.numeric(system(paste0('cat ',tmp_stat),intern = TRUE))
    Library_sizes[i] <- tot_reads
    matrix_norm[,i] <- matrix[,i]/tot_reads*1e6
  }
  names(Library_sizes) <- colnames(matrix)
  NewList <- list("norm_cpm_M" = matrix_norm, "lib_size" = Library_sizes)
  return(NewList)
}

#################################################
# load data
#################################################
setwd("/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/coverage_peaks")
list_files <- list.files(path = ".",pattern = '.union_bg4_peak')

Matrix_cov <- matrix(ncol = length(list_files), nrow = 20713)
for (i in 1:length(list_files)) {
  temp <- read.table(list_files[i],stringsAsFactors = F)
  Matrix_cov[,i] <- temp$V4
  peaks_coordinates <- temp[,c(1,2,3)]
}
colnames(Matrix_cov) <-list_files
colnames(Matrix_cov) <- gsub('.markduplicates.union_bg4_peaks','',list_files)

rownames(Matrix_cov) <- paste(peaks_coordinates$V1,peaks_coordinates$V2,peaks_coordinates$V3,sep="_")

# specific samples need to be removed
colnames(Matrix_cov)[]

Matrix_cov_selected <- Matrix_cov[,-c(33,37,39,40)]


stat5_files <- list.files(path = ".", pattern = "stat5")
norm_Matrix_cov_selected <- norm_sign_to_tot_map_reads(Matrix_cov_selected,stat5_files)


## == clustering raw data == 

# # clustering - row data
pdf('hclust_raw_counts.pdf')
h_Matrix_cov_selected <- hclust(dist(t(Matrix_cov_selected)),method="ward.D2")
p <- ggdendrogram(h_Matrix_cov_selected,rotate=T) +  ggtitle('clustering - raw counts')
print(p)
dev.off()

## == clustering - data norm to total number of reads == 
pdf('hclust_RPM.pdf')
h_norm_Matrix_cov_selected <- hclust(dist(t(norm_Matrix_cov_selected$norm_cpm_M)),method="ward.D2")
p <- ggdendrogram(h_norm_Matrix_cov_selected,rotate=T) +  ggtitle('clustering - RPM')
print(p)  
dev.off()

## == clustering - data norm to total number of reads -- NO INPUT == 
pdf('hclust_RPM_no_input.pdf')
h_norm_Matrix_cov_selected_noInput <- hclust(dist(t(norm_Matrix_cov_selected$norm_cpm_M[,-grep('input',colnames(norm_Matrix_cov_selected$norm_cpm_M))])),method="ward.D2")
p <- ggdendrogram(h_norm_Matrix_cov_selected_noInput,rotate=T) +  ggtitle('clustering - RPM - no input')
print(p)  
dev.off()

## == clustering - data norm to total number of reads -- NO INPUT == 
pdf('hclust_RPM_no_input_excluding_low10percent.pdf')
A <- norm_Matrix_cov_selected$norm_cpm_M[,-grep('input',colnames(norm_Matrix_cov_selected$norm_cpm_M))]
A_sum <- rowSums(A)
hist(A_sum)
summary(A_sum)
quantile(A_sum, 0.1)
A_selected <- A[A_sum<quantile(A_sum, 0.1),]
h_A_selected <- hclust(dist(t(A_selected)),method="ward.D2")
p <- ggdendrogram(h_A_selected,rotate=T) +  ggtitle('clustering - RPM - no input')
plot(p)
dev.off()

pdf('hclust_RPM_no_input_excluding_NOlow10percent_NOTop10percent.pdf')
A <- norm_Matrix_cov_selected$norm_cpm_M[,-grep('input',colnames(norm_Matrix_cov_selected$norm_cpm_M))]
A_sum <- rowSums(A)
hist(A_sum)
summary(A_sum)
quantile(A_sum, 0.1)
A_selected <- A[A_sum>quantile(A_sum, 0.1) & A_sum<quantile(A_sum, 0.9) ,]
#A_selected <- A[A_sum>quantile(A_sum, 0.5) ,]
h_A_selected <- hclust(dist(t(A_selected)),method="ward.D2")
p <- ggdendrogram(h_A_selected,rotate=T) +  ggtitle('clustering - RPM - no input')
plot(p)
dev.off()


#################################
# Differential analysis
#################################

# in the differential analysis we will account for:
# bio rep [4 different ones]
# technical rep [3 tech rep]
# condition [3 different cell lines]

# create variables
#A_selected <- A[A_sum<quantile(A_sum, 0.1),]
#Counts_to_use_all <- A_selected
Counts_to_use_all <- Matrix_cov_selected[,-grep('input',colnames(norm_Matrix_cov_selected$norm_cpm_M))]
colnames(Counts_to_use_all)
conditions <- c(rep('CNCC',9),rep('ESC',9),rep('NSC',9))
conditions
# "CNCC" "CNCC" "CNCC" "CNCC" "CNCC" "CNCC" "CNCC" "CNCC" "CNCC" "ESC"  "ESC"  "ESC"  "ESC"  "ESC"  "ESC"  "ESC"  "ESC"  "ESC"  "NSC"  "NSC"  "NSC"  "NSC" 
# [23] "NSC"  "NSC"  "NSC"  "NSC"  "NSC" 
#      
conditions <- factor(conditions) # there are 4 factors

conditions <- relevel(conditions,ref="ESC") # relevel as reference base the condition no_TPL_1hr_no_dm6


bio_rep <- rep(c(rep(1,3), rep(2,3), rep(3,3)),3)
bio_rep <- factor(bio_rep)
bio_rep
#[1] 1 1 1 2 2 2 3 3 3 1 1 1 2 2 2 3 3 3 1 1 1 2 2 2 3 3 3
#Levels: 1 2 3


tec_rep <- rep(1:3,9)
tec_rep <- factor(tec_rep)
# lib_size: norm_Matrix_cov_selected$lib_size[-grep('input',colnames(norm_Matrix_cov_selected$norm_cpm_M))]
curr_lib_size <- norm_Matrix_cov_selected$lib_size[-grep('input',colnames(norm_Matrix_cov_selected$norm_cpm_M))]
design_blocking <- model.matrix(~ bio_rep + tec_rep + conditions)
design <- model.matrix(~conditions)
colnames(design_blocking)
# [1]  "(Intercept)"    "bio_rep2"       "bio_rep3"       "tec_rep2"       "tec_rep3"       "conditionsCNCC" "conditionsNSC" 

# create EdgeR object
library(edgeR)
y <- DGEList(counts = Counts_to_use_all, group= conditions)
stopifnot(rownames(y$samples) == names(curr_lib_size))
y$samples$lib.size<- curr_lib_size

y <- calcNormFactors(y, method= 'none')
y <- estimateDisp(y,design_blocking)
y <- estimateGLMCommonDisp(y,design_blocking)
y <- estimateGLMTagwiseDisp(y,design_blocking)


# check MDS plot
pdf('MDS_edger_mdsplot_stem_cells.pdf', 7, 7)
 par(las= 1, cex=1.2)
col_map <- conditions
levels(col_map) <- 1:length(levels(col_map))
col_map <- as.numeric(conditions)
temp_rainbow <- rainbow(3)
col_map_rainbow <- temp_rainbow[col_map]
plotMDS(y,
        #pch=c(0,3,0,3,0,3),
        #xlim=c(-3,3), ylim= c(-0.5,0.5), 
        xlim=c(-1.5,1.5), ylim= c(-1,0.7), 
        #xaxp= c(-40, 40, 8),  yaxp= c(-10,10, 8),
        labels= paste(y$samples$group,bio_rep, sep= '--'),
        col=col_map_rainbow,
        main= 'MDS plot ', cex= 0.7)
grid()
dev.off()

pdf('MDS_edger_mdsplot_stem_cells_symb.pdf', 7, 7)
 par(las= 1, cex=1.2)
col_map <- conditions
levels(col_map) <- 1:length(levels(col_map))
col_map <- as.numeric(conditions)
temp_rainbow <- rainbow(3)
col_map_rainbow <- temp_rainbow[col_map]
plotMDS(y,
        pch=c(rep(0,9),rep(1,9),rep(3,9)),
        #xlim=c(-3,3), ylim= c(-0.5,0.5), 
        xlim=c(-1.5,1.5), ylim= c(-1,0.7), 
        #xaxp= c(-40, 40, 8),  yaxp= c(-10,10, 8),
        #labels= paste(y$samples$group,bio_rep, sep= '--'),
        col=col_map_rainbow,
        main= 'MDS plot ', cex= 0.7)

grid()
dev.off()

pdf('PCA2_edger_mdsplot_stem_cells.pdf', 9, 7)
pcaResult<-prcomp(t(y$counts))
plot(pcaResult$x,
     main= 'Principal components, counts',
     xlab="PC1",
     ylab="PC2",
     #xlab= sprintf('PC1 (sd: %s%%)', round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev)))),
     #ylab= sprintf('PC2 (sd: %s%%)', round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev)))),
     type= 'n')
text(x= pcaResult$x[,1], y= pcaResult$x[,2], labels= paste(y$samples$group,bio_rep, sep= '--'), cex= 0.9, 
     #col= c(rep('#CCFF00FF',5),rep('#00FF66FF',5),rep('#FF0000FF',5)))
     col = col_map_rainbow)
dev.off()


# fit the generalized linear model (quasi-likelihood negative binomial generalized log-linear model)
fit <- glmQLFit(y, design_blocking)
colnames(fit)
# "(Intercept)"    "bio_rep2"       "bio_rep3"       "tec_rep2"       "tec_rep3"       "conditionsCNCC" "conditionsNSC" 

## ==  CNCC vs ESC 
N_condition <- 6

colnames(fit)[N_condition]

label_case <-  paste0(colnames(fit)[N_condition],"_vs_ESC")

qlf.CNCC_vs_ESC <- glmQLFTest(fit, coef=N_condition) 

de_pairing_CNCC_vs_ESC <- topTags(qlf.CNCC_vs_ESC, n= nrow(y$counts))$table

write.table(de_pairing_CNCC_vs_ESC,file = paste0('EdgeRtable_differential_binding_analysis_',label_case,'.csv'), col.names = NA,sep = ",",quote = F)

res_CNCC_vs_ESC <- decideTestsDGE(qlf.CNCC_vs_ESC,adjust.method="BH", p.value = 0.05,lfc=0)

summary(res_CNCC_vs_ESC)


tss_res_CNCC_vs_ESC <- export_bed_and_TSS_in_bed_regions(tss_file,'./',label_case,res_CNCC_vs_ESC)

pdf(paste0('MD_MA_plot_differential_binding_analysis_',label_case,'.pdf'),8,7)
plotMD(qlf.CNCC_vs_ESC,status = res_CNCC_vs_ESC,values =c(1,0,-1),col=c("red","black","blue"))
title(paste0("\n\nDEup:",length(which(res_CNCC_vs_ESC==1)),
             " - DEdown: ",length(which(res_CNCC_vs_ESC==-1)),
             " - notDE: ",length(which(res_CNCC_vs_ESC==0))))
dev.off()

## ==  NSC vs ESC 
N_condition <- 7

colnames(fit)[N_condition]

label_case <-  paste0(colnames(fit)[N_condition],"_vs_ESC")

qlf.NSC_vs_ESC <- glmQLFTest(fit, coef=N_condition) 

de_pairing_NSC_vs_ESC <- topTags(qlf.NSC_vs_ESC, n= nrow(y$counts))$table

write.table(de_pairing_NSC_vs_ESC,file = paste0('EdgeRtable_differential_binding_analysis_',label_case,'.csv'), col.names = NA,sep = ",",quote = F)

res_NSC_vs_ESC <- decideTestsDGE(qlf.NSC_vs_ESC,adjust.method="BH", p.value = 0.1,lfc=0)

summary(res_NSC_vs_ESC)

tss_res_NSC_vs_ESC <- export_bed_and_TSS_in_bed_regions(tss_file,'./',label_case,res_NSC_vs_ESC)

pdf(paste0('MD_MA_plot_differential_binding_analysis_',label_case,'.pdf'),8,7)
plotMD(qlf.NSC_vs_ESC,status = res_NSC_vs_ESC,values =c(1,0,-1),col=c("red","black","blue"))
title(paste0("\n\nDEup:",length(which(res_NSC_vs_ESC==1)),
             " - DEdown: ",length(which(res_NSC_vs_ESC==-1)),
             " - notDE: ",length(which(res_NSC_vs_ESC==0))))
dev.off()


label_case <-  "NSC_vs_CNCC"

qlf.NSC_vs_CNCC <- glmQLFTest(fit, contrast=c(0,0,0,0,0,-1,1))

de_pairing_NSC_vs_CNCC <- topTags(qlf.NSC_vs_CNCC, n= nrow(y$counts))$table

write.table(de_pairing_NSC_vs_CNCC,file = paste0('EdgeRtable_differential_binding_analysis_',label_case,'.csv'), col.names = NA,sep = ",",quote = F)

res_NSC_vs_CNCC <- decideTestsDGE(qlf.NSC_vs_CNCC,adjust.method="BH", p.value = 0.1,lfc=0)

summary(res_NSC_vs_CNCC)

tss_res_NSC_vs_CNCC <- export_bed_and_TSS_in_bed_regions(tss_file,'./',label_case,res_NSC_vs_CNCC)

pdf(paste0('MD_MA_plot_differential_binding_analysis_',label_case,'.pdf'),8,7)
plotMD(qlf.NSC_vs_CNCC,status = res_NSC_vs_CNCC,values =c(1,0,-1),col=c("red","black","blue"))
title(paste0("\n\nDEup:",length(which(res_NSC_vs_CNCC==1)),
             " - DEdown: ",length(which(res_NSC_vs_CNCC==-1)),
             " - notDE: ",length(which(res_NSC_vs_CNCC==0))))
dev.off()

DBA_outcome_Nregions <- cbind(summary(res_NSC_vs_ESC),
                              summary(res_CNCC_vs_ESC),
                             summary(res_NSC_vs_CNCC))

DBA_outcome_N_TSS <- cbind(tss_res_NSC_vs_ESC$placeholder,
                           tss_res_CNCC_vs_ESC$placeholder,
                           tss_res_NSC_vs_CNCC$placeholder)

colnames(DBA_outcome_N_TSS) <- colnames(DBA_outcome_Nregions)
rownames(DBA_outcome_N_TSS) <- tss_res_NSC_vs_ESC$type

write.table(DBA_outcome_N_TSS,file='stem_cells_DBA_outcome_N_TSS.csv',col.names = NA,quote = F, sep = ",")

write.table(DBA_outcome_Nregions,file='stem_cells_DBA_outcome_Nregions.csv',col.names = NA,quote = F, sep = ",")
```

Fold enrichment of regions output of the DBA analysis
-----------------------------------------------------

``` bash
workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed

genomic_annotation=annotations_genomic_regions_gencode_v28.anno.bed
repetitive_annotations=repeatMasker_hg38.elements_no.anno.bed

ESC_anno=Rada_Iglesias_lifted_hg38.bed
CNCC_anno=CNCC_annotations_hg38_no_.bed
NEC_anno=NEC_GSE24447_annotations_hg38_no_.bed
NSC_anno=encodeH1_histone_released_H3histMark_hg38_no_.bed

#annotations derived from HiC paper Freire_Pritchett
FP_CRUs_anno=PCRUs_clusters.anno.hg38_no_.bed
FP_PIR_anno=PIR.anno.hg38_no_.bed
FP_PIR_bait_anno=PIR_bait.anno.hg38_no_.bed
FP_PIR_FALSE_FALSE_anno=PIR.FALSE_FALSE.anno.hg38_no_.bed
FP_PIR_FALSE_TRUE_anno=PIR.FALSE_TRUE.anno.hg38_no_.bed
FP_PIR_TRUE_FALSE_anno=PIR.TRUE_FALSE.anno.hg38_no_.bed

annotation_list1=($genomic_annotation $repetitive_annotations $ESC_anno $NSC_anno $CNCC_anno $FP_CRUs_anno $FP_PIR_anno $FP_PIR_bait_anno $FP_PIR_FALSE_FALSE_anno $FP_PIR_FALSE_TRUE_anno $FP_PIR_TRUE_FALSE_anno)


cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/DBA_output_regions
path_anno='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells'
for file in *bed
do
  echo $file
  for anno in ${annotation_list1[@]}
  do
    echo $anno
    cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$path_anno/${anno} --segments=${file} -S ${file%%.bed}_${anno%%.bed}.dat -t 4"
    echo ${file%%.bed}_${anno%%.bed}.dat
    echo $cmd_gat
    sbatch --mem 6G --wrap "$cmd_gat"
    echo " "
  done
done


  
```

Generate list of genes that include specifically BG4-regions in promoter or in genebody minus exon1
---------------------------------------------------------------------------------------------------

``` bash

# for *differentiallyUP.bed, differentiallyDOWN.bed and differentiallyNOTDE.bed
# intersect with TSS gencode.v28.TSS.slop1000.bed
# intersect with TSS gencode.v28.TSS.slop1000.bed
# intersect with gencode.v28.genebody_minus_exon1.bed
tss_1000=/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS.slop1000.bed
tss_genebody_minus_exon1=/Users/simeon01/Documents/genomes/hg38/gencode.v28.genebody_minus_exon1.bed

cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/coverage_peaks
all_bed_regions_up=`ls *differentiallyUP.bed`
all_bed_regions_down=`ls *differentiallyDOWN.bed`
all_bed_regions_notde=`ls *differentiallyNOTDE.bed`
for file in ${all_bed_regions_up[@]}
do
  intersectBed -a $tss_1000 -b $file -F 1 -wa -wb | sort | uniq > ${file%%.bed}.BG4region_in_promoter.bed
  intersectBed -a $tss_genebody_minus_exon1 -b $file -F 1 -wa -wb | sort | uniq > ${file%%.bed}.BG4region_in_genebody_minus_exon1.bed
done

for file in ${all_bed_regions_down[@]}
do
  intersectBed -a $tss_1000 -b $file -F 1 -wa -wb | sort | uniq > ${file%%.bed}.BG4region_in_promoter.bed
  intersectBed -a $tss_genebody_minus_exon1 -b $file -F 1 -wa -wb | sort | uniq > ${file%%.bed}.BG4region_in_genebody_minus_exon1.bed
done

for file in ${all_bed_regions_notde[@]}
do
  intersectBed -a $tss_1000 -b $file -F 1 -wa -wb | sort | uniq > ${file%%.bed}.BG4region_in_promoter.bed
  intersectBed -a $tss_genebody_minus_exon1 -b $file -F 1 -wa -wb | sort | uniq > ${file%%.bed}.BG4region_in_genebody_minus_exon1.bed
done
```

Additional comparison of static peaks with promoter and bivalent promoters
--------------------------------------------------------------------------

The analysis is done in 2 steps: 1. annotate promoter for bivalent marks in the 3 different cell types - keep results as a bed file where first 3 fields are the coordinates of the promoters 2. overlap promoter and BG4 peaks and keep both informations into the same bed files 3. import everything in R and plot sizes

``` bash
# attempt1
####################  .  ####################  .  ####################  .  ####################  .  ####################

promoter_1kb=/scratcha/sblab/simeon01/reference_genomes/hg38/gencode.v28.TSS.slop1000.bed
Bival_Promoters_ESC=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/Additional_annotations/Rada_Iglesias_lifted_hg38_bivalentPromoters_curated.bed
Bival_Promoters_CNCC=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/Additional_annotations/CNCC_hg38_bivalentPromoters.bed
Bival_Promoters_NSC=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/Additional_annotations/NSC_hg38_bivalentPromoters.bed

cd /scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells
grep H3K4me3 Rada_Iglesias_lifted_hg38.bed | sortBed -i - | mergeBed -i - > Rada_Iglesias_lifted_hg38.H3K4me3.bed
grep H3K27me3 Rada_Iglesias_lifted_hg38.bed | sortBed -i - | mergeBed -i -  > Rada_Iglesias_lifted_hg38.H3K27me3.bed
grep K4me3 CNCC_annotations_hg38_no_.bed | sortBed -i - | mergeBed -i -  > CNCC_annotations_hg38_no_.H3K4me3.bed
grep K27me3 CNCC_annotations_hg38_no_.bed | sortBed -i - | mergeBed -i - > CNCC_annotations_hg38_no_.H3K27me3.bed
grep H3K4me3 encodeH1_histone_released_H3histMark_hg38_no_.bed | sortBed -i - | mergeBed -i - > encodeH1_histone_released_H3histMark_hg38_no_.H3K4me3.bed
grep H3K27me3 encodeH1_histone_released_H3histMark_hg38_no_.bed | sortBed -i - | mergeBed -i - > encodeH1_histone_released_H3histMark_hg38_no_.H3K27me3.bed

# annotate promoters for K4, K27 and bivalent

sbatch --mem 12G --wrap " bedtools annotate -i $promoter_1kb -files $Bival_Promoters_ESC Rada_Iglesias_lifted_hg38.H3K4me3.bed Rada_Iglesias_lifted_hg38.H3K27me3.bed $Bival_Promoters_CNCC CNCC_annotations_hg38_no_.H3K4me3.bed CNCC_annotations_hg38_no_.H3K27me3.bed $Bival_Promoters_NSC encodeH1_histone_released_H3histMark_hg38_no_.H3K4me3.bed encodeH1_histone_released_H3histMark_hg38_no_.H3K27me3.bed -names ESC_bival ESC_H3K4me3 ESC_H3K27me3 CNCC_bival CNCC_H3K4me3 CNCC_H3K27me3 NSC_bival NSC_H3K4me3 NSC_H3K27me3  -counts > Promoters_annotated_For_epigenetic_marks.bed"

####################  .  ####################  .  ####################  .  ####################  .  ####################
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2

# common between 3 sets
multiIntersectBed -i *bio2out3.bed | awk '{if ($4==3) print $0}' | sortBed -i - | mergeBed -i - | sortBed -i -  > ESC_CNCC_NSC.multi3_common.bed

cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2
rm temp1 temp2
intersectBed -a NSC.bio2out3.bed -b CNCC.bio2out3.bed -wa > temp1
intersectBed -b NSC.bio2out3.bed -a CNCC.bio2out3.bed -wa > temp2
cat temp1 temp2 |sortBed -i - | mergeBed -i - | sortBed -i - > NSC_CNCC_intermediate1.bed # union of intersectio of NSC and CNCC

rm temp1 temp2
intersectBed -a NSC.bio2out3.bed -b ESC.bio2out3.bed -wa > temp1
intersectBed -b NSC.bio2out3.bed -a ESC.bio2out3.bed -wa > temp2
cat temp1 temp2 |sortBed -i - | mergeBed -i - | sortBed -i - > NSC_ESC_intermediate2.bed # union of intersectio of NSC and ESC

rm temp1 temp2
intersectBed -a ESC.bio2out3.bed -b CNCC.bio2out3.bed -wa > temp1
intersectBed -b ESC.bio2out3.bed -a CNCC.bio2out3.bed -wa > temp2
cat temp1 temp2 |sortBed -i -  | mergeBed -i - | sortBed -i - > ESC_CNCC_intermediate3.bed # union of intersectio of CNCC and ESC


 

# specific to each cell: cell line 1 minus merge of the other 2 (so one peak has to be exclusive to one case only)
cat NSC.bio2out3.bed CNCC.bio2out3.bed | sortBed -i - |mergeBed -i - | intersectBed -a ESC.bio2out3.bed -b - -v > ESC_specific.20200421.bed
cat ESC.bio2out3.bed CNCC.bio2out3.bed | sortBed -i - |mergeBed -i - | intersectBed -a NSC.bio2out3.bed -b - -v > NSC_specific.20200421.bed
cat NSC.bio2out3.bed ESC.bio2out3.bed | sortBed -i - |mergeBed -i - | intersectBed -a CNCC.bio2out3.bed -b - -v > CNCC_specific.20200421.bed

# attempt2 - APRIL MAY 2020 -- CNCC and NPC (like NSC) have been reprocessed from scratch!!!
####################  .  ####################  .  ####################  .  ####################  .  ####################

path_annotations=/scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/final_reference_sets_stem_cells_published_data
cd $path_annotations
ls
cd $path_annotations
promoter_1kb=/scratcha/sblab/simeon01/reference_genomes/hg38/gencode.v28.TSS.slop1000.bed
#CNCC_bival_regions.H3K27me3_SRR2096422.H3K4me3_human1.multi2.intersect.bed
#CNCC_H3K27me3_SRR2096422_vs_SRR2096456_peaks.bed
#CNCC_H3K4me3_human1.multi2.bed
#ESC_bival_regions.Rada_Iglesias_lifted_hg38.H3K27me3.H3K4me3.intersect.bed
#NSC_bival_regions.renlab.H3K27me3.H3K4me3.intersect.bed
#NSC_H3K27me3_NPC_renlab.bed
#NSC_H3K4me3_NPC_renlab.bed
#Rada_Iglesias_lifted_hg38.H3K27me3.bed
#Rada_Iglesias_lifted_hg38.H3K4me3.bed
# after updating annotations on bivalent promoters (April 2020), re-run this inluding
#-ESC_k27
#-ESC_K27_minus_bival
#-ESC_k4
#-ESC_K4_minus_bival
#-ESC_bival
ESC_bival=$path_annotations/ESC_bival_regions.Rada_Iglesias_lifted_hg38.H3K27me3.H3K4me3.intersect.bed
ESC_K27=$path_annotations/Rada_Iglesias_lifted_hg38.H3K27me3.bed
ESC_K4=$path_annotations/Rada_Iglesias_lifted_hg38.H3K4me3.bed
ESC_G4=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC.bio2out3.bed
ESC_G4_specific=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC_specific.20200421.bed
ESC_NSC_G4_union_of_intersection=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC_specific.20200421.bed

ESC_annotated_promoters=$path_annotations/ESC_annotated_promoters.Rada_Iglesias_lifted_hg38.bed

#bedtools annotate -i $promoter_1kb -files $ESC_K4 $ESC_K27 $ESC_bival -names ESC_H3K4me3 ESC_H3K27me3 ESC_bival -counts > $ESC_annotated_promoters
#-CNCC_k27
#-CNCC_k4
#-CNCC_K27_minus_bival
#-CNCC_K4_minus_bival
#-CNCC_bival
CNCC_bival=$path_annotations/CNCC_bival_regions.H3K27me3_SRR2096422.broad.H3K4me3_human1.multi2.intersect.bed #CNCC_bival_regions.H3K27me3_SRR2096422.H3K4me3_human1.multi2.intersect.bed # before it was base on narrow Peaks
CNCC_K27=$path_annotations/CNCC_H3K27me3_SRR2096422_vs_SRR2096456_peaks.broad.bed #CNCC_H3K27me3_SRR2096422_vs_SRR2096456_peaks.bed
CNCC_K4=$path_annotations/CNCC_H3K4me3_human1.multi2.bed
CNCC_G4=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/CNCC.bio2out3.bed
CNCC_G4_specific=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/CNCC_specific.20200421.bed
CNCC_annotated_promoters=$path_annotations/CNCC_annotated_promoters.CNCC_reprocessed.bed
#bedtools annotate -i $promoter_1kb -files $CNCC_K4 $CNCC_K27 $CNCC_bival -names CNCC_H3K4me3 CNCC_H3K27me3 CNCC_bival -counts > $CNCC_annotated_promoters

#-NSC_k27
#-NSC_k4
#-NSC_K27_minus_bival
#-NSC_K4_minus_bival
#-NSC_bival
NSC_bival=$path_annotations/NSC_bival_regions.renlab.H3K27me3.broad.H3K4me3.intersect.bed #NSC_bival_regions.renlab.H3K27me3.H3K4me3.intersect.bed
NSC_K27=$path_annotations/NSC_H3K27me3_NPC_renlab.broad.bed #NSC_H3K27me3_NPC_renlab.bed
NSC_K4=$path_annotations/NSC_H3K4me3_NPC_renlab.bed
NSC_G4=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/NSC.bio2out3.bed
NSC_G4_specific=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/NSC_specific.20200421.bed
NSC_annotated_promoters=$path_annotations/NSC_annotated_promoters.NPC_renlab.bed
NSC_annotated_promoters_080720=$path_annotations/NSC_annotated_promoters.NPC_renlab.080720.bed

bedtools annotate -i $promoter_1kb -files $NSC_K4 $NSC_K27 $NSC_bival -names NSC_H3K4me3 NSC_H3K27me3 NSC_bival -counts > $NSC_annotated_promoters
bedtools annotate -i $promoter_1kb -files $NSC_K4 $NSC_K27 $NSC_bival -names NSC_H3K4me3 NSC_H3K27me3 NSC_bival -counts > $NSC_annotated_promoters_080720

ESC_CNCC_NSC_annotated_promoters_G4=$path_annotations/ESC_CNCC_NSC_annotated_promoters_G4.bed
### annotations of promoters with G4 signal
ESC_CNCC_NSC_G4_common_all=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC_CNCC_NSC.multi3_common.bed

####################  .  ####################  .  ####################  .  ####################  .  ####################

# use also the intermediate cases (intermediates)
# ESC_NSC=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/NSC_ESC_intermediate2.bed
# ESC_CNCC=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC_CNCC_intermediate3.bed
# NSC_CNCC=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/NSC_CNCC_intermediate1.bed

#common_ESC_vs_CNCC.strict.merged.bed  common_ESC_vs_NSC.strict.merged.bed  common_NSC_vs_CNCC.strict.merged.bed
ESC_NSC=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/bed_intersections_common/common_ESC_vs_NSC.strict.merged.bed
ESC_CNCC=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/bed_intersections_common/common_ESC_vs_CNCC.strict.merged.bed
NSC_CNCC=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/bed_intersections_common/common_NSC_vs_CNCC.strict.merged.bed



bedtools annotate -i $promoter_1kb -files $ESC_G4 $CNCC_G4 $NSC_G4 $ESC_G4_specific $CNCC_G4_specific $NSC_G4_specific $ESC_CNCC_NSC_G4_common_all $ESC_NSC $ESC_CNCC $NSC_CNCC -names ESC_G4 CNCC_G4 NSC_G4 ESC_G4_specific CNCC_G4_specific NSC_G4_specific ESC_CNCC_NSC_G4_common common_ESC_vs_NSC.strict common_ESC_vs_CNCC.strict common_NSC_vs_CNCC.strict -counts > $ESC_CNCC_NSC_annotated_promoters_G4


## functional annotation for genes that are in proximity of the BG4 peaks by cell line
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Comparison_promoters_BG4_static_peaks
output_folder='/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Comparison_promoters_BG4_static_peaks'
bed_scenario1=`ls /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/*bed`
sbatch --mem 12G --wrap "bedtools annotate -i $promoter_1kb -files /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/*bed -names CNCC_BG4 NSC_BG4 ESC_BG4 -counts > $output_folder/Promoters_annotated_For_BG4peaks.bed"

# save to files those peaks that overlap with promoters
intersectBed -a /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC.bio2out3.bed -b $promoter_1kb -wa -wb > $output_folder/ESC_overlap_promoter.bed
intersectBed -a /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/NSC.bio2out3.bed -b $promoter_1kb -wa -wb > $output_folder/NSC_overlap_promoter.bed
intersectBed -a /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/CNCC.bio2out3.bed -b $promoter_1kb -wa -wb > $output_folder/CNCC_overlap_promoter.bed

# peaks that overlap promoter - check if they overlap with other histone modifications (K4_, K27, bivalent)

cp /scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/Promoters_annotated_For_epigenetic_marks.bed $output_folder

# include CpG methylated sites

CpG_annot=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/CpG_island/pnas.1613300114.sd02_hESC_eexported_rearranged.hg38.annotations.bed

grep 'CGI_m' $CpG_annot  > /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/CpG_island/3300114.sd02_hESC_eexported_rearranged.hg38.annotations.CGI_m.bed
grep 'CGI_um' $CpG_annot  > /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/CpG_island/3300114.sd02_hESC_eexported_rearranged.hg38.annotations.CGI_um.bed

CpG_unmeth=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/CpG_island/3300114.sd02_hESC_eexported_rearranged.hg38.annotations.CGI_um.bed
CpG_meth=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/CpG_island/3300114.sd02_hESC_eexported_rearranged.hg38.annotations.CGI_m.bed
#annotate promoters with CpG met and unmethylated
path_annotations=/scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/final_reference_sets_stem_cells_published_data
cd $path_annotations
output_folder='/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Comparison_promoters_BG4_static_peaks'
promoter_1kb=/scratcha/sblab/simeon01/reference_genomes/hg38/gencode.v28.TSS.slop1000.bed
sbatch --mem 12G --wrap "bedtools annotate -i $promoter_1kb -files $CpG_unmeth $CpG_meth -names ESC_CpG_unmeth ESC_CpG_meth -counts > $path_annotations/ESC_CpG_met_annotated_promoters.bed"


## ATAC-SEQ
# annotate promoters with atac seq data
cd /scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output
ls /scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output*.multi2.sorted.bed
CNCC_ATAC.multi2.sorted.bed  ESC_ATAC.multi2.sorted.bed  NSC_ATAC.multi2.sorted.bed

ESC_atac=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output/ESC_ATAC.multi2.sorted.bed
NSC_atac=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output/NSC_ATAC.multi2.sorted.bed
CNCC_atac=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output/CNCC_ATAC.multi2.sorted.bed

# common
# common between 3 sets
multiIntersectBed -i *.multi2.sorted.bed | awk '{if ($4==3) print $0}' | sortBed -i - | mergeBed -i - | sortBed -i -  > ESC_CNCC_NSC.atac.multi3_common.bed

ESC_NSC_CNCC_atac_common=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output/ESC_CNCC_NSC.atac.multi3_common.bed

# pariwise intersect.strict
sbatch --mem 6G --wrap "intersectBed -a $ESC_atac -b $NSC_atac | sortBed -i - | mergeBed -i -> common_ESC_vs_NSC.atac.strict.merged.bed"
sbatch --mem 6G --wrap "intersectBed -a $ESC_atac -b $CNCC_atac| sortBed -i - | mergeBed -i -> common_ESC_vs_CNCC.atac.strict.merged.bed"
sbatch --mem 6G --wrap "intersectBed -a $NSC_atac -b $CNCC_atac| sortBed -i - | mergeBed -i -> common_NSC_vs_CNCC.atac.strict.merged.bed"
  

#annotate promoters with atac-seq data
path_annotations=/scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/final_reference_sets_stem_cells_published_data
cd $path_annotations
promoter_1kb=/scratcha/sblab/simeon01/reference_genomes/hg38/gencode.v28.TSS.slop1000.bed
ESC_NSC_CNCC_atac_common=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output/ESC_CNCC_NSC.atac.multi3_common.bed
ESC_atac=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output/ESC_ATAC.multi2.sorted.bed
NSC_atac=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output/NSC_ATAC.multi2.sorted.bed
CNCC_atac=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output/CNCC_ATAC.multi2.sorted.bed
common_ESC_vs_NSC_atac_strict=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output/common_ESC_vs_NSC.atac.strict.merged.bed
common_ESC_vs_CNCC_atac_strict=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output/common_ESC_vs_CNCC.atac.strict.merged.bed
common_NSC_vs_CNCC_atac_strict=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output/common_NSC_vs_CNCC.atac.strict.merged.bed

sbatch --mem 12G --wrap "bedtools annotate -i $promoter_1kb -files $ESC_atac $CNCC_atac $NSC_atac $ESC_NSC_CNCC_atac_common $common_ESC_vs_NSC_atac_strict $common_ESC_vs_CNCC_atac_strict $common_NSC_vs_CNCC_atac_strict -names ESC_atac CNCC_atac NSC_atac ESC_NSC_CNCC_atac_common common_ESC_vs_NSC_atac_strict common_ESC_vs_CNCC_atac_strict common_NSC_vs_CNCC_atac_strict -counts > $path_annotations/ESC_CNCC_NSC_annotated_promoters_atac_seq.bed"



# 
```

generate table in R with annotated promoters obtained from attemp2
------------------------------------------------------------------

``` bash
cd /scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/final_reference_sets_stem_cells_published_data
ls *annotated_promoters*.bed
#   58382 CNCC_annotated_promoters.CNCC_reprocessed.bed
#   58382 ESC_annotated_promoters.Rada_Iglesias_lifted_hg38.bed
#   58382 ESC_CNCC_NSC_annotated_promoters_G4.bed
#   58382 NSC_annotated_promoters.NPC_renlab.bed

# on local machine
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2
#mkdir 20200421_annotated_promoters_using_histones_G4
cd 20200421_annotated_promoters_using_histones_G4
scp simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/final_reference_sets_stem_cells_published_data/*annotated_promoters*.bed .

scp simeon01@10.20.236.34:/scratcha/sblab/simeon01/reference_genomes/hg38/gencode.v28.annotation.gene_bed6_1kbaroundTSS.bed .
```

Use dplyR to merge tables and produce counts

``` r
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/20200421_annotated_promoters_using_histones_G4')
library(dplyr)
library(plyr)
library(networkD3)
library(tidyr)
# load tables

# ESC_epig_annot
ESC_epig_annot <- read.delim('ESC_annotated_promoters.Rada_Iglesias_lifted_hg38.bed',sep = "\t", header = T, stringsAsFactors = F)
# CNCC_epig_annot
CNCC_epig_annot <- read.delim('CNCC_annotated_promoters.CNCC_reprocessed.bed',sep = "\t", header = T, stringsAsFactors = F)
# NSC_epig_annot
NSC_epig_annot <- read.delim('NSC_annotated_promoters.NPC_renlab.bed',sep = "\t", header = T, stringsAsFactors = F)

#G4_annotations
G4_annotations <- read.delim('ESC_CNCC_NSC_annotated_promoters_G4.bed',sep = "\t", header = T, stringsAsFactors = F)

#atac annotations
atac_annotations <- read.delim('ESC_CNCC_NSC_annotated_promoters_atac_seq.bed',sep = "\t", header = T, stringsAsFactors = F)

#ESC methylation annotations
ESC_CpG_met <- read.delim('ESC_CpG_met_annotated_promoters.bed',sep = "\t", header = T, stringsAsFactors = F)

colnames(ESC_epig_annot)[1:6] <- colnames(CNCC_epig_annot)[1:6] <- colnames(NSC_epig_annot)[1:6] <-colnames(G4_annotations)[1:6] <- colnames(atac_annotations)[1:6] <- colnames(ESC_CpG_met)[1:6] <- c('chr_name','start_pos_prom','end_pos_prom','ens_id','feat1','strand')


# annotation of typology of transcripts (proein_coding, others...)
annot_trascripts <- read.delim('gencode.v28.annotation.gene_bed6_1kbaroundTSS.bed',sep = "\t", header = F, stringsAsFactors = F)
annot_trascripts <- annot_trascripts[,c(7,8)]
colnames(annot_trascripts) <- c('ens_id','transcript_type')
annot_trascripts$transcript_type <- gsub(";","",annot_trascripts$transcript_type)
head(annot_trascripts)

TPM_file <- '/Users/simeon01/Documents/Katie/2018_09_03/RNA_seq/TPM_hg38_all_bio_rep_meanTechReps.csv'
TPM_rnaseq <- read.delim(TPM_file, stringsAsFactors = F,header = T,row.names=NULL, sep = ",")
colnames(TPM_rnaseq)[c(1,2,3)] <- c('ens_symb','ens_id','symb')

#TPM_median_expression <-  TPM_rnaseq %>% group_by(ens_id) %>% mutate(median_TPM_ESC = median(H9_hESC_NA_p10,H9_hESC_NA_p11,H9_hESC_NA_p12,H9_hESC_NA_p3,H9_hESC_NA_p2,na.rm = T)) %>% mutate(median_TPM_NSC = median(H9_NSC_NA_p4.2,H9_NSC_NA_p4.3,H9_NSC_NA_p4.4,H9_NSC_NA_p4.1,na.rm =T)) %>% mutate(median_TPM_NEC = median(H9_neurosphere_NA_p5,H9_neurosphere_NA_p6,H9_neurosphere_NA_p7,H9_neurosphere_NA_p8,H9_neurosphere_NA_p9,na.rm =T)) %>% mutate(median_TPM_CNCC = median(H9_CNCC_A1_p15,H9_CNCC_A2_p16,H9_CNCC_B4_p17,H9_CNCC_A3_p18,H9_CNCC_B3_p19,H9_CNCC_4_p20,H9_CNCC_A1_p21,H9_CNCC_A5_p22,na.rm =T))

library(matrixStats)
TPM_median_expression <- TPM_rnaseq
TPM_median_expression$median_TPM_ESC <- rowMedians(cbind(TPM_rnaseq$H9_hESC_NA_p10,
                                              TPM_rnaseq$H9_hESC_NA_p11,
                                          TPM_rnaseq$H9_hESC_NA_p12,
                                          TPM_rnaseq$H9_hESC_NA_p3,
                                          TPM_rnaseq$H9_hESC_NA_p2),na.rm = T)

TPM_median_expression$median_TPM_NSC <- rowMedians(cbind(TPM_rnaseq$H9_NSC_NA_p4.2,
                                              TPM_rnaseq$H9_NSC_NA_p4.3,
                                          TPM_rnaseq$H9_NSC_NA_p4.4,
                                          TPM_rnaseq$H9_NSC_NA_p4.1,
                                          TPM_rnaseq$H9_hESC_NA_p2),na.rm = T)

TPM_median_expression$median_TPM_NEC <- rowMedians(cbind(TPM_rnaseq$H9_neurosphere_NA_p5,
                                          TPM_rnaseq$H9_neurosphere_NA_p6,
                                          TPM_rnaseq$H9_neurosphere_NA_p7,
                                          TPM_rnaseq$H9_neurosphere_NA_p8,
                                          TPM_rnaseq$H9_neurosphere_NA_p9),na.rm = T)

TPM_median_expression$median_TPM_CNCC <- rowMedians(cbind(TPM_rnaseq$H9_CNCC_A1_p15,
                                          TPM_rnaseq$H9_CNCC_A2_p16,
                                          TPM_rnaseq$H9_CNCC_A3_p18,
                                          TPM_rnaseq$H9_CNCC_B3_p19,
                                          TPM_rnaseq$H9_CNCC_4_p20,
                                          TPM_rnaseq$H9_CNCC_A1_p21,
                                          TPM_rnaseq$H9_CNCC_A5_p22),na.rm = T)

TPM_median_expression$ens_id <- paste0(TPM_median_expression$ens_id,";")

All_anno <- join_all(list(TPM_median_expression,ESC_epig_annot,CNCC_epig_annot,NSC_epig_annot,G4_annotations,atac_annotations,ESC_CpG_met,annot_trascripts), by='ens_id', type='left')
dim(All_anno)

All_anno <- All_anno %>% setNames(make.names(names(.), unique = TRUE)) %>% 
    select(-matches("*\\.[1-9]+$"))
colnames(All_anno)


# create specific new features for lose bival
All_anno <- dplyr::mutate(All_anno,ESC_relax_biv= ifelse(ESC_H3K4me3 >= 1 & ESC_H3K27me3 >= 1 &  ESC_bival ==0, 1,0))
All_anno <- dplyr::mutate(All_anno,CNCC_relax_biv= ifelse(CNCC_H3K4me3 >= 1 & CNCC_H3K27me3 >= 1 &  CNCC_bival ==0, 1,0))
All_anno <- dplyr::mutate(All_anno,NSC_relax_biv= ifelse(NSC_H3K4me3 >= 1 & NSC_H3K27me3 >= 1 &  NSC_bival ==0, 1,0))


All_anno <- dplyr::mutate(All_anno,ESC_status = case_when(
                          (ESC_H3K4me3 == 0 & ESC_bival == 0 &  ESC_relax_biv == 0 & ESC_H3K27me3 == 0) ~ 'Unmarked',
                          (ESC_H3K4me3 >= 1 & ESC_bival == 0 &  ESC_relax_biv == 0 & ESC_H3K27me3 == 0) ~ 'H3K4me3_only',
                          (ESC_H3K4me3 == 0 & ESC_bival == 0 &  ESC_relax_biv == 0 & ESC_H3K27me3 >= 1) ~ 'H3K27me3_only',
                          (ESC_bival >= 1) ~ 'bival',
                          (ESC_relax_biv >= 1 ) ~ 'relax_bival'))


All_anno <- dplyr::mutate(All_anno,CNCC_status = case_when(
                          (CNCC_H3K4me3 == 0 & CNCC_bival == 0 &  CNCC_relax_biv == 0 & CNCC_H3K27me3 == 0) ~ 'Unmarked',
                          (CNCC_H3K4me3 >= 1 & CNCC_bival == 0 &  CNCC_relax_biv == 0 & CNCC_H3K27me3 == 0) ~ 'H3K4me3_only',
                          (CNCC_H3K4me3 == 0 & CNCC_bival == 0 &  CNCC_relax_biv == 0 & CNCC_H3K27me3 >= 1) ~ 'H3K27me3_only',
                          (CNCC_bival >= 1) ~ 'bival',
                          (CNCC_relax_biv >= 1 ) ~ 'relax_bival'))

All_anno <- dplyr::mutate(All_anno,NSC_status = case_when(
                          (NSC_H3K4me3 == 0 & NSC_bival == 0 &  NSC_relax_biv == 0 & NSC_H3K27me3 == 0) ~ 'Unmarked',
                          (NSC_H3K4me3 >= 1 & NSC_bival == 0 &  NSC_relax_biv == 0 & NSC_H3K27me3 == 0) ~ 'H3K4me3_only',
                          (NSC_H3K4me3 == 0 & NSC_bival == 0 &  NSC_relax_biv == 0 & NSC_H3K27me3 >= 1) ~ 'H3K27me3_only',
                          (NSC_bival >= 1) ~ 'bival',
                          (NSC_relax_biv >= 1 ) ~ 'relax_bival'))

All_anno %>% filter(CNCC_G4==1 & ESC_G4==0) %>% tally()
dim(All_anno)
#[1] 58381    57

save(All_anno,file="All_anno.RData")
write.table(All_anno,file="All_anno_updated.txt", quote = F, sep= "\t",col.names = NA)
assign('OLD_all_anno', get(load('OLD_')))
# ### ### ### # ### ### ### # ### ### ### # ### ### ### # ### ### ### # ### ### ###
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/20200421_annotated_promoters_using_histones_G4')

load("All_anno.RData")
library(dplyr)
library(plyr)
library(networkD3)

# ||    **                     FILTER                       **       ||
# || ########## !!!!! HERE IS THE FILTERING CRITERIA !!! ########### ||
# \/                                                                 \/
# T <- All_anno %>% filter(CNCC_G4>=1 & ESC_G4>=1) %>% select(ens_id,ESC_status,CNCC_status) 
# T <- All_anno %>% filter(common_ESC_vs_CNCC.strict>=1) %>% select(ens_id,ESC_status,CNCC_status)
 
T <- All_anno %>% filter(ESC_G4>=1 & transcript_type=="protein_coding") %>% select(ens_id,ESC_status,CNCC_status)
TT <- All_anno %>%  select(ens_id,ESC_CpG_meth,ESC_CpG_unmeth)
# T <- All_anno %>% filter(ESC_G4>=1) %>% select(ens_id,ESC_status,CNCC_status) 
# /\                                                                 /\
# || ########## ########### ############## ###########  ############ ||

print(dim(T))

links <- plyr::count(T[,-1])

colnames(links) <- c("source","target","value")
links$source <- paste0(links$source,'_1')
links$target <- paste0(links$target,'_2')
print(links)

nodes <- data.frame(
  name=c(as.character(links$source), 
  as.character(links$target)) %>% unique())

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

p <- sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE)

p


## select only one mark per time and stratify by G4

D <- All_anno %>% filter(CNCC_status =="H3K4me3_only") %>% select(ens_id, ESC_status, CNCC_status, ESC_G4, CNCC_G4)
 
#create new column with groups based on G4peakstatus
D$ESC_G4_status = ifelse(D$ESC_G4 >0,"G4pos", "G4neg")
D$CNCC_G4_status = ifelse(D$CNCC_G4 >0,"G4pos", "G4neg")   
 
D$start <- paste(D$ESC_status, D$ESC_G4_status, sep='_')
D$finish <- paste(D$CNCC_status, D$CNCC_G4_status, sep='_')
 
T <- D %>%  select(ens_id,start,finish)
 
# ### ### ### # ### ### ### # ### ### ### # ### ### ### # ### ### ### # ### ### ###


All_anno2 %>% 
    dplyr::mutate(n = 1) %>% 
    spread(CNCC_status, n, fill=0) %>% 
    select(CNCC_status,ESC_status) %>% 
    {crossprod(as.matrix(.))} %>% 
    replace(lower.tri(., diag=T), NA) %>%
    reshape2::melt(na.rm=T) %>%
    unite('Pair', c('Var1', 'Var2'), sep=", ")
# create nodes dataframe
regions <- unique(unique(All_anno$ESC_status))
nodes <- data.frame(node = c(1:5), 
                    chr_status = regions)
#create links dataframe
results <- dplyr::left_join(results,nodes, by="chr_status")
results <- merge(results, nodes, by.x = "Region", by.y = "name")
results <- merge(results, nodes, by.x = "result", by.y = "name")
links <- results[ , c("node.x", "node.y", "vote")]
colnames(links) <- c("source", "target", "value")

# check various tally

All_anno %>% dplyr::group_by('ens_id') %>% filter(CNCC_G4==1 & (CNCC_relax_biv ==1|CNCC_H3K4me3==1 |CNCC_H3K27me3 ==1)) %>% tally()
All_anno %>% dplyr::group_by('ens_id') %>% filter(CNCC_G4_specific==1 & ESC_G4_specific==1) %>% tally()
All_anno %>% dplyr::group_by('ens_id') %>% filter(CNCC_G4==1 & ESC_G4==0 & NSC_G4==0) %>% tally()
All_anno %>% dplyr::group_by('ens_id') %>% filter(CNCC_G4==0 & ESC_G4==1 & NSC_G4==0) %>% tally()
All_anno %>% dplyr::group_by('ens_id') %>% filter(CNCC_G4==0 & ESC_G4==0 & NSC_G4==1) %>% tally()


temp_focus_on_common <-  All_anno %>% dplyr::group_by('ens_id') %>% dplyr::filter(ESC_CNCC_G4_common==1) 

temp_focus_on_common2 <- temp_focus_on_common %>% dplyr::select(ESC_H3K4me3,ESC_H3K27me3,ESC_bival,ESC_relax_biv,ESC_CNCC_G4_common,ESC_G4,CNCC_H3K4me3,CNCC_H3K27me3,CNCC_bival) 

temp_focus_on_common2[temp_focus_on_common2>1] <-  1
```

analysis on peak sizes
----------------------

-   We explore the general peak size across the 3 cell lines

-   We explore also how the peak size is distrubted when we focus on peaks overlapping promoters of genes known to be involved in specific functional categories (based on GO)

-   we look at the regions that are common across cell lines and trace how size changes if they are \*\* maked by a specific mark (and stay the same) \*\* not marked \*\* belong to contact regions in HiC data (from ESC)

paper that explored maker size and relationshipt with other features.

<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6249108/> <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6085762/> <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4137894/>

G4 regions overlapping promoters
================================

Here we focus first on regions that are common across CNCC, ESC and NSC

``` bash

cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2

# cell line BG4 peaks consensus are bio2out3.bed
# regions present in all 3 also with 100bp proximity
multiIntersectBed -i *bio2out3.bed | awk '{if ($4==3) print $0}'| sortBed -i - | mergeBed -i - > common_regions_ESC_NSC_CNCC.bed

multiIntersectBed -i *bio2out3.bed | awk '{if ($4==3) print $0}'| sortBed -i - | mergeBed -i - -d 1000 > common_regions_ESC_NSC_CNCC.merge1000.bed

# for each bio2out3.bed intersect with common regions
for file in *bio2out3.bed
do
  intersectBed -a $file -b common_regions_ESC_NSC_CNCC.bed -wa -wb > ${file%%.bed}.intersect_common_regions_ESC_NSC_CNCC.bed
  intersectBed -a $file -b common_regions_ESC_NSC_CNCC.merge1000.bed -wa -wb > ${file%%.bed}.intersect_common_regions_ESC_NSC_CNCC.merge1000.bed
done  


cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2 
intersectBed -a ESC.bio2out3.intersect_common_regions_ESC_NSC_CNCC.bed -b /Users/simeon01/Documents/genomes/hg38/gencode.v28.annotation.gene_bed6_1kbaroundTSS.bed -wa -wb > ESC.bio2out3.intersect_common_regions_ESC_NSC_CNCC.G4_overlap_prom.bed 


annotateBed -i ESC.bio2out3.intersect_common_regions_ESC_NSC_CNCC.bed -files /Users/simeon01/Documents/genomes/hg38/gencode.v28.annotation.gene_bed6_1kbaroundTSS.bed -names G4_in_1kbProm -counts > ESC.bio2out3.intersect_common_regions_ESC_NSC_CNCC.G4_peaks_anotated_by_promoter.bed

annotateBed -i ESC.bio2out3.intersect_common_regions_ESC_NSC_CNCC.bed -files /Users/simeon01/Documents/OQs/OQ_plus_hits.lifted_hg19_to_hg38.Minus_Strand_bed6.bed /Users/simeon01/Documents/OQs/OQ_minus_hits.lifted_hg19_to_hg38.Plus_Strand_bed6.bed -names Minus_Strand_bed6 Plus_Strand_bed6  -counts > ESC.bio2out3.intersect_common_regions_ESC_NSC_CNCC.G4_peaks_anotated_by_OQsStrand.bed


# additional question - 
# annotate common peaks by k27, k4 and bivalent
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2 
ESC_K27=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/rada_iglesia_ESC_macs2_peaks/Rada_Iglesias_lifted_hg38.H3K27me3.bed
ESC_K4=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/rada_iglesia_ESC_macs2_peaks/Rada_Iglesias_lifted_hg38.H3K4me3.bed
ESC_bival=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/bivalent_promoters/ESC_bival_regions.Rada_Iglesias_lifted_hg38.H3K27me3.H3K4me3.intersect.bed

CNCC_K27=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/prescott_CNCC_macs2_peaks/H3K27me3_SRR2096422_vs_SRR2096456_peaks.broadPeak
CNCC_K4=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/prescott_CNCC_macs2_peaks/H3K4me3_human1.multi2.bed
CNCC_bival=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/bivalent_promoters/CNCC_bival_regions.H3K27me3_SRR2096422.broad.H3K4me3_human1.multi2.intersect.bed

NSC_K27=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/xie_NSC_macs2_peaks/renlab.H3K27me3.broad.NPC.intersect.bed
NSC_K4=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/xie_NSC_macs2_peaks/renlab.H3K4me3.NPC.intersect.bed 
NSC_bival=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/bivalent_promoters/NSC_bival_regions.renlab.H3K27me3.broad.H3K4me3.intersect.bed
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2 
cut -f 4,5,6 ESC.bio2out3.intersect_common_regions_ESC_NSC_CNCC.bed | annotateBed -i - -files $ESC_K27 $ESC_K4 $ESC_bival $CNCC_K27 $CNCC_K4 $CNCC_bival $NSC_K27 $NSC_K4 $NSC_bival -names ESC_K27 ESC_K4 ESC_bival CNCC_K27 CNCC_K4 CNCC_bival NSC_K27 NSC_K4 NSC_bival -counts > intersect_common_regions_ESC_NSC_CNCC.G4_peaks_anotated_by_epigenetic_marks.bed

# intersect individual peaks with TSS and then with epigenetic marks
for file in *bio2out3.bed
  do
  intersectBed -a $file -b /Users/simeon01/Documents/genomes/hg38/gencode.v28.annotation.gene_bed6_1kbaroundTSS.bed -wa -wb | annotateBed -i - -files common_regions_ESC_NSC_CNCC.bed $ESC_K27 $ESC_K4 $ESC_bival $CNCC_K27 $CNCC_K4 $CNCC_bival $NSC_K27 $NSC_K4 $NSC_bival -names common_BG4_regions ESC_K27 ESC_K4 ESC_bival CNCC_K27 CNCC_K4 CNCC_bival NSC_K27 NSC_K4 NSC_bival -counts > ${file%%.bed}.G4_peaks_overlapping_prom_anotated_by_epigenetic_marks.bed
  done
```

Do something similar for regions that are specific to each cell lines

``` bash


# additional question - 
# annotate common peaks by k27, k4 and bivalent
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2 
ESC_K27=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/rada_iglesia_ESC_macs2_peaks/Rada_Iglesias_lifted_hg38.H3K27me3.bed
ESC_K4=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/rada_iglesia_ESC_macs2_peaks/Rada_Iglesias_lifted_hg38.H3K4me3.bed
ESC_bival=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/bivalent_promoters/ESC_bival_regions.Rada_Iglesias_lifted_hg38.H3K27me3.H3K4me3.intersect.bed

CNCC_K27=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/prescott_CNCC_macs2_peaks/H3K27me3_SRR2096422_vs_SRR2096456_peaks.broadPeak
CNCC_K4=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/prescott_CNCC_macs2_peaks/H3K4me3_human1.multi2.bed
CNCC_bival=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/bivalent_promoters/CNCC_bival_regions.H3K27me3_SRR2096422.broad.H3K4me3_human1.multi2.intersect.bed

NSC_K27=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/xie_NSC_macs2_peaks/renlab.H3K27me3.broad.NPC.intersect.bed
NSC_K4=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/xie_NSC_macs2_peaks/renlab.H3K4me3.NPC.intersect.bed 
NSC_bival=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/bivalent_promoters/NSC_bival_regions.renlab.H3K27me3.broad.H3K4me3.intersect.bed
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2 

# intersect individual peaks with TSS and then with epigenetic marks


for file in *_specific.20200421.bed
do
  
  annotateBed -i $file -files /Users/simeon01/Documents/genomes/hg38/gencode.v28.annotation.gene_bed6_1kbaroundTSS.bed -names G4_in_1kbProm -counts > ${file%%.bed}.G4_peaks_anotated_by_promoter.bed
  
  annotateBed -i $file -files /Users/simeon01/Documents/OQs/OQ_plus_hits.lifted_hg19_to_hg38.Minus_Strand_bed6.bed /Users/simeon01/Documents/OQs/OQ_minus_hits.lifted_hg19_to_hg38.Plus_Strand_bed6.bed -names Minus_Strand_bed6 Plus_Strand_bed6  -counts > ${file%%.bed}.G4_peaks_anotated_by_OQsStrand.bed
  
  annotateBed -i $file -files $ESC_K27 $ESC_K4 $ESC_bival $CNCC_K27 $CNCC_K4 $CNCC_bival $NSC_K27 $NSC_K4 $NSC_bival -names ESC_K27 ESC_K4 ESC_bival CNCC_K27 CNCC_K4 CNCC_bival NSC_K27 NSC_K4 NSC_bival -counts > ${file%%.bed}.G4_peaks_overlapping_prom_anotated_by_epigenetic_marks.bed
  
done
```

Overlap the full list of cell peaks with the respective H3K4me3 marks (in the respective cell line)

``` bash
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2 

ESC_K4=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/rada_iglesia_ESC_macs2_peaks/Rada_Iglesias_lifted_hg38.H3K4me3.bed
CNCC_K4=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/prescott_CNCC_macs2_peaks/H3K4me3_human1.multi2.bed
NSC_K4=/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/bw/xie_NSC_macs2_peaks/renlab.H3K4me3.NPC.intersect.bed 

intersectBed -a ESC.bio2out3.bed -b $ESC_K4 -wa -wb > ESC.bio2out3.intersect.ESC_K4.bed
intersectBed -a CNCC.bio2out3.bed -b $CNCC_K4 -wa -wb > CNCC.bio2out3.intersect.CNCC_K4.bed
intersectBed -a NSC.bio2out3.bed -b $NSC_K4 -wa -wb > NSC.bio2out3.intersect.NSC_K4.bed
```

Analysis of peak size on G4 regions in common between ESC CNCC and NSC
----------------------------------------------------------------------

``` r
rm(list = ls())
library(tidyverse)
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2')
ESC_common <- read.table('ESC.bio2out3.intersect_common_regions_ESC_NSC_CNCC.bed', stringsAsFactors = F)
CNCC_common <- read.table('CNCC.bio2out3.intersect_common_regions_ESC_NSC_CNCC.bed', stringsAsFactors = F)
NSC_common <- read.table('NSC.bio2out3.intersect_common_regions_ESC_NSC_CNCC.bed', stringsAsFactors = F)
colnames(ESC_common) <- c('ESC_chr','ESC_start','ESC_end','common_chr','common_start','common_end')
colnames(CNCC_common) <- c('CNCC_chr','CNCC_start','CNCC_end','common_chr','common_start','common_end')
colnames(NSC_common) <- c('NSC_chr','NSC_start','NSC_end','common_chr','common_start','common_end')

promoter_with_commonG4 <- read.delim('ESC.bio2out3.intersect_common_regions_ESC_NSC_CNCC.G4_peaks_anotated_by_promoter.bed', stringsAsFactors = F, header = T)
G4_common_oqs_strand <- read.delim('ESC.bio2out3.intersect_common_regions_ESC_NSC_CNCC.G4_peaks_anotated_by_OQsStrand.bed', stringsAsFactors = F, header = T)
colnames(promoter_with_commonG4)[1:6] <- c('ESC_chr','ESC_start','ESC_end','common_chr','common_start','common_end')
colnames(G4_common_oqs_strand)[1:6] <- c('ESC_chr','ESC_start','ESC_end','common_chr','common_start','common_end')

TSS_1base <- read.delim('/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS.minus1base.bed',stringsAsFactors = F,header = F)
colnames(TSS_1base) <- c('TSS_chr','TSS_start','TSS_end','ens_id','feature','TSS_strand')

Common_regions_with_epigenetic_marks <- read.delim('intersect_common_regions_ESC_NSC_CNCC.G4_peaks_anotated_by_epigenetic_marks.bed',stringsAsFactors = F,header = T)
colnames(Common_regions_with_epigenetic_marks)[1:3] <- c('common_chr','common_start','common_end')

Promoter_overlapping_G4 <- read.delim('ESC.bio2out3.intersect_common_regions_ESC_NSC_CNCC.G4_overlap_prom.bed',stringsAsFactors = F, header = F)
colnames(Promoter_overlapping_G4)[1:6] <- c('ESC_chr','ESC_start','ESC_end','common_chr','common_start','common_end')
colnames(Promoter_overlapping_G4)[7:ncol(Promoter_overlapping_G4)] <- c('Prom_chr','Prom_start','Prom_end','prom_feature','prom_strand','prom_strand_repeat','ens_id','transcript_type')
Promoter_overlapping_G4$common_id <- paste(Promoter_overlapping_G4$common_chr,Promoter_overlapping_G4$common_start,Promoter_overlapping_G4$common_end,sep = "_")
Promoter_overlapping_G4 <- Promoter_overlapping_G4 %>% dplyr::select(common_id,transcript_type,prom_feature ,prom_strand,ens_id)


ESC_common$common_id <- paste(ESC_common$common_chr,ESC_common$common_start,ESC_common$common_end,sep = "_")
CNCC_common$common_id <- paste(CNCC_common$common_chr,CNCC_common$common_start,CNCC_common$common_end,sep = "_")
NSC_common$common_id <- paste(NSC_common$common_chr,NSC_common$common_start,NSC_common$common_end,sep = "_")
promoter_with_commonG4$common_id <- paste(promoter_with_commonG4$common_chr,promoter_with_commonG4$common_start,promoter_with_commonG4$common_end,sep = "_")
G4_common_oqs_strand$common_id <- paste(G4_common_oqs_strand$common_chr,G4_common_oqs_strand$common_start,G4_common_oqs_strand$common_end,sep = "_")
Common_regions_with_epigenetic_marks$common_id <- paste(Common_regions_with_epigenetic_marks$common_chr,Common_regions_with_epigenetic_marks$common_start,Common_regions_with_epigenetic_marks$common_end,sep = "_")

ESC_common <- ESC_common  %>% dplyr::mutate( ESC_size = ESC_end-ESC_start) %>% dplyr::select(common_id,ESC_size)
NSC_common <- NSC_common  %>% dplyr::mutate( NSC_size = NSC_end-NSC_start) %>% dplyr::select(common_id,NSC_size)
CNCC_common <- CNCC_common  %>% dplyr::mutate( CNCC_size = CNCC_end-CNCC_start)%>% dplyr::select(common_id,CNCC_size)

#ESC_common <- ESC_common %>% dplyr::select(common_id,ESC_size)
#NSC_common <- NSC_common %>% dplyr::select(common_id,NSC_size)
#CNCC_common <- CNCC_common %>% dplyr::select(common_id,CNCC_size)

All_regions_joined <- dplyr::left_join(ESC_common,CNCC_common,by="common_id") %>% dplyr::left_join(NSC_common, by="common_id") %>% dplyr::left_join(promoter_with_commonG4, by="common_id") %>% dplyr::left_join(G4_common_oqs_strand, by="common_id") %>% dplyr::left_join(Promoter_overlapping_G4,by="common_id") %>% dplyr::left_join(Common_regions_with_epigenetic_marks,by="common_id")


#export bed files with coordinates of the center for BG+ and BG-
# and export bed files with coordinates of the genes having a TSS in the common regions 
# all these will be used to generate TSS plots
All_regions_joined <- All_regions_joined %>% dplyr::mutate(G4strand = case_when(
  (Minus_Strand_bed6 == 0 & Plus_Strand_bed6 >= 1) ~ 'Plus',
  (Minus_Strand_bed6 >= 1 & Plus_Strand_bed6 == 0) ~ 'Minus',
  (Minus_Strand_bed6 >= 1 & Plus_Strand_bed6 >= 1) ~ 'Both_strand',
  (Minus_Strand_bed6 == 0 & Plus_Strand_bed6 == 0) ~ 'No_strand'))

All_regions_joined <- All_regions_joined %>% dplyr::mutate(G4strand_symb = case_when(
  (Minus_Strand_bed6 == 0 & Plus_Strand_bed6 >= 1) ~ '+',
  (Minus_Strand_bed6 >= 1 & Plus_Strand_bed6 == 0) ~ '-',
  (Minus_Strand_bed6 >= 1 & Plus_Strand_bed6 >= 1) ~ 'Both_strand',
  (Minus_Strand_bed6 == 0 & Plus_Strand_bed6 == 0) ~ 'No_strand'))

All_regions_joined  <- All_regions_joined %>% mutate(common_G4region_center_coord = round((common_start.x+common_end.x)/2))

All_regions_joined_with_TSScoord <- dplyr::left_join(All_regions_joined,TSS_1base, by = "ens_id")

# export to bed format center of regions 
Center_comon_regions <- All_regions_joined %>% dplyr::filter(G4strand=="Plus"|G4strand== "Minus") %>% dplyr::group_by(common_id) %>% dplyr::select(common_chr.x,common_G4region_center_coord ,common_G4region_center_coord,G4strand_symb)

#write.table(file = 'G4_regions_common_to_ESC_CNCC_NSC_with_strand.bed',sprintf("%s\t%i\t%i\t%s\t.\t%s\n",Center_comon_regions$common_chr.x,Center_comon_regions$common_G4region_center_coord,Center_comon_regions$common_G4region_center_coord,Center_comon_regions$common_id,Center_comon_regions$G4strand_symb), quote = F, col.names = F, row.names = F)
#system("bedtools slop -i G4_regions_common_to_ESC_CNCC_NSC_with_strand.bed -g /Users/simeon01/Documents/genomes/hg38/hg38.sorted.genome -l 0 -r 1 > G4_regions_common_to_ESC_CNCC_NSC_with_strand.slop1base.bed", intern = T)


#export to bed format the TSS overlapping G4 regions
TSS_Center_comon_regions <- All_regions_joined_with_TSScoord %>% dplyr::filter(All_regions_joined$G4_in_1kbProm >0) %>% dplyr::group_by(ens_id) %>% dplyr::select(TSS_chr,TSS_start,TSS_end,ens_id,feature,TSS_strand)

#write.table(file="TSS_overlapping_G4_regions_common_to_ESC_CNCC_NSC_with_strand.bed",TSS_Center_comon_regions, sep = "\t",quote = F, col.names = F,row.names = F)


#focusing simply to the peaks that are in common between cells lines apply KS.test and plot also CDF

All_regions_joined_to_melt <- All_regions_joined %>% dplyr::select(common_id,ESC_size,CNCC_size,NSC_size)

Filtered_regions_joined_to_melt <- All_regions_joined %>% dplyr::filter(G4_in_1kbProm==1) %>% dplyr::select(common_id,ESC_size,CNCC_size,NSC_size)

Filtered2_regions_joined_to_melt <- All_regions_joined %>% dplyr::filter(G4_in_1kbProm==0) %>% dplyr::select(common_id,ESC_size,CNCC_size,NSC_size)

Filtered_ESCbival_regions_joined_to_melt <- All_regions_joined %>% dplyr::filter(ESC_bival>=1 & CNCC_bival>=1) %>% dplyr::select(common_id,ESC_size,CNCC_size,NSC_size)
Filtered_ESCK4_regions_joined_to_melt <- All_regions_joined %>% dplyr::filter(ESC_K4>=1) %>% dplyr::select(common_id,ESC_size,CNCC_size,NSC_size)
Filtered_ESCK27_regions_joined_to_melt <- All_regions_joined %>% dplyr::filter(ESC_K27==0) %>% dplyr::select(common_id,ESC_size,CNCC_size,NSC_size)

library(reshape2)
All_regions_joined_melted <- melt(All_regions_joined_to_melt)
ggplot(All_regions_joined_melted,aes(value, color = variable )) + stat_ecdf() +xlim(0,1000) + ggtitle('all G4 common regions')

Filtered_regions_joined_melted <- melt(Filtered_regions_joined_to_melt)
ggplot(Filtered_regions_joined_melted,aes(value, color = variable )) + stat_ecdf()+xlim(0,1000) + ggtitle('all G4 common regions -  overlapping prom')

Filtered2_regions_joined_melted <- melt(Filtered2_regions_joined_to_melt)
ggplot(Filtered2_regions_joined_melted,aes(value, color = variable )) + stat_ecdf()+xlim(0,1000) + ggtitle('all G4 common regions - NOT overlapping prom')

Filtered_ESCbival_regions_joined_melted <- melt(Filtered_ESCbival_regions_joined_to_melt)
ggplot(Filtered_ESCbival_regions_joined_melted,aes(value, color = variable )) + stat_ecdf()+xlim(0,1000) + ggtitle('all G4 common regions ESbival')

Filtered_ESCK4_regions_joined_melted <- melt(Filtered_ESCK4_regions_joined_to_melt)
ggplot(Filtered_ESCK4_regions_joined_melted,aes(value, color = variable )) + stat_ecdf()+xlim(0,1000) + ggtitle('all G4 common regions ESCK4')

Filtered_ESCK27_regions_joined_melted <- melt(Filtered_ESCK27_regions_joined_to_melt)
ggplot(Filtered_ESCK27_regions_joined_melted,aes(value, color = variable )) + stat_ecdf()+xlim(0,1000) + ggtitle('all G4 common regions ESCK27')
```

density plots on ChIP-files using the regions exported before
-------------------------------------------------------------

``` bash
scp TSS_overlapping_G4_regions_common_to_ESC_CNCC_NSC_with_strand.bed simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/reference_bed/

scp G4_regions_common_to_ESC_CNCC_NSC_with_strand.bed simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/reference_bed/
scp G4_regions_common_to_ESC_CNCC_NSC_with_strand.slop1base.bed simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/reference_bed/
cd /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/
mkdir density_at_common_regions

ref1=/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/reference_bed/TSS_overlapping_G4_regions_common_to_ESC_CNCC_NSC_with_strand.bed
ref2=/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/reference_bed/G4_regions_common_to_ESC_CNCC_NSC_with_strand.slop1base.bed

cd /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW
for file in *.rpm.bw
do
  echo $file
  #sbatch --mem 8G --wrap="computeMatrix reference-point -S $file -R $ref1 -bs 10 -a 1000 -b 1000 -o /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.TSS_overlapping_G4_regions_common.gz && plotProfile --matrixFile /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.TSS_overlapping_G4_regions_common.gz -o /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.TSS_overlapping_G4_regions_common.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.TSS_overlapping_G4_regions_common.gz -o /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.TSS_overlapping_G4_regions_common.heatmap.png --averageType median"
  echo "**"
  #sbatch --mem 8G --wrap="computeMatrix reference-point -S $file -R $ref2 -bs 10 -a 1000 -b 1000 -o /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.G4_regions_common.gz && plotProfile --matrixFile /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.G4_regions_common.gz -o /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.G4_regions_common.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.G4_regions_common.gz -o /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.G4_regions_common.heatmap.png --averageType median"
  echo " ==== "
  
  
  #sbatch --mem 8G --wrap="computeMatrix reference-point -S $file -R $ref1 -bs 5 -a 500 -b 500 -o /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.TSS_overlapping_G4_regions_common.500bp_5bp.gz"
  echo "**"
  #sbatch --mem 8G --wrap="computeMatrix reference-point -S $file -R $ref2 -bs 5 -a 500 -b 500 -o /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.G4_regions_common.500bp_5bp.gz"
  
  sbatch --mem 8G --wrap="computeMatrix reference-point -S $file -R $ref1 -bs 2 -a 300 -b 300 -o /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.TSS_overlapping_G4_regions_common.300bp_2bp.gz"
  echo "**"
  sbatch --mem 8G --wrap="computeMatrix reference-point -S $file -R $ref2 -bs 2 -a 300 -b 300 -o /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/${file%%.bw}.G4_regions_common.300bp_2bp.gz"
done


#locally
# copy back the results to combine segnal by cell type
/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_at_common_G4_regions
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_at_common_G4_regions
scp simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_at_common_regions/* .
```

locally combine densities from computeMatrix
--------------------------------------------

comparison of peaks -- analysis based on old annotations
--------------------------------------------------------

``` r
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/Comparison_promoters_BG4_static_peaks')
Prom_annotation_BG4 <- read.table('Promoters_annotated_For_BG4peaks.bed',stringsAsFactors = F, sep = "\t")
Prom_annotation_epigenetics <- read.table('Promoters_annotated_For_epigenetic_marks.bed', stringsAsFactors = F, sep = "\t", header = T)
colnames(Prom_annotation_BG4) <- c('chr','start','end','ens_id','field','strand','CNCC_BG4','NSC_BG4','ESC_BG4')
colnames(Prom_annotation_epigenetics) <- c('chr','start','end','ens_id','field','strand','ESC_bival','ESC_H3K4me3','ESC_H3K27me3','CNCC_bival','CNCC_H3K4me3','CNCC_H3K27me3','NSC_bival','NSC_H3K4me3','NSC_H3K27me3')

# new annotations
load("/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/20200421_annotated_promoters_using_histones_G4/All_anno.RData")
#All_anno$ens_id <- gsub(';','',All_anno$ens_id)

ESC_peaks <- read.table('ESC_overlap_promoter.bed', stringsAsFactors = F, sep = "\t")
NSC_peaks <- read.table('NSC_overlap_promoter.bed', stringsAsFactors = F, sep = "\t")
CNCC_peaks <- read.table('CNCC_overlap_promoter.bed', stringsAsFactors = F, sep = "\t")

ESC_peaks <- data.frame(ESC_peaks)
NSC_peaks <- data.frame(NSC_peaks)
CNCC_peaks <- data.frame(CNCC_peaks)
colnames(ESC_peaks) <- colnames(NSC_peaks) <- colnames(CNCC_peaks) <- c('chr','peak_start','peak_end','chrP','startP','endP','ens_id','fieldP','strandP')

ESC_peaks$peak_size <- ESC_peaks$peak_end-ESC_peaks$peak_start
ESC_peaks$peak_id <- paste(ESC_peaks$chr,ESC_peaks$peak_start,ESC_peaks$peak_end,sep="_")

NSC_peaks$peak_size <- NSC_peaks$peak_end-NSC_peaks$peak_start
NSC_peaks$peak_id <- paste(NSC_peaks$chr,NSC_peaks$peak_start,NSC_peaks$peak_end,sep="_")

CNCC_peaks$peak_size <- CNCC_peaks$peak_end-CNCC_peaks$peak_start
CNCC_peaks$peak_id <- paste(CNCC_peaks$chr,CNCC_peaks$peak_start,CNCC_peaks$peak_end,sep="_")


# ESC_table <- left_join(ESC_peaks,Prom_annotation_epigenetics,by="ens_id")
# NSC_table <- left_join(NSC_peaks,Prom_annotation_epigenetics,by="ens_id")
# CNCC_table <- left_join(CNCC_peaks,Prom_annotation_epigenetics,by="ens_id")

ESC_table <- left_join(ESC_peaks,All_anno,by="ens_id")
NSC_table <- left_join(NSC_peaks,All_anno,by="ens_id")
CNCC_table <- left_join(CNCC_peaks,All_anno,by="ens_id")

#write.table(ESC_table,file = 'ESC_table_with_prom_annotations.csv', col.names = NA, sep=",",quote = F)
#write.table(NSC_table,file = 'NSC_table_with_prom_annotations.csv', col.names = NA, sep=",",quote = F)
#write.table(CNCC_table,file = 'CNCC_table_with_prom_annotations.csv', col.names = NA, sep=",",quote = F)

#ESC_table1 <- ESC_table %>% dplyr::select(ens_id,peak_id,peak_size,ESC_bival,ESC_H3K4me3,ESC_H3K27me3,CNCC_bival,CNCC_H3K4me3,CNCC_H3K27me3,NSC_bival,NSC_H3K4me3,NSC_H3K27me3) 
ESC_table1 <- ESC_table %>% dplyr::select(ens_id,peak_id,peak_size,ESC_bival,ESC_H3K4me3,ESC_H3K27me3,CNCC_bival,CNCC_H3K4me3,CNCC_H3K27me3,NSC_bival,NSC_H3K4me3,NSC_H3K27me3) 
ESC_table1_extended <- mutate(ESC_table1,cell_epi_type= ifelse(ESC_bival >0 ,'ESC_bival',ifelse(ESC_H3K4me3>0,'ESC_K4',ifelse(ESC_H3K27me3>0,'ESC_K27','NA'))))

ESC_table1_extended <- mutate(ESC_table1,cell_epi_type= ifelse(ESC_bival >0 ,'ESC_bival',ifelse(ESC_H3K4me3>0,'ESC_K4',ifelse(ESC_H3K27me3>0,'ESC_K27','NA'))))

data_temp <- dim(ESC_table)

 
data_frame_bival <- data.frame(p_sizes = c(ESC_table$peak_size[which(ESC_table$ESC>0)],
                                           NSC_table$peak_size[which(NSC_table$NSC_bival>0)],
                                           CNCC_table$peak_size[which(CNCC_table$CNCC_bival>0)]),
                               cell_type = c(rep('ESC',length(which(ESC_table$ESC_bival>0))),
                                             rep('NSC',length(which(NSC_table$NSC_bival>0))),
                                             rep('CNCC',length(which(CNCC_table$CNCC_bival>0)))))

data_frame_K4 <- data.frame(p_sizes = c(ESC_table$peak_size[which(ESC_table$ESC_H3K4me3>0)],
                                           NSC_table$peak_size[which(NSC_table$NSC_H3K4me3>0)],
                                           CNCC_table$peak_size[which(CNCC_table$CNCC_H3K4me3>0)]),
                               cell_type = c(rep('ESC',length(which(ESC_table$ESC_H3K4me3>0))),
                                             rep('NSC',length(which(NSC_table$NSC_H3K4me3>0))),
                                             rep('CNCC',length(which(CNCC_table$CNCC_H3K4me3>0)))))

data_frame_K27 <- data.frame(p_sizes = c(ESC_table$peak_size[which(ESC_table$ESC_H3K27me3>0)],
                                           NSC_table$peak_size[which(NSC_table$NSC_H3K27me3>0)],
                                           CNCC_table$peak_size[which(CNCC_table$CNCC_H3K27me3>0)]),
                               cell_type = c(rep('ESC',length(which(ESC_table$ESC_H3K27me3>0))),
                                             rep('NSC',length(which(NSC_table$NSC_H3K27me3>0))),
                                             rep('CNCC',length(which(CNCC_table$CNCC_H3K27me3>0)))))


data_frame_K27 %>% group_by(cell_type) %>%summarise_each(funs(median)) %>% gather(Var, Val, -vs) %>%
    ggplot(., aes(x=Var, y=Val, fill=vs))+
    geom_bar(stat='identity', position='dodge')
ESC_table1_extended <- mutate(ESC_table1_extended,NSC_epi_type= ifelse(NSC_bival >0 ,'NSC_bival',ifelse(NSC_H3K4me3>0,'NSC_K4',ifelse(NSC_H3K27me3>0,'NSC_K27','NA'))))
ESC_table1_extended <- mutate(ESC_table1_extended,NSC_epi_type= ifelse(NSC_bival >0 ,'NSC_bival',ifelse(NSC_H3K4me3>0,'NSC_K4',ifelse(NSC_H3K27me3>0,'NSC_K27','NA'))))

pdf('BG4_peakSize_at_bivalent_cell_specific.pdf')
p <- data_frame_bival %>% ggplot(aes(y=p_sizes,x=cell_type )) +geom_boxplot() +ylim(0,500) + ggtitle("BG4 peaks at bivalent promoters")
print(p)
dev.off()

pdf('BG4_peakSize_at_K4_cell_specific.pdf')
p <- data_frame_K4 %>% ggplot(aes(y=p_sizes,x=cell_type )) +geom_boxplot() +ylim(0,500) + ggtitle("BG4 peaks at H3K4me3 promoters")
print(p)
dev.off()

pdf('BG4_peakSize_at_K27_cell_specific.pdf')
p <- data_frame_K27 %>% ggplot(aes(y=p_sizes,x=cell_type )) +geom_boxplot() +ylim(0,500) + ggtitle("BG4 peaks at H3K27me3 promoters")
print(p)
dev.off()
```

``` bash
# for the confirmed peaks generate list of transcripts that have TSS in BG4 peaks
# tss-bed file
tss=/scratcha/sblab/simeon01/reference_genomes/hg38/gencode.v28.TSS.bed
path_analysis=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
scenario=(scenario1_multi2 scenario2_5out9 scenario3_multi2_extended_1kb)
for case in ${scenario[@]}
  do
  cd $path_analysis/$case
  esc_bed=`ls ES*.bed`
  nsc_bed=`ls NSC*.bed`
  cncc_bed=`ls CNCC*.bed`
  mkdir gene_lists
  #intersect with TSS and extract list of genes 
  intersectBed -a $tss -b $esc_bed -wa | sort | uniq | cut -f 4 | sed 's/["\;]//g' |sed  's/\..*//g' > ./gene_lists/ESC.tss_in_BG4peaks.txt
  intersectBed -a $tss -b $esc_bed -wa | sort | uniq | cut -f 4 | sed 's/["\;]//g' > ./gene_lists/ESC.noEnsVer.tss_in_BG4peaks.txt
  intersectBed -a $tss -b $nsc_bed -wa | sort | uniq | cut -f 4 | sed 's/["\;]//g' |sed  's/\..*//g' > ./gene_lists/NSC.tss_in_BG4peaks.txt
  intersectBed -a $tss -b $nsc_bed -wa | sort | uniq | cut -f 4 | sed 's/["\;]//g' > ./gene_lists/NSC.noEnsVer.tss_in_BG4peaks.txt
  intersectBed -a $tss -b $cncc_bed -wa | sort | uniq | cut -f 4 | sed 's/["\;]//g' |sed  's/\..*//g' > ./gene_lists/CNCC.tss_in_BG4peaks.txt
  intersectBed -a $tss -b $cncc_bed -wa | sort | uniq | cut -f 4 | sed 's/["\;]//g' > ./gene_lists/CNCC.noEnsVer.tss_in_BG4peaks.txt
done
```

run topGO gene enrichment analysis
==================================

``` r
require(DOSE)
require(clusterProfiler)
source('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/GOenrichment_topgo_2lists.R')
#cut -f 4 gencode.v28.TSS.bed |sed 's/["\;]//g' |sed  's/\..*//g' | sort | uniq > gencode.v28.TSS.geneEns.txt
all_genes <- read.table('/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS.geneEns.txt',stringsAsFactors = F)
path_list_files <- "/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/gene_lists/"
list_ens_files <- list.files(path_list_files,pattern = 'BG4peaks.txt')

for (i in 1:length(list_ens_files)) {
  print(list_ens_files[i])
  list_ensembl <- read.table(paste0(path_list_files,list_ens_files[i]), stringsAsFactors = F)
  GO_BP<- GOenrichment_topgo(list_ensembl$V1)
  write.table(GO_BP$filteredRed,file=paste0(path_list_files,gsub('.txt','.topGO.BP.txt',list_ens_files[i])), col.names = NA, sep='\t',quote = F)
  write.table(GO_BP$allRes_BP,file=paste0(path_list_files,gsub('.txt','.all.res.topGO.BP.txt',list_ens_files[i])), col.names = NA, sep='\t',quote = F)
  david<- enrichDAVID(gene = list_ensembl$V1, idType="ENSEMBL_GENE_ID", annotation="GOTERM_BP_FAT",
                                               david.user="angela.simeone@cruk.cam.ac.uk")
  pdf(paste0(path_list_files,gsub('.txt','.all.res.DAVID.bp_fat.pdf',list_ens_files[i])))
    barplot(david)
  dev.off()
  rm(list_ensembl)
  
}

for (i in 1:length(list_ens_files)) {
  print(list_ens_files[i])
  rm(david,list_ensembl)
  list_ensembl <- read.table(paste0(path_list_files,list_ens_files[i]), stringsAsFactors = F)
  david<- enrichDAVID(gene = list_ensembl$V1, idType="ENSEMBL_GENE_ID", annotation="GOTERM_BP_FAT",
                                               david.user="angela.simeone@cruk.cam.ac.uk")
  pdf(paste0(path_list_files,gsub('.txt','.all.res.DAVID.bp_fat.pdf',list_ens_files[i])), 20,10)
      barplot(david,showCategory=50)
  dev.off()
}

# path_list_files <- "/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/coverage_peaks/"
# list_ens_files <- list.files(path_list_files,pattern = 'TSS.ENS.txt')
path_list_files <- "/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/gene_lists/"
list_ens_files <- list.files(path_list_files,pattern = 'BG4peaks.txt')

for (i in 1:length(list_ens_files)) {
  print(list_ens_files[i])
  rm(david, david_all,list_ensembl)
  list_ensembl <- read.table(paste0(path_list_files,list_ens_files[i]), stringsAsFactors = F)
  #david<- enrichDAVID(gene = list_ensembl$V1, idType="ENSEMBL_GENE_ID", annotation="GOTERM_BP_FAT",david.user="angela.simeone@cruk.cam.ac.uk")
  david_all<- enrichDAVID(gene = list_ensembl$V1, idType="ENSEMBL_GENE_ID", annotation="GOTERM_BP_ALL",david.user="angela.simeone@cruk.cam.ac.uk")
  # pdf(paste0(path_list_files,gsub('.txt','.all.res.DAVID.bp_fat.pdf',list_ens_files[i])), 20,10)
  #     bb <- barplot(david,showCategory=50)
  #     print(bb)
  # dev.off()
  
  pdf(paste0(path_list_files,gsub('.txt','.all.res.DAVID.bp_all.pdf',list_ens_files[i])), 20,10)
      bb <- barplot(david_all,showCategory=50)
      print(bb)
  dev.off()
}
```

analysis on specific categories generated before
------------------------------------------------

``` bash

cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2
out_dir=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/selected_function_categories
mkdir $out_dir
for ref in `ls /scratcha/sblab/simeon01/reference_genomes/biomart/*geneSymbl_promoter.bed`
do
  for file in *bio2out3.bed
  do
    tmp=${ref%%_geneSymbl_promoter.bed}
    file_n=${tmp##/scratcha/sblab/simeon01/reference_genomes/biomart/}
    echo $file_n
    sbatch -o %j.tmp.out -e %j.tmp.err --mem 2000 --wrap "intersectBed -b $ref -a $file -wa | sort | uniq > $out_dir/${file%%.bio2out3.bed}_G4confirmedPeak.${file_n}.bed"
  done
done

promoter=/scratcha/sblab/simeon01/reference_genomes/hg38/TSS_1kb_around_strand_genes_UCSC_hg38_igenomes.merged.filtered.bed
for file in *bio2out3.bed
do
#save peaks overlapping any promoter
intersectBed -b $promoter -a $file -wa | sort | uniq > $out_dir/${file%%.bed}_G4confirmedPeak.promoters.bed
done

Bival_Promoters_ESC=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/Additional_annotations/Rada_Iglesias_lifted_hg38_bivalentPromoters_curated.bed
Bival_Promoters_CNCC=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/Additional_annotations/CNCC_hg38_bivalentPromoters.bed
Bival_Promoters_NSC=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/Additional_annotations/NSC_hg38_bivalentPromoters.bed

intersectBed -b $Bival_Promoters_ESC -a ESC.bio2out3.bed -wa | sort | uniq > $out_dir/ESC.bio2out3_G4confirmedPeak.ESC_bival_prom.bed
intersectBed -b $Bival_Promoters_NSC -a NSC.bio2out3.bed -wa | sort | uniq > $out_dir/NSC.bio2out3_G4confirmedPeak.NSC_bival_prom.bed
intersectBed -b $Bival_Promoters_CNCC -a CNCC.bio2out3.bed -wa | sort | uniq > $out_dir/CNCC.bio2out3_G4confirmedPeak.CNCC_bival_prom.bed
```

check peak sizes for various functional categories
--------------------------------------------------

``` r
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/selected_function_categories')
list_peak_files <- list.files(".",'*.bed')

myfiles2 = lapply(list_peak_files, read.table)
names(myfiles2) <- list_peak_files

Size <- c()
flag <- c()
for (i in 1:length(myfiles2)){
  temp_size <- myfiles2[[i]]$V3-myfiles2[[i]]$V2
  case_flag <- rep.int(names(myfiles2[i]),length(myfiles2[[i]]$V1))
  Size <- c(Size,temp_size)
  flag <- c(flag,case_flag)
}

G4_peaks <- data.frame(size=Size[order(flag)],cond=flag[order(flag)], stringsAsFactors = F)

library(ggplot2)

p <- ggplot(G4_peaks, aes(x=flag, y=size)) + 
  geom_boxplot()
# Rotate the box plot
p + coord_flip()
sts <- boxplot.stats(G4_peaks$size)$stats
p1 = p + coord_cartesian(ylim = c(sts*1.05,sts/1.05))
pdf('boxplot_peaks_sizes_various_functional_groups.pdf',12,10)
print(p1+ coord_flip())
dev.off()
pdf('boxplot_peaks_sizes_various_functional_groups_zoom_in.pdf',12,10)
p1+ coord_flip() +ylim(0,400)
dev.off()

G4_peaks_diff <- G4_peaks %>% filter(stringr::str_detect(flag, 'differentiation'))
p_G4_peaks_diff <- ggplot(G4_peaks_diff, aes(x=cond, y=size)) + 
  geom_boxplot() 
p_G4_peaks_diff + coord_flip()
sts <- boxplot.stats(G4_peaks_diff$size)$stats
p1_G4_peaks_diff = p_G4_peaks_diff + coord_cartesian(ylim = c(sts*1.05,sts/1.05))
pdf('G4_peaks_size_by_cell_type.differentiation.zoom_in.pdf')
p1_G4_peaks_diff+ coord_flip()+ ylim(0,300)
dev.off()

G4_peaks_promoters <- G4_peaks %>% filter(stringr::str_detect(flag, 'promoters'))
G4_peaks_biv <- G4_peaks %>% filter(stringr::str_detect(flag, 'bival'))
p_G4_peaks_biv <- ggplot(G4_peaks_biv, aes(x=cond, y=size)) + 
  geom_boxplot() 
p_G4_peaks_biv + coord_flip()
sts <- boxplot.stats(G4_peaks_biv$size)$stats
p1_G4_peaks_biv = p_G4_peaks_biv + coord_cartesian(ylim = c(sts*1.05,sts/1.05))
pdf('G4_peaks_size_by_cell_type.bival.zoom_in.pdf')
p1_G4_peaks_biv + coord_flip()+ ylim(0,300)
dev.off()

G4_peaks_cellCyce <- G4_peaks %>% filter(stringr::str_detect(flag, 'cellCycle'))
p_G4_peaks_cellCycle <- ggplot(G4_peaks_cellCyce, aes(x=cond, y=size)) + 
  geom_boxplot() 
p_G4_peaks_cellCycle + coord_flip()
sts <- boxplot.stats(G4_peaks_cellCyce$size)$stats
p1_G4_peaks_cellCycle = p_G4_peaks_cellCycle + coord_cartesian(ylim = c(sts*1.05,sts/1.05))
p1_G4_peaks_cellCycle + coord_flip()+ ylim(0,300)


G4_peaks_splicing <- G4_peaks %>% filter(stringr::str_detect(flag, 'splicing'))
p_G4_peaks_splicing <- ggplot(G4_peaks_splicing, aes(x=cond, y=size)) + 
  geom_boxplot() 
p_G4_peaks_splicing + coord_flip()
sts <- boxplot.stats(G4_peaks_splicing$size)$stats
p1_G4_peaks_splicing = p_G4_peaks_splicing + coord_cartesian(ylim = c(sts*1.05,sts/1.05))
p1_G4_peaks_splicing + coord_flip() + ylim(0,300)

pdf('N_G4_peaks_overlapping_bival_prom.zoom_in.pdf')
G4_peaks_biv %>% ggplot(aes(cond)) + geom_bar(colour="black") + coord_flip()
dev.off()

pdf('N_G4_peaks_overlapping_promoters.zoom_in.pdf')
G4_peaks_promoters %>% ggplot(aes(cond)) + geom_bar(colour="black") + coord_flip()
dev.off()
```

Additional external data to download and process (peak detection) and fold enrichment analysis (chromatin architecture elements)
--------------------------------------------------------------------------------------------------------------------------------

``` bash
mkdir /mnt/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/GEO_series_GSE105028

#IgG
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816613/suppl/GSM2816613%5FIgG%5FNT%2Ebw

#RAD21
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE105nnn/GSE105028/suppl/GSE105028%5FRad21%5FNT%2Ebw

# atac
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE105nnn/GSE105028/suppl/GSE105028%5FATAC%2DSeq%2Ebw

# H3K27ac
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE105nnn/GSE105028/suppl/GSE105028%5FH3K27ac%5FNT%2Ebw

#H3K27me3
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816652/suppl/GSM2816652%5FH3K27me3%5FNT%2Ebw

#H3K4me1
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816654/suppl/GSM2816654%5FH3K4me1%5FNT%2Ebw

#H3K4me3
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816656/suppl/GSM2816656%5FH3K4me3%5FNT%2Ebw

#OCT4
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816629/suppl/GSM2816629%5FOCT4%5FNT%2Ebw

#NANOG
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816625/suppl/GSM2816625%5FNANOG%5FNT%2Ebw

#KLF4
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816627/suppl/GSM2816627%5FKLF4%5FNT%2Ebw

#ZNF143
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816621/suppl/GSM2816621%5FZNF143%5FNT%2Ebw

#pol2
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816623/suppl/GSM2816623%5FPolII%5FNT%2Ebw

#RING1B
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816631/suppl/GSM2816631%5FRING1B%5FNT%2Ebw

#cFos
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816633/suppl/GSM2816633%5FcFos%5FNT%2Ebw

#Brg1
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816635/suppl/GSM2816635%5FBrg1%5FNT%2Ebw

#HSF1
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816637/suppl/GSM2816637%5FHSF1%5FNT%2Ebw

#WAPL
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816640/suppl/GSM2816640%5FWAPL%5FNT%2Ebw

#NIPBL
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2816nnn/GSM2816642/suppl/GSM2816642%5FNIPBL%5FNT%2Ebw

## convert bigwig to bedgraph
cd  /mnt/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/GEO_series_GSE105028

for file in *bw
do 
sbatch --mem 12G --wrap "/Users/simeon01/applications/bigWigToBedGraph $file ${file%%.bw}.bedgraph"
done

cd /mnt/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/GEO_series_GSE105028
mkdir macs2_output
for file in *bedgraph
do
  sbatch --mem 12G --wrap "macs2 bdgpeakcall -i $file -c 0.7 -l 300 -g 100 -o ./macs2_output/${file%%.bedgraph}_narrowPeaks.bed"
done


# reerun macs2 peaks calls with different thredhols for histone modifications - transcription factors - architecture proteeins
 cd /mnt/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/GEO_series_GSE105028
for file in *H3K*bedgraph
do
  sbatch --mem 12G --wrap "macs2 bdgpeakcall -i $file -c 0.5 -l 300 -g 100 -o ./macs2_output/${file%%.bedgraph}_narrowPeaks.bed"
done

for file in `ls *bedgraph | grep -v H3K`
do
  sbatch --mem 12G --wrap "macs2 bdgpeakcall -i $file -c 0.6 -l 300 -g 100 -o ./macs2_output/${file%%.bedgraph}_narrowPeaks.bed"
done


# for bed in 150 proximity merge bed files
for file in *.bed
do
  sbatch --mem 4G --wrap "sortBed -i $file | mergeBed -i - -d 150 | sortBed -i - > ${file%%.bed}.merged_d150.bed"
done

# for all merged files create annotations bed files (intervals and name of the annotation)

temp_annotation=chromatin_architecture_elements.annotation.bed
touch $temp_annotation
for file in *merged_d150.bed
do
label_temp=${file%%_[nN]*}
echo $label_temp
label=${label_temp#GS*_}
echo $label
echo " == "
awk -v a="$label" '{print $1"\t"$2"\t"$3"\t"a}' $file >> $temp_annotation
done
# note:
# GSM2816625 is NANOG
# GSM2816642 is NIPBL

temp_annotation_nexus=chromatin_architecture_elements.ChIPNexus.annotation.bed
touch $temp_annotation
#/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/GEO_series_GSE105028/ChIPNexus
awk '{print $0"\tCTCF"}' GSM2816644_CTCF_ChIPNexus_peaks.bed >>$temp_annotation_nexus
awk '{print $0"\tRAD21"}' GSM2816646_Rad21_ChIPNexus_peaks.bed >>$temp_annotation_nexus


# run gat analysis on this newly generated annotation file containing the chromatin architecture elements


cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
katie_path=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed
scenario=(scenario1_multi2 scenario2_5out9 scenario3_multi2_extended_1kb)
annotation_path=/mnt/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/GEO_series_GSE105028/macs2_output
annotation_chrom_architecture=(chromatin_architecture_elements.annotation.bed)
for case in ${scenario[@]}
do
  cd $katie_path/$case
  pwd
  mkdir Chromatine_architecture_enrichments_$case
  # run gat analysis on full bed files 
  for file in *bed
  do
    echo $file
    for anno in ${annotation_chrom_architecture[@]}
    do
      echo $anno
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$annotation_path/${anno} --segments=${file} -S ./Chromatine_architecture_enrichments_$case/${file%%.bed}_${anno%%.bed}.txt -t 4"
      echo ${file%%.bed}_${anno%%.bed}.txt
      echo $cmd_gat
      sbatch --mem 6G --wrap "$cmd_gat"
      echo " "
      
    #sbatch --mem 6G --wrap "$cmd_gat"
    done
  done
  #zip -r CpG_fold_enrichments_$case.zip CpG_fold_enrichments_$case
done

## ChIP-nexus
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
katie_path=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed
scenario=(scenario1_multi2 scenario2_5out9 scenario3_multi2_extended_1kb)
annotation_path=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/GEO_series_GSE105028/ChIPNexus
annotation_chrom_architecture=(chromatin_architecture_elements.ChIPNexus.annotation.bed)
for case in ${scenario[@]}
do
  cd $katie_path/$case
  pwd
  mkdir Chromatine_architecture_enrichments_$case
  # run gat analysis on full bed files 
  for file in *bed
  do
    echo $file
    for anno in ${annotation_chrom_architecture[@]}
    do
      echo $anno
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$annotation_path/${anno} --segments=${file} -S ./Chromatine_architecture_enrichments_$case/${file%%.bed}_${anno%%.bed}.txt -t 4"
      echo ${file%%.bed}_${anno%%.bed}.txt
      echo $cmd_gat
      sbatch --mem 6G --wrap "$cmd_gat"
      echo " "
      
    #sbatch --mem 6G --wrap "$cmd_gat"
    done
  done
  #zip -r CpG_fold_enrichments_$case.zip CpG_fold_enrichments_$case
done


zip -r Chromatine_architecture_enrichments_scenario1_multi2.zip Chromatine_architecture_enrichments_scenario1_multi2

bedtools makewindows -g hg38_selected.sorted.genome -w 500 > hg38_selected.500bp_windows.bed
```

Check fold enrichment of BG4 peaks at CpG methylation sites
-----------------------------------------------------------

``` bash

tail -n +2 pnas.1613300114.sd02_hESC_eexported.csv | head -n 2 > pnas.1613300114.sd02_hESC_eexported_rearranged.csv

cut -d "," -f 1,2,3,18 pnas.1613300114.sd02_hESC_eexported_rearranged.csv | sed 's/\,/\t/g' | tail -n +2 > pnas.1613300114.sd02_hESC_eexported_rearranged.bed

#lifover to hg38
chain_hg19_to_hg38=/scratcha/sblab/simeon01/reference_genomes/liftover_chain/hg19ToHg38.over.chain.gz
/Users/simeon01/applications/liftOver pnas.1613300114.sd02_hESC_eexported_rearranged.bed $chain_hg19_to_hg38 pnas.1613300114.sd02_hESC_eexported_rearranged.hg38.bed  pnas.1613300114.sd02_hESC_eexported_rearranged.unmapped

#CpG island only
awk '{print $1"\t"$2"\t"$3"\tCpG"}' hg38_GRC38_ucsc_CpG.bed > hg38_GRC38_ucsc_CpG_annotations.bed
```

Run enrichment over all methylation sites

``` bash

cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
katie_path=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed
path_anno='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells'
scenario=(scenario1_multi2 scenario2_5out9 scenario3_multi2_extended_1kb)
annotation_path=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/CpG_island
annotation_CpG=(pnas.1613300114.sd02_hESC_eexported_rearranged.hg38.bed hg38_GRC38_ucsc_CpG_annotations.bed)
for case in ${scenario[@]}
do
  cd $katie_path/$case
  pwd
  mkdir CpG_fold_enrichments_$case
  # run gat analysis on full bed files 
  for file in *bed
  do
    echo $file
    for anno in ${annotation_CpG[@]}
    do
      echo $anno
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$annotation_path/${anno} --segments=${file} -S ./CpG_fold_enrichments_$case/${file%%.bed}_${anno%%.bed}.txt -t 4"
      echo ${file%%.bed}_${anno%%.bed}.txt
      echo $cmd_gat
      sbatch --mem 6G --wrap "$cmd_gat"
      echo " "
      
    #sbatch --mem 6G --wrap "$cmd_gat"
    done
  done
  #zip -r CpG_fold_enrichments_$case.zip CpG_fold_enrichments_$case
done

for case in ${scenario[@]}
do
  cd $katie_path/$case
  zip -r CpG_fold_enrichments_$case.zip CpG_fold_enrichments_$case
done


## check total number of overlaps



path_anno='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells'
scenario=(scenario1_multi2 scenario2_5out9 scenario3_multi2_extended_1kb)
annotation_path=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/CpG_island
annotation_CpG=(pnas.1613300114.sd02_hESC_eexported_rearranged.hg38.annotations.bed hg38_GRC38_ucsc_CpG_annotations.bed)

for case in ${scenario[@]}
do
  cd $katie_path/$case
  echo "** ==================================================== ** "
  pwd
  for file in *bed
  do
    wc -l $file
    for anno in ${annotation_CpG[@]}
    do
      echo $anno
      intersectBed -a $file -b $annotation_path/$anno -wa | sort | uniq | wc -l
      echo "==="
    done
  done
done
```

BG4 coverage at promoters
-------------------------

``` bash

promoters=/scratcha/sblab/simeon01/reference_genomes/hg38/gencode.v28.TSS_plus_1kb.bed
Bam_repA=/scratcha/sblab/simeon01/Data/190829_Katie_stemcells/SLX-18496/aligned #.hg38.sort.markduplicates.bam
Bam_repB=/scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4/aligned/merged # markduplicates.bam
Bam_repC=/scratcha/sblab/simeon01/Data/20190917_SLX-18355_Katie_StemCells_Rep2_3xBG4/SLX-18355/aligned #.hg38.sort.markduplicates.bam

#./peaks_coverages
#./promoter_coverages

#loop over scenario in order to run coverages -  for the moment I run it only for scenario1
bam_folders=($Bam_repA $Bam_repB $Bam_repC)
scenario=(scenario1_multi2)
path_analysis=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
for case in ${scenario[@]}
  do
  cd $path_analysis/$case
  
  out_dir=$path_analysis/$case/bed_merge_all/promoter_coverage_peaks
  mkdir $out_dir
  echo $bed_all_peaks
  for folder_bam in ${bam_folders[@]}
    do
    cd $folder_bam
    pwd
    for bam in `ls *.markduplicates.bam`
      do
      echo $bam
      pwd
      echo sbatch --mem 8G -o %j.log.out -e %j.log.err --wrap "bedtools coverage -a $promoters -b $bam  -counts > $out_dir/${bam%%.bam}.bg4_at_promoters.bedgraph"
      sbatch --mem 8G -o %j.log.out -e %j.log.err --wrap "bedtools coverage -a $promoters -b $bam  -counts > $out_dir/${bam%%.bam}.bg4_at_promoters.bedgraph"
      echo "   "
      done
  
    done
  
  done
```

Differential binding analysis at promoter
-----------------------------------------

``` r
## == fuctions  and packages ==
library(dendextend)
library(dplyr)
library(ggdendro)
library(ggplot2)
library(biomaRt)
library(fgsea)


pathways.hallmark <- gmtPathways("/Users/simeon01/Documents/genomes/hg38/c2.cp.v7.0.entrez.gmt.txt")
go_anno <- gmtPathways("/Users/simeon01/Documents/genomes/hg38/c5.all.v7.0.entrez.gmt.txt")
short_sequence_motifs <- gmtPathways("/Users/simeon01/Documents/genomes/hg38/c3.all.v7.0.entrez.gmt.txt")
oncogenic_signatures <- gmtPathways("/Users/simeon01/Documents/genomes/hg38/c6.all.v7.0.entrez.gmt.txt")
immunologic_signatures <- gmtPathways("/Users/simeon01/Documents/genomes/hg38/c7.all.v7.0.entrez.gmt.txt")

tss_file <- '/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS.bed'
source('/Users/simeon01/Documents/Karen/20190613_Karen_continue_pol2_G4_maps/G4_signal_DRB/function_to_export_tables.R')
source('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/promoter_coverage_peaks/GOenrichment_topgo_all.R')
norm_sign_to_tot_map_reads <- function(matrix,list_stat5_mapped_reads){
  Library_sizes <- c()
  matrix_norm <- matrix
  for (i in 1:dim(matrix)[2]){
    print(i)
    tmp_cond <- colnames(matrix)[i]
    tmp_cond <- gsub('.union_bg4_peaks','',tmp_cond)
    print(tmp_cond)
    tmp_stat <- grep(substr(tmp_cond,1,20),list_stat5_mapped_reads,value = T)
    tot_reads <- as.numeric(system(paste0('cat ',tmp_stat),intern = TRUE))
    Library_sizes[i] <- tot_reads
    matrix_norm[,i] <- matrix[,i]/tot_reads*1e6
  }
  names(Library_sizes) <- colnames(matrix)
  NewList <- list("norm_cpm_M" = matrix_norm, "lib_size" = Library_sizes)
  return(NewList)
}

#################################################
# load data
#################################################
setwd("/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/promoter_coverage_peaks")
list_files <- list.files(path = ".",pattern = 'bg4_at_promoters.bedgraph')

Matrix_cov <- matrix(ncol = length(list_files), nrow = 58381)
for (i in 1:length(list_files)) {
  temp <- read.table(list_files[i],stringsAsFactors = F)
  Matrix_cov[,i] <- temp$V8
  peaks_coordinates <- temp[,c(1,2,3)]
}
colnames(Matrix_cov) <-list_files
rownames(Matrix_cov) <-temp$V4
colnames(Matrix_cov) <- gsub('.markduplicates','',list_files)


# specific samples need to be removed
colnames(Matrix_cov)[]

#Matrix_cov_selected <- Matrix_cov[,-c(33,37,39,40)]

setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/coverage_peaks')
stat5_files <- list.files(path = ".", pattern = "stat5")
#
norm_Matrix_cov_selected <- norm_sign_to_tot_map_reads(Matrix_cov,stat5_files)
setwd("/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/promoter_coverage_peaks")

## == clustering raw data == 

# # clustering - row data
pdf('hclust_raw_counts.promoter.pdf')
h_Matrix_cov <- hclust(dist(t(Matrix_cov)),method="ward.D2")
p <- ggdendrogram(h_Matrix_cov,rotate=T) +  ggtitle('clustering - raw counts')
print(p)
dev.off()

## == clustering - data norm to total number of reads == 
pdf('hclust_RPM.promoter.pdf')
h_norm_Matrix_cov_selected <- hclust(dist(t(norm_Matrix_cov_selected$norm_cpm_M)),method="ward.D2")
p <- ggdendrogram(h_norm_Matrix_cov_selected,rotate=T) +  ggtitle('clustering - RPM')
print(p)  
dev.off()

## == clustering - data norm to total number of reads -- NO INPUT == 
pdf('hclust_RPM_no_input.promoter.pdf')
h_norm_Matrix_cov_selected_noInput <- hclust(dist(t(norm_Matrix_cov_selected$norm_cpm_M[,-grep('input',colnames(norm_Matrix_cov_selected$norm_cpm_M))])),method="ward.D2")
p <- ggdendrogram(h_norm_Matrix_cov_selected_noInput,rotate=T) +  ggtitle('clustering - RPM - no input')
print(p)  
dev.off()

## == clustering - data norm to total number of reads -- NO INPUT == 
pdf('hclust_RPM_no_input_excluding_low10percent.promoter.pdf')
A <- norm_Matrix_cov_selected$norm_cpm_M[,-grep('input',colnames(norm_Matrix_cov_selected$norm_cpm_M))]
A_sum <- rowSums(A)
hist(A_sum)
summary(A_sum)
quantile(A_sum, 0.1)
A_selected <- A[A_sum<quantile(A_sum, 0.1),]
h_A_selected <- hclust(dist(t(A_selected)),method="ward.D2")
p <- ggdendrogram(h_A_selected,rotate=T) +  ggtitle('clustering - RPM - no input')
plot(p)
dev.off()

pdf('hclust_RPM_no_input_excluding_NOlow10percent_NOTop10percent.pdf')
A <- norm_Matrix_cov_selected$norm_cpm_M[,-grep('input',colnames(norm_Matrix_cov_selected$norm_cpm_M))]
A_sum <- rowSums(A)
hist(A_sum)
summary(A_sum)
quantile(A_sum, 0.1)
A_selected <- A[A_sum>quantile(A_sum, 0.1) & A_sum<quantile(A_sum, 0.9) ,]
#A_selected <- A[A_sum>quantile(A_sum, 0.5) ,]
h_A_selected <- hclust(dist(t(A_selected)),method="ward.D2")
p <- ggdendrogram(h_A_selected,rotate=T) +  ggtitle('clustering - RPM - no input')
plot(p)
dev.off()

## == visualize as heatmap == 
#set.seed(123)
#Heatmap(mat)
#Heatmap(A[A_sum>quantile(A_sum, 0.1) & A_sum<quantile(A_sum, 0.9) ,],show_row_names =F)

#################################
# Differential analysis
#################################

# in the differential analysis we will account for:
# bio rep [4 different ones]
# technical rep [3 tech rep]
# condition [3 different cell lines]

# create variables
#A_selected <- A[A_sum<quantile(A_sum, 0.1),]
#Counts_to_use_all <- A_selected
colnames(Matrix_cov)[-c(c(33,37),grep('input',colnames(Matrix_cov)))]
Matrix_cov_selected <- Matrix_cov[,-c(c(33,37),grep('input',colnames(Matrix_cov)))]
Counts_to_use_all <- Matrix_cov_selected
colnames(Counts_to_use_all)
conditions <- c(rep('CNCC',9),rep('ESC',9),rep('NSC',9))
conditions
# "CNCC" "CNCC" "CNCC" "CNCC" "CNCC" "CNCC" "CNCC" "CNCC" "CNCC" "ESC"  "ESC"  "ESC"  "ESC"  "ESC"  "ESC"  "ESC"  "ESC"  "ESC"  "NSC"  "NSC"  "NSC"  "NSC" 
# [23] "NSC"  "NSC"  "NSC"  "NSC"  "NSC" 
#      
conditions <- factor(conditions) # there are 4 factors

conditions <- relevel(conditions,ref="ESC") # relevel as reference base the condition no_TPL_1hr_no_dm6


bio_rep <- rep(c(rep(1,3), rep(2,3), rep(3,3)),3)
bio_rep <- factor(bio_rep)
bio_rep
#[1] 1 1 1 2 2 2 3 3 3 1 1 1 2 2 2 3 3 3 1 1 1 2 2 2 3 3 3
#Levels: 1 2 3


tec_rep <- rep(1:3,9)
tec_rep <- factor(tec_rep)
# lib_size: norm_Matrix_cov_selected$lib_size[-grep('input',colnames(norm_Matrix_cov_selected$norm_cpm_M))]
curr_lib_size <- norm_Matrix_cov_selected$lib_size[-c(c(33,37),grep('input',colnames(Matrix_cov)))]
design_blocking <- model.matrix(~ bio_rep + tec_rep + conditions)
design <- model.matrix(~conditions)
colnames(design_blocking)
# [1]  "(Intercept)"    "bio_rep2"       "bio_rep3"       "tec_rep2"       "tec_rep3"       "conditionsCNCC" "conditionsNSC" 

# create EdgeR object
library(edgeR)
y <- DGEList(counts = Counts_to_use_all, group= conditions)
stopifnot(rownames(y$samples) == names(curr_lib_size))
y$samples$lib.size<- curr_lib_size

y <- calcNormFactors(y, method= 'none')
y <- estimateDisp(y,design_blocking)
y <- estimateGLMCommonDisp(y,design_blocking)
y <- estimateGLMTagwiseDisp(y,design_blocking)


# check MDS plot
pdf('MDS_edger_mdsplot_stem_cells.promoter.pdf', 7, 7)
par(las= 1, cex=1.2)
col_map <- conditions
levels(col_map) <- 1:length(levels(col_map))
col_map <- as.numeric(conditions)
temp_rainbow <- rainbow(3)
col_map_rainbow <- temp_rainbow[col_map]
plotMDS(y,
        #pch=c(0,3,0,3,0,3),
        #xlim=c(-3,3), ylim= c(-0.5,0.5), 
        xlim=c(-1.5,1.5), ylim= c(-1,0.7), 
        #xaxp= c(-40, 40, 8),  yaxp= c(-10,10, 8),
        labels= paste(y$samples$group,bio_rep, sep= '--'),
        col=col_map_rainbow,
        main= 'MDS plot ', cex= 0.7)
grid()
dev.off()

pdf('MDS_edger_mdsplot_stem_cells_symb.promoter.pdf', 7, 7)
 par(las= 1, cex=1.2)
col_map <- conditions
levels(col_map) <- 1:length(levels(col_map))
col_map <- as.numeric(conditions)
temp_rainbow <- rainbow(3)
col_map_rainbow <- temp_rainbow[col_map]
plotMDS(y,
        pch=c(rep(0,9),rep(1,9),rep(3,9)),
        #xlim=c(-3,3), ylim= c(-0.5,0.5), 
        xlim=c(-1.5,1.5), ylim= c(-1,0.7), 
        #xaxp= c(-40, 40, 8),  yaxp= c(-10,10, 8),
        #labels= paste(y$samples$group,bio_rep, sep= '--'),
        col=col_map_rainbow,
        main= 'MDS plot ', cex= 0.7)

grid()
dev.off()

pdf('PCA2_edger_mdsplot_stem_cells.promoter.pdf', 9, 7)
pcaResult<-prcomp(t(y$counts))
plot(pcaResult$x,
     main= 'Principal components, counts',
     xlab="PC1",
     ylab="PC2",
     #xlab= sprintf('PC1 (sd: %s%%)', round(100 * (pcaResult$sdev[1] / sum(pcaResult$sdev)))),
     #ylab= sprintf('PC2 (sd: %s%%)', round(100 * (pcaResult$sdev[2] / sum(pcaResult$sdev)))),
     type= 'n')
text(x= pcaResult$x[,1], y= pcaResult$x[,2], labels= paste(y$samples$group,bio_rep, sep= '--'), cex= 0.9, 
     #col= c(rep('#CCFF00FF',5),rep('#00FF66FF',5),rep('#FF0000FF',5)))
     col = col_map_rainbow)
dev.off()


# fit the generalized linear model (quasi-likelihood negative binomial generalized log-linear model)
fit <- glmQLFit(y, design_blocking)
colnames(fit)
source('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/promoter_coverage_peaks/perform_annotation_and_enrichment_analysis.R')
# "(Intercept)"    "bio_rep2"       "bio_rep3"       "tec_rep2"       "tec_rep3"       "conditionsCNCC" "conditionsNSC" 

## ==  CNCC vs ESC 
N_condition <- 6

colnames(fit)[N_condition]

label_case <-  paste0(colnames(fit)[N_condition],"_vs_ESC")

qlf.CNCC_vs_ESC <- glmQLFTest(fit, coef=N_condition) 

de_pairing_CNCC_vs_ESC <- topTags(qlf.CNCC_vs_ESC, n= nrow(y$counts))$table

de_pairing_CNCC_vs_ESC$ens_id <- gsub('\\..*','',rownames(de_pairing_CNCC_vs_ESC))

write.table(de_pairing_CNCC_vs_ESC,file = paste0('EdgeRtable_differential_binding_analysis_',label_case,'.csv'), col.names = NA,sep = ",",quote = F)


res_CNCC_vs_ESC01 <- decideTestsDGE(qlf.CNCC_vs_ESC,adjust.method="BH", p.value = 0.1,lfc=0)
summary(res_CNCC_vs_ESC01)

res_CNCC_vs_ESC <- decideTestsDGE(qlf.CNCC_vs_ESC,adjust.method="BH", p.value = 0.05,lfc=0)

summary(res_CNCC_vs_ESC)

Table_GSEA_CNCC_vs_ESC <- data.frame( ensembl_gene_id= gsub('\\..*','',rownames(de_pairing_CNCC_vs_ESC)),
                                      metric=-log10(de_pairing_CNCC_vs_ESC$PValue)/sign(de_pairing_CNCC_vs_ESC$logFC))


#write.table(Table_GSEA_CNCC_vs_ESC,file='Table_GSEA.CNCC_vs_ESC.ESCspecific.rnk',quote = F, sep = "\t",col.names = F,row.names = F)

Res_fgsea_CNCC_vs_ESC <- perform_annotation_and_enrichment_analysis(Table_GSEA_CNCC_vs_ESC)

explort_fsga_tables_to_txt(label_case,Res_fgsea_CNCC_vs_ESC)
# plotEnrichment(go_anno[["GO_EPITHELIUM_DEVELOPMENT"]],
#                 Res_fgsea_CNCC_vs_ESC$ranks) + labs(title="GO_EPITHELIUM_DEVELOPMENT")

pdf(paste0('MD_MA_plot_differential_binding_analysis_',label_case,'.pdf'),8,7)
plotMD(qlf.CNCC_vs_ESC,status = res_CNCC_vs_ESC,values =c(1,0,-1),col=c("red","black","blue"))
title(paste0("\n\nDEup:",length(which(res_CNCC_vs_ESC==1)),
             " - DEdown: ",length(which(res_CNCC_vs_ESC==-1)),
             " - notDE: ",length(which(res_CNCC_vs_ESC==0))))
dev.off()

## ==  NSC vs ESC 
N_condition <- 7

colnames(fit)[N_condition]

label_case <-  paste0(colnames(fit)[N_condition],"_vs_ESC")

qlf.NSC_vs_ESC <- glmQLFTest(fit, coef=N_condition) 

de_pairing_NSC_vs_ESC <- topTags(qlf.NSC_vs_ESC, n= nrow(y$counts))$table

write.table(de_pairing_NSC_vs_ESC,file = paste0('EdgeRtable_differential_binding_analysis_',label_case,'.csv'), col.names = NA,sep = ",",quote = F)

res_NSC_vs_ESC <- decideTestsDGE(qlf.NSC_vs_ESC,adjust.method="BH", p.value = 0.05,lfc=0)

summary(res_NSC_vs_ESC)

Table_GSEA_NSC_vs_ESC <- data.frame( ensembl_gene_id= gsub('\\..*','',rownames(de_pairing_NSC_vs_ESC)),
                                      metric=-log10(de_pairing_NSC_vs_ESC$PValue)/sign(de_pairing_NSC_vs_ESC$logFC))


#write.table(Table_GSEA_CNCC_vs_ESC,file='Table_GSEA.CNCC_vs_ESC.ESCspecific.rnk',quote = F, sep = "\t",col.names = F,row.names = F)

Res_fgsea_NSC_vs_ESC <- perform_annotation_and_enrichment_analysis(Table_GSEA_NSC_vs_ESC)

explort_fsga_tables_to_txt(label_case,Res_fgsea_NSC_vs_ESC)
# plotEnrichment(go_anno[["GO_EPITHELIUM_DEVELOPMENT"]],
#                 Res_fgsea_CNCC_vs_ESC$ranks) + labs(title="GO_EPITHELIUM_DEVELOPMENT")

#tss_res_NSC_vs_ESC <- export_bed_and_TSS_in_bed_regions(tss_file,'./',label_case,res_NSC_vs_ESC)

pdf(paste0('MD_MA_plot_differential_binding_analysis_',label_case,'.pdf'),8,7)
plotMD(qlf.NSC_vs_ESC,status = res_NSC_vs_ESC,values =c(1,0,-1),col=c("red","black","blue"))
title(paste0("\n\nDEup:",length(which(res_NSC_vs_ESC==1)),
             " - DEdown: ",length(which(res_NSC_vs_ESC==-1)),
             " - notDE: ",length(which(res_NSC_vs_ESC==0))))
dev.off()

## == case NSC vs CNCC
label_case <-  "NSC_vs_CNCC"

qlf.NSC_vs_CNCC <- glmQLFTest(fit, contrast=c(0,0,0,0,0,-1,1))

de_pairing_NSC_vs_CNCC <- topTags(qlf.NSC_vs_CNCC, n= nrow(y$counts))$table

write.table(de_pairing_NSC_vs_CNCC,file = paste0('EdgeRtable_differential_binding_analysis_',label_case,'.csv'), col.names = NA,sep = ",",quote = F)

res_NSC_vs_CNCC <- decideTestsDGE(qlf.NSC_vs_CNCC,adjust.method="BH", p.value = 0.05,lfc=0)

summary(res_NSC_vs_CNCC)

Table_GSEA_NSC_vs_CNCC <- data.frame( ensembl_gene_id= gsub('\\..*','',rownames(de_pairing_NSC_vs_CNCC)),
                                      metric=-log10(de_pairing_NSC_vs_CNCC$PValue)/sign(de_pairing_NSC_vs_CNCC$logFC))


#write.table(Table_GSEA_CNCC_vs_ESC,file='Table_GSEA.CNCC_vs_ESC.ESCspecific.rnk',quote = F, sep = "\t",col.names = F,row.names = F)

Res_fgsea_NSC_vs_CNCC <- perform_annotation_and_enrichment_analysis(Table_GSEA_NSC_vs_CNCC)

explort_fsga_tables_to_txt(label_case,Res_fgsea_NSC_vs_CNCC)
#tss_res_NSC_vs_CNCC <- export_bed_and_TSS_in_bed_regions(tss_file,'./',label_case,res_NSC_vs_CNCC)

pdf(paste0('MD_MA_plot_differential_binding_analysis_',label_case,'.pdf'),8,7)
plotMD(qlf.NSC_vs_CNCC,status = res_NSC_vs_CNCC,values =c(1,0,-1),col=c("red","black","blue"))
title(paste0("\n\nDEup:",length(which(res_NSC_vs_CNCC==1)),
             " - DEdown: ",length(which(res_NSC_vs_CNCC==-1)),
             " - notDE: ",length(which(res_NSC_vs_CNCC==0))))
dev.off()


## ----  save a summary of the findings ---- 
DBA_outcome_Nregions <- cbind(summary(res_NSC_vs_ESC),
                              summary(res_CNCC_vs_ESC),
                             summary(res_NSC_vs_CNCC))

write.table(DBA_outcome_Nregions,file='stem_cells_DBA_outcome_Npromoters.csv',col.names = NA,quote = F, sep = ",")


## import now digital info about RNA-seq and about expression levels (DE  lists and actual average TPM in cell line) -- cluster data. 
## == import promoter bed file
TSS_coordinates <-read.table('/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS.bed', stringsAsFactors = F)

## == import epigentic features == 
path_anno <- '/Users/simeon01/Documents/Katie/Annotations_Cells_Angela/promoter_based_anno'
setwd(path_anno)
# import promoter info: 1- promoter overlap broad peak
faire_ESC <- read.delim('ESC_promoter_annot_H9ESC_faireseq.final.bed', stringsAsFactors = F,header = T,row.names=NULL)
head(faire_ESC)
colnames(faire_ESC)[7] <- 'ens'

epi_ESC <- read.delim('ESC_promoter_annot_H9ESC_rada_igles.final.bed', stringsAsFactors = F,header = T,row.names=NULL)
head(epi_ESC)
colnames(epi_ESC)[7] <- 'ens'

atac_CNCC <- read.delim('CNCC_promoter_annot_atac.final.bed', stringsAsFactors = F,header = T,row.names=NULL)
head(atac_CNCC)
colnames(atac_CNCC)[7] <- 'ens'

epi_CNCC <- read.delim('CNCC_promoter_annot_encodeH1derived.final.bed', stringsAsFactors = F,header = T,row.names=NULL)
head(epi_CNCC)
colnames(epi_CNCC)[7] <- 'ens'

epi_NEC <- read.delim('NEC_promoter_annot.bed', stringsAsFactors = F,header = T,row.names=NULL)
head(epi_NEC)
colnames(epi_NEC)[7] <- 'ens'

epi_NSC<- read.delim('NSC_promoter_annot_encodeH1derived.final.bed', sep = "\t",stringsAsFactors = F,header = T,row.names=NULL)
head(epi_NSC)
colnames(epi_NSC)[7] <- 'ens'


# combine together all annotations (binary)
ESC_anno_complete <- left_join(faire_ESC,epi_ESC,by = 'ens') %>% dplyr::select(c(ens,FAIRE_H9ESC,contains("hg38")))
NSC_anno_complete <- epi_NSC %>% dplyr::select(ens,contains("NSC"))
CNCC_anno_complete <- left_join(atac_CNCC,epi_CNCC, by='ens') %>% dplyr::select(ens,contains("CNCC"))


# import transcriptomic data
TPM_file <- '/Users/simeon01/Documents/Katie/2018_09_03/RNA_seq/TPM_hg38_all_bio_rep_meanTechReps.csv'
out_path <- '/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/promoter_coverage_peaks'

# import output of differential expression analysis keeping attention to the fact that the labels could be reverted
# DEA_outcome_NSC_vs_ESC <- read.table('/Users/simeon01/Documents/Katie/2018_09_03/RNA_seq/summary_list_decide_test_H9_NSC_vs_H9_hESC.csv',sep = ',') # these labels are reverted!!! 1: ESC, -1 NSC
# DEA_outcome_NSC_vs_CNCC <- read.table('/Users/simeon01/Documents/Katie/2018_09_03/RNA_seq/summary_list_decide_test_H9_NSC_vs_H9_CNCC.csv',sep = ',')
# DEA_outcome_CNCC_vs_ESC <- read.table('/Users/simeon01/Documents/Katie/2018_09_03/RNA_seq/summary_list_decide_test_H9_CNCC_vs_H9_ESC.csv',sep = ',')
DEA_outcome_NSC_vs_ESC <- read.table('/Users/simeon01/Documents/Katie/2018_09_03/RNA_seq/updated2_summary_list_decide_test_H9_NSC_vs_H9_hESC.csv',sep = ',') # 
DEA_outcome_NSC_vs_CNCC <- read.table('/Users/simeon01/Documents/Katie/2018_09_03/RNA_seq/updated2_summary_list_decide_test_H9_NSC_vs_H9_CNCC.csv',sep = ',')
DEA_outcome_CNCC_vs_ESC <- read.table('/Users/simeon01/Documents/Katie/2018_09_03/RNA_seq/updated2_summary_list_decide_test_H9_CNCC_vs_H9_ESC.csv',sep = ',')



colnames(DEA_outcome_NSC_vs_ESC) <- c('ens_id','DEA_NSC_vs_ESC')
colnames(DEA_outcome_NSC_vs_CNCC) <- c('ens_id','DEA_NSC_vs_CNCC') 
colnames(DEA_outcome_CNCC_vs_ESC) <- c('ens_id','DEA_CNCC_vs_ESC')

# additional filtering on the DEA genes
# remove row 1
DEA_outcome_NSC_vs_ESC <- DEA_outcome_NSC_vs_ESC[-1,]
DEA_outcome_NSC_vs_CNCC <- DEA_outcome_NSC_vs_CNCC[-1,]
DEA_outcome_CNCC_vs_ESC <- DEA_outcome_CNCC_vs_ESC[-1,]
# remove gene symbols
DEA_outcome_NSC_vs_ESC$ens_id <- gsub('\\..*','',DEA_outcome_NSC_vs_ESC$ens_id)
DEA_outcome_NSC_vs_CNCC$ens_id <- gsub('\\..*','',DEA_outcome_NSC_vs_CNCC$ens_id)
DEA_outcome_CNCC_vs_ESC$ens_id <- gsub('\\..*','',DEA_outcome_CNCC_vs_ESC$ens_id)

# select only those different from 0
DEA_outcome_NSC_vs_ESC <- DEA_outcome_NSC_vs_ESC[which(DEA_outcome_NSC_vs_ESC$DEA_NSC_vs_ESC!=0),]
DEA_outcome_NSC_vs_CNCC <- DEA_outcome_NSC_vs_CNCC[which(DEA_outcome_NSC_vs_CNCC$DEA_NSC_vs_CNCC!=0),]
DEA_outcome_CNCC_vs_ESC <- DEA_outcome_CNCC_vs_ESC[which(DEA_outcome_CNCC_vs_ESC$DEA_CNCC_vs_ESC!=0),]

TPM_rnaseq <- read.delim(TPM_file, stringsAsFactors = F,header = T,row.names=NULL, sep = ",")
colnames(TPM_rnaseq)[c(1,2,3)] <- c('ens_symb','ens','symb')
head(TPM_rnaseq)
setwd(out_path)

TPM_median_expression <-  TPM_rnaseq %>% group_by(ens) %>% mutate(median_TPM_ESC = median(H9_hESC_NA_p10,H9_hESC_NA_p11,H9_hESC_NA_p12,H9_hESC_NA_p3,H9_hESC_NA_p2,na.rm = T)) %>% mutate(median_TPM_NSC = median(H9_NSC_NA_p4.2,H9_NSC_NA_p4.3,H9_NSC_NA_p4.4,H9_NSC_NA_p4.1,na.rm =T)) %>% mutate(median_TPM_NEC = median(H9_neurosphere_NA_p5,H9_neurosphere_NA_p6,H9_neurosphere_NA_p7,H9_neurosphere_NA_p8,H9_neurosphere_NA_p9,na.rm =T)) %>% mutate(median_TPM_CNCC = median(H9_CNCC_A1_p15,H9_CNCC_A2_p16,H9_CNCC_B4_p17,H9_CNCC_A3_p18,H9_CNCC_B3_p19,H9_CNCC_4_p20,H9_CNCC_A1_p21,H9_CNCC_A5_p22,na.rm =T))

TPM_median_expression_standard_0_1_temp <- as.matrix(sapply(TPM_median_expression[,4:34], as.numeric))
TPM_median_expression_standard_0_1 <- TPM_median_expression_standard_0_1_temp
for (i in (1:dim(TPM_median_expression_standard_0_1_temp)[2])) {
  temp <- TPM_median_expression_standard_0_1_temp[,i]
  TPM_median_expression_standard_0_1[,i]<- (temp-min(temp))/(max(temp)-min(temp))
}

TPM_median_expression_standard_0_1 <- as.data.frame(TPM_median_expression_standard_0_1)
TPM_median_expression_standard_0_1$ens <- TPM_median_expression$ens

res_CNCC_vs_ESC_df <- data.frame(ens = rownames(res_CNCC_vs_ESC), DE_CNCC_vs_ESC = unname(res_CNCC_vs_ESC))
res_NSC_vs_ESC_df <- data.frame(ens = rownames(res_NSC_vs_ESC), DE_NSC_vs_ESC = unname(res_NSC_vs_ESC))
res_NSC_vs_CNCC_df <- data.frame(ens = rownames(res_NSC_vs_CNCC), DE_NSC_vs_CNCC = unname(res_NSC_vs_CNCC))

res_CNCC_vs_ESC_df$ens_id <- gsub('\\..*','',res_CNCC_vs_ESC_df$ens)
res_NSC_vs_ESC_df$ens_id <- gsub('\\..*','',res_NSC_vs_ESC_df$ens)
res_NSC_vs_CNCC_df$ens_id <- gsub('\\..*','',res_NSC_vs_CNCC_df$ens)


## ==  comparison between DBA and DEA outcome == 
# create a single summary matrix
DA_DBA_NSC_vs_ESC <- left_join(DEA_outcome_NSC_vs_ESC,res_NSC_vs_ESC_df,by="ens_id")
DA_DBA_CNCC_vs_ESC <- left_join(DEA_outcome_CNCC_vs_ESC,res_CNCC_vs_ESC_df,by="ens_id")
DA_DBA_NSC_vs_CNCC <- left_join(DEA_outcome_NSC_vs_CNCC,res_NSC_vs_CNCC_df,by="ens_id")

DA_DBA_NSC_vs_ESC <-DA_DBA_NSC_vs_ESC %>% mutate(combined=paste0(DEA_NSC_vs_ESC,DE_NSC_vs_ESC))
DA_DBA_CNCC_vs_ESC <- DA_DBA_CNCC_vs_ESC %>% mutate(combined=paste0(DEA_CNCC_vs_ESC,DE_CNCC_vs_ESC))
DA_DBA_NSC_vs_CNCC <- DA_DBA_NSC_vs_CNCC %>% mutate(combined=paste0(DEA_NSC_vs_CNCC,DE_NSC_vs_CNCC))

# write to file results of pairwise analysis
write.table(DA_DBA_NSC_vs_ESC,file='Gene_tables_DBA_DE_NSC_vs_ESC.txt',sep= "\t",quote = F,col.names = T,row.names = F)
write.table(DA_DBA_CNCC_vs_ESC,file='Gene_tables_DBA_DE_CNCC_vs_ESC.txt',sep= "\t",quote = F,col.names = T,row.names = F)
write.table(DA_DBA_NSC_vs_CNCC,file='Gene_tables_DBA_DE_NSC_vs_CNCC.txt',sep= "\t",quote = F,col.names = T,row.names = F)

# what is the relationthip between differential expression and differential binding
table(paste0(DA_DBA_NSC_vs_ESC$DEA_NSC_vs_ESC,DA_DBA_NSC_vs_ESC$DE_NSC_vs_ESC))
# report number of occurrances for each of the possible cases:
#1) -1 -1 : differentially expressed in denominator, differentially bound in denominator
#2) -1 0 : differentially expressed in denominator, not differentially bound
#3) -1 1 : differentially expressed in denominator, differentially bound in numerator
#4) 1 -1 : differentially expressed in numerator, differentially bound in denominator
#5) 1 0 : differentially expressed in numerator, not differentially bound
#6)  1 1 : differentially expressed in numerator, differentially bound in numerator

write.table(table(paste0(DA_DBA_NSC_vs_ESC$DEA_NSC_vs_ESC,DA_DBA_NSC_vs_ESC$DE_NSC_vs_ESC)),file = 'PieChart_NSC_vs_ESC_counts.txt', col.names = NA,sep = "\t",quote = F)
table(paste0(DA_DBA_NSC_vs_ESC$DEA_NSC_vs_ESC,DA_DBA_NSC_vs_ESC$DE_NSC_vs_ESC))/dim(DA_DBA_NSC_vs_ESC)[1] *100
#      -1-1        -10        -11        1-1         10         11 
# 55.8161351 11.4446529  0.1407129 19.1838649 12.1013133  1.3133208 
pdf('PieChart_NSC_ESC.pdf')
pie_NSC_ESC <- as.data.frame(table(paste0(DA_DBA_NSC_vs_ESC$DEA_NSC_vs_ESC,DA_DBA_NSC_vs_ESC$DE_NSC_vs_ESC))/dim(DA_DBA_NSC_vs_ESC)[1] *100)
colnames(pie_NSC_ESC)[1] <- "class"
ggplot(pie_NSC_ESC, aes(x = "", y = Freq, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
   coord_polar("y", start = 0)+geom_text(aes(label = round(Freq,1)), position = position_stack(vjust = 0.5)) +
  ggtitle("NSC_vs_ESC") + 
  theme_void() 
dev.off()

write.table(table(paste0(DA_DBA_CNCC_vs_ESC$DEA_CNCC_vs_ESC,DA_DBA_CNCC_vs_ESC$DE_CNCC_vs_ESC)),file = 'PieChart_CNCC_vs_ESC_counts.txt', col.names = NA,sep = "\t",quote = F)
table(paste0(DA_DBA_CNCC_vs_ESC$DEA_CNCC_vs_ESC,DA_DBA_CNCC_vs_ESC$DE_CNCC_vs_ESC))/dim(DA_DBA_CNCC_vs_ESC)[1] *100
#      -1-1        -10        -11        1-1         10         11 
# 52.5563258 18.1542461  0.3466205 10.3119584 16.8977470  1.7331023 
pdf('PieChart_CNCC_ESC.pdf')
pie_CNCC_vs_ESC <- as.data.frame(table(paste0(DA_DBA_CNCC_vs_ESC$DEA_CNCC_vs_ESC,DA_DBA_CNCC_vs_ESC$DE_CNCC_vs_ESC))/dim(DA_DBA_CNCC_vs_ESC)[1] *100)
colnames(pie_CNCC_vs_ESC)[1] <- "class"
ggplot(pie_CNCC_vs_ESC, aes(x = "", y = Freq, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+geom_text(aes(label = round(Freq,1)), position = position_stack(vjust = 0.5)) +
  ggtitle("CNCC_vs_ESC") + 
  theme_void() 
dev.off()

write.table(table(paste0(DA_DBA_NSC_vs_CNCC$DEA_NSC_vs_CNCC,DA_DBA_NSC_vs_CNCC$DE_NSC_vs_CNCC)),file = 'PieChart_NSC_CNCC_counts.txt', col.names= NA,sep = "\t",quote = F)
table(paste0(DA_DBA_NSC_vs_CNCC$DEA_NSC_vs_CNCC,DA_DBA_NSC_vs_CNCC$DE_NSC_vs_CNCC))/dim(DA_DBA_NSC_vs_CNCC)[1] *100
#      -1-1        -10        -11        1-1         10         11 
# 10.7586207 16.7356322  0.6896552 12.1379310 57.7931034  1.8850575 
pdf('PieChart_NSC_CNCC.pdf')
pie_NSC_vs_CNCC <- as.data.frame(table(paste0(DA_DBA_NSC_vs_CNCC$DEA_NSC_vs_CNCC,DA_DBA_NSC_vs_CNCC$DE_NSC_vs_CNCC))/dim(DA_DBA_NSC_vs_CNCC)[1] *100)
colnames(pie_NSC_vs_CNCC)[1] <- "class"
ggplot(pie_NSC_vs_CNCC, aes(x = "", y = Freq, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+geom_text(aes(label = round(Freq)), position = position_stack(vjust = 0.5)) +
  ggtitle("NSC_vs_CNCC") + 
  theme_void() 
dev.off()

# check enrichment for each of the different classes :
classes_of_cases <- unique(DA_DBA_NSC_vs_CNCC$combined)
library(clusterProfiler)
for (i in 1:length(classes_of_cases)) {
  curr_case <- classes_of_cases[i]
  print(length(which(DA_DBA_NSC_vs_ESC$combined==curr_case)))

  curr_ensembl <- DA_DBA_NSC_vs_ESC %>% dplyr::group_by(ens_id) %>% dplyr::filter(combined==curr_case) %>% dplyr::select(ens_id)
  GO_BP<- GOenrichment_topgo_all(curr_ensembl$ens_id)
  # david <- enrichDAVID(gene = curr_ensembl$ens_id,
  #                      universe = DA_DBA_NSC_vs_ESC$ens_id,
  #                      idType="ENSEMBL_GENE_ID", annotation="GOTERM_BP_ALL", pvalueCutoff = 1.1,
  #                      pAdjustMethod = "bonferroni",
  #                      minGSSize = 1,david.user="angela.simeone@cruk.cam.ac.uk")
  # barplot(david)
  write.table(GO_BP$allRes_BP,file = paste0('Combined_DA_DBA_NSC_vs_ESC_topGO_',as.character(curr_case),'.txt'),sep = "\t",col.names = NA)

}

## == combine DBA, DEA and epigenetic annotations
# combine DBA_DEA_ESC_anno_NSC_anno
combined_DBA_DEA_ESC_anno <- left_join(DA_DBA_NSC_vs_ESC,ESC_anno_complete,by="ens")
combined_DBA_DEA_NSC_anno <- left_join(DA_DBA_NSC_vs_ESC,NSC_anno_complete ,by="ens")
combined_DBA_DEA_ESC_anno_NSC_anno <- left_join(combined_DBA_DEA_ESC_anno,combined_DBA_DEA_NSC_anno,by="ens")

DBA_DEA_ESC_refined <- combined_DBA_DEA_ESC_anno %>% group_by(ens) %>% dplyr::mutate(ESC_K4me3_binary = if_else(hg38.ESC_H3K4me3_calls >0 ,1,0)) %>% dplyr::mutate(ESC_K27me3_binary = if_else(hg38.ESC_H3K27me3_calls >0 ,1,0)) %>% dplyr::mutate(ESC_bivalent_binary = if_else(hg38.ESC_H3K4me3_calls_overlap_H3K27me3_calls.AngelaBivalent >0 ,1,0))
# force to zero those rows that have the binary all equal to one otherwise the frequencies get altered
DBA_DEA_ESC_refined <- DBA_DEA_ESC_refined %>% group_by(ens) %>% mutate(ESC_bivalent_binary=replace(ESC_bivalent_binary,ESC_K27me3_binary==1 &ESC_K4me3_binary ==1 ,1)) %>%mutate(ESC_K27me3_binary=replace(ESC_K27me3_binary,ESC_bivalent_binary==1,0)) %>% mutate(ESC_K4me3_binary=replace(ESC_K4me3_binary,ESC_bivalent_binary==1,0)) #%>% dplyr::filter(ESC_bivalent_binary==1) %>% dplyr::select(combined,ESC_K27me3_binary,ESC_K4me3_binary,ESC_bivalent_binary)

DBA_DEA_NSC_refined <- combined_DBA_DEA_NSC_anno %>% group_by(ens) %>% dplyr::mutate(NSC_K4me3_binary = if_else(NSC_H3K4me3 >0 ,1,0)) %>% dplyr::mutate(NSC_K27me3_binary = if_else(NSC_H3K27me3 >0 ,1,0)) %>% dplyr::mutate(NSC_bivalent_binary = if_else(NSC_H3K4me3_overlap_H3K27me3 >0 ,1,0))

# force to zero those rows that have the binary all equal to one otherwise the frequencies get altered
DBA_DEA_NSC_refined <- DBA_DEA_NSC_refined %>% group_by(ens) %>% mutate(NSC_bivalent_binary=replace(NSC_bivalent_binary,NSC_K27me3_binary==1 &NSC_K4me3_binary ==1 ,1)) %>% mutate(NSC_K27me3_binary=replace(NSC_K27me3_binary,NSC_bivalent_binary==1,0)) %>% mutate(NSC_K4me3_binary=replace(NSC_K4me3_binary,NSC_bivalent_binary==1,0)) #%>% dplyr::filter(ESC_bivalent_binary==1) %>% dplyr::select(combined,ESC_K27me3_binary,ESC_K4me3_binary,ESC_bivalent_



#count frequencies grouping by the different cases
DBA_DEA_ESC_refined_frequencies <- DBA_DEA_ESC_refined %>%  group_by(combined,ESC_K27me3_binary,ESC_K4me3_binary,ESC_bivalent_binary) %>% summarise(count = n(), proportion = count / nrow(.))
write.table(DBA_DEA_ESC_refined_frequencies,file='combined_DBA_DEA_ESC_anno_counts_and_freq.txt',sep = "\t",quote = F, col.names = NA)

#count frequencies grouping by the different cases

DBA_DEA_NSC_refined_frequencies <- DBA_DEA_NSC_refined %>%  group_by(combined,NSC_K27me3_binary,NSC_K4me3_binary,NSC_bivalent_binary) %>% summarise(count = n(), proportion = count / nrow(.))
write.table(DBA_DEA_NSC_refined_frequencies,file='combined_DBA_DEA_NSC_anno_counts_and_freq.txt',sep = "\t",quote = F, col.names = NA)

for (i in 1:length(unique(unique(combined_DBA_DEA_ESC_anno$combined)))) {
  print(i)
  spec_case <- unique(combined_DBA_DEA_ESC_anno$combined)[i]
  pdf(paste0('combined_DBA_DEA_ESC_anno',spec_case,'.pdf'),8,8)
  panel <- DBA_DEA_ESC_refined_frequencies %>% dplyr::filter(combined==spec_case) %>% ggplot(aes(x=paste(ESC_K4me3_binary,ESC_K27me3_binary,ESC_bivalent_binary,sep = "=="),y=proportion)) + geom_bar(stat="identity") + ggtitle(spec_case)
  print(panel)
dev.off()
}


for (i in 1:length(unique(unique(combined_DBA_DEA_NSC_anno$combined)))) {
  print(i)
  spec_case <- unique(combined_DBA_DEA_NSC_anno$combined)[i]
  pdf(paste0('combined_DBA_DEA_NSC_anno',spec_case,'.pdf'),8,8)
  panel <- DBA_DEA_NSC_refined_frequencies %>% dplyr::filter(combined==spec_case) %>% ggplot(aes(x=paste(NSC_K4me3_binary,NSC_K27me3_binary,NSC_bivalent_binary,sep = "=="),y=proportion)) + geom_bar(stat="identity", fill='azure4') + ggtitle(spec_case)
  print(panel)
dev.off()
  
}


## === === ===
## ==import BG4 rpm (signal normalized to total library size)
norm_G4 <- norm_Matrix_cov_selected$norm_cpm_M[,-grep('input',colnames(norm_Matrix_cov_selected$norm_cpm_M))]

norm_G4_standardized_0_1<- norm_G4
for (i in (1:dim(norm_G4)[2])) {
  temp <- norm_G4[,i]
  norm_G4_standardized_0_1[,i]<- (temp-min(temp))/(max(temp)-min(temp))
}

bg4_promoter_values_df <- as.data.frame(cbind(rownames(norm_G4_standardized_0_1),norm_G4_standardized_0_1))
colnames(bg4_promoter_values_df)[1] <- 'ens'
colnames(bg4_promoter_values_df) <- gsub('_SLX*.*','',colnames(bg4_promoter_values_df))
colnames(bg4_promoter_values_df) <- gsub('.all.*.*','',colnames(bg4_promoter_values_df))
colnames(bg4_promoter_values_df) <- gsub('.hg38.sort*.*','',colnames(bg4_promoter_values_df))

merged_NSC_vs_ESC <- left_join(res_NSC_vs_ESC_df,TPM_median_expression,by="ens")
merged_CNCC_vs_ESC <- left_join(res_CNCC_vs_ESC_df,TPM_median_expression,by="ens")
merged_NSC_vs_CNCC <- left_join(res_NSC_vs_CNCC_df,TPM_median_expression,by="ens")

merged_NSC_vs_ESC_bg4 <- left_join(merged_NSC_vs_ESC,bg4_promoter_values_df,by="ens")



columns_to_select <- c(grep('DE_NSC_vs_ESC',colnames(merged_NSC_vs_ESC_bg4)),
                       grep('ESC_Rep',colnames(merged_NSC_vs_ESC_bg4)),
                       grep('H9_hES',colnames(merged_NSC_vs_ESC_bg4)),
                       grep('NSC_Rep',colnames(merged_NSC_vs_ESC_bg4)),
                       grep('H9_NSC',colnames(merged_NSC_vs_ESC_bg4)))
columns_to_select <- columns_to_select[-c(22,25)] #  remove NSC_Rep3_3xBG4_ChIP2 and NSC_Rep3b_3xBG

M_temp <- as.matrix(sapply(merged_NSC_vs_ESC_bg4[,columns_to_select], as.numeric))
rownames(M_temp) <- merged_NSC_vs_ESC_bg4$ens
ind_sorted_DBA_groups <- order(as.numeric((M_temp[,1])))
M_temp <- M_temp[ind_sorted_DBA_groups,]

M_temp[M_temp>=50]<- 50
#M_temp <- scale(M_temp)
library(ComplexHeatmap)
hp_list_pam <- Heatmap(M_temp[,-1],
                       name = "promoter_BG4_expression",
                       cluster_columns = FALSE,
                       cluster_rows =F ,
                       show_row_names = FALSE) +
  Heatmap(M_temp[,1],col = c('blue','black','red'),name = "DBA_group", width = unit(5, "mm"),show_row_names = FALSE)


pdf('global_heatmap_G4_individual_rep.pdf',width = 10, height = 12)
draw(hp_list_pam)
dev.off()

# cluster within each group
M_temp_minus1 <- M_temp[which(M_temp[,1]== -1),-1]
M_temp_0 <- M_temp[which(M_temp[,1]== 0),-1]
M_temp_plus1 <- M_temp[which(M_temp[,1]== -1),-1]

# for each group PAM cluster them 
library("cluster")
set.seed(1)
#pamx_minus1 <- pam(M_temp_minus1, 4)
#pamx_0 <- pam(M_temp_0, 4)
#pamx_plus1 <- pam(M_temp_plus1, 4)

#sort data for visuzalization
mat_after_PAM_minus1 <-rbind(M_temp_minus1[rownames(M_temp_minus1)%in%names(which(pamx_minus1$clustering==1)),],
                      M_temp_minus1[rownames(M_temp_minus1)%in%names(which(pamx_minus1$clustering==2)),],
                      M_temp_minus1[rownames(M_temp_minus1)%in%names(which(pamx_minus1$clustering==3)),],
                      M_temp_minus1[rownames(M_temp_minus1)%in%names(which(pamx_minus1$clustering==4)),])
cluster_M_after_PAM_minus1 <- c(rep(1,length(which(pamx_minus1$clustering==1))),
                         rep(2,length(which(pamx_minus1$clustering==2))),
                         rep(3,length(which(pamx_minus1$clustering==3))),
                         rep(4,length(which(pamx_minus1$clustering==4))))

hp_list_pam <- Heatmap(mat_after_PAM_minus1,name = "promoter", 
                       cluster_columns = FALSE,cluster_rows =F,
                       show_row_names = FALSE, use_raster = TRUE) +
  Heatmap(cluster_M_after_PAM_minus1,col = rainbow(4),name = "cluster", width = unit(5, "mm"),show_row_names = FALSE)
pdf('test_mat_PAM_minus1.pdf',width = 6,height = 10)
draw(hp_list_pam)
dev.off()

#sort data for visuzalization
mat_after_PAM_plus1 <-rbind(M_temp_plus1[rownames(M_temp_plus1)%in%names(which(pamx_plus1$clustering==1)),],
                      M_temp_plus1[rownames(M_temp_plus1)%in%names(which(pamx_plus1$clustering==2)),],
                      M_temp_plus1[rownames(M_temp_plus1)%in%names(which(pamx_plus1$clustering==3)),],
                      M_temp_plus1[rownames(M_temp_plus1)%in%names(which(pamx_plus1$clustering==4)),])
cluster_M_after_PAM_plus1 <- c(rep(1,length(which(pamx_plus1$clustering==1))),
                         rep(2,length(which(pamx_plus1$clustering==2))),
                         rep(3,length(which(pamx_plus1$clustering==3))),
                         rep(4,length(which(pamx_plus1$clustering==4))))

hp_list_pam_plus1<- Heatmap(mat_after_PAM_plus1,name = "promoter", 
                       cluster_columns = FALSE,cluster_rows =F,
                       show_row_names = FALSE, use_raster = TRUE) +
  Heatmap(cluster_M_after_PAM_plus1,col = rainbow(4),name = "cluster", width = unit(5, "mm"),show_row_names = FALSE)
pdf('test_mat_PAM_plus1.pdf',width = 6,height = 10)
draw(hp_list_pam_plus1)
dev.off()

# cluster entire matrix without stratification by DBA outcome



mat_after_PAM_0 <-rbind(M_temp_0[rownames(M_temp_0)%in%names(which(pamx_0$clustering==1)),],
                      M_temp_0[rownames(M_temp_0)%in%names(which(pamx_0$clustering==2)),],
                      M_temp_0[rownames(M_temp_0)%in%names(which(pamx_0$clustering==3)),],
                      M_temp_0[rownames(M_temp_0)%in%names(which(pamx_0$clustering==4)),])
cluster_M_after_PAM_0 <- c(rep(1,length(which(pamx_0$clustering==1))),
                         rep(2,length(which(pamx_0$clustering==2))),
                         rep(3,length(which(pamx_0$clustering==3))),
                         rep(4,length(which(pamx_0$clustering==4))))

hp_list_pam_0 <- Heatmap(mat_after_PAM_0,name = "promoter", 
                       cluster_columns = FALSE,cluster_rows =F,
                       show_row_names = FALSE, use_raster = TRUE) +
  Heatmap(cluster_M_after_PAM_0,col = rainbow(4),name = "cluster", width = unit(5, "mm"),show_row_names = FALSE)
pdf('test_mat_PAM_0.pdf',width = 6,height = 10)
draw(hp_list_pam_0)
dev.off()

#pamx_fullMatrix <- pam(M_temp,4)




## == explore relationthip between BG4 static data and expression

# load gene id for genes that have a peak in promoter 
ESC_BG4peak_in_promoter <- read.table('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/gene_lists/ESC.noEnsVer.tss_in_BG4peaks.txt',stringsAsFactors = F)
NSC_BG4peak_in_promoter <- read.table('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/gene_lists/NSC.noEnsVer.tss_in_BG4peaks.txt',stringsAsFactors = F)
CNCC_BG4peak_in_promoter <- read.table('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/gene_lists/CNCC.noEnsVer.tss_in_BG4peaks.txt',stringsAsFactors = F)

ESC_BG4peak_in_promoter_df <- data.frame(ens=ESC_BG4peak_in_promoter$V1,
                                         BG4_ESC=rep(1,length(ESC_BG4peak_in_promoter$V1)))
NSC_BG4peak_in_promoter_df <- data.frame(ens=NSC_BG4peak_in_promoter$V1,
                                         BG4_NSC=rep(1,length(NSC_BG4peak_in_promoter$V1)))
CNCC_BG4peak_in_promoter_df <- data.frame(ens=CNCC_BG4peak_in_promoter$V1,
                                         BG4_CNCC=rep(1,length(CNCC_BG4peak_in_promoter$V1)))



# merge this lists with expressions
TPM_BG4peak_in_promoter <- TPM_median_expression %>% dplyr::select(ens,contains("median_"))%>% left_join(ESC_BG4peak_in_promoter_df,by="ens") %>% left_join(NSC_BG4peak_in_promoter_df,by="ens") %>% left_join(CNCC_BG4peak_in_promoter_df,by="ens")

#stratify gene expression by not_detec/low/med/high expression levels
ranges_expression_ESC <- c(1,quantile(TPM_BG4peak_in_promoter$median_TPM_ESC[TPM_BG4peak_in_promoter$median_TPM_ESC>1],c(0.3,0.6)))
ranges_expression_NSC <- c(1,quantile(TPM_BG4peak_in_promoter$median_TPM_NSC[TPM_BG4peak_in_promoter$median_TPM_NSC>1],c(0.3,0.6)))
ranges_expression_CNCC <- c(1,quantile(TPM_BG4peak_in_promoter$median_TPM_CNCC[TPM_BG4peak_in_promoter$median_TPM_CNCC>1],c(0.3,0.6)))

TPM_BG4peak_in_promoter <- TPM_BG4peak_in_promoter %>% mutate(ESC_strat_expression = if_else(median_TPM_ESC <ranges_expression_ESC[1]  | is.na(median_TPM_ESC),'NotExpr',if_else(median_TPM_ESC>ranges_expression_ESC[1] & median_TPM_ESC<ranges_expression_ESC[2],'Low',if_else(median_TPM_ESC>ranges_expression_ESC[2] & median_TPM_ESC<ranges_expression_ESC[3],'Med',if_else(median_TPM_ESC>ranges_expression_ESC[3],'High','NA')))))

TPM_BG4peak_in_promoter <- TPM_BG4peak_in_promoter %>% mutate(NSC_strat_expression = if_else(median_TPM_NSC <ranges_expression_NSC[1]  | is.na(median_TPM_NSC),'NotExpr',if_else(median_TPM_NSC>ranges_expression_NSC[1] & median_TPM_NSC<ranges_expression_NSC[2],'Low',if_else(median_TPM_NSC>ranges_expression_NSC[2] & median_TPM_NSC<ranges_expression_NSC[3],'Med',if_else(median_TPM_NSC>ranges_expression_NSC[3],'High','NA')))))

TPM_BG4peak_in_promoter <- TPM_BG4peak_in_promoter %>% mutate(CNCC_strat_expression = if_else(median_TPM_CNCC <ranges_expression_CNCC[1]  | is.na(median_TPM_CNCC),'NotExpr',if_else(median_TPM_CNCC>ranges_expression_CNCC[1] & median_TPM_CNCC<ranges_expression_CNCC[2],'Low',if_else(median_TPM_CNCC>ranges_expression_CNCC[2] & median_TPM_CNCC<ranges_expression_CNCC[3],'Med',if_else(median_TPM_CNCC>ranges_expression_CNCC[3],'High','NA')))))

TPM_BG4peak_in_promoter <- TPM_BG4peak_in_promoter %>% mutate(BG4_ESC=replace(BG4_ESC,is.na(BG4_ESC),'BG-')) %>% mutate(BG4_ESC=replace(BG4_ESC,BG4_ESC==1,'BG+'))
TPM_BG4peak_in_promoter <- TPM_BG4peak_in_promoter %>% mutate(BG4_NSC=replace(BG4_NSC,is.na(BG4_NSC),'BG-')) %>% mutate(BG4_NSC=replace(BG4_NSC,BG4_NSC==1,'BG+'))
TPM_BG4peak_in_promoter <- TPM_BG4peak_in_promoter %>% mutate(BG4_CNCC=replace(BG4_CNCC,is.na(BG4_CNCC),'BG-')) %>% mutate(BG4_CNCC=replace(BG4_CNCC,BG4_CNCC==1,'BG+'))




pdf('Boxplot_expression_levels_BG4Promoters_vs_notBG4.ESC.pdf')
counts_elements <- TPM_BG4peak_in_promoter %>% filter(median_TPM_ESC>1) %>% group_by(BG4_ESC) %>% dplyr::summarize(count = n())
cell_boxplot <- TPM_BG4peak_in_promoter %>% filter(median_TPM_ESC>1) %>% group_by(BG4_ESC) %>% ggplot(aes(y=log2(median_TPM_ESC),x=factor(BG4_ESC))) + geom_boxplot() + ylim(0,15) + ggtitle('ESC') + 
  scale_x_discrete(labels = c(paste0('G4- prom:',counts_elements$count[1]),paste0('G4+ prom:',counts_elements$count[2])))
print(cell_boxplot)
dev.off()

pdf('Boxplot_expression_levels_BG4Promoters_vs_notBG4.NSC.pdf')
counts_elements <- TPM_BG4peak_in_promoter %>% filter(median_TPM_NSC>1) %>% group_by(BG4_NSC) %>% dplyr::summarize(count = n())
cell_boxplot <- TPM_BG4peak_in_promoter %>% filter(median_TPM_NSC>1) %>% group_by(BG4_NSC) %>% ggplot(aes(y=log2(median_TPM_NSC),x=factor(BG4_NSC))) + geom_boxplot() + ylim(0,12) + ggtitle('NSC')+ 
  scale_x_discrete(labels = c(paste0('G4- prom:',counts_elements$count[1]),paste0('G4+ prom:',counts_elements$count[2])))
print(cell_boxplot)
dev.off()

pdf('Boxplot_expression_levels_BG4Promoters_vs_notBG4.CNCC.pdf')
counts_elements <- TPM_BG4peak_in_promoter %>% filter(median_TPM_CNCC>1) %>% group_by(BG4_CNCC) %>% dplyr::summarize(count = n())
counts_elements
cell_boxplot <- TPM_BG4peak_in_promoter %>% filter(median_TPM_CNCC>1) %>% group_by(BG4_CNCC) %>% ggplot(aes(y=log2(median_TPM_CNCC),x=factor(BG4_CNCC))) + geom_boxplot() + ylim(0,10) + ggtitle('CNCC')+ 
  scale_x_discrete(labels = c(paste0('G4- prom:',counts_elements$count[1]),paste0('G4+ prom:',counts_elements$count[2])))
print(cell_boxplot)
dev.off()


pdf('Boxplot_expression_BG4plus_mius_stratified_by_expression.ESC.pdf')
TPM_BG4peak_in_promoter %>%  ggplot(aes(y=log2(median_TPM_ESC),x=ESC_strat_expression,fill = BG4_ESC)) + geom_boxplot() + scale_x_discrete(limits=c("NotExpr", "Low", "Med","High"))+ ylim(0,7.5) +ggtitle('ESC - Expression of BG4+- promoters stratified by expression groups')
dev.off()
# willcoxon expression
cases_expression=c("Low", "Med","High")
wilcoxon_test_final_ESC <- tibble()
for (jj in (1:length(cases_expression))){
  wilcoxon_test <- TPM_BG4peak_in_promoter %>% filter(BG4_ESC=="BG+",median_TPM_ESC>1 &ESC_strat_expression ==cases_expression[jj]) %>% group_by(BG4_ESC) %>% summarise(p_value = wilcox.test(log2(median_TPM_ESC),log2(TPM_BG4peak_in_promoter$median_TPM_ESC[TPM_BG4peak_in_promoter$BG4_ESC=="BG-" &TPM_BG4peak_in_promoter$median_TPM_ESC>1]),alternative = "greater")$p.value)
  print(wilcoxon_test)
  wilcoxon_test$case=cases_expression[jj]
  wilcoxon_test_final_ESC <- rbind(wilcoxon_test_final_ESC,wilcoxon_test)
}
write.table(wilcoxon_test_final_ESC,'Boxplot_expression_BG4plus_mius_stratified_by_expression.ESC.stats.txt', quote = F)


pdf('Boxplot_expression_BG4plus_mius_stratified_by_expression.NSC.pdf')
TPM_BG4peak_in_promoter %>%  ggplot(aes(y=log2(median_TPM_NSC),x=NSC_strat_expression,fill = BG4_NSC)) + geom_boxplot() + scale_x_discrete(limits=c("NotExpr", "Low", "Med","High"))+ ylim(0,7) +ggtitle('NSC - Expression of BG4+- promoters stratified by expression groups')
dev.off()

wilcoxon_test_final_NSC <- tibble()
for (jj in (1:length(cases_expression))){
  wilcoxon_test <- TPM_BG4peak_in_promoter %>% filter(BG4_NSC=="BG+",median_TPM_NSC>1 &NSC_strat_expression ==cases_expression[jj]) %>% group_by(BG4_NSC) %>% summarise(p_value = wilcox.test(log2(median_TPM_NSC),log2(TPM_BG4peak_in_promoter$median_TPM_NSC[TPM_BG4peak_in_promoter$BG4_NSC=="BG-" &TPM_BG4peak_in_promoter$median_TPM_NSC>1]),alternative = "greater")$p.value)
  print(wilcoxon_test)
  wilcoxon_test$case=cases_expression[jj]
  wilcoxon_test_final_NSC <- rbind(wilcoxon_test_final_NSC,wilcoxon_test)
}
write.table(wilcoxon_test_final_NSC,'Boxplot_expression_BG4plus_mius_stratified_by_expression.NSC.stats.txt', quote=F)


pdf('Boxplot_expression_BG4plus_mius_stratified_by_expression.CNCC.pdf')
TPM_BG4peak_in_promoter %>%  ggplot(aes(y=median_TPM_CNCC,x=CNCC_strat_expression,fill = BG4_CNCC)) + geom_boxplot() + scale_x_discrete(limits=c("NotExpr", "Low", "Med","High"))+ ylim(0,30) +ggtitle('CNCC - Expression of BG4+- promoters stratified by expression groups')
dev.off()
wilcoxon_test_final_CNCC <- tibble()
for (jj in (1:length(cases_expression))){
  wilcoxon_test <- TPM_BG4peak_in_promoter %>% filter(BG4_CNCC=="BG+",median_TPM_CNCC>1 &ESC_strat_expression ==cases_expression[jj]) %>% group_by(BG4_CNCC) %>% summarise(p_value = wilcox.test(log2(median_TPM_CNCC),log2(TPM_BG4peak_in_promoter$median_TPM_CNCC[TPM_BG4peak_in_promoter$BG4_CNCC=="BG-" &TPM_BG4peak_in_promoter$median_TPM_CNCC>1]),alternative = "greater")$p.value)
  print(wilcoxon_test)
  wilcoxon_test$case=cases_expression[jj]
  wilcoxon_test_final_CNCC <- rbind(wilcoxon_test_final_CNCC,wilcoxon_test)
}
write.table(wilcoxon_test_final_CNCC,'Boxplot_expression_BG4plus_mius_stratified_by_expression.CNCC.stats.txt', quote = F)


TPM_BG4peak_in_promoter_frequencies_ESC <- TPM_BG4peak_in_promoter %>%  group_by(ESC_strat_expression,BG4_ESC) %>% summarise(count = n())
write.table(TPM_BG4peak_in_promoter_frequencies_ESC,file ='TPM_BG4peak_in_promoter_frequencies_ESC.data.txt',quote = F, sep = "\t",col.names = NA )

TPM_BG4peak_in_promoter_frequencies_NSC <- TPM_BG4peak_in_promoter %>%  group_by(NSC_strat_expression,BG4_NSC) %>% summarise(count = n())
write.table(TPM_BG4peak_in_promoter_frequencies_NSC,file ='TPM_BG4peak_in_promoter_frequencies_NSC.data.txt',quote = F, sep = "\t",col.names = NA )

TPM_BG4peak_in_promoter_frequencies_CNCC <- TPM_BG4peak_in_promoter %>%  group_by(CNCC_strat_expression,BG4_CNCC) %>% summarise(count = n())
write.table(TPM_BG4peak_in_promoter_frequencies_CNCC,file ='TPM_BG4peak_in_promoter_frequencies_CNCC.data.txt',quote = F, sep = "\t",col.names = NA )


TPM_BG4peak_in_promoter_frequencies <- TPM_BG4peak_in_promoter %>%  group_by(ESC_strat_expression) %>% summarise(count = n())
write.table(TPM_BG4peak_in_promoter_frequencies,file ='TPM_BG4peak_in_promoter_frequencies.expression_subgroups.ESC.data.txt',quote = F, sep = "\t",col.names = NA )


TPM_BG4peak_in_promoter_frequencies <- TPM_BG4peak_in_promoter %>%  group_by(NSC_strat_expression) %>% summarise(count = n())
write.table(TPM_BG4peak_in_promoter_frequencies,file ='TPM_BG4peak_in_promoter_frequencies.expression_subgroups.NSC.data.txt',quote = F, sep = "\t",col.names = NA )


TPM_BG4peak_in_promoter_frequencies <- TPM_BG4peak_in_promoter %>%  group_by(CNCC_strat_expression) %>% summarise(count = n())
write.table(TPM_BG4peak_in_promoter_frequencies,file ='TPM_BG4peak_in_promoter_frequencies.expression_subgroups.CNCC.data.txt',quote = F, sep = "\t",col.names = NA )


## ==  print to file bed intervals for TSS density plots ==
TSS_coordinates <- data.frame(TSS_coordinates)
TSS_coordinates$ens_id <- gsub('\\..*','',TSS_coordinates$V4)
TSS_coordinates <- TSS_coordinates[!duplicated(TSS_coordinates$ens_id), ]


# genes that are DE in ES, DE in NSC, notDE -- > check density plots in ESC and NSC
TSS_coordinates <- TSS_coordinates %>% dplyr::select("V1","V2","V3","ens_id","V6","V7")

TSS_DEA_NSC_vs_ESC <- left_join(DEA_outcome_NSC_vs_ESC,TSS_coordinates, by = 'ens_id') %>% dplyr::select("V1","V2","V3","ens_id","V6","V7","DEA_NSC_vs_ESC")
TSS_DEA_NSC_vs_ESC_for_setdiff <- TSS_DEA_NSC_vs_ESC %>% dplyr::select("V1","V2","V3","ens_id","V6","V7")

TSS_DEA_CNCC_vs_ESC <- left_join(DEA_outcome_CNCC_vs_ESC,TSS_coordinates, by = 'ens_id') %>% dplyr::select("V1","V2","V3","ens_id","V6","V7","DEA_CNCC_vs_ESC")
TSS_DEA_CNCC_vs_ESC_for_setdiff <- TSS_DEA_CNCC_vs_ESC %>% dplyr::select("V1","V2","V3","ens_id","V6","V7")

TSS_DEA_NSC_vs_CNCC <- left_join(DEA_outcome_NSC_vs_CNCC,TSS_coordinates, by = 'ens_id') %>% dplyr::select("V1","V2","V3","ens_id","V6","V7","DEA_NSC_vs_CNCC")
TSS_DEA_NSC_vs_CNCC_for_setdiff <- TSS_DEA_NSC_vs_CNCC %>% dplyr::select("V1","V2","V3","ens_id","V6","V7")

TSS_DEA_NSC_vs_ESC_plus1 <- TSS_DEA_NSC_vs_ESC %>% filter(DEA_NSC_vs_ESC==1)
TSS_DEA_NSC_vs_ESC_minus1 <- TSS_DEA_NSC_vs_ESC %>% filter(DEA_NSC_vs_ESC== -1)
TSS_NOTDE_NSC_vs_ESC <- dplyr::setdiff(TSS_coordinates,TSS_DEA_NSC_vs_ESC_for_setdiff)

write.table(TSS_NOTDE_NSC_vs_ESC, file='TSS_NOTDE_NSC_vs_ESC.bed',sep ="\t",quote = F,col.names = F,row.names = F)
write.table(TSS_DEA_NSC_vs_ESC_plus1, file='TSS_DEA_NSC_vs_ESC_plus1.bed',sep ="\t",quote = F,col.names = F,row.names = F)
write.table(TSS_DEA_NSC_vs_ESC_minus1, file='TSS_DEA_NSC_vs_ESC_minus1.bed',sep ="\t",quote = F,col.names = F,row.names = F)

TSS_DEA_CNCC_vs_ESC_plus1 <- TSS_DEA_CNCC_vs_ESC %>% filter(DEA_CNCC_vs_ESC==1)
TSS_DEA_CNCC_vs_ESC_minus1 <- TSS_DEA_CNCC_vs_ESC %>% filter(DEA_CNCC_vs_ESC== -1)
TSS_NOTDE_CNCC_vs_ESC <- dplyr::setdiff(TSS_coordinates,TSS_DEA_CNCC_vs_ESC_for_setdiff)

write.table(TSS_NOTDE_CNCC_vs_ESC, file='TSS_NOTDE_CNCC_vs_ESC.bed',sep ="\t",quote = F,col.names = F,row.names = F)
write.table(TSS_DEA_CNCC_vs_ESC_plus1, file='TSS_DEA_CNCC_vs_ESC_plus1.bed',sep ="\t",quote = F,col.names = F,row.names = F)
write.table(TSS_DEA_CNCC_vs_ESC_minus1, file='TSS_DEA_CNCC_vs_ESC_minus1.bed',sep ="\t",quote = F,col.names = F,row.names = F)

TSS_DEA_NSC_vs_CNCC_plus1 <- TSS_DEA_NSC_vs_CNCC %>% filter(DEA_NSC_vs_CNCC==1)
TSS_DEA_NSC_vs_CNCC_minus1 <- TSS_DEA_NSC_vs_CNCC %>% filter(DEA_NSC_vs_CNCC== -1)
TSS_NOTDE_NSC_vs_CNCC <- dplyr::setdiff(TSS_coordinates,TSS_DEA_NSC_vs_CNCC_for_setdiff)

write.table(TSS_NOTDE_NSC_vs_CNCC, file='TSS_NOTDE_NSC_vs_CNCC.bed',sep ="\t",quote = F,col.names = F,row.names = F)
write.table(TSS_DEA_NSC_vs_CNCC_plus1, file='TSS_DEA_NSC_vs_CNCC_plus1.bed',sep ="\t",quote = F,col.names = F,row.names = F)
write.table(TSS_DEA_NSC_vs_CNCC_minus1, file='TSS_DEA_NSC_vs_CNCC_minus1.bed',sep ="\t",quote = F,col.names = F,row.names = F)
```

density plots for genes that are DE
-----------------------------------

``` bash
#move the TSS files generated in R to the dedicated folder
scp TSS*bed simeon01@clust1-headnode:/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/DEA_cells/*
for file in *NOTDE*
  do
  shuf -n 5000 ${file}| sortBed -i - | slopBed -i - -l 1 -r 0 -g /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome > ${file%%.bed}.shuf.sorted.touse.bed
  done


for file in *TSS_DEA_*bed*
  do
  sortBed -i $file | slopBed -i - -l 1 -r 0 -g /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome > ${file%%.bed}.sorted.touse.bed
  done

Bam_repA=/scratcha/sblab/simeon01/Data/190829_Katie_stemcells/SLX-18496/aligned #.hg38.sort.markduplicates.bam
Bam_repB=/scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4/aligned/merged # markduplicates.bam
Bam_repC=/scratcha/sblab/simeon01/Data/20190917_SLX-18355_Katie_StemCells_Rep2_3xBG4/SLX-18355/aligned #.hg38.sort.markduplicates.bam

bed_regions_path=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/DEA_cells
cd $bed_regions_path


## == case NSC vs ESC
for file in $Bam_repA/*SC*bw
do
  filename=${file##*/}
  echo $filename
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $bed_regions_path/*NSC_vs_ESC*.touse.bed -bs 20 -a 1500 -b 1500 -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.gz && plotProfile --matrixFile $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.gz -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.gz -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.heatmap.png --averageType median"
  echo " ==== "
done
for file in $Bam_repB/*SC*bw
do
  filename=${file##*/}
  echo $filename
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $bed_regions_path/*NSC_vs_ESC*.touse.bed -bs 20 -a 1500 -b 1500 -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.gz && plotProfile --matrixFile $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.gz -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.gz -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.heatmap.png --averageType median"
  echo " ==== "
done
for file in $Bam_repC/*SC*bw
do
  filename=${file##*/}
  echo $filename
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $bed_regions_path/*NSC_vs_ESC*.touse.bed -bs 20 -a 1500 -b 1500 -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.gz && plotProfile --matrixFile $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.gz -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.png --averageType median --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.gz -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_ESC.heatmap.png --averageType median"
  echo " ==== "
done
## ---------------------------------------------------------------------------------------------------------------------------------------------------
## == CNCC and ESC
for file in `ls $Bam_repA/*bw | grep -v NSC`
do
  filename=${file##*/}
  echo $filename
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $bed_regions_path/*CNCC_vs_ESC*.touse.bed -bs 20 --missingDataAsZero -a 1500 -b 1500 -o $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.gz && plotProfile --matrixFile $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.gz -o $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESS.profile.png --averageType median --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.gz -o $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.heatmap.png --averageType mean"
  echo " ==== "
done
for file in `ls $Bam_repB/*bw | grep -v NSC`
do
  filename=${file##*/}
  echo $filename
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $bed_regions_path/*CNCC_vs_ESC*.touse.bed -bs 20 -a 1500 -b 1500 -o $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.gz && plotProfile --matrixFile $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.gz -o $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.profile.png --averageType median --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.gz -o $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.heatmap.png --averageType median"
  echo " ==== "
done
for file in `ls $Bam_repC/*bw | grep -v NSC`
do
  filename=${file##*/}
  echo $filename
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $bed_regions_path/*CNCC_vs_ESC*.touse.bed -bs 20 -a 1500 -b 1500 -o $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.gz && plotProfile --matrixFile $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.gz -o $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.profile.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.gz -o $bed_regions_path/${filename%%.bw}.DEA_CNCC_vs_ESC.heatmap.png --averageType median"
  echo " ==== "
done


## ---------------------------------------------------------------------------------------------------------------------------------------------------
## case NSC vs CNCC
for file in `ls $Bam_repA/*bw | grep -v ESC`
do
  filename=${file##*/}
  echo $filename
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $bed_regions_path/*NSC_vs_CNCC*.touse.bed -bs 20 -a 1500 -b 1500 -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.gz && plotProfile --matrixFile $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.gz -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.gz -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.heatmap.png --averageType median"
  echo " ==== "
done
for file in `ls $Bam_repB/*bw | grep -v ESC`
do
  filename=${file##*/}
  echo $filename
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $bed_regions_path/*NSC_vs_CNCC*.touse.bed -bs 20 -a 1500 -b 1500  -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.gz && plotProfile --matrixFile $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.gz -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.gz -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.heatmap.png --averageType median"
  echo " ==== "
done
for file in `ls $Bam_repC/*bw | grep -v ESC`
do
  filename=${file##*/}
  echo $filename
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $bed_regions_path/*NSC_vs_CNCC*.touse.bed -bs 20 -a 1500 -b 1500  -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.gz && plotProfile --matrixFile $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.gz -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.gz -o $bed_regions_path/${filename%%.bw}.DEA_NSC_vs_CNCC.heatmap.png --averageType median"
  echo " ==== "
done
```

Check single matrix produced - from deeptools and generate density plot within R
--------------------------------------------------------------------------------

``` r
setwd("/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/DEA_cells/newFigures/matrices_RPM_densities")

files_CNCC_vs_ESC <- list.files(pattern="CNCC_vs_ESC")
ESC_files_CNCC_vs_ESC<- grep("ESC_Rep",files_CNCC_vs_ESC, value = T)
ESC_files_CNCC_vs_ESC_no_input<- ESC_files_CNCC_vs_ESC[-grep("input",ESC_files_CNCC_vs_ESC)]

ESC_big_matrix_minus1 <- c()
ESC_big_matrix_plus1 <- c()
ESC_big_matrix_zeros <- c()

for (i in 1:length(ESC_files_CNCC_vs_ESC_no_input)){
  temp_Matrix <- read.table(ESC_files_CNCC_vs_ESC_no_input[i],stringsAsFactors = F, skip = 1)
  Matrix_values <- temp_Matrix[,7:dim(temp_Matrix)[2]]
  
  ESC_big_matrix_minus1 <- rbind(ESC_big_matrix_minus1,Matrix_values[1:1521,])
  ESC_big_matrix_plus1 <- rbind(ESC_big_matrix_plus1,Matrix_values[(1521+1):(1521+491),])
  ESC_big_matrix_zeros <- rbind(ESC_big_matrix_zeros,Matrix_values[(1521+491+1):7012,])
}


y_minus1 <-colMedians(as.matrix(ESC_big_matrix_minus1),na.rm = T)
y_plus1 <-colMedians(as.matrix(ESC_big_matrix_plus1),na.rm = T)
y_zeros <-colMedians(as.matrix(ESC_big_matrix_zeros),na.rm = T)

y_minus1s <- loess.smooth(x=seq(1:length(y_minus1)),y_minus1,span = 0.15)
y_minus1s_a <- approx(y_minus1s$x, y_minus1s$y,seq(1,150))$y

y_plus1s <- loess.smooth(x=seq(1:length(y_plus1)),y_plus1,span = 0.15)
y_plus1s_a <- approx(y_plus1s$x, y_plus1s$y,seq(1,150))$y

y_zeross <- loess.smooth(x=seq(1:length(y_zeros)),y_zeros,span = 0.15)
y_zeross_a <- approx(y_zeross$x, y_zeross$y,seq(1,150))$y
ESCbg4_DEA_CNCC_ESC <-data.frame(median_norm_BG4=c(y_minus1s_a,
                                 y_plus1s_a,
                                 y_zeross_a),
                                 TSS_x <- rep(seq(-1500,(1500-20),20),3),
                                 case_DE = c(rep('minus1',150),rep('plus1',150),rep('zero',150)))

pdf("ESCbg4_DE_CNCCvsESC.pdf")
plot <- ESCbg4_DEA_CNCC_ESC %>% ggplot(aes(x=TSS_x,y=median_norm_BG4,group = case_DE,color=case_DE)) + geom_line(color=c(rep('blue',150),rep('cyan',150),rep('yellow',150)),size=0.7) + theme_bw() + ggtitle('ESCbg4 - CNCCvsESC')
print(plot)
dev.off()

CNCC_files_CNCC_vs_ESC<- grep("CNCC_Rep",files_CNCC_vs_ESC, value = T)
CNCC_files_CNCC_vs_ESC_no_input<- CNCC_files_CNCC_vs_ESC[-grep("input",CNCC_files_CNCC_vs_ESC)]

CNCC_big_matrix_minus1 <- c()
CNCC_big_matrix_plus1 <- c()
CNCC_big_matrix_zeros <- c()

for (i in 1:length(CNCC_files_CNCC_vs_ESC_no_input)){
  temp_Matrix <- read.table(CNCC_files_CNCC_vs_ESC_no_input[i],stringsAsFactors = F, skip = 1)
  Matrix_values <- temp_Matrix[,7:dim(temp_Matrix)[2]]
  CNCC_big_matrix_minus1 <- rbind(CNCC_big_matrix_minus1,Matrix_values[1:1521,])
  CNCC_big_matrix_plus1 <- rbind(CNCC_big_matrix_plus1,Matrix_values[(1521+1):(1521+491),])
  CNCC_big_matrix_zeros <- rbind(CNCC_big_matrix_zeros,Matrix_values[(1521+491+1):7012,])
}
y_minus1 <-colMedians(as.matrix(CNCC_big_matrix_minus1),na.rm = T)
y_plus1 <-colMedians(as.matrix(CNCC_big_matrix_plus1),na.rm = T)
y_zeros <-colMedians(as.matrix(CNCC_big_matrix_zeros),na.rm = T)

y_minus1s <- loess.smooth(x=seq(1:length(y_minus1)),y_minus1,span = 0.1)
y_minus1s_a <- approx(y_minus1s$x, y_minus1s$y,seq(1,150))$y

y_plus1s <- loess.smooth(x=seq(1:length(y_plus1)),y_plus1,span = 0.1)
y_plus1s_a <- approx(y_plus1s$x, y_plus1s$y,seq(1,150))$y

y_zeross <- loess.smooth(x=seq(1:length(y_zeros)),y_zeros,span = 0.1)
y_zeross_a <- approx(y_zeross$x, y_zeross$y,seq(1,150))$y
CNCCbg4_DEA_CNCC_ESC <-data.frame(median_norm_BG4=c(y_minus1s_a,
                                 y_plus1s_a,
                                 y_zeross_a),
                                 TSS_x <- rep(seq(-1500,(1500-20),20),3),
                                 case_DE = c(rep('minus1',150),rep('plus1',150),rep('zero',150)))

pdf("CNCCbg4_DE_CNCCvsESC.pdf")
plot <- CNCCbg4_DEA_CNCC_ESC %>% ggplot(aes(x=TSS_x,y=median_norm_BG4,group = case_DE,color=case_DE)) + geom_line(color=c(rep('blue',150),rep('cyan',150),rep('yellow',150)),size=0.7) + theme_bw()+ ggtitle('CNCCbg4 - CNCCvsESC')
print(plot)
dev.off()

## -------------------------------------------------------------------------------------------------------------------------------------

files_NSC_vs_ESC <- list.files(pattern="NSC_vs_ESC")
ESC_files_NSC_vs_ESC<- grep("ESC_Rep",files_NSC_vs_ESC, value = T)
ESC_files_NSC_vs_ESC_no_input<- ESC_files_NSC_vs_ESC[-grep("input",ESC_files_NSC_vs_ESC)]

ESC_big_matrix_minus1 <- c()
ESC_big_matrix_plus1 <- c()
ESC_big_matrix_zeros <- c()

for (i in 1:length(ESC_files_NSC_vs_ESC_no_input)){
  temp_Matrix <- read.table(ESC_files_NSC_vs_ESC_no_input[i],stringsAsFactors = F, skip = 1)
  Matrix_values <- temp_Matrix[,7:dim(temp_Matrix)[2]]
  
  ESC_big_matrix_minus1 <- rbind(ESC_big_matrix_minus1,Matrix_values[1:1436,])
  ESC_big_matrix_plus1 <- rbind(ESC_big_matrix_plus1,Matrix_values[(1436+1):(1436+695),])
  ESC_big_matrix_zeros <- rbind(ESC_big_matrix_zeros,Matrix_values[(1436+695+1):7131,])
}


y_minus1 <-colMedians(as.matrix(ESC_big_matrix_minus1),na.rm = T)
y_plus1 <-colMedians(as.matrix(ESC_big_matrix_plus1),na.rm = T)
y_zeros <-colMedians(as.matrix(ESC_big_matrix_zeros),na.rm = T)

y_minus1s <- loess.smooth(x=seq(1:length(y_minus1)),y_minus1,span = 0.19)
y_minus1s_a <- approx(y_minus1s$x, y_minus1s$y,seq(1,150))$y

y_plus1s <- loess.smooth(x=seq(1:length(y_plus1)),y_plus1,span = 0.19)
y_plus1s_a <- approx(y_plus1s$x, y_plus1s$y,seq(1,150))$y

y_zeross <- loess.smooth(x=seq(1:length(y_zeros)),y_zeros,span = 0.19)
y_zeross_a <- approx(y_zeross$x, y_zeross$y,seq(1,150))$y

ESCbg4_DEA_NSC_ESC <-data.frame(median_norm_BG4=c(y_minus1s_a,
                                 y_plus1s_a,
                                 y_zeross_a),
                                 TSS_x <- rep(seq(-1500,(1500-20),20),3),
                                 case_DE = c(rep('minus1',150),rep('plus1',150),rep('zero',150)))

pdf("ESCbg4_DE_NSCvsESC.pdf")
plot <- ESCbg4_DEA_NSC_ESC %>% ggplot(aes(x=TSS_x,y=median_norm_BG4,group = case_DE,color=case_DE)) + geom_line(color=c(rep('blue',150),rep('cyan',150),rep('yellow',150)),size=0.7) + theme_bw() + ggtitle('ESCbg4 - NSCvsESC ')
print(plot)
dev.off()

NSC_files_NSC_vs_ESC<- grep("NSC_Rep",files_NSC_vs_ESC, value = T)
NSC_files_NSC_vs_ESC_no_input<- NSC_files_NSC_vs_ESC[-grep("input",NSC_files_NSC_vs_ESC)]

NSC_big_matrix_minus1 <- c()
NSC_big_matrix_plus1 <- c()
NSC_big_matrix_zeros <- c()

for (i in 1:length(NSC_files_NSC_vs_ESC_no_input)){
  temp_Matrix <- read.table(NSC_files_NSC_vs_ESC_no_input[i],stringsAsFactors = F, skip = 1)
  Matrix_values <- temp_Matrix[,7:dim(temp_Matrix)[2]]
  NSC_big_matrix_minus1 <- rbind(NSC_big_matrix_minus1,Matrix_values[1:1436,])
  NSC_big_matrix_plus1 <- rbind(NSC_big_matrix_plus1,Matrix_values[(1436+1):(1436+695),])
  NSC_big_matrix_zeros <- rbind(NSC_big_matrix_zeros,Matrix_values[(1436+695+1):7131,])
}
y_minus1 <-colMedians(as.matrix(NSC_big_matrix_minus1),na.rm = T)
y_plus1 <-colMedians(as.matrix(NSC_big_matrix_plus1),na.rm = T)
y_zeros <-colMedians(as.matrix(NSC_big_matrix_zeros),na.rm = T)

y_minus1s <- loess.smooth(x=seq(1:length(y_minus1)),y_minus1,span = 0.15)
y_minus1s_a <- approx(y_minus1s$x, y_minus1s$y,seq(1,150))$y

y_plus1s <- loess.smooth(x=seq(1:length(y_plus1)),y_plus1,span = 0.15)
y_plus1s_a <- approx(y_plus1s$x, y_plus1s$y,seq(1,150))$y

y_zeross <- loess.smooth(x=seq(1:length(y_zeros)),y_zeros,span = 0.15)
y_zeross_a <- approx(y_zeross$x, y_zeross$y,seq(1,150))$y
NSCbg4_DEA_NSC_ESC <-data.frame(median_norm_BG4=c(y_minus1s_a,
                                 y_plus1s_a,
                                 y_zeross_a),
                                 TSS_x <- rep(seq(-1500,(1500-20),20),3),
                                 case_DE = c(rep('minus1',150),rep('plus1',150),rep('zero',150)))

pdf("NSCbg4_DE_NSCvsESC.pdf")
plot <- NSCbg4_DEA_NSC_ESC %>% ggplot(aes(x=TSS_x,y=median_norm_BG4,group = case_DE,color=case_DE)) + geom_line(color=c(rep('blue',150),rep('cyan',150),rep('yellow',150)),size=0.7) + theme_bw()+ ggtitle('NSCbg4 - NSCvsESC')
print(plot)
dev.off()

# ----------------------------------------------------------------------------------------------------------

files_NSC_vs_CNCC <- list.files(pattern="NSC_vs_CNCC")
CNCC_files_NSC_vs_CNCC<- grep("CNCC_Rep",files_NSC_vs_CNCC, value = T)
CNCC_files_NSC_vs_CNCC_no_input<- CNCC_files_NSC_vs_CNCC[-grep("input",CNCC_files_NSC_vs_CNCC)]

CNCC_big_matrix_minus1 <- c()
CNCC_big_matrix_plus1 <- c()
CNCC_big_matrix_zeros <- c()

for (i in 1:length(CNCC_files_NSC_vs_CNCC_no_input)){
  temp_Matrix <- read.table(CNCC_files_NSC_vs_CNCC_no_input[i],stringsAsFactors = F, skip = 1)
  Matrix_values <- temp_Matrix[,7:dim(temp_Matrix)[2]]
  
  CNCC_big_matrix_minus1 <- rbind(CNCC_big_matrix_minus1,Matrix_values[1:728,])
  CNCC_big_matrix_plus1 <- rbind(CNCC_big_matrix_plus1,Matrix_values[(728+1):(728+894),])
  CNCC_big_matrix_zeros <- rbind(CNCC_big_matrix_zeros,Matrix_values[(728+894+1):6622,])
}


y_minus1 <-colMedians(as.matrix(CNCC_big_matrix_minus1),na.rm = T)
y_plus1 <-colMedians(as.matrix(CNCC_big_matrix_plus1),na.rm = T)
y_zeros <-colMedians(as.matrix(CNCC_big_matrix_zeros),na.rm = T)

y_minus1s <- loess.smooth(x=seq(1:length(y_minus1)),y_minus1,span = 0.2)
y_minus1s_a <- approx(y_minus1s$x, y_minus1s$y,seq(1,150))$y

y_plus1s <- loess.smooth(x=seq(1:length(y_plus1)),y_plus1,span = 0.2)
y_plus1s_a <- approx(y_plus1s$x, y_plus1s$y,seq(1,150))$y

y_zeross <- loess.smooth(x=seq(1:length(y_zeros)),y_zeros,span = 0.2)
y_zeross_a <- approx(y_zeross$x, y_zeross$y,seq(1,150))$y
CNCCbg4_DEA_NSC_CNCC <-data.frame(median_norm_BG4=c(y_minus1s_a,
                                 y_plus1s_a,
                                 y_zeross_a),
                                 TSS_x <- rep(seq(-1500,(1500-20),20),3),
                                 case_DE = c(rep('minus1',150),rep('plus1',150),rep('zero',150)))

pdf("CNCCbg4_DE_NSCvsCNCC.pdf")
plot <- CNCCbg4_DEA_NSC_CNCC %>% ggplot(aes(x=TSS_x,y=median_norm_BG4,group = case_DE,color=case_DE)) + geom_line(color=c(rep('blue',150),rep('cyan',150),rep('yellow',150)),size=0.7) + theme_bw() + ggtitle('CNCCbg4 - NSCvsCNCC ')
print(plot)
dev.off()

NSC_files_NSC_vs_CNCC<- grep("NSC_Rep",files_NSC_vs_CNCC, value = T)
NSC_files_NSC_vs_CNCC_no_input<- NSC_files_NSC_vs_CNCC[-grep("input",NSC_files_NSC_vs_CNCC)]

NSC_big_matrix_minus1 <- c()
NSC_big_matrix_plus1 <- c()
NSC_big_matrix_zeros <- c()

for (i in 1:length(NSC_files_NSC_vs_CNCC_no_input)){
  temp_Matrix <- read.table(NSC_files_NSC_vs_CNCC_no_input[i],stringsAsFactors = F, skip = 1)
  Matrix_values <- temp_Matrix[,7:dim(temp_Matrix)[2]]
  NSC_big_matrix_minus1 <- rbind(NSC_big_matrix_minus1,Matrix_values[1:728,])
  NSC_big_matrix_plus1 <- rbind(NSC_big_matrix_plus1,Matrix_values[(728+1):(728+894),])
  NSC_big_matrix_zeros <- rbind(NSC_big_matrix_zeros,Matrix_values[(728+894+1):6622,])
}
y_minus1 <-colMedians(as.matrix(NSC_big_matrix_minus1),na.rm = T)
y_plus1 <-colMedians(as.matrix(NSC_big_matrix_plus1),na.rm = T)
y_zeros <-colMedians(as.matrix(NSC_big_matrix_zeros),na.rm = T)

y_minus1s <- loess.smooth(x=seq(1:length(y_minus1)),y_minus1,span = 0.15)
y_minus1s_a <- approx(y_minus1s$x, y_minus1s$y,seq(1,150))$y

y_plus1s <- loess.smooth(x=seq(1:length(y_plus1)),y_plus1,span = 0.15)
y_plus1s_a <- approx(y_plus1s$x, y_plus1s$y,seq(1,150))$y

y_zeross <- loess.smooth(x=seq(1:length(y_zeros)),y_zeros,span = 0.15)
y_zeross_a <- approx(y_zeross$x, y_zeross$y,seq(1,150))$y

NSCbg4_DEA_NSC_CNCC <-data.frame(median_norm_BG4=c(y_minus1s_a,
                                 y_plus1s_a,
                                 y_zeross_a),
                                 TSS_x <- rep(seq(-1500,(1500-20),20),3),
                                 case_DE = c(rep('minus1',150),rep('plus1',150),rep('zero',150)))

pdf("NSCbg4_DE_NSCvsCNCC.pdf")
plot <- NSCbg4_DEA_NSC_CNCC %>% ggplot(aes(x=TSS_x,y=median_norm_BG4,group = case_DE,color=case_DE)) + geom_line(color=c(rep('blue',150),rep('cyan',150),rep('yellow',150)),size=0.7) + theme_bw()+ ggtitle('NSCbg4 - NSCvsCNCC')
print(plot)
dev.off()
```

Check relationship between cell-specific and cell common peaks and open chromatin state either in the specific cell line or in the other cell line
--------------------------------------------------------------------------------------------------------------------------------------------------

``` bash

# I run this check only for scenario1 - for either common or specific sets

#open chromatine in ES is the faire-seq
open_ESC=/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/GSM602298_ESC_FAIRE_aligned_reduced_lifted.hg38.bed

#open_chromatine in CNCC 
open_CNCC=/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/common_atac_sites_CNCC_liftOver_hg19_hg38.bed

open_chrom_list=($open_ESC $open_CNCC)

for open in ${open_chrom_list[@]}
do
  echo $open
  for file in `ls *CNCC*bed | grep -v NSC`
  do
    wc -l $file
    intersectBed -a $file -b $open -wa  | sort | uniq | wc -l
  done
  echo "===="
done

cd 
for open in ${open_chrom_list[@]}
do
  echo $open
  for file in *bed
  do
    wc -l $file
    intersectBed -a $file -b $open -wa  | sort | uniq | wc -l
  done
  echo "===="
done
```

check overlap between stem cells and other cancer cell type
-----------------------------------------------------------

``` bash
other_cell_types_folder='/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/nature_protocol_cell_lines_from_gio_archive'
ESC_bed_scenario1=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC.bio2out3.bed
multi3_stemCells=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC_NSC_CNCC_bio2out3.multi3.bed
genome_file='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.sorted.genome'

cd $other_cell_types_folder

# lift peaks from hg19 to hg38
for file in *peaks.bed; 
  do echo $file; 
  /home/simeon01/applications/liftOver ${file} /scratcha/sblab/simeon01/reference_genomes/liftover_chain/hg19ToHg38.over.chain.gz ${file%%.bed}.hg19_lifted_To_hg38.bed unmapped.bed; 
done

for file in *G4.bed ; 
  do echo $file; 
  /home/simeon01/applications/liftOver ${file} /scratcha/sblab/simeon01/reference_genomes/liftover_chain/hg19ToHg38.over.chain.gz ./peaks_lifted_to_hg38/${file%%.bed}.hg19_lifted_To_hg38.bed ./peaks_lifted_to_hg38/unmapped.bed; 
done


/home/simeon01/applications/liftOver DMSO_no_Dm6.hg19all.25M.inp.q005.all_peaks.narrowPeak.multi2.bed /scratcha/sblab/simeon01/reference_genomes/liftover_chain/hg19ToHg38.over.chain.gz DMSO_no_Dm6.hg19all.25M.inp.q005.all_peaks.narrowPeak.multi2.hg19_lifted_To_hg38.bed unmapped.bed
for file in *_lifted_To_hg38.bed
do
#wc -l $file | awk '{print $1}'
intersectBed -a $file -b $ESC_bed_scenario1 -wa | sort | uniq | wc -l
#echo "==========="
done

# sort other cancer systems
for file in *hg38.bed
do
sortBed -i $file > ${file%%.bed}.sorted.bed
done

# compute core of regions conserved in all systems previoously studied
multiIntersectBed -i *common*lifted*sorted.bed | awk '{if($4>=3) print $0 }' | sortBed -i - | mergeBed -i - > common_5_cell_systems.hg38.bed
wc -l common_5_cell_systems.hg38.bed
#5944

#intervene pairsise
mkdir intervene_output
sbatch --mem 12G --wrap "intervene pairwise -i ./*lifted*sorted.bed  $ESC_bed_scenario1  --type genomic  --compute frac --htype pie --filenames --output ./intervene_output"


#check overlap common_5_cell_systems with stem cells
intersectBed -a common_5_cell_systems.hg38.bed -b $ESC_bed_scenario1 -wa | sort |uniq | wc -l
#5197

# specific peaks per cell line
for file in *common*lifted*sorted.bed
do
intersectBed -a $file -b common_5_cell_systems.hg38.bed -v | sortBed -i - | mergeBed -i - > ${file%%.sorted.bed}.sorted.specific.bed
done

for file in *common*specific.bed
do
intersectBed -a $file -b $ESC_bed_scenario1 -wa | sort |uniq | wc -l
done
#1682
#307
#1851
#3008

wc -l *common*specific.bed
#6039
#894
#3373
#8159

intersectBed -a common_5_cell_systems.hg38.bed -b /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC_NSC_CNCC_bio2out3.multi3.bed -wa | sort | uniq | wc -l
2654

wc -l /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC_NSC_CNCC_bio2out3.multi3.bed
#3426

intersectBed -b common_5_cell_systems.hg38.bed -a /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC_NSC_CNCC_bio2out3.multi3.bed -wa | sort | uniq | wc -l
#2721

for file in *_lifted_To_hg38.bed; do intersectBed -a $file -b $OQs_hg38 -wa | sort | uniq | wc -l; done

for file in *_lifted_To_hg38.bed; 
  do 
  slopBed -b 500 -i $ESC_bed_scenario1 -g $genome_file | sortBed -i - | intersectBed -a $file -b - -wa | sort | uniq | wc -l; 
  done

ls $other_cell_types_folder/*_lifted_To_hg38.bed


for file in *_lifted_To_hg38.bed
do
#wc -l $file | awk '{print $1}'
intersectBed -a $file -b $multi3_stemCells -wa | sort | uniq | wc -l
#echo "==========="
done

# include cell-specific (ESC-specific, NSC-specific and CNCC-specific) and cell-common peaks and re-run overlap

# considert overlap of all cancers and check overlap of that with ESC-all, ESC-specific and ESC-NSC-CNCC common
```

check occurance of BG4 peaks in promoter of genes from specific gene lists -curated by Katie
--------------------------------------------------------------------------------------------

``` r
library(tidyverse)
path_gene_lists <- '/Users/simeon01/Documents/Katie/Gene_lists_for_Angela/'

biomart_mapping_file <- read.delim('/Users/simeon01/Documents/genomes/hg38/mart_export_ensemb_GRCh38.p13_ensembl_gene_id_entrezgene_id.txt',sep = ",",stringsAsFactors = F,header = T) # 1st column is ensembl_id, 2nd column is entrez_id

Big_Table <- as_tibble(biomart_mapping_file)

txt_gene_lists <- list.files(path_gene_lists,pattern = '.txt')

for (i in 1:length(txt_gene_lists)){
  
  temp_file <- read.table(paste0(path_gene_lists,txt_gene_lists[i]),stringsAsFactors = F)
  names(temp_file) <- 'entrezgene_id'
  temp_file$case <- rep(1,length(temp_file$entrezgene_id))
  colnames(temp_file)[2] <- gsub('.txt','',txt_gene_lists[i])
  #left join with Big_Table
  #Big_Table <- head(Big_Table)
  Big_Table <- Big_Table %>%dplyr::left_join(temp_file, by='entrezgene_id')
}

#Big_Table contains info about all of them 

#Big_Table <- Big_Table %>% left_join(BG4_DBA,by=ensembl_gene_id)
#colnames(Big_Table)[1] <- "ens"

TPM_BG4peak_in_promoter$ensembl_gene_id <- gsub('\\..*','',TPM_BG4peak_in_promoter$ens)
TPM_BG4peak_in_promoter_joinedTo_Big_Table <- left_join(TPM_BG4peak_in_promoter,Big_Table,by = "ensembl_gene_id")

length(which(TPM_BG4peak_in_promoter_joinedTo_Big_Table$Chang_2011_HK_genes==1))

#check now frequency of BG+ and BG4-  for various gene lists 
# barplot of fraction of genes in each group BG4+ and BG4-

TPM_BG4peak_in_promoter_joinedTo_Big_Table %>%group_by('ens') %>% dplyr::select(c('median_TPM_ESC','BG4_ESC'), contains('Xie_2013_all_lineage'))
temp_Matrix_Xie_2013 <- TPM_BG4peak_in_promoter_joinedTo_Big_Table %>% dplyr::select(c('ens','median_TPM_ESC'),contains('BG4_'),contains('Xie_2013')) %>% dplyr::filter_at(vars(matches("Xie_2013")), all_vars(!is.na(.)))

temp_Matrix_Chang <- TPM_BG4peak_in_promoter_joinedTo_Big_Table %>% dplyr::select(c('ens','median_TPM_ESC'),contains('BG4_'),contains('Chang_2011')) %>% dplyr::filter(!is.na(Chang_2011_HK_genes) | !is.na(Chang_2011_TS_genes))

temp_Matrix_Chang %>% ggplot(aes(x=Chang_2011_HK_genes, group=BG4_ESC)) + geom_bar()

temp_Matrix_Chang %>%filter(Chang_2011_HK_genes==1) %>% group_by(BG4_ESC) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% ggplot(aes(x=BG4_ESC,group=BG4_ESC,fill = BG4_ESC)) + geom_bar(aes(y=freq), stat="identity") + ggtitle('Chang_2011_HK_genes')

temp_Matrix_Chang %>%filter(Chang_2011_TS_genes==1) %>% group_by(BG4_ESC) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% ggplot(aes(x=BG4_ESC,group=BG4_ESC,fill = BG4_ESC)) + geom_bar(aes(y=freq), stat="identity") + ggtitle('Chang_2011_TS_genes')

temp_Matrix_Chang %>%filter(Chang_2011_HK_genes==1) %>% group_by(BG4_CNCC) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% ggplot(aes(x=BG4_CNCC,group=BG4_CNCC,fill = BG4_CNCC)) + geom_bar(aes(y=freq), stat="identity") + ggtitle('Chang_2011_HK_genes')

temp_Matrix_Chang %>%filter(Chang_2011_TS_genes==1) %>% group_by(BG4_CNCC) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% ggplot(aes(x=BG4_CNCC,group=BG4_CNCC,fill = BG4_CNCC)) + geom_bar(aes(y=freq), stat="identity") + ggtitle('Chang_2011_TS_genes')

temp_Matrix_Chang %>%filter(Chang_2011_HK_genes==1) %>% group_by(BG4_NSC) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% ggplot(aes(x=BG4_NSC,group=BG4_NSC,fill = BG4_NSC)) + geom_bar(aes(y=freq), stat="identity") + ggtitle('Chang_2011_HK_genes')

temp_Matrix_Chang %>%filter(Chang_2011_TS_genes==1) %>% group_by(BG4_NSC) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% ggplot(aes(x=BG4_NSC,group=BG4_NSC,fill = BG4_NSC)) + geom_bar(aes(y=freq), stat="identity") + ggtitle('Chang_2011_TS_genes')


temp_Matrix_Xie_2013 %>% filter(Xie_2013_all_lineage_gene_H1_specific==1) %>% group_by(BG4_ESC) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% ggplot(aes(x=BG4_ESC,group=BG4_ESC,fill = BG4_ESC)) + geom_bar(aes(y=freq), stat="identity") + ggtitle('Xie_2013_all_lineage_gene_H1_specific')

temp_Matrix_Xie_2013 %>% filter(Xie_2013_all_lineage_gene_ME_enrich_noH1==1) %>% group_by(BG4_ESC) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% ggplot(aes(x=BG4_ESC,group=BG4_ESC,fill = BG4_ESC)) + geom_bar(aes(y=freq), stat="identity") + ggtitle('Xie_2013_all_lineage_gene_ME_enrich_noH1')


temp_Matrix_Xie_2013 %>% filter(Xie_2013_all_lineage_gene_H1_enrich==1) %>% group_by(BG4_ESC) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% ggplot(aes(x=BG4_ESC,group=BG4_ESC,fill = BG4_ESC)) + geom_bar(aes(y=freq), stat="identity") + ggtitle('Xie_2013_all_lineage_gene_H1_enrich')

temp_Matrix_Xie_2013 %>% filter(Xie_2013_all_lineage_gene_ME_enrich==1) %>% group_by(BG4_ESC) %>% summarize(n=n()) %>% mutate(freq=n/sum(n)) %>% ggplot(aes(x=BG4_ESC,group=BG4_ESC,fill = BG4_ESC)) + geom_bar(aes(y=freq), stat="identity") + ggtitle('Xie_2013_all_lineage_gene_ME_enrich')
```

G4 motifs analysis
------------------

In order to do this, I have to extract fasta from the bed of interest. Each of this fasta need than to be shuffles (radom). Giovanni performed 3 shuffling over regions that have a potetial of G4 formation. In our case, we have open chromatin maps for stem cells: I will use that as a shuffling regions

``` bash
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2
mkdir fasta_G4_regions
fasta_hg38='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa'
genome_hg38='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome'
g='/Users/simeon01/Documents/genomes/hg38/hg38.sorted.genome'
open_chrom='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/GSM602298_ESC_FAIRE_aligned_reduced_lifted.hg38.bed'
for file in *bed
do
  sbatch --mem 4G --wrap "bedtools getfasta -fi $fasta_hg38 -bed $file -fo fasta_G4_regions/${file%%.bed}.fa"
done



cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2
n=(1 2 3) 
for f in *.bed*; 
do 
  for seed_n in ${n[@]} 
  do 
    sbatch --mem 4G --wrap "bedtools shuffle -incl $open_chrom -chromFirst -noOverlapping -seed $seed_n -i $f -g $genome_hg38 | sort -k1,1 -k2,2n > ${f%%.bed}.shuffle_${seed_n}.bed && bedtools getfasta -fi $fasta_hg38 -bed ${f%%.bed}.shuffle_${seed_n}.bed -fo ${f%%.bed}.shuffle_${seed_n}.fa"
    done
done

bedtools shuffle -incl SortBed_on_G4_seq_hits.bed -chromFirst -noOverlapping -seed ${n[i]} -i $f -g /Users/hansel01/data_seq/SLX-15831/OQS_count_analysis/hg19.genome3 | sort -k1,1 -k2,2n > ${f%%.bed}.shuffle_${n[i]}.bed; done; done

# create the set of peaks that are specific for each cell line
fasta_hg38='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa'
genome_hg38='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome'
g='/Users/simeon01/Documents/genomes/hg38/hg38.sorted.genome'
open_chrom='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/GSM602298_ESC_FAIRE_aligned_reduced_lifted.hg38.bed'
common_in_3_cellTypes=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC_NSC_CNCC_bio2out3.multi3.bed

for file in *bio2out3.bed
do 
  intersectBed -a $file -b $common_in_3_cellTypes -v > ${file%%.bed}.specific_to_cell.bed
done

# shuffle those specific peaks
n=(1 2 3) 
for f in *.specific_to_cell.bed; 
do 
  for seed_n in ${n[@]} 
  do 
    sbatch --mem 4G --wrap "bedtools getfasta -fi $fasta_hg38 -bed $f -fo fasta_G4_regions/${f%%.bed}.fa"
    sbatch --mem 4G --wrap "bedtools shuffle -incl $open_chrom -chromFirst -noOverlapping -seed $seed_n -i $f -g $genome_hg38 | sort -k1,1 -k2,2n > ${f%%.bed}.shuffle_${seed_n}.bed && bedtools getfasta -fi $fasta_hg38 -bed ${f%%.bed}.shuffle_${seed_n}.bed -fo fasta_G4_regions/${f%%.bed}.shuffle_${seed_n}.fa"
    done
done


## extract fasta for K562 jochen and Karen ## ==> folder: /scratcha/sblab/simeon01/Data/fasta_other_cells
cd /scratcha/sblab/simeon01/Data/fasta_other_cells
fasta_hg38='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa'
genome_hg38='/scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome'
g='/Users/simeon01/Documents/genomes/hg38/hg38.sorted.genome'
open_chrom='/scratchb/sblab/simeon01/20191119_Karen_K563_atac_perturbation/other_atac_seq_and_previous/50000_K562_wholecells_fragment_1hr_ATACseq_test01.all.hg38.merged.macs2_peaks_peaks.narrowPeak'
karen_path='/scratcha/sblab/simeon01/Data/20190613_Karen_continue_pol2_G4_maps/peaks_DRB/multi_bed'
k562_karen='chipseq_noDRB.multi5.bed'
jochen_path='/scratcha/sblab/simeon01/Data/Jochen_K562_peaks'
k562_jochen='20180108_K562_async_rep1-3.mult.5of8.hg19_lifted_hg38.bed'
out_dir='/scratcha/sblab/simeon01/Data/fasta_other_cells'

n=(1 2 3) 

sbatch --mem 4G --wrap "bedtools getfasta -fi $fasta_hg38 -bed $karen_path/$k562_karen -fo $out_dir/${k562_karen%%.bed}.fa"
for seed_n in ${n[@]} 
  do 
    sbatch --mem 4G --wrap "bedtools shuffle -incl $open_chrom -chromFirst -noOverlapping -seed $seed_n -i $karen_path/$k562_karen -g $genome_hg38 | sort -k1,1 -k2,2n > $out_dir/${k562_karen%%.bed}.shuffle_${seed_n}.bed && bedtools getfasta -fi $fasta_hg38 -bed  $out_dir/${k562_karen%%.bed}.shuffle_${seed_n}.bed -fo $out_dir/${k562_karen%%.bed}.shuffle_${seed_n}.fa"
done

sbatch --mem 4G --wrap "bedtools getfasta -fi $fasta_hg38 -bed $jochen_path/$k562_jochen -fo $out_dir/${k562_jochen%%.bed}.fa"
for seed_n in ${n[@]} 
  do 
    sbatch --mem 4G --wrap "bedtools shuffle -incl $open_chrom -chromFirst -noOverlapping -seed $seed_n -i $jochen_path/$k562_jochen -g $genome_hg38 | sort -k1,1 -k2,2n > $out_dir/${k562_jochen%%.bed}.shuffle_${seed_n}.bed && bedtools getfasta -fi $fasta_hg38 -bed $out_dir/${k562_jochen%%.bed}.shuffle_${seed_n}.bed -fo $out_dir/${k562_jochen%%.bed}.shuffle_${seed_n}.fa"
done
```

locally I run a script that evaluates those structural motifs in the fasta generated above

``` bash

cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/fasta_G4_regions
Rscript /Users/simeon01/Documents/code_G4seq_giovanni/sequence_hits_analysis_ChIP_C.R NSC.bio2out3.fa NSC.bio2out3.shuffle_1.fa NSC.bio2out3.shuffle_2.fa NSC.bio2out3.shuffle_3.fa

Rscript /Users/simeon01/Documents/code_G4seq_giovanni/sequence_hits_analysis_ChIP_C.R ESC.bio2out3.fa ESC.bio2out3.shuffle_1.fa ESC.bio2out3.shuffle_2.fa ESC.bio2out3.shuffle_3.fa

Rscript /Users/simeon01/Documents/code_G4seq_giovanni/sequence_hits_analysis_ChIP_C.R CNCC.bio2out3.fa CNCC.bio2out3.shuffle_1.fa CNCC.bio2out3.shuffle_2.fa CNCC.bio2out3.shuffle_3.fa

Rscript /Users/simeon01/Documents/code_G4seq_giovanni/sequence_hits_analysis_ChIP_C.R ESC_NSC_CNCC_bio2out3.multi3.fa ESC_NSC_CNCC_bio2out3.multi3.shuffle_1.fa ESC_NSC_CNCC_bio2out3.multi3.shuffle_2.fa ESC_NSC_CNCC_bio2out3.multi3.shuffle_3.fa

## check structural motifs also for cell-specific regions
Rscript /Users/simeon01/Documents/code_G4seq_giovanni/sequence_hits_analysis_ChIP_C.R ESC.bio2out3.specific_to_cell.fa ESC.bio2out3.specific_to_cell.shuffle_1.fa ESC.bio2out3.specific_to_cell.shuffle_2.fa ESC.bio2out3.specific_to_cell.shuffle_3.fa

Rscript /Users/simeon01/Documents/code_G4seq_giovanni/sequence_hits_analysis_ChIP_C.R NSC.bio2out3.specific_to_cell.fa NSC.bio2out3.specific_to_cell.shuffle_1.fa NSC.bio2out3.specific_to_cell.shuffle_2.fa NSC.bio2out3.specific_to_cell.shuffle_3.fa

Rscript /Users/simeon01/Documents/code_G4seq_giovanni/sequence_hits_analysis_ChIP_C.R CNCC.bio2out3.specific_to_cell.fa CNCC.bio2out3.specific_to_cell.shuffle_1.fa CNCC.bio2out3.specific_to_cell.shuffle_2.fa CNCC.bio2out3.specific_to_cell.shuffle_3.fa


# test hacat
cd 
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/fasta_G4_regions
Rscript /Users/simeon01/Documents/code_G4seq_giovanni/sequence_hits_analysis_ChIP_C.R HaCaT_common_peaks.hg19_lifted_To_hg38.sorted.fa HaCaT_common_peaks.hg19_lifted_To_hg38.sorted.shuffle_1.fa HaCaT_common_peaks.hg19_lifted_To_hg38.sorted.shuffle_2.fa HaCaT_common_peaks.hg19_lifted_To_hg38.sorted.shuffle_3.fa

# test other cell lines
Rscript /Users/simeon01/Documents/code_G4seq_giovanni/sequence_hits_analysis_ChIP_C.R 20180108_K562_async_rep1-3.mult.5of8.hg19_lifted_hg38.fa 20180108_K562_async_rep1-3.mult.5of8.hg19_lifted_hg38.shuffle_1.fa 20180108_K562_async_rep1-3.mult.5of8.hg19_lifted_hg38.shuffle_2.fa 20180108_K562_async_rep1-3.mult.5of8.hg19_lifted_hg38.shuffle_3.fa

Rscript /Users/simeon01/Documents/code_G4seq_giovanni/sequence_hits_analysis_ChIP_C.R chipseq_noDRB.multi5.fa chipseq_noDRB.multi5.shuffle_1.fa chipseq_noDRB.multi5.shuffle_2.fa chipseq_noDRB.multi5.shuffle_3.fa
```

generate plots after sequence motifs analysis (G4 motifs, loops etc...)
-----------------------------------------------------------------------

``` r
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/fasta_G4_regions')

piePlot <- function(count, categories) {
    dat <- data.frame(count = count, category = categories)
    dat$fraction <- dat$count / sum(dat$count)
    dat$ymax <- cumsum(dat$fraction)
    dat$ymin <- c(0, head(dat$ymax, n = -1))
    dat$label <- factor(paste(dat$category, dat$count), levels = paste(dat$category, dat$count))
    plot <-
        ggplot(dat, aes(
            fill = label, # fill by label not category
            ymax = ymax,
            ymin = ymin,
            xmin = 0,
            xmax = 1
        )) +
        geom_rect() +
        coord_polar(theta = "y") +
        theme(legend.position="top") + theme_void() # no need for labels anymore
    plot
}

generate_structure_plots <- function(dataframe_peaks,outuput_file,outuput_file_pie,outuput_file_fe, label_to_print,labels_legend){
  
  #pie chart
  pie_all_peaks <- ggplot(dataframe_peaks, aes(x = "", y = actual_count, fill = motif_name)) +
  geom_bar(width = 1, stat = "identity", 
           color = "white") +
  coord_polar("y", start = 0) + 
  scale_fill_manual(values=c("skyblue4","tomato3","darkseagreen3","mediumpurple1","lightskyblue","coral","azure3"),
                    label=labels_legend)+
  #geom_text(aes(label = actual_count), position = position_stack(vjust = 0.5)) +
  ggtitle(label_to_print) + 
  theme(legend.position ="rigth") +
  theme_void()
  
  bar_fold_enrichments <- ggplot(dataframe_peaks, 
       aes(x = motif_name, y = fold_change_over_random, fill = motif_name)) +
  geom_bar(stat = "identity", color = "white",
           width = 0.7) + 
  scale_fill_manual(values=c("skyblue4","tomato3","darkseagreen3","mediumpurple1","lightskyblue","coral","azure3"),
                    label=label_to_print)+
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.8) +
  ggtitle(paste0("ESC - fold_enrich")) + ylab('fold_enrichm') + 
  theme(legend.position = "none",
        text = element_text(size=12),
        #labels=labels_legend,
        axis.text.x = element_text(angle = 90))
  
  Final_plot <- grid.arrange(pie_all_peaks,bar_fold_enrichments)
  ggsave(file=outuput_file,Final_plot)
  ggsave(file=outuput_file_pie,pie_all_peaks)
  ggsave(file=outuput_file_fe,bar_fold_enrichments)
  rm(Final_plot)
}

assign('ESC_all_peaks_motifs', get(load('ESC.bio2out3.results_motifs.Rdata')))
ESC_all_peaks_motifs <- ESC_all_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
assign('NSC_all_peaks_motifs', get(load('NSC.bio2out3.results_motifs.Rdata')))
NSC_all_peaks_motifs <- NSC_all_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
assign('CNCC_all_peaks_motifs', get(load('CNCC.bio2out3.results_motifs.Rdata')))
CNCC_all_peaks_motifs <- CNCC_all_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))

assign('ESC_NSC_CNCC_common_peaks_motifs', get(load('ESC_NSC_CNCC_bio2out3.multi3.results_motifs.Rdata')))
ESC_NSC_CNCC_common_peaks_motifs <- ESC_NSC_CNCC_common_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))

assign('ESC_specific_peaks_motifs', get(load('ESC.bio2out3.specific_to_cell.results_motifs.Rdata')))
ESC_specific_peaks_motifs <- ESC_specific_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
assign('NSC_specific_peaks_motifs', get(load('NSC.bio2out3.specific_to_cell.results_motifs.Rdata')))
NSC_specific_peaks_motifs <- NSC_specific_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
assign('CNCC_specific_peaks_motifs', get(load('CNCC.bio2out3.specific_to_cell.results_motifs.Rdata')))
CNCC_specific_peaks_motifs <- CNCC_specific_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
# plot: counts_3, counts_5 - counts_3, counts_7 - counts_5, counts_long, counts_sbulges, counts_cbulges_GG, counts_other

#piePlot(HaCaT$actual_count,categories = HaCaT$motif_name)

#piePlot(round(ESC_all_peaks_motifs$actual_count/sum(ESC_all_peaks_motifs$actual_count)*100,2),categories = ESC_all_peaks_motifs$motif_name)

#piePlot(NSC_all_peaks_motifs$actual_count,categories = HaCaT$motif_name)

#piePlot(CNCC_all_peaks_motifs$actual_count,categories = HaCaT$motif_name)

labels_legend_ESC <- paste0(ESC_all_peaks_motifs$actual_count,' (',ESC_all_peaks_motifs$percentage,'%)')
labels_legend_NSC <- paste0(NSC_all_peaks_motifs$actual_count,' (',NSC_all_peaks_motifs$percentage,'%)')
labels_legend_CNCC <- paste0(CNCC_all_peaks_motifs$actual_count,' (',CNCC_all_peaks_motifs$percentage,'%)')
labels_legend_ESC_NSC_CNCC_common <- paste0(ESC_NSC_CNCC_common_peaks_motifs$actual_count,' (',ESC_NSC_CNCC_common_peaks_motifs$percentage,'%)')
labels_legend_ESC_specific <- paste0(ESC_specific_peaks_motifs$actual_count,' (',ESC_specific_peaks_motifs$percentage,'%)')
labels_legend_NSC_specific <- paste0(NSC_specific_peaks_motifs$actual_count,' (',NSC_specific_peaks_motifs$percentage,'%)')
labels_legend_CNCC_specific <- paste0(CNCC_specific_peaks_motifs$actual_count,' (',CNCC_specific_peaks_motifs$percentage,'%)')


generate_structure_plots(ESC_all_peaks_motifs,
                         'ESC_all_peaks_motifs.pdf',
                         'ESC_all_peaks_motifs_pie_counts.pdf',
                         'ESC_all_peaks_motifs_fold_enrich.pdf',
                         'ESC_all_peaks',labels_legend_ESC)

generate_structure_plots(NSC_all_peaks_motifs,
                         'NSC_all_peaks_motifs.pdf',
                         'NSC_all_peaks_motifs_pie_counts.pdf',
                         'NSC_all_peaks_motifs_fold_enrich.pdf',
                         'NSC_all_peaks',labels_legend_NSC)

generate_structure_plots(CNCC_all_peaks_motifs,
                         'CNCC_all_peaks_motifs.pdf',
                         'CNCC_all_peaks_motifs_pie_counts.pdf',
                         'CNCC_all_peaks_motifs_fold_enrich.pdf',
                         'CNCC_all_peaks',labels_legend_CNCC)

generate_structure_plots(ESC_NSC_CNCC_common_peaks_motifs,
                         'ESC_NSC_CNCC_common_peaks_motifs.pdf',
                         'ESC_NSC_CNCC_common_peaks_motifs_pie_counts.pdf',
                         'ESC_NSC_CNCC_common_peaks_motifs_fold_enrich.pdf',
                         'ESC_NSC_CNCC_common_peaks',labels_legend_ESC_NSC_CNCC_common)

generate_structure_plots(ESC_specific_peaks_motifs,
                         'ESC_specific_peaks_motifs.pdf',
                         'ESC_specific_peaks_motifs_pie_counts.pdf',
                         'ESC_specific_peaks_motifs_fold_enrich.pdf',
                         'ESC_specific_peaks',labels_legend_ESC_specific)

generate_structure_plots(NSC_specific_peaks_motifs,
                         'NSC_specific_peaks_motifs.pdf',
                         'NSC_specific_peaks_motifs_pie_counts.pdf',
                         'NSC_specific_peaks_motifs_fold_enrich.pdf',
                         'NSC_specific_peaks',labels_legend_NSC_specific)

generate_structure_plots(CNCC_specific_peaks_motifs,
                         'CNCC_specific_peaks_motifs.pdf',
                         'CNCC_specific_peaks_motifs_pie_counts.pdf',
                         'CNCC_specific_peaks_motifs_fold_enrich.pdf',
                         'CNCC_specific_peaks',labels_legend_CNCC_specific)

#other cell lines
#HaCaT case
assign('HaCaT', get(load('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/fasta_G4_regions/HaCaT_common_peaks.hg19_lifted_To_hg38.sorted.results_motifs.Rdata')))
HaCaT <- HaCaT %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
labels_legend_HaCaT <- paste0(HaCaT$actual_count,' (',HaCaT$percentage,'%)')

generate_structure_plots(HaCaT,
                         'HaCaT_peaks_motifs.pdf',
                         'HaCaT_peaks_motifs_pie_counts.pdf',
                         'HaCaT_peaks_motifs_fold_enrich.pdf',
                         'HaCaT_peaks',labels_legend_HaCaT)


#K562 cell line
assign('K562_karen', get(load('chipseq_noDRB.multi5.results_motifs.Rdata')))
K562_karen <- K562_karen %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
labels_legend_K562_karen <- paste0(K562_karen$actual_count,' (',K562_karen$percentage,'%)')
generate_structure_plots(K562_karen,
                         'K562_karen_peaks_motifs.pdf',
                         'K562_karen_peaks_motifs_pie_counts.pdf',
                         'K562_karen_peaks_motifs_fold_enrich.pdf',
                         'K562_karen_peaks',labels_legend_K562_karen)

assign('K562_jochen', get(load('20180108_K562_async_rep1-3.mult.5of8.hg19_lifted_hg38.results_motifs.Rdata')))
K562_jochen <- K562_jochen %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
labels_legend_K562_jochen <- paste0(K562_jochen$actual_count,' (',K562_jochen$percentage,'%)')
generate_structure_plots(K562_jochen,
                         'K562_jochen_peaks_motifs.pdf',
                         'K562_jochen_peaks_motifs_pie_counts.pdf',
                         'K562_jochen_peaks_motifs_fold_enrich.pdf',
                         'K562_jochen_peaks',labels_legend_K562_jochen)
```

density plots at Katie's curated gene lists and at COSMIC gene lists
--------------------------------------------------------------------

In those casees I compare them to the GW distribution

``` bash

#locally
cd /Users/simeon01/Documents/Katie/Gene_lists_for_Angela
scp *TSScoord.bed simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Katie_curated_gene_list_and_cosmic_lists/

# on cluster
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Katie_curated_gene_list_and_cosmic_lists
output_dir=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Katie_curated_gene_list_and_cosmic_lists

#full list of TSS
file_gw_tss=/scratcha/sblab/simeon01/reference_genomes/hg38/gencode.v28.TSS.bed
# extract 10000 case from genome-widee
shuf -n 10000 ${file_gw_tss}| sortBed -i - | slopBed -i - -l 1 -r 0 -g /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome > $output_dir/gencode.v28.TSS.shuf10000.sorted.touse.bed

#preepare bed files to use for density plots
for file in *TSScoord*bed*
  do
  sortBed -i $file | slopBed -i - -l 1 -r 0 -g /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome > ${file%%.bed}.sorted.touse.bed
  done

#path with bed files  
Bam_repA=/scratcha/sblab/simeon01/Data/190829_Katie_stemcells/SLX-18496/aligned #.hg38.sort.markduplicates.bam
Bam_repB=/scratcha/sblab/simeon01/Data/20190913_SLX18495_Katie_StemCell_Rep3_3xBG4/aligned/merged # markduplicates.bam
Bam_repC=/scratcha/sblab/simeon01/Data/20190917_SLX-18355_Katie_StemCells_Rep2_3xBG4/SLX-18355/aligned #.hg38.sort.markduplicates.bam

for file in $Bam_repA/*bw
do
  filename=${file##*/}
  echo $filename
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $output_dir/*cosmic*.touse.bed $output_dir/gencode.v28.TSS.shuf10000.sorted.touse.bed   -bs 20 -a 1500 -b 1500 -o $output_dir/${filename%%.bw}.cosmic_vs_gw.gz && plotProfile --matrixFile $output_dir/${filename%%.bw}.cosmic_vs_gw.gz -o $output_dir/${filename%%.bw}.cosmic_vs_gw.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $output_dir/${filename%%.bw}.cosmic_vs_gw.gz -o $output_dir/${filename%%.bw}.cosmic_vs_gw.heatmap.png --averageType median"
  echo " ==== "
done

for file in $Bam_repA/*bw
do
  filename=${file##*/}
  echo $filename
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $output_dir/*cosmic*.touse.bed $output_dir/gencode.v28.TSS.shuf10000.sorted.touse.bed   -bs 20 -a 1500 -b 1500 -o $output_dir/${filename%%.bw}.cosmic_vs_gw.gz && plotProfile --matrixFile $output_dir/${filename%%.bw}.cosmic_vs_gw.gz -o $output_dir/${filename%%.bw}.cosmic_vs_gw.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $output_dir/${filename%%.bw}.cosmic_vs_gw.gz -o $output_dir/${filename%%.bw}.cosmic_vs_gw.heatmap.png --averageType median"
  echo " ==== "
done

for file in $Bam_repA/*bw
do
  filename=${file##*/}
  echo $filename
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $output_dir/*Chang_2011*.touse.bed $output_dir/gencode.v28.TSS.shuf10000.sorted.touse.bed   -bs 20 -a 1500 -b 1500 -o $output_dir/${filename%%.bw}.Chang_2011_vs_gw.gz && plotProfile --matrixFile $output_dir/${filename%%.bw}.Chang_2011_vs_gw.gz -o $output_dir/${filename%%.bw}.Chang_2011.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $output_dir/${filename%%.bw}.Chang_2011.gz -o $output_dir/${filename%%.bw}.cosmic_vs_gw.heatmap.png --averageType median"
  
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $output_dir/*Xie_2013_all_lineage_gene_H1*.touse.bed $output_dir/gencode.v28.TSS.shuf10000.sorted.touse.bed   -bs 20 -a 1500 -b 1500 -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_H1.gz && plotProfile --matrixFile $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_H1.gz -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_H1.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_H1.gz -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_H1.heatmap.png --averageType median"
  
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $output_dir/*Xie_2013_all_lineage_gene_ME*.touse.bed $output_dir/gencode.v28.TSS.shuf10000.sorted.touse.bed   -bs 20 -a 1500 -b 1500 -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_ME.gz && plotProfile --matrixFile $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_ME.gz -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_ME.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_ME.gz -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_ME.heatmap.png --averageType median"
  
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $output_dir/*Xie_2013_all_lineage_gene_MSC_enrich*.touse.bed $output_dir/gencode.v28.TSS.shuf10000.sorted.touse.bed   -bs 20 -a 1500 -b 1500 -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_MSC_enrich.gz && plotProfile --matrixFile $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_MSC_enrich.gz -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_MSC_enrich.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_MSC_enrich.gz -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_MSC_enrich.heatmap.png --averageType median"
  
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $output_dir/*Xie_2013_all_lineage_gene_noIMR90*.touse.bed $output_dir/gencode.v28.TSS.shuf10000.sorted.touse.bed   -bs 20 -a 1500 -b 1500 -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_noIMR90.gz && plotProfile --matrixFile $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_noIMR90.gz -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_noIMR90.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_noIMR90.gz -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_noIMR90.heatmap.png --averageType median"
  
  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $output_dir/*Xie_2013_all_lineage_gene_NPC_enrich*.touse.bed $output_dir/gencode.v28.TSS.shuf10000.sorted.touse.bed   -bs 20 -a 1500 -b 1500 -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_NPC_enrich.gz && plotProfile --matrixFile $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_NPC_enrich.gz -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_NPC_enrich.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_NPC_enrich.gz -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_NPC_enrich.heatmap.png --averageType median"
  

  sbatch --mem 12G --wrap="computeMatrix reference-point -S $file -R $output_dir/*Xie_2013_all_lineage_gene_TBL_enrich*.touse.bed $output_dir/gencode.v28.TSS.shuf10000.sorted.touse.bed   -bs 20 -a 1500 -b 1500 -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_TBL_enrich.gz && plotProfile --matrixFile $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_TBL_enrich.gz -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_TBL_enrich.png --averageType mean --plotWidth 18 --plotHeight 10 --yMin 0 && plotHeatmap --matrixFile $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_TBL_enrich.gz -o $output_dir/${filename%%.bw}.Xie_2013_all_lineage_gene_TBL_enrich.heatmap.png --averageType median"
  echo " ==== "
  
  
done
```

Check single matrix produced - from deeptools and generate density plot within R
--------------------------------------------------------------------------------

``` r
setwd("/Users/simeon01/Documents/Katie/Gene_lists_for_Angela/matrices_densities")
list_comparison <- c('cosmic_vs_gw','Xie_2013_all_lineage_gene_H1',
                    'Xie_2013_all_lineage_gene_ME',
                    'Xie_2013_all_lineage_gene_MSC_enric',
                    'Xie_2013_all_lineage_gene_noIMR90',
                    'Xie_2013_all_lineage_gene_NPC_enrich',
                    'Xie_2013_all_lineage_gene_TBL_enrich')
list_cells <-c('ESC','NSC','CNCC')
for (i in (1:length(list_cells)){
  for (j in (1:length(list_comparison))){
    files_interest <- grep(list_comparison[j],list.files(pattern=c(list_cells[i])), value = T)
    big_matrix_minus1 <- c()
    big_matrix_plus1 <- c()
    big_matrix_zeros <- c()
    for (k in 1:length(files_interest)){
      temp_Matrix <- read.table(files_interest[k],stringsAsFactors = F, skip = 1)
      Matrix_values <- temp_Matrix[,7:dim(temp_Matrix)[2]]
      temp_boundaties <- read.table(files_interest[k],stringsAsFactors = F,nrows=1)
  
  ESC_big_matrix_minus1 <- rbind(ESC_big_matrix_minus1,Matrix_values[1:1521,])
  ESC_big_matrix_plus1 <- rbind(ESC_big_matrix_plus1,Matrix_values[(1521+1):(1521+491),])
  ESC_big_matrix_zeros <- rbind(ESC_big_matrix_zeros,Matrix_values[(1521+491+1):7012,])
}


y_minus1 <-colMedians(as.matrix(ESC_big_matrix_minus1),na.rm = T)
y_plus1 <-colMedians(as.matrix(ESC_big_matrix_plus1),na.rm = T)
y_zeros <-colMedians(as.matrix(ESC_big_matrix_zeros),na.rm = T)

y_minus1s <- loess.smooth(x=seq(1:length(y_minus1)),y_minus1,span = 0.15)
y_minus1s_a <- approx(y_minus1s$x, y_minus1s$y,seq(1,150))$y

y_plus1s <- loess.smooth(x=seq(1:length(y_plus1)),y_plus1,span = 0.15)
y_plus1s_a <- approx(y_plus1s$x, y_plus1s$y,seq(1,150))$y

y_zeross <- loess.smooth(x=seq(1:length(y_zeros)),y_zeros,span = 0.15)
y_zeross_a <- approx(y_zeross$x, y_zeross$y,seq(1,150))$y
ESCbg4_DEA_CNCC_ESC <-data.frame(median_norm_BG4=c(y_minus1s_a,
                                 y_plus1s_a,
                                 y_zeross_a),
                                 TSS_x <- rep(seq(-1500,(1500-20),20),3),
                                 case_DE = c(rep('minus1',150),rep('plus1',150),rep('zero',150)))

pdf("ESCbg4_DE_CNCCvsESC.pdf")
plot <- ESCbg4_DEA_CNCC_ESC %>% ggplot(aes(x=TSS_x,y=median_norm_BG4,group = case_DE,color=case_DE)) + geom_line(color=c(rep('blue',150),rep('cyan',150),rep('yellow',150)),size=0.7) + theme_bw() + ggtitle('ESCbg4 - CNCCvsESC')
print(plot)
dev.off()
```

additional annotation for stem cells - provided by Katies
---------------------------------------------------------

Various bed files have been copied at:

``` bash

cd /scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/additional_annotations_various_stem_cells

# for hiC data generate additional individual files for each of the A and B compartment
awk '{print "chr"$1"\t"$2"\t"$3}' HiC_matrix_Compartments_coordinated.bed > HiC_matrix_Compartments_coordinated_chr.bed
tail -n +2 Selected_cells_HiC_compartments_H1_NPC.txt | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' > ESC_HiC_data_Schmitt_2016_hg19.bed
tail -n +2 Selected_cells_HiC_compartments_H1_NPC.txt | awk '{print "chr"$1"\t"$2"\t"$3"\t"$5}' > NPC_HiC_data_Schmitt_2016_hg19.bed
for file in *HiC_data_Schmitt_2016_hg19.bed
do
  echo  $file
  awk '{print $1"\t"$2"\t"$3"\t"$4}' $file | grep -w "A" > ${file/2016_hg19.bed/Acomp_hg19.bed}
  awk '{print $1"\t"$2"\t"$3"\t"$4}' $file | grep -w "B" > ${file/2016_hg19.bed/Bcomp_hg19.bed}

done

#for tad borders extract the single base that limit the tad itself
for file in *TAD_Schmitt_2016_hg19.bed
do
awk '{print $1"\t"$2"\t"$2+1}' $file >  temp1.tad.bed
awk '{print $1"\t"$3"\t"$3+1}' $file >  temp2.tad.bed

sortBed -i $file | awk '{print $1"\t"$2"\t"$3"\tTAD"}' > ${file/TAD_Schmitt_2016_hg19.bed/TAD_Schmitt_2016_final_hg19.bed}

cat temp1.tad.bed temp2.tad.bed | sortBed -i - | awk '{print $0"\tTADborder1bp"}'> ${file/TAD_Schmitt_2016_hg19.bed/TADborder_Schmitt_2016_hg19.bed}
rm temp1.tad.bed temp2.tad.bed
done

#lifover chains
chain_hg19_to_hg38=/scratcha/sblab/simeon01/reference_genomes/liftover_chain/hg19ToHg38.over.chain.gz
chain_hg18_to_hg38=/scratcha/sblab/simeon01/reference_genomes/liftover_chain/hg18ToHg38.over.chain.gz

# lifover from hg19 to hg38
/Users/simeon01/applications/liftOver pnas.1613300114.sd02_hESC_eexported_rearranged.bed $chain_hg19_to_hg38 pnas.1613300114.sd02_hESC_eexported_rearranged.hg38.bed  pnas.1613300114.sd02_hESC_eexported_rearranged.unmapped

for file in *hg19.bed
do
/Users/simeon01/applications/liftOver $file $chain_hg19_to_hg38 ${file/hg19.bed/hg19_lifted_hg38.hg38.bed} temp.unmapped
done



# liftover from hg18 to hg38
for file in *hg18.sorted.bed
do
# slop them by few bases
bedtools slop -i $file -g hg18.chrom.sizes -b 2 | sortBed -i - > ${file/hg18.sorted.bed/hg18.slop2bases.sorted.bed}
/Users/simeon01/applications/liftOver ${file/hg18.sorted.bed/hg18.slop2bases.sorted.bed} $chain_hg18_to_hg38 ${file/hg18.sorted.bed/hg18_lifted_hg38.hg38.bed} temp.unmapped
done


final_dir=/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/additional_annotations_various_stem_cells/Annotations_for_GAT

mkdir $final_dir

# copy all files of interest to the destination folder
#ESC_HiC_data_Schmitt_2016_hg19_lifted_hg38.hg38.bed
#ESC_HiC_data_Schmitt_Acomp_hg19_lifted_hg38.hg38.bed
#ESC_HiC_data_Schmitt_Bcomp_hg19_lifted_hg38.hg38.bed
#ESC_HiC_data_TADborder_Schmitt_2016_hg19_lifted_hg38.hg38.bed
#ESC_HiC_data_TAD_Schmitt_2016_hg19_lifted_hg38.hg38.bed
#ESC_NPC.enhancerAtlas.hg19_lifted_hg38.hg38.bed
#H1ESC_Xie_2013_all_ehancers_hg18_lifted_hg38.hg38.bed
#H1ESC_Xie_2013_enriched_ehancers_hg18_lifted_hg38.hg38.bed
#H1_Hniscz2013_superenhancers_hg19_lifted_hg38.hg38.bed
#H1NPC_Xie_2013_all_ehancers_hg18_lifted_hg38.hg38.bed
#H1NPC_Xie_2013_enriched_ehancers_hg18_lifted_hg38.hg38.bed
#NPC_HiC_data_Schmitt_2016_hg19_lifted_hg38.hg38.bed
#H9.enhancerAtlas.hg19_lifted_hg38.hg38.bed
#hNCC.enhancerAtlas.hg19_lifted_hg38.hg38.bed
#NPC_HiC_data_Schmitt_Acomp_hg19_lifted_hg38.hg38.bed
#NPC_HiC_data_Schmitt_Bcomp_hg19_lifted_hg38.hg38.bed
#NPC_HiC_data_TADborder_Schmitt_2016_hg19_lifted_hg38.hg38.bed
#NPC_HiC_data_TAD_Schmitt_2016_hg19_lifted_hg38.hg38.bed
#Prescott_enhancers_6000_CNCC_hg19_lifted_hg38.hg38.bed



cp ESC_HiC_data_Schmitt_Acomp_hg19_lifted_hg38.hg38.bed $final_dir
cp ESC_HiC_data_Schmitt_Bcomp_hg19_lifted_hg38.hg38.bed $final_dir
cp ESC_HiC_data_TADborder_Schmitt_2016_hg19_lifted_hg38.hg38.bed $final_dir


awk '{print $1"\t"$2"\t"$3"\tNPC.enhancerAtlas"}' ESC_NPC.enhancerAtlas.hg19_lifted_hg38.hg38.bed > $final_dir/ESC_NPC.enhancerAtlas.hg19_lifted_hg38.hg38.bed

awk '{print $1"\t"$2"\t"$3"\tH1ESC_Xie_2013_all_ehancers"}' H1ESC_Xie_2013_all_ehancers_hg18_lifted_hg38.hg38.bed > $final_dir/H1ESC_Xie_2013_all_ehancers_hg18_lifted_hg38.hg38.bed

awk '{print $1"\t"$2"\t"$3"\tH1ESC_Xie_2013_enriched_ehancer"}' H1ESC_Xie_2013_enriched_ehancers_hg18_lifted_hg38.hg38.bed > $final_dir/H1ESC_Xie_2013_enriched_ehancers_hg18_lifted_hg38.hg38.bed

awk '{print $1"\t"$2"\t"$3"\tH1_Hniscz2013_superenhancers"}' H1_Hniscz2013_superenhancers_hg19_lifted_hg38.hg38.bed > $final_dir/H1_Hniscz2013_superenhancers_hg19_lifted_hg38.hg38.bed

awk '{print $1"\t"$2"\t"$3"\tH1NPC_Xie_2013_all_ehancers"}' H1NPC_Xie_2013_all_ehancers_hg18_lifted_hg38.hg38.bed > $final_dir/H1NPC_Xie_2013_all_ehancers_hg18_lifted_hg38.hg38.bed

awk '{print $1"\t"$2"\t"$3"\tH1NPC_Xie_2013_enriched_ehancers"}' H1NPC_Xie_2013_enriched_ehancers_hg18_lifted_hg38.hg38.bed > $final_dir/H1NPC_Xie_2013_enriched_ehancers_hg18_lifted_hg38.hg38.bed

awk '{print $1"\t"$2"\t"$3"\tH9_ESC.enhancerAtlas"}' H9.enhancerAtlas.hg19_lifted_hg38.hg38.bed > $final_dir/H9.enhancerAtlas.hg19_lifted_hg38.hg38.bed

awk '{print $1"\t"$2"\t"$3"\thNCC.enhancerAtlas"}' hNCC.enhancerAtlas.hg19_lifted_hg38.hg38.bed > $final_dir/hNCC.enhancerAtlas.hg19_lifted_hg38.hg38.bed


cp NPC_HiC_data_Schmitt_Acomp_hg19_lifted_hg38.hg38.bed $final_dir/NPC_HiC_data_Schmitt_Acomp_hg19_lifted_hg38.hg38.bed

cp NPC_HiC_data_Schmitt_Bcomp_hg19_lifted_hg38.hg38.bed $final_dir/NPC_HiC_data_Schmitt_Bcomp_hg19_lifted_hg38.hg38.bed

cp NPC_HiC_data_TADborder_Schmitt_2016_hg19_lifted_hg38.hg38.bed  $final_dir/NPC_HiC_data_TADborder_Schmitt_2016_hg19_lifted_hg38.hg38.bed



awk '{print $1"\t"$2"\t"$3"\tPrescott_enhancers_6000_CNCC"}' Prescott_enhancers_6000_CNCC_hg19_lifted_hg38.hg38.bed > $final_dir/Prescott_enhancers_6000_CNCC_hg19_lifted_hg38.hg38.bed

# For each bed file (all, specific and intersection) run gat analysis over each of the annotations


# perform fold enrichmeent analysis on any bed file in each subfolder
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed
path_anno='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/additional_annotations_various_stem_cells/Annotations_for_GAT'
path_output='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/additional_annotations_various_stem_cells/Annotations_for_GAT/output_enrichments_gat'
mkdir $path_output
scenario=(scenario1_multi2)
#genomic_annotation=annotations_genomic_regions_gencode_v28.anno.bed
annotation_list1=`ls $path_anno/*.hg38.bed`
#scenario=(scenario3_multi2_extended_1kb)
for case in ${scenario[@]}
do
  cd ./$case
  # run gat analysis on full bed files 
  for file in `ls *bed| grep -v shuffle`
  do
    echo "===>> " $file "<<==="
    for anno in ${annotation_list1[@]}
    do
      echo $anno
      echo "===="
      base_annotation=`basename $anno`
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=${anno} --segments=${file} -S $path_output/${file%%.bed}_intervals_${base_annotation%%.bed}.dat -t 4"
      echo ${file%%.bed}_intervals_${base_annotation%%.bed}.dat
      echo $cmd_gat
      sbatch --mem 6G --wrap "$cmd_gat"
      echo " "
    done
  done

  #run the enrichment for subsets (common and specific ones)
  cd ./bed_intersections_common
  pwd
  for file in *bed
  do
    echo $file
    for anno in ${annotation_list1[@]}
    do
      base_annotation=`basename $anno`
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=${anno} --segments=${file} -S $path_output/${file%%.bed}_intervals_${base_annotation%%.bed}.dat -t 4"
      echo ${file%%.bed}_intervals_${base_annotation%%.bed}.dat
      pwd
      echo $cmd_gat
      sbatch --mem 6G --wrap "$cmd_gat"
      echo " "
      
    done
  done
  echo " ** == ** =="
  echo " ++++++++++ END BLOCK +++++++++++++++++"
  cd ../..
done






annotation_list1=/scratcha/sblab/simeon01/reference_genomes/hg38/annotations_genomic_regions_gencode_v28.anno2.bed
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
for case in ${scenario[@]}
do
  cd ./$case
  # run gat analysis on full bed files 
  for file in `ls *bed| grep -v shuffle`
  do
    echo "===>> " $file "<<==="
    for anno in ${annotation_list1[@]}
    do
      echo $anno
      echo "===="
      base_annotation=`basename $anno`
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=${anno} --segments=${file} -S $path_output/${file%%.bed}_intervals_${base_annotation%%.bed}.dat -t 4"
      echo ${file%%.bed}_intervals_${base_annotation%%.bed}.dat
      echo $cmd_gat
      sbatch --mem 6G --wrap "$cmd_gat"
      echo " "
    done
  done

  #run the enrichment for subsets (common and specific ones)
  cd ./bed_intersections_common
  pwd
  for file in *bed
  do
    echo $file
    for anno in ${annotation_list1[@]}
    do
      base_annotation=`basename $anno`
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=${anno} --segments=${file} -S $path_output/${file%%.bed}_intervals_${base_annotation%%.bed}.dat -t 4"
      echo ${file%%.bed}_intervals_${base_annotation%%.bed}.dat
      pwd
      echo $cmd_gat
      sbatch --mem 6G --wrap "$cmd_gat"
      echo " "
      
    done
  done
  echo " ** == ** =="
  echo " ++++++++++ END BLOCK +++++++++++++++++"
  cd ../..
done
```

Cross comparison between Differntial expression and differential binding analysis on all G4 peaks
-------------------------------------------------------------------------------------------------

Prepare data

``` bash

cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/coverage_peaks

for file in EdgeRtable_differential_binding_analysis*csv
do
  sed $'s/_/\t/g' $file | tail -n +2 | sed $'s/,/\t/g' |intersectBed -a - -b /Users/simeon01/Documents/genomes/hg38/gencode.v28.annotation.gene_bed6_1kbaroundTSS.bed -wa -wb | sed $'s/"//g'| sed $'s/;//g'  > $file.peak_in_gencode.v28.annotation.gene_bed6_1kbaroundTSS.bed
done
```

``` r
# import files with promoters with 
path_diff_analysis <- '/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/coverage_peaks'
list_output_DiffBind_peaks <- list.files(path=path_diff_analysis,pattern= 'peak_in_gencode.v28.annotation.gene_bed6_1kbaroundTSS.bed')
setwd(path_diff_analysis)
M_edgeR_differential_analysis <- read.table(list_output_DiffBind_peaks[1])[,c(1,2,3)]
M_edgeR_differential_analysis$peak_id <- paste(M_edgeR_differential_analysis[,1],M_edgeR_differential_analysis[,2],M_edgeR_differential_analysis[,3],sep="_")
for (i in (1:length(list_output_DiffBind_peaks))){
  temp <- read.table(list_output_DiffBind_peaks[i])
  colnames(temp) <- c('chr','start','end',"logFC","logCPM","F","PValue","FDR","chr_prom","start_prom","end_prom","par1","par2", "par3","Ens_id",    "typeProt")
  temp$peak_id <- paste(temp[,1],temp[,2],temp[,3],sep="_")
  M_edgeR_differential_analysis <- left_join(M_edgeR_differential_analysis,temp,by="peak_id")
}

M_diffBind_CNCC_vs_ESC <- read.table(list_output_DiffBind_peaks[1],stringsAsFactors = F)
M_diffBind_NSC_vs_ESC <- read.table(list_output_DiffBind_peaks[2],stringsAsFactors = F)
M_diffBind_NSC_vs_CNCC <- read.table(list_output_DiffBind_peaks[3],stringsAsFactors = F)
colnames(M_diffBind_CNCC_vs_ESC) <- colnames(M_diffBind_NSC_vs_ESC) <- colnames(M_diffBind_NSC_vs_CNCC) <- c('chr','start','end',"logFC","logCPM","F","PValue","FDR","chr_prom","start_prom","end_prom","par1","par2",  "par3","Ens_id",    "typeProt")
M_diffBind_CNCC_vs_ESC$ens_id <- gsub('\\..*','',M_diffBind_CNCC_vs_ESC$Ens_id)
M_diffBind_NSC_vs_ESC$ens_id <- gsub('\\..*','',M_diffBind_NSC_vs_ESC$Ens_id)
M_diffBind_NSC_vs_CNCC$ens_id <- gsub('\\..*','',M_diffBind_NSC_vs_CNCC$Ens_id)

M_diffBind_CNCC_vs_ESC_reduced <- M_diffBind_CNCC_vs_ESC %>% dplyr::group_by(ens_id) %>% dplyr::summarize(min_logFC = min(logFC),minFDR = min(FDR))  
M_diffBind_NSC_vs_ESC_reduced <- M_diffBind_NSC_vs_ESC %>% dplyr::group_by(ens_id) %>% dplyr::summarize(min_logFC = min(logFC),minFDR = min(FDR))  
M_diffBind_NSC_vs_CNCC_reduced <- M_diffBind_NSC_vs_CNCC %>% dplyr::group_by(ens_id) %>% dplyr::summarize(min_logFC = min(logFC),minFDR = min(FDR))  


# import expression from differntial expression analysis
# import output of differential expression analysis keeping attention to the fact that the labels could be reverted
DEA_outcome_NSC_vs_ESC <- read.delim('/Users/simeon01/Documents/Katie/2018_09_03/RNA_seq/summary_edgeR_H9_CNCC_vs_H9_ESC.txt',sep = '\t') # these labels are reverted!!! 1: ESC, -1 NSC
DEA_outcome_NSC_vs_CNCC <- read.delim('/Users/simeon01/Documents/Katie/2018_09_03/RNA_seq/summary_edgeR_H9_NSC_vs_H9_CNCC.txt',sep = '\t')
DEA_outcome_CNCC_vs_ESC <- read.delim('/Users/simeon01/Documents/Katie/2018_09_03/RNA_seq/summary_edgeR_H9_NSC_vs_H9_hESC.txt',sep = '\t')

#colnames(DEA_outcome_NSC_vs_ESC) <- c('X','DEA_NSC_vs_ESC')
#colnames(DEA_outcome_CNCC_vs_ESC) <- c('X','DEA_CNCC_vs_ESC')
#colnames(DEA_outcome_NSC_vs_CNCC) <- c('X','DEA_NSC_vs_CNCC') 

# additional filtering on the DEA genes
# remove row 1
DEA_outcome_NSC_vs_ESC <- DEA_outcome_NSC_vs_ESC[-1,]
DEA_outcome_NSC_vs_CNCC <- DEA_outcome_NSC_vs_CNCC[-1,]
DEA_outcome_CNCC_vs_ESC <- DEA_outcome_CNCC_vs_ESC[-1,]
# remove gene symbols
DEA_outcome_NSC_vs_ESC$ens_id <- gsub('\\..*','',DEA_outcome_NSC_vs_ESC$X)
DEA_outcome_NSC_vs_CNCC$ens_id <- gsub('\\..*','',DEA_outcome_NSC_vs_CNCC$X)
DEA_outcome_CNCC_vs_ESC$ens_id <- gsub('\\..*','',DEA_outcome_CNCC_vs_ESC$X)

DEA_outcome_NSC_vs_ESC$symb <- gsub('.*._','',DEA_outcome_NSC_vs_ESC$X)
DEA_outcome_NSC_vs_CNCC$symb <- gsub('.*._','',DEA_outcome_NSC_vs_CNCC$X)
DEA_outcome_CNCC_vs_ESC$symb <- gsub('.*._','',DEA_outcome_CNCC_vs_ESC$X)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

DEA_outcome_NSC_vs_ESC_reduced <- DEA_outcome_NSC_vs_ESC %>% dplyr::group_by(ens_id) %>% dplyr::summarise(mean_logFC=mean(logFC),mean_logCPM=mean(logCPM),comb_PValue=gm_mean(PValue),comb_FDR=gm_mean(FDR))

DEA_outcome_NSC_vs_CNCC_reduced <- DEA_outcome_NSC_vs_CNCC %>% dplyr::group_by(ens_id) %>% dplyr::summarise(mean_logFC=mean(logFC),mean_logCPM=mean(logCPM),comb_PValue=gm_mean(PValue),comb_FDR=gm_mean(FDR))

DEA_outcome_CNCC_vs_ESC_reduced <- DEA_outcome_CNCC_vs_ESC %>% dplyr::group_by(ens_id) %>% dplyr::summarise(mean_logFC=mean(logFC),mean_logCPM=mean(logCPM),comb_PValue=gm_mean(PValue),comb_FDR=gm_mean(FDR))



# left join Differential expression data and differential binding
DiffExpr_DiffBind_NSC_ESC <- left_join(DEA_outcome_NSC_vs_ESC_reduced,M_diffBind_NSC_vs_ESC_reduced,by="ens_id")

DiffExpr_DiffBind_CNCC_ESC <- left_join(DEA_outcome_CNCC_vs_ESC_reduced,M_diffBind_CNCC_vs_ESC_reduced,by="ens_id")

DiffExpr_DiffBind_NSC_vs_CNCC <- left_join(DEA_outcome_NSC_vs_CNCC_reduced,M_diffBind_NSC_vs_CNCC_reduced,by="ens_id")


DiffExpr_DiffBind_NSC_ESC <- DiffExpr_DiffBind_NSC_ESC %>% 
  mutate(DiffExpr = ifelse(mean_logFC >2 & comb_FDR < 0.05, "NSC",
                    ifelse(mean_logFC < -2 & comb_FDR < 0.05 ,"ESC","NOTDE")), 
         DiffBind = ifelse(min_logFC>1 & minFDR<=0.05,"NSCbg4",
                       ifelse(min_logFC< -1 & minFDR<=0.05,"ESCbg4","NOTDBbg4")))

DiffExpr_DiffBind_CNCC_ESC <- DiffExpr_DiffBind_CNCC_ESC %>% 
  mutate(DiffExpr = ifelse(mean_logFC >2 & comb_FDR < 0.05, "CNCC",
                    ifelse(mean_logFC < -2 & comb_FDR < 0.05 ,"ESC","NOTDE")), 
         DiffBind = ifelse(min_logFC>1 & minFDR<=0.05,"CNCCbg4",
                       ifelse(min_logFC< -1 & minFDR<=0.05,"ESCbg4","NOTDBbg4")))

DiffExpr_DiffBind_NSC_vs_CNCC <- DiffExpr_DiffBind_NSC_vs_CNCC %>% 
  mutate(DiffExpr = ifelse(mean_logFC >2 & comb_FDR < 0.05, "NSC",
                    ifelse(mean_logFC < -2 & comb_FDR < 0.05 ,"CNCC","NOTDE")), 
         DiffBind = ifelse(min_logFC>1 & minFDR<=0.05,"NSCbg4",
                       ifelse(min_logFC< -1 & minFDR<=0.05,"CNCCbg4","NOTDBbg4")))


pdf('Comparison_NSC_ESC_differentialExpression_differentialBinding_peaks_expressionlogFC_boxplot.pdf')
DiffExpr_DiffBind_NSC_ESC %>%dplyr::filter(DiffExpr != "NOTDE" & !is.na(DiffExpr) & !is.na(DiffBind)  & DiffBind!="NOTDBbg4") %>% ggplot(aes(x=DiffExpr,y=mean_logFC, fill=DiffBind )) + geom_boxplot()
dev.off()

pdf('Comparison_CNCC_ESC_differentialExpression_differentialBinding_peaks_expressionlogFC_boxplot.pdf')
DiffExpr_DiffBind_CNCC_ESC %>%dplyr::filter(DiffExpr != "NOTDE" & !is.na(DiffExpr) & !is.na(DiffBind)  & DiffBind!="NOTDBbg4") %>% ggplot(aes(x=DiffExpr,y=mean_logFC, fill=DiffBind )) + geom_boxplot()
dev.off()

pdf('Comparison_NSC_vs_CNCC_differentialExpression_differentialBinding_peaks_expressionlogFC_boxplot.pdf')
DiffExpr_DiffBind_NSC_vs_CNCC %>%dplyr::filter(DiffExpr != "NOTDE" & !is.na(DiffExpr) & !is.na(DiffBind)  & DiffBind!="NOTDBbg4") %>% ggplot(aes(x=DiffExpr,y=mean_logFC, fill=DiffBind )) + geom_boxplot()
dev.off()



# do the same including NOTDE <<<<<<<<<<<<<<<<<<<<<===================================

DiffExpr_DiffBind_NSC_ESC %>%dplyr::filter(DiffExpr != "NOTDE" & !is.na(DiffExpr) & !is.na(DiffBind)  & DiffBind!="NOTDBbg4") %>% ggplot(aes(x=DiffExpr,y=mean_logFC, fill=DiffBind )) + geom_boxplot()

DiffExpr_DiffBind_NSC_ESC$all_labels <- rep('allGenes',dim(DiffExpr_DiffBind_NSC_vs_CNCC)[1])
DiffExpr_DiffBind_CNCC_ESC$all_labels <- rep('allGenes',dim(DiffExpr_DiffBind_CNCC_ESC)[1])
DiffExpr_DiffBind_NSC_vs_CNCC$all_labels <- rep('allGenes',dim(DiffExpr_DiffBind_NSC_vs_CNCC)[1])

pdf('Comparison_NSC_ESC_differentialExpression_differentialBinding_peaks_expressionlogFC_boxplot_allgenes.pdf')
DiffExpr_DiffBind_NSC_ESC %>%dplyr::filter(!is.na(DiffExpr) & !is.na(DiffBind)  & DiffBind!="NOTDBbg4") %>% ggplot(aes(x=all_labels,y=mean_logFC, fill=DiffBind )) + geom_boxplot()+ggtitle("NSC_ESC")
dev.off()

pdf('Comparison_CNCC_ESC_differentialExpression_differentialBinding_peaks_expressionlogFC_boxplot_allgenes.pdf')
DiffExpr_DiffBind_CNCC_ESC %>%dplyr::filter(!is.na(DiffExpr) & !is.na(DiffBind)  & DiffBind!="NOTDBbg4") %>% ggplot(aes(x=all_labels,y=mean_logFC, fill=DiffBind )) + geom_boxplot() +ggtitle("CNCC_ESC")
dev.off()

pdf('Comparison__NSC_vs_CNCC_differentialExpression_differentialBinding_peaks_expressionlogFC_boxplot_allgenes.pdf')
DiffExpr_DiffBind_NSC_vs_CNCC %>%dplyr::filter(!is.na(DiffExpr) & !is.na(DiffBind)  & DiffBind!="NOTDBbg4") %>% ggplot(aes(x=all_labels,y=mean_logFC, fill=DiffBind )) + geom_boxplot() +ggtitle("NSC_vs_CNCC")
dev.off()

write.table(DiffExpr_DiffBind_NSC_ESC, file='DiffExpr_DiffBind_NSC_ESC.csv',quote = F, sep = ",", col.names = NA)
write.table(DiffExpr_DiffBind_CNCC_ESC, file='DiffExpr_DiffBind_CNCC_ESC.csv',quote = F, sep = ",", col.names = NA)
write.table(DiffExpr_DiffBind_NSC_vs_CNCC, file='DiffExpr_DiffBind_NSC_vs_CNCC.csv',quote = F, sep = ",", col.names = NA)
```

### comparison BG4 peaks and ATAC peaks

``` bash

# consensus BG4 regions
/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/*bed

# consensus ATAC regions
path_atac=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output
mkdir $path_atac/intervene_atac_only
cd $path_atac
cmd_intervene="intervene venn -i *.multi2.bed -o $path_atac/intervene_atac_only"
sbatch --mem 2G --wrap "$cmd_intervene"

mkdir $path_atac/intervene_atac_only_pariwise
cmd_intervene2="intervene pairwise -i *.multi2.bed --compute frac --htype pie -o $path_atac/intervene_atac_only_pariwise"
sbatch --mem 2G --wrap "$cmd_intervene2"

mkdir $path_atac/intervene_pariwise_atac_bg4
cmd_intervene2="intervene pairwise -i *.multi2.bed ./BG4_peaks/*bed --compute frac --htype pie -o $path_atac/intervene_atac_only_pariwise_atac_bg4"
sbatch --mem 2G --wrap "$cmd_intervene2"

path_atac=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output
cd $path_atac
mkdir comparison_atac_bg4
cells=(ESC NSC CNCC)
for cell in ${cells[@]}
do
  atac_cell=`ls ${cell}_ATAC.multi2.sorted.bed`
  bg4_cell=`ls BG4_peaks/${cell}.bg4.bio2out3.sorted.bed`
  
  # for each cell compute intersection and differeence set
  intersectBed -a $atac_cell -b $bg4_cell -wa -v | sortBed -i - > ./comparison_atac_bg4/$cell.atac.specific.bed
  intersectBed -b $atac_cell -a $bg4_cell -wa -v | sortBed -i -> ./comparison_atac_bg4/$cell.bg4.specific.bed
  intersectBed -b $atac_cell -a $bg4_cell -wa | sortBed -i -> ./comparison_atac_bg4/temp1
  intersectBed -a $atac_cell -b $bg4_cell -wa | sortBed -i -> ./comparison_atac_bg4/temp2
  cat ./comparison_atac_bg4/temp1 ./comparison_atac_bg4/temp2 | mergeBed -i - > ./comparison_atac_bg4/$cell.atac.bg4.common.bed
  rm ./comparison_atac_bg4/temp1 ./comparison_atac_bg4/temp2
  
  # save to bed the coordinates of the intersections
  
done
cells=(ESC NSC CNCC)
for cell in ${cells[@]}
do
  common_intervals=`ls comparison_atac_bg4/${cell}.atac.bg4.common.bed`
  intersectBed -a ${cell}_ATAC.multi2.sorted.bed -b BG4_peaks/${cell}.bg4.bio2out3.sorted.bed -wa -wb > comparison_atac_bg4/${cell}.atac.bg4.common.peaks_intervals.bed
done


# compute coverages across common intervals in both bg4 and atac
atac_bam=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352 #merged.sorted.bam
bg4_bam=/scratcha/sblab/simeon01/Data/20190917_all_katie_bams_stem_cells
out_dir=/scratcha/sblab/simeon01/Data/20190917_all_katie_bams_stem_cells/coverages_atac_bg4
mkdir $out_dir

#cells=(ESC NSC CNCC)
cells=(CNCC)
path_atac=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output
cd $path_atac
for cell in ${cells[@]}
do
  common_regions=comparison_atac_bg4/$cell.atac.bg4.common.bed
  for bam in `ls $atac_bam/${cell}*.merged.sorted.bam`
  do
    bam_file=`basename $bam`
    sbatch --mem 8G --wrap "bedtools coverage -a $path_atac/$common_regions -b $bam  -counts > $out_dir/${bam_file%%.bam}.atac.coverages.txt"
    echo " "
  done
  echo "++++++"
  
  for bam in `ls $bg4_bam/${cell}*.markduplicates.bam`
  do
    bam_file=`basename $bam`
    #sbatch --mem 8G --wrap "bedtools coverage -a $common_regions -b $bam  -counts > $out_dir/${bam_file%%.bam}.bg4.coverages.txt"
    echo " "
  done
  echo "======"
done

# for atac_seq and BG4-chipseq find regions that are common across all cells

multiIntersectBed -i *_ATAC.multi2.sorted.bed | awk '{ if ($4>=3) print $0}'| sortBed -i - | mergeBed -i - > ATAC_ESC_CNCC_NSC.common.bed
path_atac=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output
cd $path_atac
atac_bam=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352 #merged.sorted.bam
bg4_bam=/scratcha/sblab/simeon01/Data/20190917_all_katie_bams_stem_cells
out_dir=/scratcha/sblab/simeon01/Data/20190917_all_katie_bams_stem_cells/coverages_atac_bg4_conservedATAC
mkdir $out_dir
cells=(ESC NSC CNCC)
for cell in ${cells[@]}
do
  common_regions=comparison_atac_bg4/$cell.atac.bg4.common.bed
  for bam in `ls $atac_bam/${cell}*.merged.sorted.bam`
  do
    bam_file=`basename $bam`
    sbatch --mem 8G --wrap "bedtools coverage -a ATAC_ESC_CNCC_NSC.common.bed -b $bam  -counts > $out_dir/${bam_file%%.bam}.atac.coverages.txt"
    echo " "
  done
  echo "++++++"
  
  for bam in `ls $bg4_bam/${cell}*.markduplicates.bam`
  do
    bam_file=`basename $bam`
    #sbatch --mem 8G --wrap "bedtools coverage -a ATAC_ESC_CNCC_NSC.common.bed -b $bam  -counts > $out_dir/${bam_file%%.bam}.bg4.coverages.txt"
    echo " "
  done
  echo "======"
done

# do the same with only open chromatine sites that are mantained across all and 
path_atac=/scratcha/sblab/simeon01/Data/20200128_Katie_Atac_seq/merged_SLX-18379_SLX-18352/macs2_output
cd $path_atac
intersectBed -a ATAC_ESC_CNCC_NSC.common.bed -b /scratcha/sblab/simeon01/reference_genomes/hg38/gencode.v28.TSS.bed -wa > ATAC_ESC_CNCC_NSC.common.gencode.v28.TSS.bed
out_dir=/scratcha/sblab/simeon01/Data/20190917_all_katie_bams_stem_cells/coverages_atac_bg4_conservedATAC_TSSpositive
mkdir $out_dir
cells=(ESC NSC CNCC)
for cell in ${cells[@]}
do
  common_regions=comparison_atac_bg4/$cell.atac.bg4.common.bed
  for bam in `ls $atac_bam/${cell}*.merged.sorted.bam`
  do
    bam_file=`basename $bam`
    sbatch --mem 8G --wrap "bedtools coverage -a ATAC_ESC_CNCC_NSC.common.bed -b $bam  -counts > $out_dir/${bam_file%%.bam}.atac.coverages.txt"
    echo " "
  done
  echo "++++++"
  
  for bam in `ls $bg4_bam/${cell}*.markduplicates.bam`
  do
    bam_file=`basename $bam`
    sbatch --mem 8G --wrap "bedtools coverage -a ATAC_ESC_CNCC_NSC.common.bed -b $bam  -counts > $out_dir/${bam_file%%.bam}.bg4.coverages.txt"
    echo " "
  done
  echo "======"
done
```

local folder /Users/simeon01/Documents/Katie/Atac\_seq/analysis\_after\_second\_sequencing
==========================================================================================

``` r
setwd("/Users/simeon01/Documents/Katie/Atac_seq/analysis_after_second_sequencing")

ESC_case <- read.table('ESC.atac.bg4.common.peaks_intervals.bed',stringsAsFactors = F)
NSC_case <- read.table('NSC.atac.bg4.common.peaks_intervals.bed',stringsAsFactors = F)
CNCC_case <- read.table('CNCC.atac.bg4.common.peaks_intervals.bed',stringsAsFactors = F)

ESC_case$atac_size <- ESC_case$V3-ESC_case$V2
ESC_case$bg4_size <- ESC_case$V6-ESC_case$V5
plot(ESC_case$atac_size,ESC_case$bg4_size)

NSC_case$atac_size <- NSC_case$V3-NSC_case$V2
NSC_case$bg4_size <- NSC_case$V6-NSC_case$V5

CNCC_case$atac_size <- CNCC_case$V3-CNCC_case$V2
CNCC_case$bg4_size <- CNCC_case$V6-CNCC_case$V5

#setwd("/scratcha/sblab/simeon01/Data/20190917_all_katie_bams_stem_cells/coverages_atac_bg4")
setwd("/Users/simeon01/Documents/Katie/Atac_seq/analysis_after_second_sequencing/coverages_atac_bg4_conservedATAC")
#coverages
ESC_atac_list <- list.files(path='.',pattern='ESC.*ATAC.*txt')

ESC_bg4_list <- list.files(path='.',pattern='ESC.*BG4.*txt')
ESC_bg4_list <- ESC_bg4_list[-grep('input',ESC_bg4_list)]

NSC_atac_list <- list.files(path='.',pattern='NSC.*ATAC.*txt')

NSC_bg4_list <- list.files(path='.',pattern='NSC.*BG4.*txt')
NSC_bg4_list <- NSC_bg4_list[-grep('input',NSC_bg4_list)]

CNCC_atac_list <- list.files(path='.',pattern='CNCC.*ATAC.*txt')

CNCC_bg4_list <- list.files(path='.',pattern='CNCC.*BG4.*txt')
CNCC_bg4_list <- CNCC_bg4_list[-grep('input',CNCC_bg4_list)]


stat5_list <- list.files(path='/Users/simeon01/Documents/Katie/Atac_seq/analysis_after_second_sequencing/coverages_atac_bg4',pattern='stat5')

import_coverages_from_files <- function(file_list){
  temp=read.table(file_list[1],stringsAsFactors = F)[,4]
  M_initial <- matrix(ncol=length(file_list),nrow=length(temp))
  for (i in (1:length(file_list))){
    temp1 <- read.table(file_list[i],stringsAsFactors = F)
    M_initial[,i] <- temp1$V4
  
  }
  colnames(M_initial) <- file_list
  return(M_initial)
  
}
norm_sign_to_tot_map_reads <- function(matrix,list_stat5_mapped_reads){
  Library_sizes <- c()
  matrix_norm <- matrix
  for (i in 1:dim(matrix)[2]){
    print(i)
    tmp_cond <- colnames(matrix)[i]
    tmp_cond <- gsub('.atac.coverages.txt|.markduplicates.bg4.coverages.txt','',tmp_cond)
    print(tmp_cond)
    tmp_stat <- grep(substr(tmp_cond,1,20),list_stat5_mapped_reads,value = T)
    tot_reads <- as.numeric(system(paste0('cat ',tmp_stat),intern = TRUE))
    Library_sizes[i] <- tot_reads
    matrix_norm[,i] <- matrix[,i]/tot_reads*1e6
  }
  return(matrix_norm)
}

setwd('/Users/simeon01/Documents/Katie/Atac_seq/analysis_after_second_sequencing/coverages_atac_bg4_conservedATAC/')

## ESC ##

ESC_atac_M <- import_coverages_from_files(ESC_atac_list)

ESC_atac_M_norm <- norm_sign_to_tot_map_reads(ESC_atac_M,stat5_list)

ESC_bg4_M <- import_coverages_from_files(ESC_bg4_list)

ESC_bg4_M_norm <- norm_sign_to_tot_map_reads(ESC_bg4_M,stat5_list)

ESC_atac_M_norm_avg <- rowMeans(ESC_atac_M_norm)
ESC_bg4_M_norm_avg <- rowMeans(ESC_bg4_M_norm)
 
## NSC ##

NSC_atac_M <- import_coverages_from_files(NSC_atac_list)

NSC_atac_M_norm <- norm_sign_to_tot_map_reads(NSC_atac_M,stat5_list)

NSC_bg4_M <- import_coverages_from_files(NSC_bg4_list)

NSC_bg4_M_norm <- norm_sign_to_tot_map_reads(NSC_bg4_M,stat5_list)

NSC_atac_M_norm_avg <- rowMeans(NSC_atac_M_norm)
NSC_bg4_M_norm_avg <- rowMeans(NSC_bg4_M_norm)

## CNCC ##

CNCC_atac_M <- import_coverages_from_files(CNCC_atac_list)

CNCC_atac_M_norm <- norm_sign_to_tot_map_reads(CNCC_atac_M,stat5_list)

CNCC_bg4_M <- import_coverages_from_files(CNCC_bg4_list)

CNCC_bg4_M_norm <- norm_sign_to_tot_map_reads(CNCC_bg4_M,stat5_list)

CNCC_atac_M_norm_avg <- rowMeans(CNCC_atac_M_norm)
CNCC_bg4_M_norm_avg <- rowMeans(CNCC_bg4_M_norm)


boxplot(ESC_case$atac_size,ESC_case$bg4_size,NSC_case$atac_size,NSC_case$bg4_size,CNCC_case$atac_size,CNCC_case$bg4_size, outline = F)

plot(log2(ESC_atac_M_norm_avg),log2(ESC_bg4_M_norm_avg))
plot(log2(NSC_atac_M_norm_avg),log2(NSC_bg4_M_norm_avg))
plot(log2(CNCC_atac_M_norm_avg),log2(CNCC_bg4_M_norm_avg))

M <- cbind(ESC_atac_M_norm_avg,ESC_bg4_M_norm_avg,NSC_atac_M_norm_avg,NSC_bg4_M_norm_avg,CNCC_atac_M_norm_avg,CNCC_bg4_M_norm_avg)
M_atac <- cbind(ESC_atac_M_norm_avg,CNCC_atac_M_norm_avg,NSC_atac_M_norm_avg)
M_bg4 <- cbind(ESC_bg4_M_norm_avg,CNCC_bg4_M_norm_avg,NSC_bg4_M_norm_avg)
library(ComplexHeatmap)
Heatmap(cor(M),cluster_rows = F, cluster_columns = F)
```

Fold enrichments of G4 peaks at different promoters as annotated in All\_anno.RData
-----------------------------------------------------------------------------------

``` bash
library(tidyverse)
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/20200421_annotated_promoters_using_histones_G4')
load("All_anno.RData")

All_anno_prom_sel <- All_anno %>% select(ens_id,ESC_status,CNCC_status,NSC_status)


promoters_annotations <- read.delim('/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS.slop1000.bed',sep ="\t",stringsAsFactors = F, header = F)

colnames(promoters_annotations) <- c("chr","start","end","ens_id","feat","strand")

All_anno_prom_sel_prom_coord <- dplyr::left_join(All_anno_prom_sel,promoters_annotations,by="ens_id")




## steps: select only promoters with 1 in each case (ESC, NSC and CNCC)

ESC_anno_promoters <- All_anno_prom_sel_prom_coord %>% select(chr,start,end,ens_id,feat,strand,ESC_status)

CNCC_anno_promoters <- All_anno_prom_sel_prom_coord %>% select(chr,start,end,ens_id,feat,strand,CNCC_status)

NSC_anno_promoters <- All_anno_prom_sel_prom_coord %>% select(chr,start,end,ens_id,feat,strand,NSC_status)

write.table(ESC_anno_promoters,file ='ESC_anno_promoters_for_foldEnrichments.bed',sep = "\t",quote = F,col.names = F,row.names = F )

write.table(NSC_anno_promoters,file ='NSC_anno_promoters_for_foldEnrichments.bed',sep = "\t",quote = F,col.names = F,row.names = F )

write.table(CNCC_anno_promoters,file ='CNCC_anno_promoters_for_foldEnrichments.bed',sep = "\t",quote = F,col.names = F,row.names = F )
```

Go to the folder with those annotations and sort the bed files

``` bash
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/20200421_annotated_promoters_using_histones_G4
for file in *for_foldEnrichments.bed
do
  bedtools sort -i $file | cut -f 1,2,3,7 > ${file%%.bed}.sorted.bed
done

scp *_foldEnrichments.sorted.bed simeon01@10.20.236.34:/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells
```

Move annotation files to cluster and run enrichment analysis

``` bash



# perform fold enrichmeent analysis on any bed file in each subfolder
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019
path_anno='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells'
scenario=(scenario1_multi2 scenario2_5out9 scenario3_multi2_extended_1kb)
#genomic_annotation=annotations_genomic_regions_gencode_v28.anno.bed
workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed

ESC_prom=ESC_anno_promoters_for_foldEnrichments.sorted.bed
CNCC_prom=CNCC_anno_promoters_for_foldEnrichments.sorted.bed
NSC_prom=NSC_anno_promoters_for_foldEnrichments.sorted.bed
annotation_list1=($ESC_prom $CNCC_prom $NSC_prom)
scenario=(scenario1_multi2)

for case in ${scenario[@]}
do
  cd ./$case
  # run gat analysis on full bed files 
  for file in `ls *bed | grep -v "shuffle"`
  do
    echo $file
    for anno in ${annotation_list1[@]}
    do
      echo $anno
      pwd
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$path_anno/${anno} --segments=${file} -S ${file%%.bed}_${anno%%.bed}.dat -t 4"
      echo ${file%%.bed}_${anno%%.bed}.dat
      echo $cmd_gat
      sbatch --mem 6G --wrap "$cmd_gat"
      echo " "
    done
  done

  #run the enrichment for subsets (common and specific ones)
  cd ./bed_intersections_common
  pwd
  for file in *bed
  do
    echo $file
    for anno in ${annotation_list1[@]}
    do
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$path_anno/${anno} --segments=${file} -S ${file%%.bed}_${anno%%.bed}.dat -t 4"
      echo ${file%%.bed}.${anno%%.bed}.dat
      pwd
      echo $cmd_gat
      sbatch --mem 6G --wrap "$cmd_gat"
      echo " "
    done
  done
  echo " ** == ** =="
  cd ../..
done



for case in ${scenario[@]}
do
  cd ./$case
  # run gat analysis on full bed files 
  for file in `ls CNCC_specific.20200421.bed`
  do
    echo $file
    for anno in ${annotation_list1[@]}
    do
      echo $anno
      pwd
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$path_anno/${anno} --segments=${file} -S ${file%%.bed}_${anno%%.bed}.dat -t 4"
      echo ${file%%.bed}_${anno%%.bed}.dat
      echo $cmd_gat
      sbatch --mem 6G --wrap "$cmd_gat"
      echo " "
    done
  done
done
```

annotate BG4 regions based on OQs info
--------------------------------------

``` bash

# for all G4 regions (BG4 peaks by cell type) extract bed files with center and imputed strand (based on OQs)
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2
path_oqs=/scratcha/sblab/simeon01/reference_genomes/OQs
for file in **bio2out3.bed  
do
  sbatch --mem 2G --wrap "annotateBed -i $file -files $path_oqs/OQ_plus_hits.lifted_hg19_to_hg38.Minus_Strand_bed6.bed $path_oqs/OQ_minus_hits.lifted_hg19_to_hg38.Plus_Strand_bed6.bed -names Minus_Strand_bed6 Plus_Strand_bed6  -counts > ${file%%.bed}.G4_peaks_anotated_by_OQsStrand.bed"
done


for file in *.G4_peaks_anotated_by_OQsStrand.bed
do
  awk '{if ($4>=0 && $5==0) {print $1"\t"$2+int(($3-$2)/2+0.49)-1 "\t" $2+int(($3-$2)/2+0.49)"\t.\t.\t-"} else if ($4==0 && $5>=0) {print $1"\t"$2+int(($3-$2)/2+0.49)"\t"$2 +int(($3-$2)/2+0.49)+1 "\t.\t.\t+"}}' $file > ${file%%.G4_peaks_anotated_by_OQsStrand.bed}.center.OQsStrand.bed
done


# for regions that are common and specific, assign a sign (based on OQs) and use the coordinates of the center
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/bed_intersections_common
path_oqs=/scratcha/sblab/simeon01/reference_genomes/OQs
for file in *C.bed
do
  sbatch --mem 2G --wrap "annotateBed -i $file -files $path_oqs/OQ_plus_hits.lifted_hg19_to_hg38.Minus_Strand_bed6.bed $path_oqs/OQ_minus_hits.lifted_hg19_to_hg38.Plus_Strand_bed6.bed -names Minus_Strand_bed6 Plus_Strand_bed6  -counts > ${file%%.bed}.G4_peaks_anotated_by_OQsStrand.bed"
done


for file in *.G4_peaks_anotated_by_OQsStrand.bed
do
  awk '{if ($4>=0 && $5==0) {print $1"\t"$2+int(($3-$2)/2+0.49)-1 "\t" $2+int(($3-$2)/2+0.49)"\t.\t.\t-"} else if ($4==0 && $5>=0) {print $1"\t"$2+int(($3-$2)/2+0.49)"\t"$2 +int(($3-$2)/2+0.49)+1 "\t.\t.\t+"}}' $file > ${file%%.G4_peaks_anotated_by_OQsStrand.bed}.center.OQsStrand.bed
done

#
```

intervene pairwise -i ESC\_specific.20200421.bed CNCC\_specific.20200421.bed NSC\_specific.20200421.bed ESC\_CNCC\_NSC.multi3\_common.bed /scratcha/sblab/simeon01/Data/20200128\_Katie\_Atac\_seq/merged\_SLX-18379\_SLX-18352/macs2\_output/\*multi2.bed --type genomic --compute frac --htype pie --filenames --output ./intervene\_bg4\_atac

density plots of BG4 chip-seq
-----------------------------

``` bash

cells=(CNCC NSC)
folder_out_result=/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_second_BG4_seq_at_G4_ppromoters
mkdir $folder_out_result

for cell in ${cells[@]}
do
cd /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW
  for bwfile in ${cell}*bw
  do
    # DAUTHER CELL 1
    sbatch --mem 6G --wrap="computeMatrix reference-point -S ${bwfile} -R ${cell}_bed_lists/*bed  -a 1000 -b 1000 --outFileName $folder_out_result/${bwfile}_${cell}.matrix.gz --missingDataAsZero && 
    plotProfile --matrixFile $folder_out_result/${bwfile}_${cell}.matrix.gz -o $folder_out_result/${bwfile}_${cell}.matrix.pdf --averageType median --plotWidth 18 --plotHeight 10 &&
    plotHeatmap --matrixFile $folder_out_result/${bwfile}_${cell}.matrix.gz -o $folder_out_result/${bwfile}_${cell}.matrix_heatmap.pdf"
    echo "========================"
  done
  for bwfile in ESC*bw
  do
    # ESC
    sbatch --mem 6G --wrap="computeMatrix reference-point -S ${bwfile} -R ${cell}_bed_lists/*bed  -a 1000 -b 1000 --outFileName $folder_out_result/${bwfile}_${cell}.matrix.gz --missingDataAsZero && 
    plotProfile --matrixFile $folder_out_result/${bwfile}_${cell}.matrix.gz -o $folder_out_result/${bwfile}_${cell}.matrix.pdf --averageType median --plotWidth 18 --plotHeight 10 &&
    plotHeatmap --matrixFile $folder_out_result/${bwfile}_${cell}.matrix.gz -o $folder_out_result/${bwfile}_${cell}.matrix_heatmap.pdf"
    echo "========================"
  done
  echo "**********"
done


#local machine
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2

scp simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_second_BG4_seq_at_G4_ppromoters 
```

``` r
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_plots_second_BG4_seq_at_G4_ppromoters')
#NSC_ATAC*matrix*gz
#ESC_ATAC*matrix*gz
#CNCC_ATAC*matrix*gz

ending_label_CNCC="_CNCC.matrix.gz"
ending_label_NSC="_NSC.matrix.gz"
ending_label_ESC="_ESC.matrix.gz"

case_label_ESC='ESC_'

case_label_CNCC='CNCC_'

case_label_NSC='NSC_'


average_densities <- function(case_label,ending_label){
  matrix_files <- list.files(".",paste0(case_label,".*",ending_label,"$"))
  matrix_files
  # extract params for initialize matrix
  header_temp <- readLines(gzfile(matrix_files[1]),n=1)
  
  A <- unlist(strsplit(header_temp, split=',\"'))
  A <- gsub('\"','',A)
  A <- gsub('[','',A, fixed = T)
  A <- gsub(']','',A, fixed = T)
  N_rows <- as.numeric(tail(unlist(strsplit(grep('group_boundaries',A, value =T),',|:')),1))
  groups_ind <- as.numeric(tail(unlist(strsplit(grep('group_boundaries',A, value =T),',|:')),-2)) # to exclude first 2
  N_col <- as.numeric(tail(unlist(strsplit(grep('sample_boundaries',A, value =T),',|:')),1))

  temp1_sum <- matrix(0,ncol=N_col,nrow=N_rows)
  for (i in 1:length(matrix_files)){
    header <- readLines(gzfile(matrix_files[i]),n=1)
    temp1 <- read.table(gzfile(matrix_files[i]),skip = 1, stringsAsFactors = F) 
    temp1_sum <- temp1_sum +  temp1[,7:dim(temp1)[2]]
    print(i)
  }
  temp1_avg <- temp1_sum/length(matrix_files)
  temp1_avg_complete <- cbind(temp1[,1:6],temp1_avg)
  file_name <- paste0(case_label,ending_label,'_average.matrix')
  fileConn<-file(file_name)
  writeLines(header,fileConn)
  write.table(temp1_avg_complete, file = file_name, append = T,sep = "\t",quote = F, col.names = F, row.names = F)
  flag='done'
  close(fileConn)
  list_sizes <- list("N_rows"=N_rows,"groups_ind"=groups_ind,"N_col"=N_col)
  return(list_sizes)
}  


#average_densities(case_label_ESC,ending_label_ESC)
average_densities(case_label_ESC,ending_label_CNCC)
average_densities(case_label_ESC,ending_label_NSC)


average_densities(case_label_CNCC,ending_label_CNCC)
average_densities(case_label_NSC,ending_label_NSC) 

# cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_plots_second_atac_seq_at_bivalent_promoters

scp *average.matrix.gz simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_second_BG4_seq_at_G4_ppromoters
```

generate Heatmap and plot with deeptools

``` bash
for file in *average*.matrix.gz
do
  sbatch --mem 2G --wrap="plotHeatmap --matrixFile $file  --averageTypeSummaryPlot median -o ${file%%.matrix.gz}.matrix_heatmap.pdf --yMax 1 --zMax 1"
  echo "========================"
done


# locally 
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_plots_second_BG4_seq_at_G4_ppromoters
scp simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_second_BG4_seq_at_G4_ppromoters/*average.matrix_heatmap.pdf . 
```

Compare Stem cells peaks to other cancers bed files (BG4-ChIP-seq)
------------------------------------------------------------------

===================================================================== For this analysis I consider the peaks for - K562 (from Karen, from Jochen) - HeLa cells (YuQi) - HaCat cells (Robert) - U2OS

-   IMR90
-   Breast cancers models
-   "Normal cancer model"
-   HEK cells (YuQi) =====================================================================

after Katie's email 2 July Excel sheet containing number and % overlap of stem cell G4 map with cancer cell lines. We discussed it might be possible to report K562, HaCAT, (two cancer) NHEK, IMR90 (two normal) in the paper as three of these have been already published (and maybe U20S)-I will speak to David?

After this email i have to include - U2OS - NHEK

``` bash
 cd /scratcha/sblab/simeon01/Data/mlmod_d/nature_gen
intersectBed -a GSE76688_HEKnp_Lonza_1472015_BG4.1e4_peaks.narrowPeak -b GSE76688_HEKnp_Lonza_1572015_BG4.1e4_peaks.narrowPeak | sort | uniq | sortBed -i -  | cut -f 1,2,3 > GSE76688_HEKnp_Lonza.BG4.second.multi2.bed
file=GSE76688_HEKnp_Lonza.BG4.second.multi2.bed
chain_hg19_to_hg38=/scratcha/sblab/simeon01/reference_genomes/liftover_chain/hg19ToHg38.over.chain.gz
/Users/simeon01/applications/liftOver $file $chain_hg19_to_hg38 ${file%%.bed}.hg38.bed  ${file%%.bed}.unmapped
```

``` bash
 
 K562_karen=
K562_jochen=


cd /scratcha/sblab/simeon01/Data/20190523_Angela_metrics/nature_prot/peaks/
for file in *common_peaks.bed
do
chain_hg19_to_hg38=/scratcha/sblab/simeon01/reference_genomes/liftover_chain/hg19ToHg38.over.chain.gz
/Users/simeon01/applications/liftOver $file $chain_hg19_to_hg38 ${file%%.bed}.hg38.bed  ${file%%.bed}.unmapped
done

#lifover consensus
chain_hg19_to_hg38=/scratcha/sblab/simeon01/reference_genomes/liftover_chain/hg19ToHg38.over.chain.gz
for file in *consensus.multi2.bed
do
  echo $file
  echo ${file%%.bed}.hg38.bed
  /Users/simeon01/applications/liftOver $file $chain_hg19_to_hg38 ${file%%.bed}.hg38.bed  ${file%%.bed}unmapped
  #/Users/simeon01/applications/liftOver $file $chain_hg19_to_hg38 ${file%%.bed}.hg38.bed
done
U2OS_rob_natureProt=/scratcha/sblab/simeon01/Data/20190523_Angela_metrics/nature_prot/peaks/U2OS_common_peaks.hg38.bed
Nhek_rob_natureGen=/scratcha/sblab/simeon01/Data/mlmod_d/nature_gen/GSE76688_HEKnp_Lonza.BG4.second.multi2.hg38.bed
K562_rob=/scratcha/sblab/simeon01/Data/20190523_Angela_metrics/nature_prot/peaks/K562_common_peaks.hg38.bed
IMR90_rob=/scratcha/sblab/simeon01/Data/20190523_Angela_metrics/nature_prot/peaks/IMR90_common_peaks.hg38.bed
HaCaT_rob=/scratcha/sblab/simeon01/Data/20190523_Angela_metrics/nature_prot/peaks/HaCaT_common_peaks.hg38.bed

HeLa_yuqi=/scratcha/sblab/simeon01/Data/20200205_yuqi_hl/SLX-18996/trimmed/aligned/macs2_individual_rep/multi2_bed/HeLa.20200421.bio2.sorted.bed
HEK_yuqi=/scratcha/sblab/simeon01/Data/20200207_yuqi_HEK293/SLX-19000/trimmed/aligned/macs2_individual_rep/multi2_bed/Hek.bio2.sorted.bed
HEK_dhaval=/scratcha/sblab/simeon01/Data/20200701_dhaval_old_HEK293/merged_runs/macs2_individual_rep/consensus_narrow_peaks/HEK293_NTC_no_ChIP_.multi2.bed
MCF7=/scratcha/sblab/simeon01/Data/mlmod_d/mcf7/MCF7.default_peaks.hg19.bio2out3.hg38.sorted.bed
MCF10=/scratcha/sblab/simeon01/Data/mlmod_d/mcf10/MCF10A.default_peaks.hg19.bio2out3.hg38.sorted.bed
MDAMB231=/scratcha/sblab/simeon01/Data/mlmod_d/MDAMB231/MDAMB231.default_peaks.hg19.bio2out3.hg38.sorted.bed

K562_karen=/scratcha/sblab/simeon01/Data/karen_K562_noDRB_peaks/chipseq_noDRB.multi5.sorted.bed
K562_jochen=/scratcha/sblab/simeon01/Data/Jochen_K562_peaks/20180108_K562_async_rep1-3.mult.5of8.hg19_lifted_hg38.sorted.bed
oqs=/scratcha/sblab/simeon01/reference_genomes/OQs/OQ_hits.lifted_hg19_to_hg38_no_.bed

mkdir /scratcha/sblab/simeon01/Data/20190523_Angela_metrics/cells_references_hg38
cd /scratcha/sblab/simeon01/Data/20190523_Angela_metrics/cells_references_hg38
# list of cells: all_cells_genomic enrichments : U2OS_rob_natureProt Nhek_rob_natureGen K562_robnatureProt IMR90_robnatureProt HaCaT_robnatureProt HeLa_yuqi HEK_yuqi HEK_dhaval MCF7 MCF10 MDAMB231 K562_karen K562_jochen oqs
# U2OS_rob_natureProt
# Nhek_rob_natureGen
# K562_robnatureProt
# IMR90_robnatureProt
# HaCaT_robnatureProt
# HeLa_yuqi
# HEK_yuqi
# HEK_dhaval
# MCF7
# MCF10
# MDAMB231
# K562_karen
# K562_jochen
all_cells_files=($U2OS_rob_natureProt $Nhek_rob_natureGen $K562_rob $IMR90_rob $HaCaT_rob $HeLa_yuqi $HEK_yuqi $HEK_dhaval $MCF7 $MCF10 $MDAMB231 $K562_karen $K562_jochen $oqs)

for cell in ${all_cells_files[@]}
do
cp $cell /scratcha/sblab/simeon01/Data/20190523_Angela_metrics/cells_references_hg38
done

cd /scratcha/sblab/simeon01/Data/20190523_Angela_metrics/cells_references_hg38

mv 20180108_K562_async_rep1-3.mult.5of8.hg19_lifted_hg38.sorted.bed K562_Jochen.cell.20180108_K562_async_rep1-3.mult.5of8.hg19_lifted_hg38.sorted.bed
mv chipseq_noDRB.multi5.sorted.bed K562_Karen.cell.chipseq_noDRB.multi5.sorted.bed
mv GSE76688_HEKnp_Lonza.BG4.second.multi2.hg38.bed NHek_Rob.cell.GSE76688_HEKnp_Lonza.BG4.second.multi2.hg38.bed
mv HaCaT_common_peaks.hg38.bed HacCaT_Rob.cell.HaCaT_common_peaks.hg38.bed
mv HEK293_NTC_no_ChIP_.multi2.bed HEK293_Dhaval.cell.HEK293_NTC_no_ChIP_.multi2.bed
mv Hek.bio2.sorted.bed HEK293_Yuqi.cell.Hek.bio2.sorted.bed
mv HeLa.20200421.bio2.sorted.bed HeLa_Yuqi.cell.HeLa.20200421.bio2.sorted.bed
mv IMR90_common_peaks.hg38.bed IMR90_Rob.cell.IMR90_common_peaks.hg38.bed
mv K562_common_peaks.hg38.bed K562_Rob.cell.K562_common_peaks.hg38.bed
mv MCF10A.default_peaks.hg19.bio2out3.hg38.sorted.bed MCF10A_Winnie.cell.MCF10A.default_peaks.hg19.bio2out3.hg38.sorted.bed
mv MCF7.default_peaks.hg19.bio2out3.hg38.sorted.bed MCF7_Winnie.cell.MCF7.default_peaks.hg19.bio2out3.hg38.sorted.bed
mv MDAMB231.default_peaks.hg19.bio2out3.hg38.sorted.bed MDAMB231_Winnie.cell.MDAMB231.default_peaks.hg19.bio2out3.hg38.sorted.bed
mv OQ_hits.lifted_hg19_to_hg38_no_.bed OQ_hits.cell.OQ_hits.lifted_hg19_to_hg38_no_.bed
mv U2OS_common_peaks.hg38.bed U2OS_Rob.cell.U2OS_common_peaks.hg38.bed


# create an annotation file for the cells
touch cells_mix_annotation.anno
for file in *bed
do
echo $file
label=${file%%.cell*}
echo $label
awk -v var=$label '{print $0"\t"var}' $file >>cells_mix_annotation.anno

done



## compute annotations overlap and compute enrichments in cells_mix_annotation.anno

## compute overlap with intervene

bed_to_use1=`ls  /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/*bed | grep -v shuffle | grep -v intermediate`
bed_to_use2=`ls  /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/bed_intersections_common/*C.bed | grep -v shuffle | grep -v intermediate`
out_dir1=/scratcha/sblab/simeon01/Data/20190523_Angela_metrics/cells_references_hg38/intervene_compare_stem_cells_to_all_other_cells
out_dir2=/scratcha/sblab/simeon01/Data/20190523_Angela_metrics/cells_references_hg38/intervene_compare_stem_cells_to_all_other_cells2
mkdir $out_dir1
mkdir $out_dir2
cd /scratcha/sblab/simeon01/Data/20190523_Angela_metrics/cells_references_hg38/stem_cells_files

sbatch --time 02:00:00 --mem 6G --wrap "intervene pairwise -i *.bed ../*cell*.bed --type genomic --compute count --htype pie --output $out_dir1"
sbatch --time 02:00:00 --mem 6G --wrap "intervene pairwise -i *.bed ../*cell*.bed --type genomic --compute count --htype pie --output $out_dir2"

## -------------------------------------------------------------------------

multiIntersectBed -i /scratcha/sblab/simeon01/Data/20190523_Angela_metrics/cells_references_hg38/stem_cells_files/ESC_CNCC_NSC.multi3_common.bed HacCaT_Rob.cell.HaCaT_common_peaks.hg38.bed K562_Rob.cell.K562_common_peaks.hg38.bed U2OS_Rob.cell.U2OS_common_peaks.hg38.bed | awk '{ if ($4>=4) print $0}' |sort | uniq | sortBed -i - >  /scratcha/sblab/simeon01/Data/20190523_Angela_metrics/cells_references_hg38/stem_cells_regions_in_common_with_rob_nat_pro/regions_multi3.all_cancer.bed


cat NHek_Rob.cell.GSE76688_HEKnp_Lonza.BG4.second.multi2.hg38.sorted.bed IMR90_Rob.cell.IMR90_common_peaks.hg38.sorted.bed MCF10A_Winnie.cell.MCF10A.default_peaks.hg19.bio2out3.hg38.sorted.bed | sortBed -i - | mergeBed -i - > normal.cell.bed

intersectBed -a /scratcha/sblab/simeon01/Data/20190523_Angela_metrics/cells_references_hg38/stem_cells_regions_in_common_with_rob_nat_pro/regions_multi3.all_cancer.bed -b /scratcha/sblab/simeon01/Data/20190523_Angela_metrics/cells_references_hg38/normal.cell.bed -wa -v > /scratcha/sblab/simeon01/Data/20190523_Angela_metrics/cells_references_hg38/stem_cells_regions_in_common_with_rob_nat_pro/regions_multi3.all_cancer.specific_to_cancer.bed


## 
#https://bioconductor.org/packages/devel/bioc/vignettes/annotatr/inst/doc/annotatr-vignette.html

all_cells=($U2OS_rob $IMR90_rob $HaCaT_rob $HeLa_yuqi $MCF7 $MCF10 $MDAMB231)
all_cells=($K562_karen $K562_jochen)


for cell in ${all_cells[@]}
do
  sortBed -i $cell > ${cell%%.bed}.sorted.bed
done


cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2


bed_to_use=`ls *bed | grep -v shuffle | grep -v intermediate`
for file in ${bed_to_use[@]}
do
  echo $file
  out_dir=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/comparison_cells_lines/comparison_to_other_cell_lines_${file%%.bed}
  mkdir $out_dir
  sbatch --mem 6G --wrap "intervene pairwise -i $file $U2OS_rob $IMR90_rob $HaCaT_rob $HeLa_yuqi $MCF7 $MCF10 $MDAMB231 $K562_karen $K562_jochen $oqs --type genomic --compute count --htype pie --output $out_dir"
  #sbatch --mem 6G --wrap "intervene upset -i $file $U2OS_rob $IMR90_rob $HaCaT_rob $HeLa_yuqi $MCF7 $MCF10 $MDAMB231 $K562_karen $K562_jochen $oqs --output $out_dir"
done


#other aspect of comparing to cancer: check what pathways are in the common sets
 intersectBed -a  /scratcha/sblab/simeon01/reference_genomes/hg38/gencode.v28.TSS.slop1000.bed -b ESC_CNCC_NSC.multi3_common.bed -wa | sort | uniq | cut -f 4 | sed 's/"//g' | sed 's/\..*$//g' > ESC_CNCC_NSC.multi3_common.TSS.list
```

Compare G4 with other embyo and germ cells.

``` bash
cd /scratcha/sblab/simeon01/Data/Human_Sperm_GSE124718/new_download
guzip *zip

chain_hg19_to_hg38=/scratcha/sblab/simeon01/reference_genomes/liftover_chain/hg19ToHg38.over.chain.gz
for file in *Peak
do
  echo $file
  echo ${file}.hg38.bed
  cut -f 1,2,3 $file  > temp.bed && /Users/simeon01/applications/liftOver temp.bed $chain_hg19_to_hg38 ${file}.hg38.bed  ${file}unmapped
done

ESC_all=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC.bio2out3.bed
ESC_spec=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC.bio2out3.specific_to_cell.bed
ESC_common=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/ESC_CNCC_NSC.multi3_common.bed
ESC_K4=/scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/final_reference_sets_stem_cells_published_data/Rada_Iglesias_lifted_hg38.H3K4me3.bed

cd /scratcha/sblab/simeon01/Data/Human_Sperm_GSE124718/New_download

mkdir intervene_develop
intervene pairwise -i *H3K4me3*hg38.bed $ESC_all $ESC_spec $ESC_common $ESC_K4 --type genomic --compute frac --htype pie --output ./intervene_develop
```

``` r
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/comparison_cells_lines')
```

create Bed reference sets for density plots
-------------------------------------------

Here i generated the bed intervals of TSSs stratyfing by promoter types where types are defined based on histone modification status

``` bash
library(tidyverse)
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/20200421_annotated_promoters_using_histones_G4')
load("All_anno.RData")

All_anno <- All_anno %>% mutate(ESC_G4_status=case_when(ESC_G4 == 0 ~ "G4-",ESC_G4> 0 ~ "G4+")) %>% mutate(NSC_G4_status=case_when(NSC_G4 == 0 ~ "G4-",NSC_G4> 0 ~ "G4+")) %>% mutate(CNCC_G4_status=case_when(CNCC_G4 == 0 ~ "G4-",CNCC_G4> 0 ~ "G4+"))


All_anno_prom_sel <- All_anno %>% select(ens_id,ESC_status,CNCC_status,NSC_status,ESC_G4_status,NSC_G4_status,CNCC_G4_status)

TSS_1base <- read.delim('/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS.minus1base.bed',stringsAsFactors = F,header = F)
colnames(TSS_1base) <- c('TSS_chr','TSS_start','TSS_end','ens_id','feature','TSS_strand')



All_anno_prom_sel_prom_coord_1base <- dplyr::left_join(All_anno_prom_sel,TSS_1base,by="ens_id")


## steps: select promoter type and generate a new field that contains G4 and promoter type by cell type then print it to bed file


ESC_anno_promoters_G4plus <- All_anno_prom_sel_prom_coord_1base %>% filter(ESC_G4_status == "G4+") %>% select(TSS_chr,TSS_start,TSS_end,ens_id,feature,TSS_strand,ESC_status,ESC_G4_status)
ESC_anno_promoters_G4minus <- All_anno_prom_sel_prom_coord_1base %>% filter(ESC_G4_status == "G4-") %>% select(TSS_chr,TSS_start,TSS_end,ens_id,feature,TSS_strand,ESC_status,ESC_G4_status)


CNCC_anno_promoters_G4plus <- All_anno_prom_sel_prom_coord_1base %>% filter(CNCC_G4_status == "G4+") %>% select(TSS_chr,TSS_start,TSS_end,ens_id,feature,TSS_strand,ESC_status,ESC_G4_status)
CNCC_anno_promoters_G4minus <- All_anno_prom_sel_prom_coord_1base %>% filter(CNCC_G4_status == "G4-") %>% select(TSS_chr,TSS_start,TSS_end,ens_id,feature,TSS_strand,ESC_status,ESC_G4_status)


NSC_anno_promoters_G4plus <- All_anno_prom_sel_prom_coord_1base %>% filter(NSC_G4_status == "G4+") %>% select(TSS_chr,TSS_start,TSS_end,ens_id,feature,TSS_strand,ESC_status,ESC_G4_status)
NSC_anno_promoters_G4minus <- All_anno_prom_sel_prom_coord_1base %>% filter(NSC_G4_status == "G4-") %>% select(TSS_chr,TSS_start,TSS_end,ens_id,feature,TSS_strand,ESC_status,ESC_G4_status)


write.table(ESC_anno_promoters_G4plus,file ='ESC_anno_promoters_G4plus.prom_coord.bed',sep = "\t",quote = F,col.names = F,row.names = F )
write.table(ESC_anno_promoters_G4minus,file ='ESC_anno_promoters_G4minus.prom_coord.bed',sep = "\t",quote = F,col.names = F,row.names = F )

write.table(NSC_anno_promoters_G4plus,file ='NSC_anno_promoters_G4plus.prom_coord.bed',sep = "\t",quote = F,col.names = F,row.names = F )
write.table(NSC_anno_promoters_G4minus,file ='NSC_anno_promoters_G4minus.prom_coord.bed',sep = "\t",quote = F,col.names = F,row.names = F )

write.table(CNCC_anno_promoters_G4plus,file ='CNCC_anno_promoters_G4plus.prom_coord.bed',sep = "\t",quote = F,col.names = F,row.names = F )
write.table(CNCC_anno_promoters_G4minus,file ='CNCC_anno_promoters_G4minus.prom_coord.bed',sep = "\t",quote = F,col.names = F,row.names = F )
```

Copy files to folder on cluster

``` bash
 
 scp *anno_promoters_G4*.prom_coord.bed simeon01@10.20.236.34:/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells
```

density plots

``` bash
## density plots of BG4 chip-seq
```

``` bash
cd /scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells
mkdir promoter_annotated_by_histones
cp **_anno_promoters_G4*.prom_coord.bed  ./promoter_annotated_by_histones
# for each annotation file generate and saved as
# split annotations in independent subfile
##
for file in *_anno_promoters_G4*.prom_coord.bed 
do
     sbatch --mem 1G --wrap "awk '{print >> \$7FILENAME; close(\$7)}' $file"
done
rm ESC*bed
rm NSC*bed
rm CNCC*bed


cells=(ESC CNCC NSC)
folder_out_result=/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_BG4_seq_at_promoter_stratified_by_histones
path_reference_bed=/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/promoter_annotated_by_histones
mkdir $folder_out_result

for cell in ${cells[@]}
do
cd /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW
  for bwfile in ${cell}*bw
  do
    sbatch --mem 6G --wrap="computeMatrix reference-point -S ${bwfile} -R $path_reference_bed/*${cell}_anno_promoters_G4*plus.prom_coord.bed   -a 1000 -b 1000 --outFileName $folder_out_result/${bwfile}_${cell}_G4plus.matrix.gz --missingDataAsZero && 
    plotProfile --matrixFile $folder_out_result/${bwfile}_${cell}_G4plus.matrix.gz -o $folder_out_result/${bwfile}_${cell}_G4plus.matrix.pdf --averageType median --plotWidth 18 --plotHeight 10 &&
    plotHeatmap --matrixFile $folder_out_result/${bwfile}_${cell}_G4plus.matrix.gz -o $folder_out_result/${bwfile}_${cell}_G4plus.matrix_heatmap.pdf"
    echo "========================"
    
    sbatch --mem 6G --wrap="computeMatrix reference-point -S ${bwfile} -R $path_reference_bed/*${cell}_anno_promoters_G4*minus.prom_coord.bed   -a 1000 -b 1000 --outFileName $folder_out_result/${bwfile}_${cell}_G4minus.matrix.gz --missingDataAsZero && 
    plotProfile --matrixFile $folder_out_result/${bwfile}_${cell}_G4minus.matrix.gz -o $folder_out_result/${bwfile}_${cell}_G4minus.matrix.pdf --averageType median --plotWidth 18 --plotHeight 10 &&
    plotHeatmap --matrixFile $folder_out_result/${bwfile}_${cell}_G4minus.matrix.gz -o $folder_out_result/${bwfile}_${cell}_G4minus.matrix_heatmap.pdf"
    echo "========================"
  done
done  
```

``` bash

#copy locally the outcome 
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2
scp -r  simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_BG4_seq_at_promoter_stratified_by_histones .
```

Average profiles with R script/function locally

``` r
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_plots_BG4_seq_at_promoter_stratified_by_histones')
#NSC_ATAC*matrix*gz
#ESC_ATAC*matrix*gz
#CNCC_ATAC*matrix*gz

ending_label_ESCminus="_ESC_G4minus.matrix.gz"
ending_label_ESCplus="_ESC_G4plus.matrix.gz"

ending_label_CNCCminus="_CNCC_G4minus.matrix.gz"
ending_label_CNCCplus="_CNCC_G4plus.matrix.gz"

ending_label_NSCminus="_NSC_G4minus.matrix.gz"
ending_label_NSCplus="_NSC_G4plus.matrix.gz"

case_label_ESC='ESC_'

case_label_CNCC='CNCC_'

case_label_NSC='NSC_'


average_densities <- function(case_label,ending_label){
  matrix_files <- list.files(".",paste0(case_label,".*",ending_label,"$"))
  matrix_files
  # extract params for initialize matrix
  header_temp <- readLines(gzfile(matrix_files[1]),n=1)
  
  A <- unlist(strsplit(header_temp, split=',\"'))
  A <- gsub('\"','',A)
  A <- gsub('[','',A, fixed = T)
  A <- gsub(']','',A, fixed = T)
  N_rows <- as.numeric(tail(unlist(strsplit(grep('group_boundaries',A, value =T),',|:')),1))
  groups_ind <- as.numeric(tail(unlist(strsplit(grep('group_boundaries',A, value =T),',|:')),-2)) # to exclude first 2
  N_col <- as.numeric(tail(unlist(strsplit(grep('sample_boundaries',A, value =T),',|:')),1))

  temp1_sum <- matrix(0,ncol=N_col,nrow=N_rows)
  for (i in 1:length(matrix_files)){
    header <- readLines(gzfile(matrix_files[i]),n=1)
    temp1 <- read.table(gzfile(matrix_files[i]),skip = 1, stringsAsFactors = F) 
    temp1_sum <- temp1_sum +  temp1[,7:dim(temp1)[2]]
    print(i)
  }
  temp1_avg <- temp1_sum/length(matrix_files)
  temp1_avg_complete <- cbind(temp1[,1:6],temp1_avg)
  file_name <- paste0(case_label,ending_label,'_average.matrix')
  fileConn<-file(file_name)
  writeLines(header,fileConn)
  write.table(temp1_avg_complete, file = file_name, append = T,sep = "\t",quote = F, col.names = F, row.names = F)
  flag='done'
  close(fileConn)
  list_sizes <- list("N_rows"=N_rows,"groups_ind"=groups_ind,"N_col"=N_col)
  return(list_sizes)
}  


average_densities(case_label_ESC,ending_label_ESCminus)
average_densities(case_label_ESC,ending_label_ESCplus)

average_densities(case_label_CNCC,ending_label_CNCCminus)
average_densities(case_label_CNCC,ending_label_CNCCplus)

average_densities(case_label_NSC,ending_label_NSCminus)
average_densities(case_label_NSC,ending_label_NSCplus)
```

copy back to cluster and run deeptools heatmap

``` bash

cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_plots_BG4_seq_at_promoter_stratified_by_histones
gzip *average.matrix
scp *_average.matrix.gz simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_BG4_seq_at_promoter_stratified_by_histones
```

generate Heatmap and plot with deeptools

``` bash
#on cluster
cd /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_BG4_seq_at_promoter_stratified_by_histones
for file in *average*.matrix.gz
do
  sbatch --mem 2G --wrap="plotHeatmap --matrixFile $file  --averageTypeSummaryPlot median -o ${file%%.matrix.gz}.matrix_heatmap.pdf --yMax 0.2"
  echo "========================"
done


# locally 
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_plots_BG4_seq_at_promoter_stratified_by_histones

scp simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_BG4_seq_at_promoter_stratified_by_histones/*average.matrix_heatmap.pdf . 
```

heatmap by promoter histone type

``` bash

prom_type=(bival H3K27me3_only H3K4me3_only relax_bival)
cells=(ESC CNCC NSC)
folder_out_result=/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_BG4_seq_at_promoter_stratified_by_histones2
path_reference_bed=/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/promoter_annotated_by_histones
mkdir $folder_out_result

for prom in ${prom_type[@]}
  do
  
  for cell in ${cells[@]}
  do
  cd /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW
    for bwfile in ${cell}*bw
    do
      sbatch --mem 6G --wrap="computeMatrix reference-point -S ${bwfile} -R $path_reference_bed/${prom}*${cell}_anno_promoters_G4*.prom_coord.bed   -a 1000 -b 1000 --outFileName $folder_out_result/${bwfile}_${prom}_${cell}_G4.matrix.gz --missingDataAsZero && 
      plotProfile --matrixFile $folder_out_result/${bwfile}_${prom}_${cell}_G4.matrix.gz -o $folder_out_result/${bwfile}_${prom}_${cell}_G4.matrix.pdf --averageType median --plotWidth 18 --plotHeight 10 &&
      plotHeatmap --matrixFile $folder_out_result/${bwfile}_${prom}_${cell}_G4.matrix.gz -o $folder_out_result/${bwfile}_${cell}_G4plus.matrix_heatmap.pdf"
      echo "========================"
      
    done
  done  
done
```

copy locally

``` bash
#copy  the outcome  locally
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2
scp -r  simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_BG4_seq_at_promoter_stratified_by_histones2 .
```

``` r
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_plots_BG4_seq_at_promoter_stratified_by_histones2')
#NSC_ATAC*matrix*gz
#ESC_ATAC*matrix*gz
#CNCC_ATAC*matrix*gz



ending_label_ESC_K27<-"_H3K27me3_only_ESC_G4.matrix.gz"
ending_label_ESC_K4<-"_H3K4me3_only_ESC_G4.matrix.gz"
ending_label_ESC_bival<-"w_bival_ESC_G4.matrix.gz"
ending_label_ESC_relax_bival<-"w_relax_bival_ESC_G4.matrix.gz"

ending_label_CNCC_K27<-"_H3K27me3_only_CNCC_G4.matrix.gz"
ending_label_CNCC_K4<-"_H3K4me3_only_CNCC_G4.matrix.gz"
ending_label_CNCC_bival<-"w_bival_CNCC_G4.matrix.gz"
ending_label_CNCC_relax_bival<-"_relax_bival_CNCC_G4.matrix.gz"

ending_label_NSC_K27<-"_H3K27me3_only_NSC_G4.matrix.gz"
ending_label_NSC_K4<-"_H3K4me3_only_NSC_G4.matrix.gz"
ending_label_NSC_bival<-"w_bival_NSC_G4.matrix.gz"
ending_label_NSC_relax_bival<-"w_relax_bival_NSC_G4.matrix.gz"


case_label_ESC='ESC_'

case_label_CNCC='CNCC_'

case_label_NSC='NSC_'

source('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/average_deeptools_densities.R')


average_densities(case_label_ESC,ending_label_ESC_K27)
average_densities(case_label_ESC,ending_label_ESC_K4)
average_densities(case_label_ESC,ending_label_ESC_bival)
average_densities(case_label_ESC,ending_label_ESC_relax_bival)

average_densities(case_label_CNCC,ending_label_CNCC_K27)
average_densities(case_label_CNCC,ending_label_CNCC_K4)
average_densities(case_label_CNCC,ending_label_CNCC_bival)
average_densities(case_label_CNCC,ending_label_CNCC_relax_bival)

average_densities(case_label_NSC,ending_label_NSC_K27)
average_densities(case_label_NSC,ending_label_NSC_K4)
average_densities(case_label_NSC,ending_label_NSC_bival)
average_densities(case_label_NSC,ending_label_NSC_relax_bival)
```

``` bash

cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_plots_BG4_seq_at_promoter_stratified_by_histones2
gzip *average.matrix
scp *_average.matrix.gz simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_BG4_seq_at_promoter_stratified_by_histones2
```

on the cluster

``` bash

cd /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_BG4_seq_at_promoter_stratified_by_histones2
for file in *average*.matrix.gz
do
  sbatch --mem 2G --wrap="plotHeatmap --matrixFile $file  --averageTypeSummaryPlot median -o ${file%%.matrix.gz}.matrix_heatmap.pdf" # --yMax 0.2
  echo "========================"
done
```

locally

``` bash
scp simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_BG4_seq_at_promoter_stratified_by_histones2/*average.matrix_heatmap.pdf . 
```

### Desity plots obtained by looking at histone modification ChIP-seq and stratyfing by promoter types

``` bash
bw_hist_ESC=/scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/rada_iglesia_ESC/trimmed/aligned/
bw_hist_NSC=/scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/xie_NSC/trimmed/aligned/
bw_hist_CNCC=/scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/prescott_CNCC/trimmed/aligned/



prom_type=(bival H3K27me3_only H3K4me3_only relax_bival)
cells=(ESC CNCC NSC)
folder_out_result=/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_histoneMarks_at_promoter_stratified_by_histones_and_G4
path_reference_bed=/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/promoter_annotated_by_histones
mkdir $folder_out_result

for prom in ${prom_type[@]}
  do
  
  for cell in ${cells[@]}
  do
  path="bw_hist_$cell"
  cd "${!path}"
  
    for bwfile in ${cell}*bw
    do
      sbatch --mem 6G --wrap="computeMatrix reference-point -S ${bwfile} -R $path_reference_bed/${prom}*${cell}_anno_promoters_G4*.prom_coord.bed   -a 1000 -b 1000 --outFileName $folder_out_result/${bwfile}_${prom}_${cell}_G4.matrix.gz --missingDataAsZero && 
      plotProfile --matrixFile $folder_out_result/${bwfile}_${prom}_${cell}_G4.matrix.gz -o $folder_out_result/${bwfile}_${prom}_${cell}_G4.matrix.pdf --averageType median --plotWidth 18 --plotHeight 10 &&
      plotHeatmap --matrixFile $folder_out_result/${bwfile}_${prom}_${cell}_G4.matrix.gz -o $folder_out_result/${bwfile}_${prom}_${cell}_G4.matrix_heatmap.pdf"
      echo "========================"
      
    done
  done  
done
```

    ## mkdir: /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW: No such file or directory
    ## bash: line 18: cd: /scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/rada_iglesia_ESC/trimmed/aligned/: No such file or directory
    ## bash: line 22: sbatch: command not found
    ## ========================
    ## bash: line 18: cd: /scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/prescott_CNCC/trimmed/aligned/: No such file or directory
    ## bash: line 22: sbatch: command not found
    ## ========================
    ## bash: line 18: cd: /scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/xie_NSC/trimmed/aligned/: No such file or directory
    ## bash: line 22: sbatch: command not found
    ## ========================
    ## bash: line 18: cd: /scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/rada_iglesia_ESC/trimmed/aligned/: No such file or directory
    ## bash: line 22: sbatch: command not found
    ## ========================
    ## bash: line 18: cd: /scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/prescott_CNCC/trimmed/aligned/: No such file or directory
    ## bash: line 22: sbatch: command not found
    ## ========================
    ## bash: line 18: cd: /scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/xie_NSC/trimmed/aligned/: No such file or directory
    ## bash: line 22: sbatch: command not found
    ## ========================
    ## bash: line 18: cd: /scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/rada_iglesia_ESC/trimmed/aligned/: No such file or directory
    ## bash: line 22: sbatch: command not found
    ## ========================
    ## bash: line 18: cd: /scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/prescott_CNCC/trimmed/aligned/: No such file or directory
    ## bash: line 22: sbatch: command not found
    ## ========================
    ## bash: line 18: cd: /scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/xie_NSC/trimmed/aligned/: No such file or directory
    ## bash: line 22: sbatch: command not found
    ## ========================
    ## bash: line 18: cd: /scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/rada_iglesia_ESC/trimmed/aligned/: No such file or directory
    ## bash: line 22: sbatch: command not found
    ## ========================
    ## bash: line 18: cd: /scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/prescott_CNCC/trimmed/aligned/: No such file or directory
    ## bash: line 22: sbatch: command not found
    ## ========================
    ## bash: line 18: cd: /scratcha/sblab/simeon01/Data/20200326_extern_stem_cells_data/xie_NSC/trimmed/aligned/: No such file or directory
    ## bash: line 22: sbatch: command not found
    ## ========================

copy from cluster to local computer

``` bash

copy locally
```

``` bash
#copy  the outcome  locally
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2
scp -r  simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_histoneMarks_at_promoter_stratified_by_histones_and_G4 .
```

combine replicates into a single map

``` r
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_plots_histoneMarks_at_promoter_stratified_by_histones_and_G4')
#NSC_ATAC*matrix*gz
#ESC_ATAC*matrix*gz
#CNCC_ATAC*matrix*gz





ending_label_CNCC_K27<-"_H3K27me3_only_CNCC_G4.matrix.gz"
ending_label_CNCC_K4<-"_H3K4me3_only_CNCC_G4.matrix.gz"
ending_label_CNCC_bival<-"w_bival_CNCC_G4.matrix.gz"
ending_label_CNCC_relax_bival<-"_relax_bival_CNCC_G4.matrix.gz"

ending_label_NSC_K27<-"_H3K27me3_only_NSC_G4.matrix.gz"
ending_label_NSC_K4<-"_H3K4me3_only_NSC_G4.matrix.gz"
ending_label_NSC_bival<-"w_bival_NSC_G4.matrix.gz"
ending_label_NSC_relax_bival<-"w_relax_bival_NSC_G4.matrix.gz"


case_label_NSC1='NSC_H3K4me3'
case_label_NSC2='NSC_H3K27me3'

case_label_CNCC='CNCC_H3K4me3_'

source('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/average_deeptools_densities.R')

average_densities(case_label_CNCC,ending_label_CNCC_K27)
average_densities(case_label_CNCC,ending_label_CNCC_K4)
average_densities(case_label_CNCC,ending_label_CNCC_bival)
average_densities(case_label_CNCC,ending_label_CNCC_relax_bival)

average_densities(case_label_NSC1,ending_label_NSC_K27)
average_densities(case_label_NSC1,ending_label_NSC_K4)
average_densities(case_label_NSC1,ending_label_NSC_bival)
average_densities(case_label_NSC1,ending_label_NSC_relax_bival)

average_densities(case_label_NSC2,ending_label_NSC_K27)
average_densities(case_label_NSC2,ending_label_NSC_K4)
average_densities(case_label_NSC2,ending_label_NSC_bival)
average_densities(case_label_NSC2,ending_label_NSC_relax_bival)
```

``` bash

/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_plots_histoneMarks_at_promoter_stratified_by_histones_and_G4
gzip *average.matrix
scp *_average.matrix.gz simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_histoneMarks_at_promoter_stratified_by_histones_and_G4
```

on the cluster

``` bash

cd /scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_histoneMarks_at_promoter_stratified_by_histones_and_G4
for file in *average*.matrix.gz
do
  sbatch --mem 2G --wrap="plotHeatmap --matrixFile $file  --averageTypeSummaryPlot median -o ${file%%.matrix.gz}.matrix_heatmap.pdf" # --yMax 0.2
  echo "========================"
done
```

locally

``` bash
cd /Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_plots_histoneMarks_at_promoter_stratified_by_histones_and_G4
scp simeon01@10.20.236.34:/scratcha/sblab/simeon01/Data/20200512_KatieStemCells_BW/density_plots_histoneMarks_at_promoter_stratified_by_histones_and_G4/*average.matrix_heatmap.pdf . 
```

NCP promoter capture Hi-C data for NPC
--------------------------------------

``` bash
mkdir /scratcha/sblab/simeon01/Data/20072020_npc_hic
cd /scratcha/sblab/simeon01/Data/20072020_npc_hic

# data have been downloaded from geo
GSE86189_npc.po.all.txt
GSE86189_npc.pp.all.txt

# select first and secon column and remove header
cut -f 1 GSE86189_npc.pp.all.txt | tail -n +2 > temp1_pp.txt
cut -f 2 GSE86189_npc.pp.all.txt | tail -n +2 > temp2_pp.txt

# cat together the pp and substitute . with \t
cat *_pp.txt | sed 's/\./\t/g' | sort | uniq | sortBed -i - > GSE86189_npc.pp.all.bed
rm *temp*

# do the same for po file
cut -f 1 GSE86189_npc.po.all.txt | tail -n +2 > temp1_po.txt
cut -f 2 GSE86189_npc.po.all.txt | tail -n +2 > temp2_po.txt

# cat together the pp and substitute . with \t
cat *_po.txt | sed 's/\./\t/g' | sort | uniq | sortBed -i - > GSE86189_npc.po.all.bed
rm *temp*

# cat together also pp and po
cat GSE86189_npc.pp.all.bed GSE86189_npc.po.all.bed | sort | uniq | sortBed - i - > GSE86189.pp_and_po.bed


# create an annotation file
touch GSE86189_HiC_NPC.pp.po.anno
awk '{print $1"\t"$2"\t"$3"\tGSE86189_npc.pp"}' GSE86189_npc.pp.all.bed >> GSE86189_HiC_NPC.pp.po.anno
awk '{print $1"\t"$2"\t"$3"\tGSE86189_npc.po"}' GSE86189_npc.po.all.bed >> GSE86189_HiC_NPC.pp.po.anno 

touch GSE86189_HiC_NPC.pp_and_po.anno
awk '{print $1"\t"$2"\t"$3"\tGSE86189_npc.pp_and_po"}' GSE86189.pp_and_po.bed >> GSE86189_HiC_NPC.pp_and_po.anno

cp /scratcha/sblab/simeon01/Data/20072020_npc_hic/*anno /scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells
```

curate H3K27ac annotations for the 3 cell lines
-----------------------------------------------

Extract H3K27ac annotations and subtract promoter annotations

``` bash

path_anno='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells'

ESC_anno=Rada_Iglesias_lifted_hg38.bed
CNCC_anno=CNCC_annotations_hg38_no_.bed
NEC_anno=NEC_GSE24447_annotations_hg38_no_.bed
NSC_anno=encodeH1_histone_released_H3histMark_hg38_no_.bed
prom=/scratcha/sblab/simeon01/reference_genomes/hg38/gencode.v28.TSS.slop1000.bed

annotation_list1=($ESC_anno $CNCC_anno $NSC_anno)
cd $path_anno
touch H3K27ac_hg38.ESC.CNCC_NSC.anno
grep "H3K27ac" $ESC_anno | intersectBed -a - -b $prom -v | awk '{print $1"\t"$2"\t"$3"\tH3K27ac_ESC"}' >>H3K27ac_hg38.ESC.CNCC_NSC.anno 

grep "K27ac" $CNCC_anno | intersectBed -a - -b $prom  -v | awk '{print $1"\t"$2"\t"$3"\tH3K27ac_CNCC"}' >>H3K27ac_hg38.ESC.CNCC_NSC.anno 

grep "H3K27ac" $NSC_anno | intersectBed -a - -b $prom -v | awk '{print $1"\t"$2"\t"$3"\tH3K27ac_NSC"}' >>H3K27ac_hg38.ESC.CNCC_NSC.anno 
```

Revisit H1ESC\_Xie\_2013 enhancers and expand them to 200bp

``` bash
cd /scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/additional_annotations_various_stem_cells
ls *Xie_2013*_hg38.hg38.bed  

touch Xie_2013_enhancers_H1ESC_NSC.anno

for file in *Xie_2013*_hg38.hg38.bed
do
  label=${file%%_hg18*}
  awk -v a="$label" '{print $1"\t"int(($3+$2)/2+0.49)"\t"int(($3+$2)/2+0.49)"\t"a}' $file > temp
  #slop bed with -b 100
  slopBed -i temp -g /scratcha/sblab/simeon01/reference_genomes/hg38/hg38.genome -b 100 >> Xie_2013_enhancers_H1ESC_NSC.anno
  rm temp  
done

cp Xie_2013_enhancers_H1ESC_NSC.anno /scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells/
```

produce fold enrichment and proportion of regions overlapping

``` bash


#scenario=(scenario1_multi2 scenario2_5out9 scenario3_multi2_extended_1kb)
scenario=(scenario1_multi2)
workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed
#genomic_annotation=genomic_regions_UCSC_hg38_igenomes.bed6.bed

H3K27ac_cells=H3K27ac_hg38.ESC.CNCC_NSC.anno
Xie_2013_enhancers=Xie_2013_enhancers_H1ESC_NSC.anno

NPC_promoter_capture_po_pp=GSE86189_HiC_NPC.pp.po.anno
NPC_promoter_capture_po_and_pp=GSE86189_HiC_NPC.pp_and_po.anno


annotation_list1=($H3K27ac_cells $Xie_2013_enhancers $NPC_promoter_capture_po_pp $NPC_promoter_capture_po_and_pp)

#/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/bed_intersections_common/*C.bed
#annotation_list1=($genomic_annotation)

cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2
bed_to_use=`ls *bed | grep -v shuffle | grep -v intermediate`
mkdir /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Fold_enrichment_scenario_additional_July20
path_anno='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells'
for file in ${bed_to_use[@]}
do
  echo $file
  for anno in ${annotation_list1[@]}
    do
      echo $anno
        
      # collect the features in the annotation file
      res_feature_annot_list=${file%%.bed}_intervals_${anno%%.bed}.overlaps.tot_${n_feat}_featNames
      feat=`cut -f 4 $path_anno/$anno | sort | uniq`
      cut -f 4 $path_anno/$anno | sort | uniq >$res_feature_annot_list 
      n_feat=`cut -f 4 $path_anno/$anno | sort | uniq | wc -l`
      res_feature_annot=${file%%.bed}_intervals_${anno%%.bed}.overlaps.tot_${n_feat}_feat
      res_feature_annot_fin=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Fold_enrichment_scenario_additional_July20/${file%%.bed}_intervals_${anno%%.bed}.overlaps.tot_${n_feat}_feat.tab
      # for each feature in the annotation count number of overlaps with bed (how many intervals in bed are in feature (unique)
      for f in ${feat[@]}
        do
          echo $f
          cat $path_anno/$anno  | awk -v var=$f '$4 ~ var' > feature_annot_tmp.bed
          intersectBed -a $file -b feature_annot_tmp.bed -wa | sort | uniq | wc -l >> $res_feature_annot
        done
        # reformat results into a single file
        paste $res_feature_annot_list $res_feature_annot >> temp_res_feat
        wc -l $file  | awk '{print $2"\t"$1}'| cat - temp_res_feat > $res_feature_annot_fin
        rm $res_feature_annot_list
        rm $res_feature_annot
        rm temp_res_feat
        rm feature_annot_tmp.bed
    done
done





#########################
cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$path_anno/${anno} --segments=${file} -S /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Fold_enrichment_scenario1/${file%%.bed}_${anno%%.bed}.dat -t 4"
    echo ${file%%.bed}_${anno%%.bed}.dat
    echo $cmd_gat
    sbatch --mem 6G --wrap "$cmd_gat"

#############
```

``` bash

scenario=(scenario1_multi2)
workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed
#genomic_annotation=genomic_regions_UCSC_hg38_igenomes.bed6.bed

H3K27ac_cells=H3K27ac_hg38.ESC.CNCC_NSC.anno
Xie_2013_enhancers=Xie_2013_enhancers_H1ESC_NSC.anno

NPC_promoter_capture_po_pp=GSE86189_HiC_NPC.pp.po.anno
NPC_promoter_capture_po_and_pp=GSE86189_HiC_NPC.pp_and_po.anno

workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed

annotation_list1=($H3K27ac_cells $Xie_2013_enhancers $NPC_promoter_capture_po_pp $NPC_promoter_capture_po_and_pp)

#cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2
#bed_to_use=`ls *.bed | grep -v shuffle | grep -v intermediate`

cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/bed_intersections_common
bed_to_use=`ls *C.bed`
#

path_anno='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells'
for file in ${bed_to_use[@]}
do
  echo $file
  for anno in ${annotation_list1[@]}
    do
      echo $anno
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$path_anno/${anno} --segments=${file} -S /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Fold_enrichment_scenario_additional_July20/${file%%.bed}_${anno%%.bed}.dat -t 4"
    echo ${file%%.bed}_${anno%%.bed}.dat
    echo $cmd_gat
    sbatch --time 01:00:00 --mem 6G --wrap "$cmd_gat"
    done
done
```

enrichements at PAX6 NSC
------------------------

``` bash

cd /scratcha/sblab/simeon01/Data/21072020_pax6ChIP_NEC
awk '{print $1"\t"$2"\t"$3"\tpax6ChIP_NEC"}' Pax6.NEC.default_peaks.top15795.qval25.bed > Pax6.NEC.default_peaks.top15795.qval25.anno

cp Pax6.NEC.default_peaks.top15795.qval25.anno /scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells

cd /scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells

#script_for_count_feat_overlaps_pax6nec_pariwise_interaction.sh

# prepare file for enrichment run#

sbatch --time 1:00:00 script_for_count_feat_overlaps_pax6nec_pariwise_interaction.sh

cat script_for_count_feat_overlaps_pax6nec_pariwise_interaction.sh
#!/bin/bash
#SBATCH --mem 4G


scenario=(scenario1_multi2)
workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed
#genomic_annotation=genomic_regions_UCSC_hg38_igenomes.bed6.bed

PAS6_NEC=Pax6.NEC.default_peaks.top15795.qval25.anno

annotation_list1=($PAS6_NEC)

#/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/bed_intersections_common/*C.bed
#annotation_list1=($genomic_annotation)

#bed_to_use=`ls *bed | grep -v shuffle | grep -v intermediate`
cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/bed_intersections_common
bed_to_use=`ls *C.bed`

path_anno='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells'
for file in ${bed_to_use[@]}
do
  echo $file
  for anno in ${annotation_list1[@]}
    do
      echo $anno
        
      # collect the features in the annotation file
      res_feature_annot_list=${file%%.bed}_intervals_${anno%%.bed}.overlaps.tot_${n_feat}_featNames
      feat=`cut -f 4 $path_anno/$anno | sort | uniq`
      cut -f 4 $path_anno/$anno | sort | uniq >$res_feature_annot_list 
      n_feat=`cut -f 4 $path_anno/$anno | sort | uniq | wc -l`
      res_feature_annot=${file%%.bed}_intervals_${anno%%.bed}.overlaps.tot_${n_feat}_feat
      res_feature_annot_fin=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Fold_enrichment_scenario_additional_PAX6Nec/${file%%.bed}_intervals_${anno%%.bed}.overlaps.tot_${n_feat}_feat.tab
      # for each feature in the annotation count number of overlaps with bed (how many intervals in bed are in feature (unique)
      for f in ${feat[@]}
        do
          echo $f
          cat $path_anno/$anno  | awk -v var=$f '$4 ~ var' > feature_annot_tmp.bed
          intersectBed -a $file -b feature_annot_tmp.bed -wa | sort | uniq | wc -l >> $res_feature_annot
        done
        # reformat results into a single file
        paste $res_feature_annot_list $res_feature_annot >> temp_res_feat
        wc -l $file  | awk '{print $2"\t"$1}'| cat - temp_res_feat > $res_feature_annot_fin
        rm $res_feature_annot_list
        rm $res_feature_annot
        rm temp_res_feat
        rm feature_annot_tmp.bed
    done
done

cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/bed_intersections_common
bed_to_use=`ls *C.bed`

path_anno='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells'
for file in ${bed_to_use[@]}
do
  echo $file
  for anno in ${annotation_list1[@]}
    do
      echo $anno
        
      # collect the features in the annotation file
      res_feature_annot_list=${file%%.bed}_intervals_${anno%%.bed}.overlaps.tot_${n_feat}_featNames
      feat=`cut -f 4 $path_anno/$anno | sort | uniq`
      cut -f 4 $path_anno/$anno | sort | uniq >$res_feature_annot_list 
      n_feat=`cut -f 4 $path_anno/$anno | sort | uniq | wc -l`
      res_feature_annot=${file%%.bed}_intervals_${anno%%.bed}.overlaps.tot_${n_feat}_feat
      res_feature_annot_fin=/scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Fold_enrichment_scenario_additional_PAX6Nec/${file%%.bed}_intervals_${anno%%.bed}.overlaps.tot_${n_feat}_feat.tab
      # for each feature in the annotation count number of overlaps with bed (how many intervals in bed are in feature (unique)
      for f in ${feat[@]}
        do
          echo $f
          cat $path_anno/$anno  | awk -v var=$f '$4 ~ var' > feature_annot_tmp.bed
          intersectBed -a $file -b feature_annot_tmp.bed -wa | sort | uniq | wc -l >> $res_feature_annot
        done
        # reformat results into a single file
        paste $res_feature_annot_list $res_feature_annot >> temp_res_feat
        wc -l $file  | awk '{print $2"\t"$1}'| cat - temp_res_feat > $res_feature_annot_fin
        rm $res_feature_annot_list
        rm $res_feature_annot
        rm temp_res_feat
        rm feature_annot_tmp.bed
    done
done
```

run gat enrichments

``` bash


scenario=(scenario1_multi2)
workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed
#genomic_annotation=genomic_regions_UCSC_hg38_igenomes.bed6.bed
PAS6_NEC=Pax6.NEC.default_peaks.top15795.qval25.anno

annotation_list1=($PAS6_NEC)

cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2
bed_to_use=`ls *.bed | grep -v shuffle | grep -v intermediate`


path_anno='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells'
for file in ${bed_to_use[@]}
do
  echo $file
  for anno in ${annotation_list1[@]}
    do
      echo $anno
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$path_anno/${anno} --segments=${file} -S /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Fold_enrichment_scenario_additional_PAX6Nec/${file%%.bed}_${anno%%.bed}.dat -t 4"
    echo ${file%%.bed}_${anno%%.bed}.dat
    echo $cmd_gat
    sbatch --time 01:00:00 --mem 6G --wrap "$cmd_gat"
    done
done

cd /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/bed_intersections_common
bed_to_use=`ls *C.bed`
#

path_anno='/scratcha/sblab/simeon01/reference_genomes/epigenetic_annotations_stem_cells'
for file in ${bed_to_use[@]}
do
  echo $file
  for anno in ${annotation_list1[@]}
    do
      echo $anno
      cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$path_anno/${anno} --segments=${file} -S /scratcha/sblab/simeon01/Data/Katie_high_level_analysis_Oct2019/scenario1_multi2/Fold_enrichment_scenario_additional_PAX6Nec/${file%%.bed}_${anno%%.bed}.dat -t 4"
    echo ${file%%.bed}_${anno%%.bed}.dat
    echo $cmd_gat
    sbatch --time 01:00:00 --mem 6G --wrap "$cmd_gat"
    done
done
```

Peak sizes - all and G4 regions in common among cells
-----------------------------------------------------

Peak size R script are the following and they also integrate external annotation in case we want to screen specific functional categories.

[R script to assess G4 wideness based on G4 density information](/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/density_at_common_G4_regions/Plot_combine_densities.R)

[R script to assess G4 wideness based on G4 peak size and/or GO annotations](/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/temp1_G4_all_regions_GO.R)

[R script to assess G4 peak size at regions also marked by other histone modifications](/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/temp1_G4_and_epigenetic.R)

[R script to create bed coordinates of promoters stratified by histone modifications and G4 presence/absence (bival, H3K4me3...etc)](/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/20200421_annotated_promoters_using_histones_G4/BG4_based_promoter_annotations_for_atacDensityPlots_at_TSS.R)

[Rmd to generate density plots at promoters stratified by histone modifications and G4 presence/absence (bival, H3K4me3...etc)](/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/20200421_annotated_promoters_using_histones_G4/Analysis_G4_location_in_respect_to_histone_marks.Rmd)
