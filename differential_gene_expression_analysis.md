Differential Gene Expression Analysis
================
Angela Simeone

Differential Analysis have been performed using the expected counts obtained after alignment (rsem using bowtie2) to transcriptome. We used counts obtained in `genes.results` files. <br /> <br />

``` bash

# PROCESSING KATIE'S RNA SEQ AND COMPARE IT TO EXTERNAL DATA
# NOTE: ALL DATA  need to be ALIGNED AT THE SAME WAY TO THE SAME REFERENCE , USING THE SAME ANNOTATIONS FILE

 ==================================================================================================================
#
# LOAD PACKAGES and FUNCTIONS
#
# ==================================================================================================================
library('edgeR')
library('ggplot2')
library('gridExtra')
library('corrplot')
library('ComplexHeatmap')
require('openxlsx')


# FUNCTION TO IMPORT RSEM results - customized for RSEM - gene.results file
extract_counts_TPM <- function(rsem_lists_to_import,condition_names_to_assign){
  Table_expected_counts <- c()
  Table_TPM <- c()
  for (i in 1:length(rsem_lists_to_import))
  {
    temp_read <- read.table(rsem_lists_to_import[i], stringsAsFactors = F, header = T)
    Table_expected_counts <- cbind(Table_expected_counts,temp_read$expected_count)
    Table_TPM <- cbind(Table_TPM,temp_read$TPM)
  }
  colnames(Table_expected_counts) <- condition_names_to_assign
  colnames(Table_TPM) <- condition_names_to_assign
  rownames(Table_expected_counts) <- temp_read$gene_id
  rownames(Table_TPM) <- temp_read$gene_id
  # store effected length
  gene_eff_lengths <- temp_read$effective_length
  names(gene_eff_lengths) <- temp_read$gene_id
  rsem_tables = list('Table_expected_counts' = Table_expected_counts,'Table_TPM' = Table_TPM, 'genes_effective_lengths' = gene_eff_lengths)
  return(rsem_tables)
}

# FUNCTION THAT CAN BE CALLED TO ITERATE ALL COMBINATIONS OF DIFFERENTIAL ANALYSIS
DE_analysis_routine <- function(y_DEGL,fit_of_interest, contrast_of_interest,produce_GO,flag_for_name_to_use) {

  
  All_DE_list_genes <- c()
  #DE_analysis_routine(y,fit,my.contrast)
  y_DEGL <- y
  fit_of_interest <- fit
  contrast_of_interest <- my.contrast
  
  for (i in 1:length(colnames(contrast_of_interest))) {
    print(i)
    print(length(colnames(contrast_of_interest)))
    
    qlf <- glmQLFTest(fit_of_interest, contrast=contrast_of_interest[,i])
    
    de_pairing_table_curr_case <- topTags(qlf, n= nrow(y_DEGL$counts))$table
    
    de_fdr05_lfg2 = de_pairing_table_curr_case[which(de_pairing_table_curr_case$FDR <=0.05 & abs(de_pairing_table_curr_case$logFC)>=2),]
    
    #decideTestsDGE(qlf,adjust.method="BH", p.value = 0.05,lfc=0)
    
    list_genes_up <- rownames(de_fdr05_lfg2[(which(de_fdr05_lfg2$logFC>0)),])
    
    list_genes_up_ensembl <- matrix(unlist(strsplit(list_genes_up,split = "\\.[0-9]{1,2}_", perl = TRUE)),length(list_genes_up), byrow = T)
    
    list_genes_down <- rownames(de_fdr05_lfg2[(which(de_fdr05_lfg2$logFC<0)),])
    
    list_genes_down_ensembl <- matrix(unlist(strsplit(list_genes_down,split = "\\.[0-9]{1,2}_", perl = TRUE)),length(list_genes_down), byrow = T)
    
    output_decide_test <- decideTests(qlf,lfc = 2)
    
    summary_report <- summary(output_decide_test) # this are default --> adjust.method = "BH", p.value = 0.05
    
    All_DE_list_genes <- c(All_DE_list_genes,union(list_genes_up,list_genes_down))
    
    Summary_current_case = list('summary' = summary_report,
                                'list_genes_UP' = list_genes_up,
                                'list_genes_UP_ensembl' = list_genes_up_ensembl,
                                'list_genes_DOWN' = list_genes_down,
                                'list_genes_DOWN_ensembl' = list_genes_down_ensembl)
    names(Summary_current_case) = c('summary','list_genes_UP','list_genes_UP_ensembl',
                                    'list_genes_DOWN','list_genes_DOWN_ensembl')
    
    EdgeR_results <- list(de_pairing_table_curr_case,de_fdr05_lfg2)
    names(EdgeR_results) <- c('full_table','FDR05_LFG2')
    
    write.table(output_decide_test, file = paste0(flag_for_name_to_use,"_summary_list_decide_test_",colnames(my.contrast)[i],".csv"),sep = ",",quote = F, col.names = NA)
    write.xlsx(Summary_current_case, file = paste0(flag_for_name_to_use,"_summary_list_",colnames(my.contrast)[i],".xlsx"))
    write.xlsx(EdgeR_results,rowNames=T,rowNames = T,file = paste0(flag_for_name_to_use,"_summary_edgeR_",colnames(my.contrast)[i],".xlsx" ))
    
    
    pdf(paste0(flag_for_name_to_use,"_",colnames(my.contrast)[i],'.MD_plot.pdf'))
    plotMD(qlf, status=output_decide_test,col=c("blue","black","red"))
    title(paste0("\n\n UP:",length(list_genes_up)," - DOWN: ",length(list_genes_down)))
    dev.off()
    
    
  }
  return(All_DE_list_genes)
}

# ==================================================================================================================
#
# LOAD DATA
#
# ==================================================================================================================

# localPath
path_data <- '/Users/simeon01/Documents/Katie/2018_09_03/RNA_seq'

setwd(path_data)

# customized analysis for Katies RNAseq
# read katie description of the samples
katie_description <- read.csv("~/Documents/Katie/2018_09_03/RNA_seq/20180703_stemcell_RNA_labels.csv", stringsAsFactors = F)
katie_description_df <- as.data.frame(katie_description)

# read filenames
genes_results_file <- list.files(path_data,pattern = 'genes.results')

# generate list of condition names
condition_names <- gsub(".genes.results", "", genes_results_file)

# generate table describing experiments
table_experiment_info <- matrix(nrow = length(condition_names),ncol = 5)
colnames(table_experiment_info) <- c('sample','sample_num','cell','cell_type','origin') 
for (i in 1:length(condition_names)){
  temp_info <- strsplit(condition_names[i],'_')  
  table_experiment_info[i,1:length(temp_info[[1]])] <- temp_info[[1]]
}

# create a 6th and 7th columns with lablel - we use this for identify replicates  
# and create factors to use in the differential expression analysis
table_experiment_info <- cbind(table_experiment_info, paste(table_experiment_info[,3],table_experiment_info[,4],table_experiment_info[,5], sep = "_"))
table_experiment_info <- cbind(table_experiment_info, paste(table_experiment_info[,3],table_experiment_info[,4], sep = "_"))
table_experiment_info_df <- as.data.frame(table_experiment_info) # data frame

#merge Katie description in order to import also additional informations
table_experiments_info_plus_katie_description <- merge(table_experiment_info_df,katie_description_df, by.x = "sample_num",by.y = "SampleID")

updated_table <- cbind(table_experiments_info_plus_katie_description[,1:7],table_experiments_info_plus_katie_description$pool)
updated_table[,6] <- paste0(table_experiments_info_plus_katie_description[,6],'_p',table_experiments_info_plus_katie_description$pool)
colnames(updated_table)[6:7] <- c("rep","cell_type")

# check number of replicates for each case
replicates_table <- table(updated_table[,6])

# import all katies data: Expected_counts and TPM
katies_rsem_data <- extract_counts_TPM(genes_results_file,condition_names)


inds_to_exclude <- c(grep("sample_24",colnames(katies_rsem_data$Table_TPM)),
                    grep("sample_27",colnames(katies_rsem_data$Table_TPM)),
                    grep("sample_29",colnames(katies_rsem_data$Table_TPM)),
                    grep("sample_30|sample_31|sample_32|sample_33|sample_34",colnames(katies_rsem_data$Table_TPM)),
                    grep("mesenchy",colnames(katies_rsem_data$Table_TPM)),
                    grep("neurosp",colnames(katies_rsem_data$Table_TPM)))


AA_sel_Table_TPM <- katies_rsem_data$Table_TPM[,-inds_to_exclude]


# corr for tech_rep_ESC
ESC_rep4_tech <- AA_sel_Table_TPM[, grep("sample_35|sample_36|sample_7|sample_8",colnames(AA_sel_Table_TPM))]
cor(ESC_rep4_tech)
# sample_35_H9_hESC sample_36_H9_hESC sample_7_H9_hESC sample_8_H9_hESC
# sample_35_H9_hESC         1.0000000         0.9918697        0.9542460        0.9754268
# sample_36_H9_hESC         0.9918697         1.0000000        0.9761753        0.9844257
# sample_7_H9_hESC          0.9542460         0.9761753        1.0000000        0.9888634
# sample_8_H9_hESC          0.9754268         0.9844257        0.9888634        1.0000000
cor(AA_sel_Table_TPM[, grep("sample_5|sample_6",colnames(AA_sel_Table_TPM))])


# ==================================================================================================================
#
# DIFFERENTIAL ANALYSIS on KATIE's DATA
#
# ==================================================================================================================
setwd(path_data)


##### in case we mananage technical replicates we consider the pool information now also store in the field n. 6
colnames(katies_rsem_data$Table_expected_counts) <- factor(updated_table[,6]) # reassign the colnames to keep track of the tecnical replicates information

# manage technical replicates
Table_expected_counts_sumTechRep <- sumTechReps(katies_rsem_data$Table_expected_counts)

# create a matrix where 1st column is the gene 'effective' length and the rest of the colkumns are the expected counts obteined after summing technical replicates
Matrix_length_counts <- as.matrix(cbind(katies_rsem_data$genes_effective_lengths,Table_expected_counts_sumTechRep))
colnames(Matrix_length_counts)[1] <- "Length"

# remove rows with only zero values & rows where the Lenght (col1) is zero -- then compute TPM
zeros_to_remove <- union(which(Matrix_length_counts[,1]==0), which(rowSums(Matrix_length_counts[,-1])==0))

# update matrix without rows with zero values
Matrix_length_counts_updated <- Matrix_length_counts[-zeros_to_remove,]

# compute TPM with calc.tpm()
TPM_sumTechRep_katie <- calc.tpm(Matrix_length_counts_updated)

# assign col names
colnames(TPM_sumTechRep_katie) <- colnames(Matrix_length_counts_updated)[-1]

# select counts to print to file for code submission
counts_to_print_for_code_ind <- grep('sample_25|sample_26|sample_28|NSC|H9_hESC',colnames(katies_rsem_data$Table_expected_counts))
counts_to_print_for_code <- katies_rsem_data$Table_expected_counts[,counts_to_print_for_code_ind]
write.table(counts_to_print_for_code,file= 'counts_individual_libs.csv',quote = F, col.names = NA)

# average now the tecnical replicates TPM and compare to the case when counts have been summed up
# use this when all experiments are used
TechReps_names_unique <- unique(colnames(katies_rsem_data$Table_expected_counts))
updated_meanTechReps_TPM = matrix(nrow = nrow(katies_rsem_data$Table_expected_counts),ncol = length(TechReps_names_unique))


for (i in 1:length(TechReps_names_unique)) {
  temp_ind <- which(colnames(katies_rsem_data$Table_expected_counts) == TechReps_names_unique[i])
  n_rep <- length(temp_ind)
  if (n_rep >1) {
    print(i)
    updated_meanTechReps_TPM[,i] <- rowMeans(katies_rsem_data$Table_TPM[,temp_ind])
  }
  else{
    updated_meanTechReps_TPM[,i] <- katies_rsem_data$Table_TPM[,temp_ind]
  }
}
colnames(updated_meanTechReps_TPM) <- TechReps_names_unique
rownames(updated_meanTechReps_TPM) <- rownames(katies_rsem_data$Table_expected_counts)

#plot and print correlation  to file
# png('Rplot_TPM_versus_meanTPM.png',width = 1640, height = 1640)
# p <- list()
# 
# for (i in 1:dim(TPM_sumTechRep_katie)[2]) {
#   # compute spearman correlation
#   Spearman_corr <- round(cor(TPM_sumTechRep_katie[,i],
#                              updated_meanTechReps_TPM[rownames(updated_meanTechReps_TPM)%in%rownames(TPM_sumTechRep_katie),i], method = 'spearman'),3)
#   print(Spearman_corr)
#   # initialize a dataframe weitht the two types of calculated TPM
#   dat <- data.frame(TPM_after_sumRep = TPM_sumTechRep_katie[,i],
#                     meanTPMreps = updated_meanTechReps_TPM[rownames(updated_meanTechReps_TPM)%in%rownames(TPM_sumTechRep_katie),i])
#   
#   p[[i]] <- ggplot(dat, aes(x=TPM_after_sumRep, y=meanTPMreps)) +
#     geom_point(shape=1) +    # Use hollow circles
#     geom_smooth(method=lm)  + 
#     ggtitle(paste0(Spearman_corr,',',colnames(TPM_sumTechRep_katie)[i]))
#   
# }
# do.call(grid.arrange,p)
# dev.off()

# update with technical replicates summed
table_experiments_after_sumTechRep <- matrix(unlist(strsplit(colnames(Table_expected_counts_sumTechRep),'_')),nrow = length(unique(updated_table[,6])), byrow = T)

# define group of cells - H9_ESC, H9_NSC, H9_CNCC etc
group_cells <- factor(paste(table_experiments_after_sumTechRep[,1],table_experiments_after_sumTechRep[,2],sep='_'))
# group_cells
# [1] H1_hESC        H9_NSC         H9_NSC         H9_NSC         H9_neurosphere H9_neurosphere H9_neurosphere H9_neurosphere H9_neurosphere H9_hESC       
# [11] H1_hESC        H9_hESC        H9_hESC        H9_mesenchymal H9_mesenchymal H9_CNCC        H9_CNCC        H9_CNCC        H9_CNCC        H9_CNCC       
# [21] H9_CNCC        H1_hESC        H9_CNCC        H9_CNCC        H9_hESC        H1_hESC        H9_hESC        H9_NSC        
# Levels: H1_hESC H9_CNCC H9_hESC H9_mesenchymal H9_neurosphere H9_NSC


inds_to_exclude <- c(grep("H9_CNCC_A1_p15",colnames(Table_expected_counts_sumTechRep)),
                     grep("H9_CNCC_A3_p18",colnames(Table_expected_counts_sumTechRep)),
                     grep("H9_CNCC_4_p20",colnames(Table_expected_counts_sumTechRep)),
                     grep("H9_CNCC_A1_p21",colnames(Table_expected_counts_sumTechRep)),
                     grep("H9_CNCC_A5_p22",colnames(Table_expected_counts_sumTechRep)))

Table_expected_counts_sumTechRep_selected_cells <- Table_expected_counts_sumTechRep[,-inds_to_exclude]
group_cells_selected_cell <- factor(paste(table_experiments_after_sumTechRep[,1],table_experiments_after_sumTechRep[,2],sep='_')[-inds_to_exclude])


# create DGEList
y <- DGEList(counts = Table_expected_counts_sumTechRep_selected_cells,group = group_cells_selected_cell)
# y <- DGEList(counts = Table_expected_counts_sumTechRep,group = group_cells)

# CPM
cpm.y <- cpm(y)

average_cmp.y_by_cell_type <- cbind(rowMeans(cpm.y[,grep('ESC',colnames(cpm.y))]),
                                    rowMeans(cpm.y[,grep('CNCC',colnames(cpm.y))]),
                                    rowMeans(cpm.y[,grep('NSC',colnames(cpm.y))]))

#filtering: selecting only genes have have average cpm > 1 in each of the group
ind_selected <- rowSums(average_cmp.y_by_cell_type >= 1) >= 1
y_orig <- y
cpm.y_orig <- cpm.y
y = y[ind_selected,] # at least 1 cell type has a average cpm > 1
cpm.y <- cpm.y[ind_selected,]

library(ComplexHeatmap)
library(circlize)
group_cells_selected_cell_df <- as.data.frame(group_cells_selected_cell)
ha = HeatmapAnnotation(cell_type = group_cells_selected_cell_df$group_cells_selected_cell , annotation_name_side = "left",
                       col = list(cell_type=c("H1_hESC" = "pink", "H9_CNCC" = "red", 
                                  "H9_hESC" = "magenta", "H9_mesenchymal"="grey",
                                  "H9_neurosphere" = "green", "H9_NSC"="blue")))
dim(cpm.y)

cpm.y_rowAvg <- rowMeans(cpm.y)
selected_cpm.y <- cpm.y[which(cpm.y_rowAvg>0 & cpm.y_rowAvg<quantile(cpm.y_rowAvg,0.90)),]
selected_cpm.y_scaled<- scale(selected_cpm.y,scale = T)
#for (i in seq(1,100,10)){
i=11
set.seed(i)
ht_list = Heatmap(selected_cpm.y_scaled, name = "scaled\nnorm_expr", row_km = 4, column_km =6,
                  col = colorRamp2(c(-1.75, 0, 1.75), c("green", "white", "red")),
                  top_annotation = ha, 
                  show_column_names = FALSE,show_row_names = FALSE, row_title = NULL, show_row_dend = FALSE,show_column_dend = T) 
pdf(paste0('clust_and_heatmap_rna-seq_cells_seed_',i,'.pdf'))
draw(ht_list)
dev.off()



# seed 11 was the choosen one ==== replot only ESC, CNCC and NSC
ind_to_select_for_heatmap <- grep('H9_hESC|H9_CNCC|H9_NSC',colnames(selected_cpm.y_scaled))

group_cells_selected_cell_df2 <- droplevels(as.factor(group_cells_selected_cell[ind_to_select_for_heatmap]))


ha2 = HeatmapAnnotation(cell_type = group_cells_selected_cell_df2, annotation_name_side = "left",
                       col = list(cell_type=c("H9_CNCC" = "deepskyblue3", 
                                              "H9_hESC" = "firebrick1",
                                              "H9_NSC"="chartreuse4")))

#for (i in seq(1,100,10)){
i=11
set.seed(i)
ht_list2 = Heatmap(selected_cpm.y_scaled[,ind_to_select_for_heatmap], name = "scaled\nnorm_expr", row_km = 4,col = colorRamp2(c(-1.75, 0, 1.75), c("blue", "white", "red")),
                  top_annotation = ha2, 
                  show_column_names = T,show_row_names = FALSE, row_title = NULL, show_row_dend = FALSE,show_column_dend = T) 
pdf(paste0('Selected_clust_and_heatmap_rna-seq_cells_seed_',i,'.pdf'))
draw(ht_list2)
dev.off()
#}

# other strategy for filtering
# ind_to_keep <- c()
# for (i in 1:length(unique(group_cells))){
#   
#   curr_group <- unique(group_cells)[i]
#   print(length(grep(curr_group,colnames(cpm.y))))
#   
#   if (length(grep(curr_group,colnames(cpm.y)))>1){ 
#     ind_to_keep_temp <- which(rowSums(cpm.y[,grep(curr_group,colnames(cpm.y))] >= 1) >= 1)
#   }
#   else if(length(grep(curr_group,colnames(cpm.y)))==1)
#   {
#     ind_to_keep_temp <- which(cpm.y[,grep(curr_group,colnames(cpm.y))] >= 1)
#   }
#   ind_to_keep <- union(ind_to_keep,ind_to_keep_temp)
# }
# ind_to_keep <- unique(ind_to_keep)
# cpm.y <- cpm.y[ind_to_keep,]

#y = y[ind_to_keep,]


desing_matrix <- model.matrix(~0+group_cells_selected_cell)
colnames(desing_matrix) <- levels(y$samples$group)

#normalization
y <- calcNormFactors(y)

y <- estimateDisp(y,desing_matrix)

y <- estimateGLMCommonDisp(y,desing_matrix)

y <- estimateGLMTagwiseDisp(y,desing_matrix)

 hc_afterNorm <- hclust(dist(t(cpm.y)),method = "ward.D2")
# pdf('hc_cpm_all_experiments_after_sumTechRep.pdf')
 plot(hc_afterNorm)
# dev.off()
# 
# pdf('edger_mdsplot_RNAseq_after_sumTechRep.pdf')
# #  par(las= 1, cex= 0.8)
# col_map <- group_cells
# levels(col_map) <- 1:length(levels(col_map))
# col_map <- as.numeric(group_cells)
# temp_rainbow <- rainbow(7)
# col_map_rainbow <- temp_rainbow[col_map]
# plotMDS(y,
#         #pch=c(0,3,0,3,0,3),
#         xlim=c(-3.5,3.8), ylim= c(-3.3,2), 
#         #xaxp= c(-40, 40, 8),  yaxp= c(-10,10, 8),
#         #labels= paste(d_esc$samples$group, colnames(d_esc), sep= '-'),
#         col=col_map_rainbow,
#         main= 'MDS plot of control and condition libraries', cex= 0.7)
# grid()
# dev.off()


# QL - Quasi Likelyhood model represenintg the sdtudy design is fitted to the data
# QL account for gene-specifi variability from both biological and technical sources. 
# In this framework, the variace of each gene count is a quadratic 
# 
fit <- glmQLFit(y, desing_matrix)

# comparisons
#qlf.3vs2 <- glmQLFTest(fit, contrast=c(0,-1,1))
my.contrast <- makeContrasts(
  # comparison of hESC
  H9_CNCC_vs_H9_hESC = H9_CNCC - H9_hESC,
  H9_NSC_vs_H9_hESC = H9_NSC - H9_hESC,
  H9_NSC_vs_H9_CNCC = H9_NSC - H9_CNCC,
  levels = desing_matrix
  #levels = c("H9_CNCC","H9_hESC","H9_NSC")
)


# DIFFERENTIAL ANALYSIS OF ALL PAIR-WISE COMPARISON
flag_GO = F # this flag is to generate or not GO enrichment with TOP GO
flag_for_name <- "updated_CNCC3samplesSel_all_others_included"
DE_genes_ALL <- DE_analysis_routine(y,fit,my.contrast,flag_GO,flag_for_name)

# y_DEGL=y
# fit_of_interest=fit,
# contrast_of_interest=my.contrast
# produce_GO=flag_GO
```
