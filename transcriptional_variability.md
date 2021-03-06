Trancritional variability based on mother-dauther G4 transmition/loss/acquisition
================
Angela Simeone

Promoters have been stratified by mother-to-daughter G4 transmittion and expression levels mother dauther have been used to fit a regression model. The residuals to that model represent the expression variability. For each of the 4 promoter groups we perform a fitting and subsequently we extract the residuals and estimate R2.

``` r
library(tidyverse)
setwd('/Users/simeon01/Documents/Katie/Sequencing_Sept19_bio_multi2_compendium/broad/scenario1_multi2/20200421_annotated_promoters_using_histones_G4')
load("All_anno.RData")


All_anno_prom_sel <- All_anno %>% select(ens_id,ESC_status,CNCC_status,NSC_status)

promoters_annotations <- read.delim('/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS.slop1000.bed',sep ="\t",stringsAsFactors = F, header = F)
colnames(promoters_annotations) <- c("chr","start","end","ens_id","feat","strand")
TSS_1base <- read.delim('/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS.minus1base.bed',stringsAsFactors = F,header = F)
colnames(TSS_1base) <- c('TSS_chr','TSS_start','TSS_end','ens_id','feature','TSS_strand')


All_anno_prom_1base <- left_join(All_anno,TSS_1base, by='ens_id')


## select G4 pos or nneg
# --> use All_anno_prom_1base
## select G4 pos or neg in the xontect of atac only
# --> use All_anno_prom_1base_atac
All_anno_prom_1base_atac <- All_anno_prom_1base %>% filter(ESC_atac>=1 & CNCC_atac>=1)

## sleect G4+K4 pos and neg
All_anno_prom_1base_G4 <- All_anno_prom_1base %>% filter(ESC_G4>=1 & CNCC_G4>=1)
All_anno_prom_1base_K4 <- All_anno_prom_1base %>% filter(ESC_H3K4me3>=1 &CNCC_H3K4me3>=1)

All_anno_prom_1base_K4andG4 <- All_anno_prom_1base %>% filter(ESC_H3K4me3>=1 &CNCC_H3K4me3>=1 & ESC_G4>=1 & CNCC_G4>=1)
All_anno_prom_1base_ATACandG4 <- All_anno_prom_1base %>% filter(ESC_atac>=1 & CNCC_atac>=1 & ESC_G4>=1 & CNCC_G4>=1)
All_anno_prom_1base_K4andATAC <- All_anno_prom_1base %>% filter(ESC_H3K4me3>=1 &CNCC_H3K4me3>=1 & ESC_atac>=1 & CNCC_atac>=1)

All_anno_prom_1base_K4_N <- All_anno_prom_1base %>% filter(ESC_H3K4me3>=1 & NSC_H3K4me3>=1)
All_anno_prom_1base_biv_to_K4 <- All_anno_prom_1base %>% filter(ESC_bival>=1 & CNCC_H3K4me3>=1)
All_anno_prom_1base_biv_to_K4_NSC <- All_anno_prom_1base %>% filter(ESC_bival>=1 &NSC_H3K4me3>=1)

## in any selection elect 2 cells lines; for these cell line what we want to select is the TPM
# select case to screen by feature name
# cell1='ESC'
# cell2='CNCC'
# cell3='NSC'

selection_matrix <- function(dataframe,query_feat,cell1,cells2,cell3){
  # query_feat='G4'
  # dataframe=All_anno_prom_1base
  # dataframe=All_anno_prom_1base
  # query_feat="G4"
  # cell1='ESC'
  # cell2='CNCC'
  # cell3='NSC'
  # 
  temp_anno <- dataframe%>% dplyr::select(-contains(cell3))
  
  # find features to filter data
  feat_to_use1 <- grep(paste0(cell1,'_',query_feat,'$'),colnames(temp_anno), value = T)
  feat_to_use2 <- grep(paste0(cell2,'_',query_feat,'$'),colnames(temp_anno), value = T)
  
  # find features to regress
  param_to_regr1 <- grep(paste0("TPM_",cell1),colnames(temp_anno), value = T)
  param_to_regr2 <- grep(paste0("TPM_",cell2),colnames(temp_anno), value = T)
  
  #select temp_anno, param_to_regr1, param_to_regr2
  temp_anno_final <- temp_anno %>% select(c(feat_to_use1,feat_to_use2,param_to_regr1,param_to_regr2))
  head(temp_anno_final)
  
  #  library(correlationfunnel) --> binarize
  temp_anno_final[feat_to_use1] <- ifelse(temp_anno_final[feat_to_use1]> 0,1,0)
  temp_anno_final[feat_to_use2] <- ifelse(temp_anno_final[feat_to_use2]> 0,1,0)
  head(temp_anno_final)
  
  # #### >>>>> !!!! select only data where TPM product is > 0   !!!!<<<<< #### 
  #temp_anno_final <- temp_anno_final %>% filter((!!sym(param_to_regr1)*!!sym(param_to_regr2))>0)
  #temp_anno_final <- temp_anno_final %>% filter((!!sym(param_to_regr1)+!!sym(param_to_regr2))>1)
  #temp_anno_final <- temp_anno_final %>% filter(!!sym(param_to_regr1)>=1 | !!sym(param_to_regr2)>=1)
  temp_anno_final <- temp_anno_final %>% filter((!!sym(param_to_regr1) + !!sym(param_to_regr2))>0)
  
  
  #create dummy variable
  temp_anno_final$combined_feat <- paste0(temp_anno_final[[feat_to_use1]],"_",temp_anno_final[[feat_to_use2]])
  head(temp_anno_final)
  
  
  # fit regression on all elements
  lm_TPM <- lm(log2(temp_anno_final[[param_to_regr1]]+0.1) ~ log2(temp_anno_final[[param_to_regr2]]+0.1),weights = rank(log2(temp_anno_final[[param_to_regr1]]+0.1)))
  
  plot(log2(temp_anno_final[[param_to_regr1]]+0.1),
       log2(temp_anno_final[[param_to_regr2]]+0.1), col= "grey50")
  abline(lm_TPM, col='orange',lwd=2)
  
  
  #paste redisuals and dummy var
  new_df <- data.frame(res=lm_TPM$residuals,case=temp_anno_final$combined_feat)
  head(new_df)
  
  hist(new_df$res)
  new_df$case <- factor(new_df$case, levels = c("0_0", "1_0","1_1","0_1"))
  gg1 <- ggplot(new_df,aes(x=case, y=res^2, fill=case)) + 
    stat_boxplot(geom = "errorbar", width = 0.2) + 
    geom_boxplot(outlier.shape=NA,notch=T) + coord_cartesian(ylim=c(0, 7.5))
  
  ## second approach fit regression with weigths
  ## other approach: loop over combined_feat, fit lm and then extract res
  res_out <- c()
  case_out <- c()
  rsquare <- c()
  lm_model_fitted=list()
  for (i in 1:length(unique(temp_anno_final$combined_feat))) {
    print(i)
    
    curr_case <- unique(temp_anno_final$combined_feat)[i]
    print(curr_case)
    a <- log2(temp_anno_final[[param_to_regr1]]+0.1)
    b <- log2(temp_anno_final[[param_to_regr2]]+0.1)
    ind_selection <- which(temp_anno_final$combined_feat==curr_case)
    length(ind_selection)
    
    lm_TPM_sel <- lm(b[ind_selection] ~ a[ind_selection],weights = rank(a[ind_selection]))
    #lm_TPM_sel <- lm(b[ind_selection] ~ a[ind_selection],weights =0.9^rank(a[ind_selection]))
    plot(a[ind_selection],b[ind_selection],col= "grey50")
    abline(lm_TPM_sel, col='orange',lwd=2)
    res_out<- c(res_out, lm_TPM_sel$residuals)
    case_out <- c(case_out, temp_anno_final$combined_feat[ind_selection])
    #print(lm_TPM_sel)
    #summary(lm_TPM_sel)
    rsquare <- c(rsquare,summary(lm_TPM_sel)$r.squared)
    lm_model_fitted[[i]] <- lm_TPM_sel
    
  }
  
  length(res_out)
  length(case_out)
  length(lm_model_fitted)
  names(rsquare) <- names(lm_model_fitted) <- unique(temp_anno_final$combined_feat)
  rsquare_df <- data.frame(rsquare_outcome=rsquare,case=unique(temp_anno_final$combined_feat))
  rsquare_df$case <- factor(rsquare_df$case,levels = c("0_0", "1_0","1_1","0_1"))
  
  new_df_second <- data.frame(res=res_out,case=case_out)
  
  new_df_second$case <- factor(new_df_second$case, levels = c("0_0", "1_0","1_1","0_1"))
  head(new_df_second)
  
  
  new_df_second <- new_df_second %>% mutate(color_by_group=case_when(case =="0_0" ~"#adadad",
                                                                     case =="1_0" ~"#1EB2AA",
                                                                     case =="1_1" ~"#956EB0",
                                                                     case =="0_1" ~"#F47E53"))
  
  
  cols = c("#adadad","#1EB2AA","#956EB0","#F47E53")
  #cols = c("#adadad50","#1EB2AA50","#956EB050","#F47E5350")
  gg2 <- ggplot(new_df_second,aes(x=case, y=res^2 )) + 
    stat_boxplot(geom = "errorbar", width = 0.2) + 
    geom_boxplot(outlier.shape=NA,width=0.5,notch=T,fill=cols) +
    #geom_boxplot(outlier.shape=NA,width=0.5)+
    coord_cartesian(ylim=c(0, 12)) +theme_bw() + theme(legend.position = "none")  
  
  gg2
  # select 3 cases to generate individual_boxplots
  new_df_second_select_11_00 <- new_df_second %>% filter(case%in% c("1_1","0_0")) 
  new_df_second_select_11_00$case2 <- factor(new_df_second_select_11_00$case, levels = c("1_1","0_0"),ordered = T)
  new_df_second_select_11_10 <- new_df_second %>% filter(case%in% c("1_1","1_0"))
  new_df_second_select_11_10$case2 <- factor(new_df_second_select_11_10$case, levels = c("1_1","0_0"),ordered = T)
  new_df_second_select_11_01 <- new_df_second %>% filter(case%in% c("1_1","0_1"))
  new_df_second_select_11_01$case2 <- factor(new_df_second_select_11_01$case, levels = c("1_1","0_0"),ordered = T)
  
  #extract ranges with function that reads a dataframe with a pair of conditions
  
  
  max_val_11_00 = 1.1 * new_df_second_select_11_00 %>%group_by(case) %>% summarize(range = quantile(res^2,0.75) + 1.5*quantile(res^2,0.75)-quantile(res^2,0.25)) %>% pull(range) %>% max()
  max_val_11_00
  
  max_val_11_10 = 1.1 * new_df_second_select_11_10 %>%group_by(case) %>% summarize(range = quantile(res^2,0.75) + 1.5*quantile(res^2,0.75)-quantile(res^2,0.25)) %>% pull(range) %>% max()
  max_val_11_10
  
  max_val_11_01 = 1.1 * new_df_second_select_11_01 %>%group_by(case) %>% summarize(range = quantile(res^2,0.75) + 1.5*quantile(res^2,0.75)-quantile(res^2,0.25)) %>% pull(range) %>% max()
  max_val_11_01
  
  max_for_all <- max(c(max_val_11_00,max_val_11_10,max_val_11_01))
  
  
  # "0_0" ~"#adadad",
  # "1_0" ~"#1EB2AA",
  # "1_1" ~"#956EB0",
  # "0_1" ~"#F47E53"))
  
  
  gg2_select_11_00 <- ggplot(new_df_second_select_11_00,aes(x=case2, y=res^2 )) + 
    stat_boxplot(geom = "errorbar", width = 0.2) + 
    geom_boxplot(outlier.shape=NA,width=0.5,notch=T,fill=c("#956EB0","#adadad")) +
    #geom_boxplot(outlier.shape=NA,width=0.5)+
    coord_cartesian(ylim=c(0, max_for_all)) +theme_bw() + theme(legend.position = "none")  
  gg2_select_11_00
  gg2_select_11_10 <- ggplot(new_df_second_select_11_10,aes(x=case2, y=res^2 )) + 
    stat_boxplot(geom = "errorbar", width = 0.2) + 
    geom_boxplot(outlier.shape=NA,width=0.5,notch=T,fill=c("#956EB0","#1EB2AA")) +
    #geom_boxplot(outlier.shape=NA,width=0.5)+
    coord_cartesian(ylim=c(0, max_for_all)) +theme_bw() + theme(legend.position = "none")  
  gg2_select_11_10
  gg2_select_11_01 <- ggplot(new_df_second_select_11_01,aes(x=case2, y=res^2 )) + 
    stat_boxplot(geom = "errorbar", width = 0.2) + 
    geom_boxplot(outlier.shape=NA,width=0.5,notch=T,fill=c("#956EB0","#F47E53")) +
    #geom_boxplot(outlier.shape=NA,width=0.5)+
    coord_cartesian(ylim=c(0, max_for_all)) +theme_bw() + theme(legend.position = "none")  
  gg2_select_11_01
  #plot_info <- ggplot_build(gg2_select_11_00)$`data`[[2]]
  
  
  
  ggplot(new_df_second,aes(x=case, y=res^2 )) + 
    stat_boxplot(geom = "errorbar", width = 0.2) + 
    geom_boxplot(outlier.shape=NA,width=0.5,notch=T,fill=cols) +
    #geom_boxplot(outlier.shape=NA,width=0.5)+
    coord_cartesian(ylim=c(0, 12)) +theme_bw() + theme(legend.position = "none")  
  
  # copute test on variance
  pval_variances <- c()
  for (i in 1:length(unique(new_df_second$case))){
    curr_case=unique(new_df_second$case)[i]
    print(curr_case)
    temp <-  var.test(new_df_second$res[which(new_df_second$case=="1_1")],
                      new_df_second$res[which(new_df_second$case==curr_case)], alternative = 'less')$p.value
    pval_variances <- c(pval_variances,temp)
    rm(temp)
  }
  pval_variances_df <- data.frame(pval_oucome=pval_variances,case=unique(new_df_second$case))
  pval_variances_df$case <- factor(pval_variances_df$case, levels = c("0_0", "1_0","1_1","0_1"))
  
  temp_anno_final$combined_feat <- factor(temp_anno_final$combined_feat,levels = c("0_0", "1_0","1_1","0_1"))
  
  # create a column with a color by group
  temp_anno_final <- temp_anno_final %>% mutate(color_by_group=case_when(combined_feat=="0_0" ~"#adadad50",
                                                                         combined_feat=="1_0" ~"#1EB2AA50",
                                                                         combined_feat=="1_1" ~"#956EB050",
                                                                         combined_feat=="0_1" ~"#F47E5350"))
  
  cols_long <- as.character(temp_anno_final$color_by_group)
  names(cols_long) <- as.character(temp_anno_final$case)
  bp <- ggplot(temp_anno_final, aes(x=log2(temp_anno_final[[param_to_regr1]]+0.1), y=log2(temp_anno_final[[param_to_regr2]]+0.1),
                                    group=combined_feat,color=color_by_group)) + 
    geom_point( size = 0.7, col =cols_long) +geom_smooth(size = 0.55, method="lm", color="#00000070", mapping=(aes(weight =rank(log2(temp_anno_final[[param_to_regr1]]+0.1)))),se=FALSE)
  
  bp
  #geom_abline(intercept = 0, col = "grey")
  
  #levels = c("0_0", "1_0","1_1","0_1")
  MyColour <- c("#adadad", "#1EB2AA", "#956EB0", "#F47E53")
  
  # scatter plots by groups
  gbp <- bp + facet_grid(~combined_feat) + xlab(paste0('log2(TPM) ',cell1)) + ylab(paste0('log2(TPM) ',cell2)) + theme_bw() + theme(legend.position = "none",strip.background = element_rect(colour="white", fill="white"))  
  gbp
  
  # scatter plot of each 3 individual group (00, 10 and 01) vs 11
  
  df_layer_1 <- temp_anno_final[which(temp_anno_final$combined_feat=="0_0"),]
  df_layer_2 <- temp_anno_final[which(temp_anno_final$combined_feat=="1_1"),]
  df_layer_1_and_2 <- temp_anno_final[which(temp_anno_final$combined_feat=="0_0" | temp_anno_final$combined_feat=="1_1"),]
  
  GG_2_groups <- ggplot(temp_anno_final,aes(x=log2(temp_anno_final[[param_to_regr1]]+0.1), y=log2(temp_anno_final[[param_to_regr2]]+0.1)), alpha=0) +
    geom_smooth(method="lm", color="gray64",mapping=(aes(weight =rank(log2(temp_anno_final[[param_to_regr1]]+0.1)))),se=FALSE) +
    geom_point(
      data=df_layer_1,
      aes(x=log2(df_layer_1[[param_to_regr1]]+0.1), y=log2(df_layer_1[[param_to_regr2]]+0.1)),
      colour="grey",
      size=9, alpha = 0.15)+
    geom_point(
      data=df_layer_2,
      aes(x=log2(df_layer_2[[param_to_regr1]]+0.1), y=log2(df_layer_2[[param_to_regr2]]+0.1)),
      colour="#956EB0",
      size=1) + theme_bw() + labs(x=paste0("log2(TPM_",cell1,")"),y=paste0("log2(TPM_",cell2,")"))
  
  # example 
  # ggplot(data = dat, aes(x, y)) +   stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200) +   scale_fill_continuous(low = "white", high = "dodgerblue4")
  
  
  # GG_2_groups <- ggplot(temp_anno_final,aes(x=log2(temp_anno_final[[param_to_regr1]]+0.1), y=log2(temp_anno_final[[param_to_regr2]]+0.1)), alpha=0) +
  #   geom_smooth(method="lm", color="gray64",mapping=(aes(weight =rank(log2(temp_anno_final[[param_to_regr1]]+0.1)))),se=FALSE) +
  #   stat_density2d(
  #     data=df_layer_1,
  #     aes(x=log2(df_layer_1[[param_to_regr1]]+0.1), y=log2(df_layer_1[[param_to_regr2]]+0.1),fill = ..density..^0.15),geom = "tile", contour = FALSE, n = 20,
  #     colour="grey") + scale_fill_continuous(low = "white", high = "grey")  + 
  #     #size=9, alpha = 0.15)+
  #   geom_point(
  #     data=df_layer_2,
  #     aes(x=log2(df_layer_2[[param_to_regr1]]+0.1), y=log2(df_layer_2[[param_to_regr2]]+0.1)),
  #     colour="#956EB0",
  #     size=1) + theme_bw() + labs(x=paste0("log2(TPM_",cell1,")"),y=paste0("log2(TPM_",cell2,")"))
  # GG_2_groups
  
  ggsave(paste0('ScatterTPM_',query_feat,'_',cell1,'_',cell2,'_group00grey_11magenta.pdf'),plot=GG_2_groups, width = 5, height = 5)
  ggsave(paste0('ScatterTPM_',query_feat,'_',cell1,'_',cell2,'_group00grey_11magenta.png'),plot=GG_2_groups, width = 5, height = 5)
  
  
  
  new_df_second_for_scatter <- new_df_second
  x_ele_tot <- matrix(ncol = 1,nrow = nrow(new_df_second))
  for (i in 1:length(unique(new_df_second$case))){
    
    curr_case=unique(new_df_second$case)[i]
    print(curr_case)
    ind_case <- which(new_df_second$case==curr_case)
    l=length(ind_case)
    print(l)
    
    x_ele_tot[ind_case] <- seq(1:l)
    
  }
  
  # new_df_second_for_scatter$x_ele <- x_ele_tot
  # bp_res <- ggplot(new_df_second_for_scatter, aes(y=res,x=x_ele,group=case)) + geom_point(color="grey", size = 1)+geom_hline(yintercept=0, linetype="dashed",color = "grey50", size=0.5)
  # gbp_res <- bp_res + facet_grid(~case) 
  # gbp_res
  
  pG4 <- ggplot(new_df_second, aes(y=res^2,x=case)) + geom_jitter(color="grey", size = 0.6,width = 0.1)+geom_hline(yintercept=0, linetype="dashed",color = "grey50", size=0.1)
  
  
  
  # pdf(paste0('Residuals_',query_feat,'_',cell1,'_',cell2,'_group00grey_11magenta.pdf'),height = 12,width = 12)
  #   xlim_max <- max(length(new_df_second$res[which(new_df_second$case=='0_0')]),length(new_df_second$res[which(new_df_second$case=='1_1')]))
  #   plot(new_df_second$res[which(new_df_second$case=='0_0')]^2, col = "grey", xlim=c(0,xlim_max),pch=19)
  #   points(new_df_second$res[which(new_df_second$case=='1_1')]^2, col = "magenta",pch=20)
  # dev.off()
  
  new_df_second_end <- new_df_second[which(new_df_second$case=='0_0'|new_df_second$case=='1_1'),]
  new_df_second_end$case <- factor(new_df_second_end$case,levels=c('1_1','0_0'))
  GG_2_groups_boxplot <- ggplot(new_df_second_end, aes(y=res^2,x=case)) +
    stat_boxplot(geom = "errorbar", width = 0.2) + 
    #geom_boxplot(outlier.shape=NA,width=0.5,notch=T)
    geom_boxplot(outlier.shape=NA,width=0.5,notch=T,fill=c('#956EB0','gray90')) +
    coord_cartesian(ylim=c(0, 12) ) + theme_bw()
  ggsave(paste0('Boxplot_Res2_',query_feat,'_',cell1,'_',cell2,'_group00grey_11magenta.pdf'),plot=GG_2_groups_boxplot, width = 3, height = 4)
  
  
  #export list of results
  list_results <- list(residual_single_lm=new_df, 
                       residual_individual_lm=new_df_second, 
                       scatterPlot=gbp,
                       boxplot_single=gg1,
                       boxplot_individual=gg2,
                       dotPlot_residuals=pG4,
                       fitted_model_single=lm_TPM,
                       tt=lm_TPM_sel,
                       pval_variances_individual = pval_variances_df,
                       rsquare_lm=rsquare_df,
                       fitted_models=lm_model_fitted,
                       boxplot_select_11_00=gg2_select_11_00,
                       boxplot_select_11_10=gg2_select_11_10,
                       boxplot_select_11_01=gg2_select_11_01)
  
  
  
  
  
  return(list_results)
}



cell1='ESC'
cell2='CNCC'
cell3='NSC'
library(ggpubr)

# G4
caseG4 <- selection_matrix(All_anno_prom_1base,"G4",cell1,cells2,cell3)
G4_C <- ggarrange( caseG4$scatterPlot, caseG4$boxplot_individual, ncol=1) 
G4_C
ggsave(plot=G4_C,filename=paste0("Scatter_boxes_different_groups_G4",cell1,cell2,".nov2020.pdf"), width = 8, height   = 8)
ggsave(plot=G4_C,filename=paste0("Scatter_boxes_different_groups_G4",cell1,cell2,".nov2020.png"), width = 8, height   = 8)
caseG4$pval_variances_individual
caseG4$fitted_models[[1]]
G4_box11_00 <- caseG4$boxplot_select_11_00
G4_box11_10 <- caseG4$boxplot_select_11_10
G4_box11_01 <- caseG4$boxplot_select_11_01
ggsave(plot=G4_box11_00,filename=paste0("boxplot_11_00_G4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)
ggsave(plot=G4_box11_10,filename=paste0("boxplot_11_10_G4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)
ggsave(plot=G4_box11_01,filename=paste0("boxplot_11_01_G4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)


write.table(cbind(caseG4$pval_variances_individual,caseG4$rsquare_lm), file="Pvalue_Rsquare_lm_G4_ESC_CNCC.csv",row.names=T, col.names = NA, quote = F,sep = ",")


# ATAC
caseAtac <- selection_matrix(All_anno_prom_1base,"atac",cell1,cells2,cell3)
atac_C <- ggarrange( caseAtac$scatterPlot, caseAtac$boxplot_individual, ncol=1)
ggsave(plot=atac_C,filename=paste0("Scatter_boxes_different_groups_atac",cell1,cell2,".nov2020.pdf"), width = 8, height   = 8)
ggsave(plot=atac_C,filename=paste0("Scatter_boxes_different_groups_atac",cell1,cell2,".nov2020.png"), width = 8, height   = 8)
caseAtac$pval_variances_individual
atac_C_box11_00 <- caseAtac$boxplot_select_11_00
atac_C_box11_10 <- caseAtac$boxplot_select_11_10
atac_C_box11_01 <- caseAtac$boxplot_select_11_01
ggsave(plot=atac_C_box11_00,filename=paste0("boxplot_11_00_atac",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)
ggsave(plot=atac_C_box11_00,filename=paste0("boxplot_11_10_atac",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)
ggsave(plot=atac_C_box11_00,filename=paste0("boxplot_11_01_atac",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)

#ggplot(caseAtac$pval_variances_individual,aes(x=case,y=-10*log(pval_oucome) )) + geom_bar(stat="identity", width = 0.3) + geom_hline(yintercept=-10*log(0.05), linetyp ="dashed",color = "red", size=0.1)
write.table(cbind(caseAtac$pval_variances_individual,caseAtac$rsquare_lm), file="Pvalue_Rsquare_lm_atac_ESC_CNCC.csv",row.names=T, col.names = NA, quote = F,sep = ",")

# H3K4me3
caseK4 <- selection_matrix(All_anno_prom_1base,"H3K4me3",cell1,cells2,cell3)
H3K4me3_C <-ggarrange( caseK4$scatterPlot, caseK4$boxplot_individual, ncol=1)
ggsave(plot=H3K4me3_C,filename=paste0("Scatter_boxes_different_groups_H3K4me3",cell1,cell2,".nov2020.pdf"), width = 8, height   = 8)
ggsave(plot=H3K4me3_C,filename=paste0("Scatter_boxes_different_groups_H3K4me3",cell1,cell2,".nov2020.png"), width = 8, height   = 8)
caseK4$pval_variances_individual
K4_C_box11_00 <- caseK4$boxplot_select_11_00
K4_C_box11_10 <- caseK4$boxplot_select_11_10
K4_C_box11_01 <- caseK4$boxplot_select_11_01
ggsave(plot=K4_C_box11_00,filename=paste0("boxplot_11_00_k4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)
ggsave(plot=K4_C_box11_10,filename=paste0("boxplot_11_10_k4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)
ggsave(plot=K4_C_box11_01,filename=paste0("boxplot_11_01_k4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)

write.table(cbind(caseK4$pval_variances_individual,caseK4$rsquare_lm), file="Pvalue_Rsquare_lm_H3K4me3_ESC_CNCC.csv",row.names=T, col.names = NA, quote = F,sep = ",")

barplot(c(caseG4$rsquare_lm$rsquare_outcome[caseG4$rsquare_lm$case=='1_1'],
          caseAtac$rsquare_lm$rsquare_outcome[caseAtac$rsquare_lm$case=='1_1'],
          caseK4$rsquare_lm$rsquare_outcome[caseK4$rsquare_lm$case=='1_1']),
        #ylim=c(0.75,0.85),
        names=c("G4","atac","K4"), ylab="R^2")

R2_ESC_CNCC <- data.frame(R2=c(caseG4$rsquare_lm$rsquare_outcome[caseG4$rsquare_lm$case=='1_1'],
                               caseAtac$rsquare_lm$rsquare_outcome[caseAtac$rsquare_lm$case=='1_1'],
                               caseK4$rsquare_lm$rsquare_outcome[caseK4$rsquare_lm$case=='1_1']),case=c("G4","atac","H3K4me3"))
R2_ESC_CNCC$case <- factor(R2_ESC_CNCC$case,levels=c("G4","atac","H3K4me3"))
ggplot(R2_ESC_CNCC,aes(x=case,y=R2)) + geom_bar(stat="identity") +
  coord_cartesian(ylim=c(0.7, 0.85) ) + theme_bw()
ggsave('R2_ESC_CNCC_G4_atac_H3K4me3.pdf', height = 4, width = 3)


### ESC to NSC ########## ========================================================================================================

cell1='ESC'
cell2='NSC'
cell3='CNCC'

# G4
caseG4_N <- selection_matrix(All_anno_prom_1base,"G4",cell1,cells2,cell3)
G4_N <- ggarrange( caseG4_N$scatterPlot, caseG4_N$boxplot_individual, ncol=1)
ggsave(plot=G4_N,filename=paste0("Scatter_boxes_different_groups_G4",cell1,cell2,".nov2020.pdf"), width = 8, height   = 8)
ggsave(plot=G4_N,filename=paste0("Scatter_boxes_different_groups_G4",cell1,cell2,".nov2020.png"), width = 8, height   = 8)
ggplot(caseG4_N$pval_variances_individual,aes(x=case,y=-10*log(pval_oucome) )) + geom_bar(stat="identity", width = 0.3) + geom_hline(yintercept=-10*log(0.05), linetyp ="dashed",color = "red", size=0.1)
G4_N_box11_00 <- G4_N$boxplot_select_11_00
G4_N_box11_10 <- G4_N$boxplot_select_11_10
G4_N_box11_01 <- G4_N$boxplot_select_11_01
ggsave(plot=G4_N_box11_00,filename=paste0("boxplot_11_00_G4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)
ggsave(plot=G4_N_box11_10,filename=paste0("boxplot_11_10_G4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)
ggsave(plot=G4_N_box11_01,filename=paste0("boxplot_11_01_G4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)
write.table(cbind(caseG4_N$pval_variances_individual,caseG4_N$rsquare_lm), file="Pvalue_Rsquare_lm_G4_ESC_NSC.csv",row.names=T, col.names = NA, quote = F,sep = ",")

# ATAC
caseAtac_N <- selection_matrix(All_anno_prom_1base,"atac",cell1,cells2,cell3)
atac_N <- ggarrange( caseAtac_N$scatterPlot, caseAtac_N$boxplot_individual, ncol=1)
ggsave(plot=atac_N,filename=paste0("Scatter_boxes_different_groups_atac",cell1,cell2,".nov2020.pdf"), width = 8, height   = 8)
ggsave(plot=atac_N,filename=paste0("Scatter_boxes_different_groups_atac",cell1,cell2,".nov2020.png"), width = 8, height   = 8)
ggplot(caseAtac_N$pval_variances_individual,aes(x=case,y=-10*log(pval_oucome) )) + geom_bar(stat="identity", width = 0.3) + geom_hline(yintercept=-10*log(0.05), linetyp ="dashed",color = "red", size=0.1)
write.table(cbind(caseAtac_N$pval_variances_individual,caseAtac_N$rsquare_lm), file="Pvalue_Rsquare_lm_atac_ESC_NSC.csv",row.names=T, col.names = NA, quote = F,sep = ",")
atac_N_box11_00 <- atac_N$boxplot_select_11_00
atac_N_box11_10 <- atac_N$boxplot_select_11_10
atac_N_box11_01 <- atac_N$boxplot_select_11_01
ggsave(plot=atac_N_box11_00,filename=paste0("boxplot_11_00_G4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)
ggsave(plot=atac_N_box11_10,filename=paste0("boxplot_11_10_G4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)
ggsave(plot=atac_N_box11_01,filename=paste0("boxplot_11_01_G4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)

# H3K4me3
caseK4_N <- selection_matrix(All_anno_prom_1base,"H3K4me3",cell1,cells2,cell3)
H3K4me3_N <- ggarrange( caseK4_N$scatterPlot, caseK4_N$boxplot_individual, ncol=1)
ggsave(plot=H3K4me3_N,filename=paste0("Scatter_boxes_different_groups_H3K4me3",cell1,cell2,".nov2020.pdf"), width = 8, height   = 8)
ggsave(plot=H3K4me3_N,filename=paste0("Scatter_boxes_different_groups_H3K4me3",cell1,cell2,".nov2020.png"), width = 8, height   = 8)
ggplot(caseK4_N$pval_variances_individual,aes(x=case,y=-10*log(pval_oucome) )) + geom_bar(stat="identity", width = 0.3) + geom_hline(yintercept=-10*log(0.05), linetyp ="dashed",color = "red", size=0.1)
k4_N_box11_00 <- caseK4_N$boxplot_select_11_00
k4_N_box11_10 <- caseK4_N$boxplot_select_11_10
k4_N_box11_01 <- caseK4_N$boxplot_select_11_01
ggsave(plot=k4_N_box11_00,filename=paste0("boxplot_11_00_G4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)
ggsave(plot=k4_N_box11_10,filename=paste0("boxplot_11_10_G4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)
ggsave(plot=k4_N_box11_01,filename=paste0("boxplot_11_01_G4",cell1,cell2,".nov2020.pdf"), width = 5, height   = 5)


write.table(cbind(caseK4_N$pval_variances_individual,caseK4_N$rsquare_lm), file="Pvalue_Rsquare_lm_H3K4me3_ESC_NSC.csv",row.names=T, col.names = NA, quote = F,sep = ",")


R2_ESC_CNCC_N <- data.frame(R2=c(caseG4_N$rsquare_lm$rsquare_outcome[caseG4_N$rsquare_lm$case=='1_1'],
                                 caseAtac_N$rsquare_lm$rsquare_outcome[caseAtac_N$rsquare_lm$case=='1_1'],
                                 caseK4_N$rsquare_lm$rsquare_outcome[caseK4_N$rsquare_lm$case=='1_1']),case=c("G4","atac","H3K4me3"))
R2_ESC_CNCC_N$case <- factor(R2_ESC_CNCC_N$case,levels=c("G4","atac","H3K4me3"))

ggplot(R2_ESC_CNCC_N,aes(x=case,y=R2)) + geom_bar(stat="identity") +
  coord_cartesian(ylim=c(0.7, 0.85) ) + theme_bw()
ggsave('R2_ESC_NSC_G4_atac_H3K4me3.pdf', height = 4, width = 3)
```
