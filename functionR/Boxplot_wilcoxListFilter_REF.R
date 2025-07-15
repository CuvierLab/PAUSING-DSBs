#PROJET H2AV 
# Cuvier'team
#Refka 2021
#askri.93.refka@gmail.com
################################################################################
## Comparison by boxplot of reads quantif upon 2 conditions / 2 features for a list of genes subset with wilcoxon test
##


#https://www.datanovia.com/en/fr/blog/comment-ajouter-des-p-values-sur-un-ggplot-groupe-avec-ggpubr/#pr%C3%A9requis
#https://www.datanovia.com/en/fr/blog/ggpubr-comment-ajouter-des-p-values-generees-ailleurs-a-un-ggplot/#pr%C3%A9paration-des-donn%C3%A9es
#https://www.datanovia.com/en/fr/blog/ggpubr-comment-ajouter-des-p-values-generees-ailleurs-a-un-ggplot/#pr%C3%A9requis
#https://www.datanovia.com/en/fr/lessons/test-de-wilcoxon-dans-r/

print('USAGE : Boxplot_wilcoxListFilter_REF(quantifWT, quantifKD, cond1 = "WT", cond2="KD", filterGNList, effMin = NULL, YLIM = c(0, 20), bxplt_color=NULL,outlierTH = 0.01, logTrans =T, outdir, readQuantif = "readQuantif", Cond = "KD", select = "", info = "")')
#library(ggpubr)
#library(rstatix)
# Required packages
#library(ggpubr)  # https://github.com/kassambara/ggpubr
#library(rstatix)  # https://github.com/kassambara/rstatix
#library("gplots")
#install.packages('textplot')
#library(textplot)
#library("ggplot2")

#> str(filterGNList)
#List of 3
# $ Random                  : chr [1:697] "FBgn0034744.1" "FBgn0262512.1" "FBgn0264075.1" "FBgn0260859.1" ...
# $ UP10_PAUSE_IND_pol2_ctrl: chr [1:697] "FBgn0013277.1" "FBgn0023083.1" "FBgn0027500.1" "FBgn0029676.1" ...
# $ DN10_PAUSE_IND_pol2_ctrl: chr [1:697] "FBgn0039067.1" "FBgn0032487.1" "FBgn0063386.1" "FBgn0262166.1" ...
# 
#> str(quantifWT)
# Named num [1:17453] 14 556 527 531 527 ...
# - attr(*, "names")= chr [1:17453] "FBgn0000003.1" "FBgn0000008.1" "FBgn0000014.1" "FBgn0000015.1" ...
#str(quantifKD)
# Named num [1:17453] 0 512 467 469 516 ...
#- attr(*, "names")= chr [1:17453] "FBgn0000003.1" "FBgn0000008.1" "FBgn0000014.1" "FBgn0000015.1" ...
#
################################################################################

#quantifWT = Q_H2AV_GB_PW_R1
#quantifKD = Q_H2AV_GB_PN_R1 
#cond1 = "Q_H2AV_GB_PW_R1"
# cond2="Q_H2AV_GB_PN_R1"
#filterGNList = LIST_UP_DN_PAUSE
#effMin =500
#test.side = "two.sided"
# SampleNorm = c(F, "NULL")
#YLIM = c(1,20)
#bxplt_color = c("#0000FF", "#FF3300")
#outlierTH = 0.01
#logTrans=T
#outdir = outfig
#readQuantif = "Q_H2AV_GB_PW_R1"
#Cond = "Q_H2AV_GB_PN_R1 "
#select = "LIST_UP_DN_PAUSE"
#info = NULL
################################################################################
Boxplot_wilcoxListFilter_REF = function(quantifWT, quantifKD, cond1 = "WT", cond2="KD", filterGNList, effMin = NULL,  SampleNorm = c(F, "NULL"), YLIM = NULL, bxplt_color=NULL, outlierTH = 0.01, logTrans =T, outdir, readQuantif = "readQuantif", Cond = "KD", select = "filterValue", info = "")
{
  # LOG transform
  if(logTrans %in% T){
        quantifWT = log2(quantifWT+1)
        quantifWT=quantifWT[!is.na(quantifWT)]
        quantifKD = log2(quantifKD+1)
        quantifKD=quantifKD[!is.na(quantifKD)]
  }
  #d <- d[!is.na(d)]
  ## remove OUTLIER
  if(outlierTH %ni% 0){
    for(i in 1:length(filterGNList)){
      QUANT_WT_fortrim = quantifWT[names(quantifWT) %in% filterGNList[[i]]]
      QUANT_KD_fortrim = quantifKD[names(quantifKD) %in% filterGNList[[i]]]
      GN_to_remove = unique(c(names(which(QUANT_WT_fortrim > quantile(QUANT_WT_fortrim, 1-outlierTH))), names(which(QUANT_WT_fortrim < quantile(QUANT_WT_fortrim, outlierTH))),
      names(which(QUANT_KD_fortrim > quantile(QUANT_KD_fortrim, 1-outlierTH))), names(which(QUANT_KD_fortrim < quantile(QUANT_KD_fortrim, outlierTH)))))
      filterGNList[[i]] = filterGNList[[i]][filterGNList[[i]] %ni% GN_to_remove]
    }
  }
  ## Norm inter sample
  if(SampleNorm[1] %in% T){
    GN_for_norm = filterGNList[[SampleNorm[2]]]
    SampleNorm_fact = median(quantifKD[GN_for_norm])/median(quantifWT[GN_for_norm])
    quantifWT = quantifWT*SampleNorm_fact
  }

  ### Normalisation des effectifs
  if(is.null(effMin)){
    effMin =  min(unlist(lapply(filterGNList, length)))
  }
  filterGNList_effMin = lapply(filterGNList, function(filterGN){
    if(length(filterGN) >= effMin){
      filterGN = sample(filterGN, effMin)
    }else{
      filterGN = filterGN
    }
  })
  ## GGplot method -> DATA FRAME
  df_toPlot = data.frame()
  if(class(quantifWT) %in% "numeric"){
    # df_toPlot = as.data.frame(cbind(as.numeric(unname(quantifWT)), rep("all_Genes", length(quantifWT)), rep(cond1, length(quantifWT))))
    # df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD)), rep("all_Genes", length(quantifKD)), rep(cond2, length(quantifKD)))))
    for(i in 1:length(filterGNList_effMin)){
      df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifWT[names(quantifWT) %in% filterGNList_effMin[[i]]])), rep(names(filterGNList_effMin)[i], length(quantifWT[names(quantifWT) %in% filterGNList_effMin[[i]]])), rep(cond1, length(quantifWT[names(quantifWT) %in% filterGNList_effMin[[i]]])))))
      df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD[names(quantifKD) %in% filterGNList_effMin[[i]]])), rep(names(filterGNList_effMin)[i], length(quantifKD[names(quantifKD) %in% filterGNList_effMin[[i]]])), rep(cond2, length(quantifKD[names(quantifKD) %in% filterGNList_effMin[[i]]])))))
    }
  }else{
    # df_toPlot = as.data.frame(cbind(as.numeric(unname(quantifWT)), rep("all_Genes", length(quantifWT)), rep(cond1, length(quantifWT))))
    # df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD)), rep("all_Genes", length(quantifKD)), rep(cond2, length(quantifKD)))))
    for(i in 1:length(filterGNList_effMin)){
      df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifWT[rownames(quantifWT) %in% filterGNList_effMin[[i]],1,drop=F])), rep(names(filterGNList_effMin)[i], length(quantifWT[rownames(quantifWT) %in% filterGNList_effMin[[i]],1,drop=F])), rep(cond1, length(quantifWT[rownames(quantifWT) %in% filterGNList_effMin[[i]],1,drop=F])))))
      df_toPlot = rbind(df_toPlot, as.data.frame(cbind(as.numeric(unname(quantifKD[rownames(quantifKD) %in% filterGNList_effMin[[i]],1,drop=F])), rep(names(filterGNList_effMin)[i], length(quantifKD[rownames(quantifKD) %in% filterGNList_effMin[[i]],1,drop=F])), rep(cond2, length(quantifKD[rownames(quantifKD) %in% filterGNList_effMin[[i]],1,drop=F])))))
    }
  }
  colnames(df_toPlot) = c("readsCount", "select", "condition")
  df_toPlot$readsCount = as.numeric(as.character(df_toPlot$readsCount))

  # bxplt_color
  if(is.null(bxplt_color)){
    bxplt_color =  c("#aaaaaa", "#5e5e5e")
  }

#####################################

 #head(df_toPlot)
#readsCount select       condition
# 9.574139 Random Q_H2AV_GB_PW_R1
# 8.763884 Random Q_H2AV_GB_PW_R1
#10.442389 Random Q_H2AV_GB_PW_R1
# 8.651634 Random Q_H2AV_GB_PW_R1
# 9.517005 Random Q_H2AV_GB_PW_R1
# 9.160054 Random Q_H2AV_GB_PW_R1


pdf(paste0(outdir, "BOXPLOT_",readQuantif ,"_reads_WTvs",Cond, "_by",select ,".pdf"), height=7, width=10)
# vecteur de position des p val dans le plot 
y_position_val= max(quantifWT)+1
y_position_val= trunc(y_position_val)
y_position_vect=rep(y_position_val,length.out= length(filterGNList))

### Paired wilcox_test
stat.test_1 <- df_toPlot %>%
group_by(select) %>%
wilcox_test(readsCount ~ condition,paired = TRUE) %>% # paire test 
adjust_pvalue(method = "bonferroni") %>%
add_significance("p")
stat.test_1
##### NO_ paired wilcox_test
stat.test_2 <- df_toPlot %>%
group_by(select) %>%
wilcox_test(readsCount ~ condition,paired = FALSE) %>% # NO paire test 
adjust_pvalue(method = "bonferroni") %>%
add_significance("p")
stat.test_2

# Créer  boxplot 1 
df_toPlot$condition <- factor(df_toPlot$condition, levels = c(cond1,cond2))
bxp_1 <- ggboxplot(df_toPlot, x = "select", y = "readsCount", color = "black", fill="condition", palette = bxplt_color,title="wilcox paire test")
# Ajoutez des p-values sur les graphiques en box plot
stat.test_1 <- stat.test_1 %>% add_xy_position(x = "select", dodge = 0.8)
# ajouter les crochets + des etoiles 
bx1=bxp_1 + stat_pvalue_manual(stat.test_1,  label = "{p.signif}", tip.length = 0,remove.bracket = FALSE, y.position=y_position_vect,x="select")
#bx1=bxp_1 + stat_pvalue_manual(stat.test_1,  label = "{p}{p.signif}",step.increase=0.001 , y.position=y_position_vect)

# Créer box plot 2 
bxp_2 <- ggboxplot(df_toPlot, x = "select", y = "readsCount", color = "black", fill="condition",  palette = bxplt_color,
title= "wilcox NO paire test")
# Ajoutez des p-values sur les graphiques en box plot
stat.test_2 <- stat.test_2 %>% add_xy_position(x = "select", dodge = 0.8)
# ajouter les crochets
bx2=bxp_2 + stat_pvalue_manual(stat.test_2,  label = "{p.signif}", tip.length = 0,remove.bracket = FALSE, y.position=y_position_vect,x="select")

print(bx1)
textplot(stat.test_1 )
print(bx2)
textplot(stat.test_2)
textplot(lapply(filterGNList, length))
textplot(lapply(filterGNList_effMin, length))


dev.off()
}
