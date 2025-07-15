#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

library(ggpubr)
library(rstatix)
library("gplots")
library(textplot)
library("ggplot2")
library(dplyr)
'%ni%' = Negate('%in%')

#####################################################################################-
#         FUNCTIONS  ----
#####################################################################################-

workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"
source(paste0(workdir,"functionR/Boxplot_wilcoxListFilter_REF.R"))

outfig = paste0(workdir,"FIGURES/BOXPLOT/BOXPLOT_RAD51/")


#####################################################################################-
#         DATA  ----
#####################################################################################-

RAD51_Nelf_KD = readRDS(paste0(workdir,"DATA/QUANTIF/Q_Rad51_GB_N_bis_R1_f.RDS"))
RAD51_WT = readRDS(paste0(workdir,"DATA/QUANTIF/Q_Rad51_GB_WT_bis_R1_f.RDS"))

#####################################################################################-

##### GENES GROUPES 

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))

Q_H3K36me3_2C4_GB_f=LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f

ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f = LIST_QUANTIF_K36$ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f


getNameList = function(Vec, topdown = "top", prct = 10){
  Vec = Vec[order(Vec, decreasing=T)]
  if(topdown %in% "top"){
    GN = names(Vec[Vec > quantile(Vec, (100-prct)/100)])
  }
  if(topdown %in% "down"){
    GN = names(Vec[Vec < quantile(Vec, (prct)/100)])
  }
  if(topdown %in% "mid"){
    tmp1 = names(Vec[Vec < quantile(Vec, (100/2-prct/2)/100)])
    tmp2 = names(Vec[Vec < quantile(Vec, (100/2-prct/2+prct)/100)])
    GN = tmp2[tmp2 %ni% tmp1]
  }
  return(GN)
}


DN5_Q_H3K36me3_2C4 = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "down", prct = 5)
UP5_Q_H3K36me3_2C4 = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 5)


#### LIST GENES GROUPS

LIST_5_UP_DN_Q_H3K36me3_2C4=list(
  DN5_Q_H3K36me3_2C4=DN5_Q_H3K36me3_2C4,
  delta_UP5_Q_H3K36me3_2C4=UP5_Q_H3K36me3_2C4
)



#####################################################################################-
#         PLOT  ----
#####################################################################################-


Boxplot_wilcoxListFilter_REF(quantifWT = RAD51_WT, quantifKD = RAD51_Nelf_KD, cond1 = "RAD51_WT", cond2="RAD51_Nelf_KD", filterGNList = LIST_5_UP_DN_Q_H3K36me3_2C4,
                             effMin =500,  SampleNorm = c(F, "NULL"),YLIM =NULL, bxplt_color = c("#285bad", "#eb3434"), outlierTH = 0.01, logTrans=T,
                             outdir = outfig, readQuantif = "RAD51_WT", Cond = "RAD51_Nelf_KD", select = "LIST_5_UP_DN_Q_H3K36me3_2C4", info = NULL)
