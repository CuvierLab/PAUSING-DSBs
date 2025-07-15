#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(seqplots)

#####################################################################################-
#         FUNCTIONS  ----
#####################################################################################-

workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"

source(paste0(workdir, "functionR/AVG_PROFILE.R"))

tmp <- paste0(workdir,"FIGURES/AVG_PROF/TMPgetPlotSetArray")

#####################################################################################-
#         DATA  ----
#####################################################################################-

### GENOME REF 

r6_ref_genes = readRDS(paste0(workdir,"DATA/r6.13/TxDb.GR.dm6.RDS"))



#### CHIPSEQ

H3K36me3_FSEA4_Luc = paste0(workdir,"DATA/H3K36me3/ChIP_H3K36me3_FSEA4_Luc.bigWig")
H3K36me3_FSEA4_Nelf = paste0(workdir,"DATA/H3K36me3/ChIP_H3K36me3_FSEA4_Nelf.bigWig")

H3K36me3_FWT3_Luc = paste0(workdir,"DATA/H3K36me3/ChIP_H3K36me3_FWT3_Luc.bigWig")
H3K36me3_FWT3_Nelf = paste0(workdir,"DATA/H3K36me3/ChIP_H3K36me3_FWT3_Nelf.bigWig")



#####################################################################################-

### GENE GROUPES 

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
Q_H3K36me3_2C4_GB_f = LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f

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

len = length(Q_H3K36me3_2C4_GB_f)
CTRL_60_65_Q_H3K36me3_2C4_GB_f = names(Q_H3K36me3_2C4_GB_f[(len*0.6) : (len*0.65)])
UP_5_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 5)
UP_25_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 25)
DN_5_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "down", prct = 5)
DN_25_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "down", prct = 25)



## GET GRANGES COORDS OF FEATURES TO PLOT AROUND
get_GR_feat = function(refGN, GNlist){
  myovlp=refGN[refGN$name %in% GNlist]
  return(myovlp)
}


GR_UP_5_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,UP_5_Q_H3K36me3_2C4_GB_f)
GR_UP_25_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,UP_25_Q_H3K36me3_2C4_GB_f)
GR_CTRL_60_65_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,CTRL_60_65_Q_H3K36me3_2C4_GB_f)
GR_DN_5_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,DN_5_Q_H3K36me3_2C4_GB_f)
GR_DN_25_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,DN_25_Q_H3K36me3_2C4_GB_f)


GR_list_toPlot_H3K36me3 = c("GR_UP_5_Q_H3K36me3_2C4_GB_f", "GR_UP_25_Q_H3K36me3_2C4_GB_f", "GR_CTRL_60_65_Q_H3K36me3_2C4_GB_f", "GR_DN_5_Q_H3K36me3_2C4_GB_f","GR_DN_25_Q_H3K36me3_2C4_GB_f")

#####################################################################################-
#         PLOT  ----
#####################################################################################-

# H3K36me3 FSEA4 - AVERAGE PLOT


pdf(paste0(workdir,"FIGURES/AVG_PROF/H3K36me3_SEA4_WT3/AVG_PROF_H3K36me3_FSEA4_H3K36me3_2C4.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H3K36me3_FSEA4_Luc, H3K36me3_FSEA4_Nelf),tmp,GR,c(0,10000),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.3, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()


# H3K36me3 WT3 - AVERAGE PLOT

pdf(paste0(workdir,"FIGURES/AVG_PROF/H3K36me3_SEA4_WT3/AVG_PROF_H3K36me3_FWT3_H3K36me3_2C4.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(H3K36me3_FWT3_Luc, H3K36me3_FWT3_Nelf),tmp,GR,c(0,10000),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.3, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()




