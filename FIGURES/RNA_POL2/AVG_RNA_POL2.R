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
seqlevels(r6_ref_genes) <- gsub("chr", "", seqlevels(r6_ref_genes))


### BIGWIG FILES

pol2_ctrl_N_RPGC = paste0(workdir, "DATA/CHIPSEQ/pol2_ctrl_N_filt_sort_RPGC.bw")
pol2_nelf_N_RPGC = paste0(workdir, "DATA/CHIPSEQ/pol2_nelf_N_filt_sort_RPGC.bw")


#####################################################################################-

### GENE GROUPES 

PAUSE_INDICE_VEC = readRDS(paste0(workdir, "DATA/LIST_FEATURES/PAUSE_INDICE_VEC.RDS"))
DELTA_PAUSE_IND_pol2_nelf_N = PAUSE_INDICE_VEC$DELTA_PAUSE_IND_pol2_nelf_N

PAUSE_IND_pol2_ctrl_N = PAUSE_INDICE_VEC$PAUSE_IND_pol2_ctrl_N

PAUSE_IND_pol2_start_TSS_WT_KD = PAUSE_INDICE_VEC$PAUSE_IND_pol2_start_TSS_WT_KD

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



len = length(PAUSE_IND_pol2_start_TSS_WT_KD)
CTRL_60_65_PAUSE_IND_pol2_start_TSS_WT_KD = names(PAUSE_IND_pol2_start_TSS_WT_KD[(len*0.6) : (len*0.65)])
UP_25_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "top", prct = 25)
UP_5_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "top", prct = 5)
UP_1_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "top", prct = 1)
DN_25_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "down", prct = 25)
DN_5_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "down", prct = 5)
DN_1_PAUSE_IND_pol2_start_TSS_WT_KD = getNameList(PAUSE_IND_pol2_start_TSS_WT_KD, topdown = "down", prct = 1)

## GET GRANGES COORDS OF FEATURES TO PLOT AROUND
get_GR_feat = function(refGN, GNlist){
  myovlp=refGN[refGN$name %in% GNlist]
  return(myovlp)
}

GR_UP_25_PAUSE_IND_pol2_start_TSS_WT_KD = get_GR_feat(r6_ref_genes,UP_25_PAUSE_IND_pol2_start_TSS_WT_KD)
GR_UP_5_PAUSE_IND_pol2_start_TSS_WT_KD = get_GR_feat(r6_ref_genes,UP_5_PAUSE_IND_pol2_start_TSS_WT_KD)
GR_UP_1_PAUSE_IND_pol2_start_TSS_WT_KD= get_GR_feat(r6_ref_genes,UP_1_PAUSE_IND_pol2_start_TSS_WT_KD)
GR_CTRL_60_65_PAUSE_IND_pol2_start_TSS_WT_KD= get_GR_feat(r6_ref_genes,CTRL_60_65_PAUSE_IND_pol2_start_TSS_WT_KD)
GR_DN_25_PAUSE_IND_pol2_start_TSS_WT_KD = get_GR_feat(r6_ref_genes,DN_25_PAUSE_IND_pol2_start_TSS_WT_KD)
GR_DN_5_PAUSE_IND_pol2_start_TSS_WT_KD = get_GR_feat(r6_ref_genes,DN_5_PAUSE_IND_pol2_start_TSS_WT_KD)
GR_DN_1_PAUSE_IND_pol2_start_TSS_WT_KD = get_GR_feat(r6_ref_genes,DN_1_PAUSE_IND_pol2_start_TSS_WT_KD)
GR_list_toPlot = c("GR_UP_5_PAUSE_IND_pol2_start_TSS_WT_KD", "GR_UP_25_PAUSE_IND_pol2_start_TSS_WT_KD","GR_UP_1_PAUSE_IND_pol2_start_TSS_WT_KD", "GR_CTRL_60_65_PAUSE_IND_pol2_start_TSS_WT_KD", "GR_DN_5_PAUSE_IND_pol2_start_TSS_WT_KD", "GR_DN_25_PAUSE_IND_pol2_start_TSS_WT_KD","GR_DN_1_PAUSE_IND_pol2_start_TSS_WT_KD")


#####################################################################################-
#         AVERAGE PLOT  ----
#####################################################################################-

# POL2 PAUSE INDICE start TSS(-200 / +300) WT-KD  

pdf(paste0(workdir,"FIGURES/AVG_PROF/RNA_POL2/AVG_RNA_POL2_PAUSE_IND_pol2_start_TSS_WT_KD.pdf"))
for(GR in GR_list_toPlot){
  seqPlotSDoutliers_scaleFact(c(pol2_ctrl_N_RPGC, pol2_nelf_N_RPGC),tmp,GR,c(0,4),c(500,500),type="af",bin=10,
                              smooth=TRUE,spar=0.3, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#594d32", "#e09d00")) 
}
dev.off()




