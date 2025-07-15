#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-


library(dplyr) # used for recover names(coverages)
library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library(seqplots.tweaked)
library("BSgenome.Dmelanogaster.UCSC.dm6")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(rtracklayer)
library(devtools)
'%ni%' = Negate('%in%')
library("BSgenome")


#####################################################################################-
#         FUNCTIONS  ----
#####################################################################################-
# Load FUNCTION

workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"

source(paste0(workdir, "functionR/AVG_PROFILE.R"))

tmp <- paste0(workdir,"FIGURES/AVG_PROF/TMPgetPlotSetArray")

#####################################################################################-
#         DATA  ----
#####################################################################################-


### LOAD BIGWIG

HypB_KD_RPGC=paste0(workdir,"DATA/CHIPSEQ/HypB_KD_RPGC_norm.bw")
LucKD_RPGC=paste0(workdir,"DATA/CHIPSEQ/LucKD_RPGC_norm.bw")



## GENENS GROUPES 
## GENE GROUPES 
LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
Q_H3K36me3_2C4_GB_f=LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f
ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f = LIST_QUANTIF_K36$ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f

PAUSE_INDICE_VEC = readRDS(paste0(workdir, "DATA/LIST_FEATURES/PAUSE_INDICE_VEC.RDS"))
PAUSE_INDICE = PAUSE_INDICE_VEC$PAUSE_IND_pol2_ctrl_N

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
UP_5_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "top", prct = 5)
DN_5_Q_H3K36me3_2C4_GB_f = getNameList(Q_H3K36me3_2C4_GB_f, topdown = "down", prct = 5)

get_GR_feat = function(refGN, GNlist){
  myovlp=refGN[refGN$name %in% GNlist]
  return(myovlp)
}

GR_UP_5_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,UP_5_Q_H3K36me3_2C4_GB_f)
GR_DN_5_Q_H3K36me3_2C4_GB_f=get_GR_feat(r6_ref_genes,DN_5_Q_H3K36me3_2C4_GB_f)


GR_list_toPlot_H3K36me3 = c("GR_UP_5_Q_H3K36me3_2C4_GB_f", "GR_DN_5_Q_H3K36me3_2C4_GB_f")



#####################################################################################-
#         PLOT  ----
#####################################################################################-

pdf(paste0(workdir,"FIGURES/AVG_PROF/H2AV/AVG_PROF_H2AV_WT_HYPB_KD_2023_BY_H3K36me3.pdf"))
for(GR in GR_list_toPlot_H3K36me3){
  seqPlotSDoutliers_scaleFact(c(LucKD_RPGC,HypB_KD_RPGC),tmp,GR,c(0,10),c(250,250),type="af",bin=10,
                              smooth=TRUE,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme=NA, colvec = c("#285bad", "#eb3434")) 
}
dev.off()
