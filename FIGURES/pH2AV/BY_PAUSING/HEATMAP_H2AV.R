#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

library(GenomicRanges)
library(ggplot2)
library(ggpubr) 
require(BiocGenerics)
library(rtracklayer)
library(gridBase)

#####################################################################################-
#         FUNCTIONS  ----
#####################################################################################-

workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"
source(paste0(workdir,"functionR/Script_HEATMAP_profile.R"))

outfig=paste0(workdir,"FIGURES/HEATMAP/H2AV/")

#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-

LIST_QUANTIF=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF.RDS"))
ZSCORE_PROFMAT_H2AV_NELF_vs_WT_R1 = LIST_QUANTIF$ZSCORE_PROFMAT_H2AV_NELF_vs_WT_R1

PAUSE_INDICE_VEC = readRDS(paste0(workdir, "DATA/LIST_FEATURES/PAUSE_INDICE_VEC.RDS"))
PAUSE_INDICE = PAUSE_INDICE_VEC$PAUSE_IND_pol2_ctrl_N

#####################################################################################-
#         PLOT  ----
#####################################################################################-

rangeheatmap = c(1:1000)

#############  ZSCORE  NELF - WT #############

pdf(paste0(outfig, "ZSCORE_PROFMAT_H2AV_NELF_vs_WT_R1","_BY_","PAUSE_INDICE",".pdf"))
heatMatrixMat(ZSCORE_PROFMAT_H2AV_NELF_vs_WT_R1[names(PAUSE_INDICE[order(PAUSE_INDICE, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_H2AV_NELF_vs_WT_R1",legend.name="PAUSE_INDICE")
dev.off()
