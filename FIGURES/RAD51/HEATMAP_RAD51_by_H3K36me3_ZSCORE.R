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

outfig=paste0(workdir,"FIGURES/HEATMAP/RAD51/")

#####################################################################################-
#          LOAD DATA  ----
#####################################################################################-

GNref = readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f = LIST_QUANTIF_K36$ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f



ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2 = readRDS(paste0(workdir,"DATA/ZSCORE/ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2.RDS"))

#####################################################################################-
#         PLOT  ----
#####################################################################################-

rangeheatmap = c(1:1000)

# ZSCORE NELF-KD vs WT

pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2_BY_ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2[names(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f[order(ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f, decreasing=T)]),rangeheatmap],winsorize=c(5,95), order = FALSE, main="ZSCORE_PROFMAT_GB_Rad51_N_vs_WT_R2",legend.name="ZSCORE_H3K36me3_2C4_H3K36me3_2N4_f")        
dev.off()




