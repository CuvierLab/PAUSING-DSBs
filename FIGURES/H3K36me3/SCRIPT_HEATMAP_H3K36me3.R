
#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

library(GenomicRanges)
library(ggplot2)
library(ggpubr) # magrittr
library(gplots)
'%ni%' = Negate('%in%')
require(BiocGenerics)
require(parallel)
library(gsubfn) # proto
library(Rsamtools)
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library(rtracklayer)
library(nucleR)
library(gridBase)

######## PATH 
workdir="/home/cperrois/work/"
outfig="/home/cperrois/work/PROJET_H2AV_2021/FIGURE/HEATMAPS/HEATMAP_H3K36me3/BY_H3K36/By_Q_H3K36m3/"

source(paste0(workdir, "functionR/Script_HEATMAP_profile.R"))

######## DATA 
GNref = readRDS(paste0(workdir, "PROJET_H2AV_2021/DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))

PAUSE_INDICE_VEC = readRDS(paste0(workdir, "DATA/LIST_FEATURES/PAUSE_INDICE_VEC.RDS"))
PAUSE_INDICE = PAUSE_INDICE_VEC$PAUSE_IND_pol2_ctrl_N

ZSCORE_PROFMAT_GB_H3K36m3_2N4_vs_2C4=readRDS("/home/cperrois/work/PROJET_H2AV_2021/DATA/LIST_FEATURES/ZSCORE/ZSCORE_PROFMAT_GB_H3K36m3_2N4_vs_2C4.RDS")

######## HEATMAP 

rangeheatmap = c(1:1000)

pdf(paste0(outfig, "ZSCORE_PROFMAT_GB_H3K36m3_2N4_vs_2C4","_BY_","PAUSE_INDICE",".pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_H3K36m3_2N4_vs_2C4[names(PAUSE_INDICE[order(PAUSE_INDICE, decreasing=T)]),rangeheatmap],RangeValue = c(-0.9,0.8),winsorize=c(5,95),order = FALSE, main="ZSCORE_PROFMAT_GB_H3K36m3_2N4_vs_2C4",legend.name="PAUSE_INDICE")        
dev.off()






