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
outfig="/home/cperrois/work/PROJET_H2AV_2021/FIGURE/HEATMAPS/HEATMAP_H2AV/BY_H3K36/BY_Q_H3K36/"

source(paste0(workdir, "functionR/Script_HEATMAP_profile.R"))

######## DATA 
GNref = readRDS(paste0(workdir, "PROJET_H2AV_2021/DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))


### ZSCORE 

ZSCORE_PROFMAT_GB_PHYPB_R1_PW_R1=readRDS("/home/cperrois/work/PROJET_H2AV_2021/DATA/LIST_FEATURES/ZSCORE/ZSCORE_PROFMAT_GB_PHYPB_R1_PW_R1.RDS")
ZSCORE_PROFMAT_GB_PHYPB_R2_PW_R2=readRDS("/home/cperrois/work/PROJET_H2AV_2021/DATA/LIST_FEATURES/ZSCORE/ZSCORE_PROFMAT_GB_PHYPB_R2_PW_R2.RDS")
ZSCORE_PROFMAT_GB_PHYPB_R1_PW_R2=readRDS("/home/cperrois/work/PROJET_H2AV_2021/DATA/LIST_FEATURES/ZSCORE/ZSCORE_PROFMAT_GB_PHYPB_R1_PW_R2.RDS")
ZSCORE_PROFMAT_GB_PHYPB_R2_PW_R1=readRDS("/home/cperrois/work/PROJET_H2AV_2021/DATA/LIST_FEATURES/ZSCORE/ZSCORE_PROFMAT_GB_PHYPB_R2_PW_R1.RDS")

LIST_QUANTIF_K36=readRDS("/home/cperrois/work/PROJET_H2AV/DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS")
Q_H3K36me3_2C4_GB_f=LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f


rangeheatmap = c(1:1000)
pdf(paste0(outfig,"ZSCORE_PROFMAT_GB_PHYPB_PW_Combin_BY_Q_H3K36me3_2C4_GB_f_upd5K.pdf"))
heatMatrixMat(ZSCORE_PROFMAT_GB_PHYPB_R1_PW_R1[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),main="ZSCORE_PROFMAT_GB_PHYPB_R1_PW_R1",legend.name="Q_H3K36me3_2C4_GB_f")        
heatMatrixMat(ZSCORE_PROFMAT_GB_PHYPB_R2_PW_R2[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),main="ZSCORE_PROFMAT_GB_PHYPB_R2_PW_R2",legend.name="Q_H3K36me3_2C4_GB_f")
heatMatrixMat(ZSCORE_PROFMAT_GB_PHYPB_R1_PW_R2[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),main="ZSCORE_PROFMAT_GB_PHYPB_R1_PW_R2",legend.name="Q_H3K36me3_2C4_GB_f")
heatMatrixMat(ZSCORE_PROFMAT_GB_PHYPB_R2_PW_R1[names(Q_H3K36me3_2C4_GB_f),rangeheatmap],winsorize=c(5,95),main="ZSCORE_PROFMAT_GB_PHYPB_R2_PW_R1",legend.name="Q_H3K36me3_2C4_GB_f")
dev.off()

