######################################################################################################################################################
######################################################################################################################################################
# Cuvier's Lab
# David Depierre - Thomas De Freitas - Refka Askri
######################################################################################################################################################
######################################################################################################################################################

#########################################################################################################################
# LIST DONE IN THIS SCRIPT
#########################################################################################################################
GENES_dm6_ORI = readRDS(paste0(workdir, "PROJET_H2AV/DATA/LIST_FEATURES/LIST_GN_ORI_dm6.RDS"))
GNref = readRDS(paste0(workdir, "PROJET_H2AV/DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))
PAUSE_INDICE_VEC = readRDS(paste0(workdir, "PROJET_H2AV/DATA/LIST_FEATURES/PAUSE_INDICE_VEC.RDS"))
LIST_QUANTIF_K36 = readRDS(paste0(workdir, "PROJET_H2AV/DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
LIST_QUANTIF = readRDS(paste0(workdir, "PROJET_H2AV/DATA/LIST_FEATURES/LIST_QUANTIF.RDS"))
LIST_H2AV_VEC = readRDS(paste0(workdir, "PROJET_H2AV/DATA/LIST_FEATURES/LIST_H2AV_VEC.RDS"))


#########################################################################################################################
# UTILS FUNCTIONS
#########################################################################################################################

matMeans <- function(X,Y){ mean(c(X,Y)) }

matReplaceNA = function(M){
  Mv = c(M)
  Mv[which(is.na(Mv))] = 0
  Mv = matrix(Mv, ncol=1)
  rownames(Mv) = rownames(M)
  return(Mv)
}

matReplaceINF = function(M){
  Mv = c(M)
  Mv[is.finite(Mv) %in% FALSE] = 0
  Mv = matrix(Mv, ncol=1)
  rownames(Mv) = rownames(M)
  return(Mv)
}

computeZscore = function(Q_KD, Q_CTRL){
  ZSCORE = (Q_KD-Q_CTRL[names(Q_KD)])/sqrt(matrix(mapply(matMeans, Q_KD, Q_CTRL), ncol=1))
  ZSCORE = matReplaceNA(ZSCORE)
  ZSCORE = matReplaceINF(ZSCORE)
  rownames(ZSCORE) = names(Q_KD)
  ZSCORE  = ZSCORE[,1]
  ZSCORE = ZSCORE[order(ZSCORE, decreasing=T),drop=F]
  return(ZSCORE)
}


#########################################################################################################################
###########################################   CHIP SEQ H3K36me3     #####################################################
# CHIP SEQ K36 QUANTIF ENTIRE GB SCALED
# LIST_QUANTIF_K36

Quantifdir = paste0(workdir, "PROJET_H2AV/DATA/QUANTIF/GB/")

# LOAD DATA

Q_H3K36me2_MNASE_RPGC_GB = readRDS(paste0(Quantifdir, "H3K36me2_MNASE_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H3K36me2_RPGC_GB = readRDS(paste0(Quantifdir, "H3K36me2_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H3K36me3_MNASE_RPGC_GB = readRDS(paste0(Quantifdir, "H3K36me3_MNASE_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H3K36me3_RPGC_GB = readRDS(paste0(Quantifdir, "H3K36me3_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H3K36me3_GB = readRDS(paste0(Quantifdir, "H3K36me3_readsCounts_GB_SCALED.RDS"))

Q_K36C_RPGC_GB = readRDS(paste0(Quantifdir, "K36C_RPGC_readsCounts_GB_SCALED.RDS"))
Q_K36N_RPGC_GB = readRDS(paste0(Quantifdir, "K36N_RPGC_readsCounts_GB_SCALED.RDS"))

Q_H3K36me2_2C3_GB = readRDS(paste0(Quantifdir, "H3K36me2_2C3_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H3K36me3_2C4_GB = readRDS(paste0(Quantifdir, "H3K36me3_2C4_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H3K36me2_2N3_GB = readRDS(paste0(Quantifdir, "H3K36me2_2N3_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H3K36me3_2N4_GB = readRDS(paste0(Quantifdir, "H3K36me3_2N4_RPGC_readsCounts_GB_SCALED.RDS"))


################################ LOGFC
## RATIO K36me2/3 3/2
LOG2RATIO_K36me3_me2 = log2(Q_H3K36me3_RPGC_GB+1)-log2(Q_H3K36me2_RPGC_GB[names(Q_H3K36me3_RPGC_GB)]+1)
LOG2RATIO_K36me2_me3 = log2(Q_H3K36me2_RPGC_GB+1)-log2(Q_H3K36me3_RPGC_GB[names(Q_H3K36me2_RPGC_GB)]+1)
LOG2RATIO_K36me3_2C4_me2_2C3 = log2(Q_H3K36me3_2C4_GB+1)-log2(Q_H3K36me2_2C3_GB[names(Q_H3K36me3_2C4_GB)]+1)
LOG2RATIO_K36me3_2N4_me2_2N3 = log2(Q_H3K36me3_2N4_GB+1)-log2(Q_H3K36me2_2N3_GB[names(Q_H3K36me3_2N4_GB)]+1)
LOG2RATIO_K36me2_2C2_me3_2C4 = log2(Q_H3K36me2_2C3_GB+1)-log2(Q_H3K36me3_2C4_GB[names(Q_H3K36me2_2C3_GB)]+1)
LOG2RATIO_K36me3_me2_MNASE = log2(Q_H3K36me3_MNASE_RPGC_GB+1)-log2(Q_H3K36me2_MNASE_RPGC_GB[names(Q_H3K36me3_MNASE_RPGC_GB)]+1)
LOG2RATIO_K36me2_me3_MNASE = log2(Q_H3K36me2_MNASE_RPGC_GB+1)-log2(Q_H3K36me3_MNASE_RPGC_GB[names(Q_H3K36me2_MNASE_RPGC_GB)]+1)

DELTA_LOG2RATIO_K36me3me2_N_WT = LOG2RATIO_K36me3_2N4_me2_2N3-LOG2RATIO_K36me3_2C4_me2_2C3


# FILTER ACTIF
LOG2RATIO_K36me3_me2_f = LOG2RATIO_K36me3_me2[names(LOG2RATIO_K36me3_me2) %in% GNref]
LOG2RATIO_K36me2_me3_f = LOG2RATIO_K36me2_me3[names(LOG2RATIO_K36me2_me3) %in% GNref]
LOG2RATIO_K36me3_2C4_me2_2C3_f = LOG2RATIO_K36me3_2C4_me2_2C3[names(LOG2RATIO_K36me3_2C4_me2_2C3) %in% GNref]
LOG2RATIO_K36me2_2C2_me3_2C4_f = LOG2RATIO_K36me2_2C2_me3_2C4[names(LOG2RATIO_K36me2_2C2_me3_2C4) %in% GNref]
LOG2RATIO_K36me3_me2_MNASE_f = LOG2RATIO_K36me3_me2_MNASE[names(LOG2RATIO_K36me3_me2_MNASE) %in% GNref]
LOG2RATIO_K36me2_me3_MNASE_f = LOG2RATIO_K36me2_me3_MNASE[names(LOG2RATIO_K36me2_me3_MNASE) %in% GNref]
LOG2RATIO_K36me3_2N4_me2_2N3_f = LOG2RATIO_K36me3_2N4_me2_2N3[names(LOG2RATIO_K36me3_2N4_me2_2N3) %in% GNref]
DELTA_LOG2RATIO_K36me3me2_N_WT_f = DELTA_LOG2RATIO_K36me3me2_N_WT[names(DELTA_LOG2RATIO_K36me3me2_N_WT) %in% GNref]


################################ ZSCORE
## ZSCORE K36me2/3 3/2
ZSCORE_K36me3_me2 = computeZscore(Q_H3K36me3_RPGC_GB, Q_H3K36me2_RPGC_GB)
ZSCORE_K36me2_me3 = computeZscore(Q_H3K36me2_RPGC_GB, Q_H3K36me3_RPGC_GB)
ZSCORE_K36me3_me2_MNASE = computeZscore(Q_H3K36me3_MNASE_RPGC_GB, Q_H3K36me2_MNASE_RPGC_GB)
ZSCORE_K36me2_me3_MNASE = computeZscore(Q_H3K36me2_MNASE_RPGC_GB, Q_H3K36me3_MNASE_RPGC_GB)
ZSCORE_K36N_K36C = computeZscore(Q_K36N_RPGC_GB, Q_K36C_RPGC_GB)
ZSCORE_H3K36me3_2N4_H3K36me3_2C4 = computeZscore(Q_H3K36me3_2N4_GB, Q_H3K36me3_2C4_GB)
ZSCORE_H3K36me2_2N3_H3K36me2_2C3 = computeZscore(Q_H3K36me2_2N3_GB, Q_H3K36me2_2C3_GB)
ZSCORE_H2AV_R3_HYPB_WT = computeZscore(Q_H2AV_hypb_R3,Q_H2AV_ctrl_R3)
ZSCORE_H3K36me3_WT3_GB_N_WT = computeZscore(Q_H3K36me3_FWT3_Nelf,Q_H3K36me3_FWT3_Luc)
ZSCORE_H3K36me3_SEA4_GB_N_WT = computeZscore(Q_H3K36me3_FSEA4_Nelf,Q_H3K36me3_FSEA4_Luc)

# FILTER ACTIF
ZSCORE_K36me3_me2_f = ZSCORE_K36me3_me2[names(ZSCORE_K36me3_me2) %in% GNref]
ZSCORE_K36me2_me3_f = ZSCORE_K36me2_me3[names(ZSCORE_K36me2_me3) %in% GNref]
ZSCORE_K36me3_me2_MNASE_f = ZSCORE_K36me3_me2_MNASE[names(ZSCORE_K36me3_me2_MNASE) %in% GNref]
ZSCORE_K36me2_me3_MNASE_f = ZSCORE_K36me2_me3_MNASE[names(ZSCORE_K36me2_me3_MNASE) %in% GNref]
ZSCORE_K36N_K36C_f = ZSCORE_K36N_K36C[names(ZSCORE_K36N_K36C) %in% GNref]
ZSCORE_H3K36me3_2N4_H3K36me3_2C4_f = ZSCORE_H3K36me3_2N4_H3K36me3_2C4[names(ZSCORE_H3K36me3_2N4_H3K36me3_2C4) %in% GNref]
ZSCORE_H3K36me2_2N3_H3K36me2_2C3_f = ZSCORE_H3K36me2_2N3_H3K36me2_2C3[names(ZSCORE_H3K36me2_2N3_H3K36me2_2C3) %in% GNref]
ZSCORE_H3K36me3_WT3_GB_N_WT = ZSCORE_H3K36me3_WT3_GB_N_WT[names(ZSCORE_H3K36me3_WT3_GB_N_WT) %in% GNref]
ZSCORE_H3K36me3_SEA4_GB_N_WT = ZSCORE_H3K36me3_SEA4_GB_N_WT[names(ZSCORE_H3K36me3_SEA4_GB_N_WT) %in% GNref]

Q_H3K36me2_MNASE_RPGC_GB_f = Q_H3K36me2_MNASE_RPGC_GB[names(Q_H3K36me2_MNASE_RPGC_GB) %in% GNref]
Q_H3K36me2_RPGC_GB_f = Q_H3K36me2_RPGC_GB[names(Q_H3K36me2_RPGC_GB) %in% GNref]
Q_H3K36me3_MNASE_RPGC_GB_f = Q_H3K36me3_MNASE_RPGC_GB[names(Q_H3K36me3_MNASE_RPGC_GB) %in% GNref]
Q_H3K36me3_RPGC_GB_f = Q_H3K36me3_RPGC_GB[names(Q_H3K36me3_RPGC_GB) %in% GNref]
Q_K36C_RPGC_GB_f = Q_K36C_RPGC_GB[names(Q_K36C_RPGC_GB) %in% GNref]
Q_K36N_RPGC_GB_f = Q_K36N_RPGC_GB[names(Q_K36N_RPGC_GB) %in% GNref]
Q_H3K36me2_2C3_GB_f = Q_H3K36me2_2C3_GB[names(Q_H3K36me2_2C3_GB) %in% GNref]
Q_H3K36me3_2C4_GB_f = Q_H3K36me3_2C4_GB[names(Q_H3K36me3_2C4_GB) %in% GNref]
Q_H3K36me2_2N3_GB_f = Q_H3K36me2_2N3_GB[names(Q_H3K36me2_2N3_GB) %in% GNref]
Q_H3K36me3_2N4_GB_f = Q_H3K36me3_2N4_GB[names(Q_H3K36me3_2N4_GB) %in% GNref]


# ORDER decreasingly

LOG2RATIO_K36me3_me2_f = LOG2RATIO_K36me3_me2_f[order(LOG2RATIO_K36me3_me2_f, decreasing=T)]
LOG2RATIO_K36me2_me3_f = LOG2RATIO_K36me2_me3_f[order(LOG2RATIO_K36me2_me3_f, decreasing=T)]
LOG2RATIO_K36me3_2C4_me2_2C3_f = LOG2RATIO_K36me3_2C4_me2_2C3_f[order(LOG2RATIO_K36me3_2C4_me2_2C3_f, decreasing=T)]
LOG2RATIO_K36me2_2C2_me3_2C4_f = LOG2RATIO_K36me2_2C2_me3_2C4_f[order(LOG2RATIO_K36me2_2C2_me3_2C4_f, decreasing=T)]
LOG2RATIO_K36me3_2N4_me2_2N3_f = LOG2RATIO_K36me3_2N4_me2_2N3_f[order(LOG2RATIO_K36me3_2N4_me2_2N3_f, decreasing=T)]
DELTA_LOG2RATIO_K36me3me2_N_WT_f = DELTA_LOG2RATIO_K36me3me2_N_WT_f[order(DELTA_LOG2RATIO_K36me3me2_N_WT_f, decreasing=T)]
ZSCORE_K36N_K36C_f = ZSCORE_K36N_K36C_f[order(ZSCORE_K36N_K36C_f, decreasing=T)]
Q_H3K36me2_RPGC_GB_f = Q_H3K36me2_RPGC_GB_f[order(Q_H3K36me2_RPGC_GB_f, decreasing=T)]
Q_H3K36me3_RPGC_GB_f = Q_H3K36me3_RPGC_GB_f[order(Q_H3K36me3_RPGC_GB_f, decreasing=T)]
Q_K36C_RPGC_GB_f = Q_K36C_RPGC_GB_f[order(Q_K36C_RPGC_GB_f, decreasing=T)]
Q_K36N_RPGC_GB_f = Q_K36N_RPGC_GB_f[order(Q_K36N_RPGC_GB_f, decreasing=T)]
Q_H3K36me2_2C3_GB_f = Q_H3K36me2_2C3_GB_f[order(Q_H3K36me2_2C3_GB_f, decreasing=T)]
Q_H3K36me3_2C4_GB_f = Q_H3K36me3_2C4_GB_f[order(Q_H3K36me3_2C4_GB_f, decreasing=T)]
Q_H3K36me2_2N3_GB_f = Q_H3K36me2_2N3_GB_f[order(Q_H3K36me2_2N3_GB_f, decreasing=T)]
Q_H3K36me3_2N4_GB_f = Q_H3K36me3_2N4_GB_f[order(Q_H3K36me3_2N4_GB_f, decreasing=T)]


# CREATE LIST

LIST_QUANTIF_K36 = list(
  LOG2RATIO_K36me3_me2_f = LOG2RATIO_K36me3_me2_f,
  LOG2RATIO_K36me2_me3_f = LOG2RATIO_K36me2_me3_f,
  LOG2RATIO_K36me3_2C4_me2_2C3_f = LOG2RATIO_K36me3_2C4_me2_2C3_f,
  LOG2RATIO_K36me2_2C2_me3_2C4_f = LOG2RATIO_K36me2_2C2_me3_2C4_f,
  LOG2RATIO_K36me3_2N4_me2_2N3_f = LOG2RATIO_K36me3_2N4_me2_2N3_f,
  DELTA_LOG2RATIO_K36me3me2_N_WT_f = DELTA_LOG2RATIO_K36me3me2_N_WT_f,
  ZSCORE_K36N_K36C_f = ZSCORE_K36N_K36C_f,
  ZSCORE_H3K36me3_2N4_H3K36me3_2C4_f = ZSCORE_H3K36me3_2N4_H3K36me3_2C4_f,
  ZSCORE_H3K36me2_2N3_H3K36me2_2C3_f = ZSCORE_H3K36me2_2N3_H3K36me2_2C3_f,
  Q_H3K36me2_RPGC_GB_f = Q_H3K36me2_RPGC_GB_f,
  Q_H3K36me3_RPGC_GB_f = Q_H3K36me3_RPGC_GB_f,
  Q_K36C_RPGC_GB_f = Q_K36C_RPGC_GB_f,
  Q_K36N_RPGC_GB_f = Q_K36N_RPGC_GB_f,
  Q_H3K36me2_2C3_GB_f = Q_H3K36me2_2C3_GB_f,
  Q_H3K36me3_2C4_GB_f = Q_H3K36me3_2C4_GB_f,
  Q_H3K36me2_2N3_GB_f = Q_H3K36me2_2N3_GB_f,
  Q_H3K36me3_2N4_GB_f = Q_H3K36me3_2N4_GB_f,
  ZSCORE_H3K36me3_WT3_GB_N_WT = ZSCORE_H3K36me3_WT3_GB_N_WT,
  ZSCORE_H3K36me3_SEA4_GB_N_WT = ZSCORE_H3K36me3_SEA4_GB_N_WT)


saveRDS(LIST_QUANTIF_K36, paste0(workdir, "PROJET_H2AV/DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
LIST_QUANTIF_K36 = readRDS(paste0(workdir, "PROJET_H2AV/DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))



#########################################################################################################################
######################################  RNA POL2 PAUSE QUANTIFICATION    ################################################
# COMPUTE PAUSE INDICE with POL2 chipseq quantification
# PAUSE INDICE = reads around TSS / reads on gene body

# LOAD DATA

pol2_ctrl_N_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "pol2_ctrl_N_RPGC_readsCounts_TSS_480_530.RDS"))
pol2_ctrl_N_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "pol2_ctrl_N_RPGC_readsCounts_TSS_550_600.RDS"))
pol2_ctrl_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "pol2_ctrl_RPGC_readsCounts_TSS_480_530.RDS"))
pol2_ctrl_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "pol2_ctrl_RPGC_readsCounts_TSS_550_600.RDS"))
pol2_hypbKD_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "pol2_hypbKD_RPGC_readsCounts_TSS_480_530.RDS"))
pol2_hypbKD_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "pol2_hypbKD_RPGC_readsCounts_TSS_550_600.RDS"))
pol2_kdm4aKD_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "pol2_kdm4aKD_RPGC_readsCounts_TSS_480_530.RDS"))
pol2_kdm4aKD_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "pol2_kdm4aKD_RPGC_readsCounts_TSS_550_600.RDS"))
pol2_mes4KD_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "pol2_mes4KD_RPGC_readsCounts_TSS_480_530.RDS"))
pol2_mes4KD_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "pol2_mes4KD_RPGC_readsCounts_TSS_550_600.RDS"))
pol2_nelf_N_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "pol2_nelf_N_RPGC_readsCounts_TSS_480_530.RDS"))
pol2_nelf_N_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "pol2_nelf_N_RPGC_readsCounts_TSS_550_600.RDS"))
pol2ser2P_ctrl_N_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "pol2ser2P_ctrl_N_RPGC_readsCounts_TSS_480_530.RDS"))
pol2ser2P_ctrl_N_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "pol2ser2P_ctrl_N_RPGC_readsCounts_TSS_550_600.RDS"))
pol2ser2P_ctrl_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "pol2ser2P_ctrl_RPGC_readsCounts_TSS_480_530.RDS"))
pol2ser2P_ctrl_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "pol2ser2P_ctrl_RPGC_readsCounts_TSS_550_600.RDS"))
pol2ser2P_hypbKD_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "pol2ser2P_hypbKD_RPGC_readsCounts_TSS_480_530.RDS"))
pol2ser2P_hypbKD_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "pol2ser2P_hypbKD_RPGC_readsCounts_TSS_550_600.RDS"))
pol2ser2P_mes4KD_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "pol2ser2P_mes4KD_RPGC_readsCounts_TSS_480_530.RDS"))
pol2ser2P_mes4KD_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "pol2ser2P_mes4KD_RPGC_readsCounts_TSS_550_600.RDS"))
pol2ser2P_nelf_N_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "pol2ser2P_nelf_N_RPGC_readsCounts_TSS_480_530.RDS"))
pol2ser2P_nelf_N_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "pol2ser2P_nelf_N_RPGC_readsCounts_TSS_550_600.RDS"))
polIIC_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "polIIC_RPGC_readsCounts_TSS_480_530.RDS"))
polIIC_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "polIIC_RPGC_readsCounts_TSS_550_600.RDS"))
polIIH_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "polIIH_RPGC_readsCounts_TSS_480_530.RDS"))
polIIH_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "polIIH_RPGC_readsCounts_TSS_550_600.RDS"))
RNAPolII_control_rep1_RPGC_readsCounts_TSS_480_530 = readRDS(paste0(Quantifdir, "RNAPolII_control_rep1_RPGC_readsCounts_TSS_480_530.RDS"))
RNAPolII_control_rep1_RPGC_readsCounts_TSS_550_600 = readRDS(paste0(Quantifdir, "RNAPolII_control_rep1_RPGC_readsCounts_TSS_550_600.RDS"))

# FILTER ACTIVE GENES
pol2_ctrl_N_RPGC_readsCounts_TSS_480_530_f = pol2_ctrl_N_RPGC_readsCounts_TSS_480_530[names(pol2_ctrl_N_RPGC_readsCounts_TSS_480_530) %in% GNref]
pol2_ctrl_N_RPGC_readsCounts_TSS_550_600_f = pol2_ctrl_N_RPGC_readsCounts_TSS_550_600[names(pol2_ctrl_N_RPGC_readsCounts_TSS_550_600) %in% GNref]
pol2_ctrl_RPGC_readsCounts_TSS_480_530_f = pol2_ctrl_RPGC_readsCounts_TSS_480_530[names(pol2_ctrl_RPGC_readsCounts_TSS_480_530) %in% GNref]
pol2_ctrl_RPGC_readsCounts_TSS_550_600_f = pol2_ctrl_RPGC_readsCounts_TSS_550_600[names(pol2_ctrl_RPGC_readsCounts_TSS_550_600) %in% GNref]
pol2_hypbKD_RPGC_readsCounts_TSS_480_530_f = pol2_hypbKD_RPGC_readsCounts_TSS_480_530[names(pol2_hypbKD_RPGC_readsCounts_TSS_480_530) %in% GNref]
pol2_hypbKD_RPGC_readsCounts_TSS_550_600_f = pol2_hypbKD_RPGC_readsCounts_TSS_550_600[names(pol2_hypbKD_RPGC_readsCounts_TSS_550_600) %in% GNref]
pol2_kdm4aKD_RPGC_readsCounts_TSS_480_530_f = pol2_kdm4aKD_RPGC_readsCounts_TSS_480_530[names(pol2_kdm4aKD_RPGC_readsCounts_TSS_480_530) %in% GNref]
pol2_kdm4aKD_RPGC_readsCounts_TSS_550_600_f = pol2_kdm4aKD_RPGC_readsCounts_TSS_550_600[names(pol2_kdm4aKD_RPGC_readsCounts_TSS_550_600) %in% GNref]
pol2_mes4KD_RPGC_readsCounts_TSS_480_530_f = pol2_mes4KD_RPGC_readsCounts_TSS_480_530[names(pol2_mes4KD_RPGC_readsCounts_TSS_480_530) %in% GNref]
pol2_mes4KD_RPGC_readsCounts_TSS_550_600_f = pol2_mes4KD_RPGC_readsCounts_TSS_550_600[names(pol2_mes4KD_RPGC_readsCounts_TSS_550_600) %in% GNref]
pol2_nelf_N_RPGC_readsCounts_TSS_480_530_f = pol2_nelf_N_RPGC_readsCounts_TSS_480_530[names(pol2_nelf_N_RPGC_readsCounts_TSS_480_530) %in% GNref]
pol2_nelf_N_RPGC_readsCounts_TSS_550_600_f = pol2_nelf_N_RPGC_readsCounts_TSS_550_600[names(pol2_nelf_N_RPGC_readsCounts_TSS_550_600) %in% GNref]
pol2ser2P_ctrl_N_RPGC_readsCounts_TSS_480_530_f = pol2ser2P_ctrl_N_RPGC_readsCounts_TSS_480_530[names(pol2ser2P_ctrl_N_RPGC_readsCounts_TSS_480_530) %in% GNref]
pol2ser2P_ctrl_N_RPGC_readsCounts_TSS_550_600_f = pol2ser2P_ctrl_N_RPGC_readsCounts_TSS_550_600[names(pol2ser2P_ctrl_N_RPGC_readsCounts_TSS_550_600) %in% GNref]
pol2ser2P_ctrl_RPGC_readsCounts_TSS_480_530_f = pol2ser2P_ctrl_RPGC_readsCounts_TSS_480_530[names(pol2ser2P_ctrl_RPGC_readsCounts_TSS_480_530) %in% GNref]
pol2ser2P_ctrl_RPGC_readsCounts_TSS_550_600_f = pol2ser2P_ctrl_RPGC_readsCounts_TSS_550_600[names(pol2ser2P_ctrl_RPGC_readsCounts_TSS_550_600) %in% GNref]
pol2ser2P_hypbKD_RPGC_readsCounts_TSS_480_530_f = pol2ser2P_hypbKD_RPGC_readsCounts_TSS_480_530[names(pol2ser2P_hypbKD_RPGC_readsCounts_TSS_480_530) %in% GNref]
pol2ser2P_hypbKD_RPGC_readsCounts_TSS_550_600_f = pol2ser2P_hypbKD_RPGC_readsCounts_TSS_550_600[names(pol2ser2P_hypbKD_RPGC_readsCounts_TSS_550_600) %in% GNref]
pol2ser2P_mes4KD_RPGC_readsCounts_TSS_480_530_f = pol2ser2P_mes4KD_RPGC_readsCounts_TSS_480_530[names(pol2ser2P_mes4KD_RPGC_readsCounts_TSS_480_530) %in% GNref]
pol2ser2P_mes4KD_RPGC_readsCounts_TSS_550_600_f = pol2ser2P_mes4KD_RPGC_readsCounts_TSS_550_600[names(pol2ser2P_mes4KD_RPGC_readsCounts_TSS_550_600) %in% GNref]
pol2ser2P_nelf_N_RPGC_readsCounts_TSS_480_530_f = pol2ser2P_nelf_N_RPGC_readsCounts_TSS_480_530[names(pol2ser2P_nelf_N_RPGC_readsCounts_TSS_480_530) %in% GNref]
pol2ser2P_nelf_N_RPGC_readsCounts_TSS_550_600_f = pol2ser2P_nelf_N_RPGC_readsCounts_TSS_550_600[names(pol2ser2P_nelf_N_RPGC_readsCounts_TSS_550_600) %in% GNref]
polIIC_RPGC_readsCounts_TSS_480_530_f = polIIC_RPGC_readsCounts_TSS_480_530[names(polIIC_RPGC_readsCounts_TSS_480_530) %in% GNref]
polIIC_RPGC_readsCounts_TSS_550_600_f = polIIC_RPGC_readsCounts_TSS_550_600[names(polIIC_RPGC_readsCounts_TSS_550_600) %in% GNref]
polIIH_RPGC_readsCounts_TSS_480_530_f = polIIH_RPGC_readsCounts_TSS_480_530[names(polIIH_RPGC_readsCounts_TSS_480_530) %in% GNref]
polIIH_RPGC_readsCounts_TSS_550_600_f = polIIH_RPGC_readsCounts_TSS_550_600[names(polIIH_RPGC_readsCounts_TSS_550_600) %in% GNref]
RNAPolII_control_rep1_RPGC_readsCounts_TSS_480_530_f = RNAPolII_control_rep1_RPGC_readsCounts_TSS_480_530[names(RNAPolII_control_rep1_RPGC_readsCounts_TSS_480_530) %in% GNref]
RNAPolII_control_rep1_RPGC_readsCounts_TSS_550_600_f = RNAPolII_control_rep1_RPGC_readsCounts_TSS_550_600[names(RNAPolII_control_rep1_RPGC_readsCounts_TSS_550_600) %in% GNref]


# COMPUTE PAUSE INDICE

PAUSE_IND_pol2_ctrl_N = log2(pol2_ctrl_N_RPGC_readsCounts_TSS_480_530_f+1)/log2(pol2_ctrl_N_RPGC_readsCounts_TSS_550_600_f+1)
PAUSE_IND_pol2_ctrl = log2(pol2_ctrl_RPGC_readsCounts_TSS_480_530_f+1)/log2(pol2_ctrl_RPGC_readsCounts_TSS_550_600_f+1)
PAUSE_IND_pol2_hypbKD = log2(pol2_hypbKD_RPGC_readsCounts_TSS_480_530_f+1)/log2(pol2_hypbKD_RPGC_readsCounts_TSS_550_600_f+1)
PAUSE_IND_pol2_kdm4aKD = log2(pol2_kdm4aKD_RPGC_readsCounts_TSS_480_530_f+1)/log2(pol2_kdm4aKD_RPGC_readsCounts_TSS_550_600_f+1)
PAUSE_IND_pol2_mes4KD = log2(pol2_mes4KD_RPGC_readsCounts_TSS_480_530_f+1)/log2(pol2_mes4KD_RPGC_readsCounts_TSS_550_600_f+1)
PAUSE_IND_pol2_nelf_N = log2(pol2_nelf_N_RPGC_readsCounts_TSS_480_530_f+1)/log2(pol2_nelf_N_RPGC_readsCounts_TSS_550_600_f+1)
PAUSE_IND_pol2ser2P_ctrl_N = log2(pol2ser2P_ctrl_N_RPGC_readsCounts_TSS_480_530_f+1)/log2(pol2ser2P_ctrl_N_RPGC_readsCounts_TSS_550_600_f+1)
PAUSE_IND_pol2ser2P_ctrl = log2(pol2ser2P_ctrl_RPGC_readsCounts_TSS_480_530_f+1)/log2(pol2ser2P_ctrl_RPGC_readsCounts_TSS_550_600_f+1)
PAUSE_IND_pol2ser2P_hypbKD = log2(pol2ser2P_hypbKD_RPGC_readsCounts_TSS_480_530_f+1)/log2(pol2ser2P_hypbKD_RPGC_readsCounts_TSS_550_600_f+1)
PAUSE_IND_pol2ser2P_mes4KD = log2(pol2ser2P_mes4KD_RPGC_readsCounts_TSS_480_530_f+1)/log2(pol2ser2P_mes4KD_RPGC_readsCounts_TSS_550_600_f+1)
PAUSE_IND_pol2ser2P_nelf_N = log2(pol2ser2P_nelf_N_RPGC_readsCounts_TSS_480_530_f+1)/log2(pol2ser2P_nelf_N_RPGC_readsCounts_TSS_550_600_f+1)
PAUSE_IND_polIIC = log2(polIIC_RPGC_readsCounts_TSS_480_530_f+1)/log2(polIIC_RPGC_readsCounts_TSS_550_600_f+1)
PAUSE_IND_polIIH = log2(polIIH_RPGC_readsCounts_TSS_480_530_f+1)/log2(polIIH_RPGC_readsCounts_TSS_550_600_f+1)
PAUSE_IND_RNAPolII_control = log2(RNAPolII_control_rep1_RPGC_readsCounts_TSS_480_530_f+1)/log2(RNAPolII_control_rep1_RPGC_readsCounts_TSS_550_600_f+1)


# REPLACE NA BY 0

PAUSE_IND_pol2_ctrl_N[is.na(PAUSE_IND_pol2_ctrl_N)] <- 0
PAUSE_IND_pol2_ctrl[is.na(PAUSE_IND_pol2_ctrl)] <- 0
PAUSE_IND_pol2_hypbKD[is.na(PAUSE_IND_pol2_hypbKD)] <- 0
PAUSE_IND_pol2_kdm4aKD[is.na(PAUSE_IND_pol2_kdm4aKD)] <- 0
PAUSE_IND_pol2_mes4KD[is.na(PAUSE_IND_pol2_mes4KD)] <- 0
PAUSE_IND_pol2_nelf_N[is.na(PAUSE_IND_pol2_nelf_N)] <- 0
PAUSE_IND_pol2ser2P_ctrl_N[is.na(PAUSE_IND_pol2ser2P_ctrl_N)] <- 0
PAUSE_IND_pol2ser2P_ctrl[is.na(PAUSE_IND_pol2ser2P_ctrl)] <- 0
PAUSE_IND_pol2ser2P_hypbKD[is.na(PAUSE_IND_pol2ser2P_hypbKD)] <- 0
PAUSE_IND_pol2ser2P_mes4KD[is.na(PAUSE_IND_pol2ser2P_mes4KD)] <- 0
PAUSE_IND_pol2ser2P_nelf_N[is.na(PAUSE_IND_pol2ser2P_nelf_N)] <- 0
PAUSE_IND_polIIC[is.na(PAUSE_IND_polIIC)] <- 0
PAUSE_IND_polIIH[is.na(PAUSE_IND_polIIH)] <- 0
PAUSE_IND_RNAPolII_control[is.na(PAUSE_IND_RNAPolII_control)] <- 0


# REPLACE INFINITE VALUES BY MAX 

PAUSE_IND_pol2_ctrl_N[is.finite(PAUSE_IND_pol2_ctrl_N)==F] <- max(PAUSE_IND_pol2_ctrl_N[is.finite(PAUSE_IND_pol2_ctrl_N)])
PAUSE_IND_pol2_ctrl[is.finite(PAUSE_IND_pol2_ctrl)==F] <- max(PAUSE_IND_pol2_ctrl[is.finite(PAUSE_IND_pol2_ctrl)])
PAUSE_IND_pol2_hypbKD[is.finite(PAUSE_IND_pol2_hypbKD)==F] <- max(PAUSE_IND_pol2_hypbKD[is.finite(PAUSE_IND_pol2_hypbKD)])
PAUSE_IND_pol2_kdm4aKD[is.finite(PAUSE_IND_pol2_kdm4aKD)==F] <- max(PAUSE_IND_pol2_kdm4aKD[is.finite(PAUSE_IND_pol2_kdm4aKD)])
PAUSE_IND_pol2_mes4KD[is.finite(PAUSE_IND_pol2_mes4KD)==F] <- max(PAUSE_IND_pol2_mes4KD[is.finite(PAUSE_IND_pol2_mes4KD)])
PAUSE_IND_pol2_nelf_N[is.finite(PAUSE_IND_pol2_nelf_N)==F] <- max(PAUSE_IND_pol2_nelf_N[is.finite(PAUSE_IND_pol2_nelf_N)])
PAUSE_IND_pol2ser2P_ctrl_N[is.finite(PAUSE_IND_pol2ser2P_ctrl_N)==F] <- max(PAUSE_IND_pol2ser2P_ctrl_N[is.finite(PAUSE_IND_pol2ser2P_ctrl_N)])
PAUSE_IND_pol2ser2P_ctrl[is.finite(PAUSE_IND_pol2ser2P_ctrl)==F] <- max(PAUSE_IND_pol2ser2P_ctrl[is.finite(PAUSE_IND_pol2ser2P_ctrl)])
PAUSE_IND_pol2ser2P_hypbKD[is.finite(PAUSE_IND_pol2ser2P_hypbKD)==F] <- max(PAUSE_IND_pol2ser2P_hypbKD[is.finite(PAUSE_IND_pol2ser2P_hypbKD)])
PAUSE_IND_pol2ser2P_mes4KD[is.finite(PAUSE_IND_pol2ser2P_mes4KD)==F] <- max(PAUSE_IND_pol2ser2P_mes4KD[is.finite(PAUSE_IND_pol2ser2P_mes4KD)])
PAUSE_IND_pol2ser2P_nelf_N[is.finite(PAUSE_IND_pol2ser2P_nelf_N)==F] <- max(PAUSE_IND_pol2ser2P_nelf_N[is.finite(PAUSE_IND_pol2ser2P_nelf_N)])
PAUSE_IND_polIIC[is.finite(PAUSE_IND_polIIC)==F] <- max(PAUSE_IND_polIIC[is.finite(PAUSE_IND_polIIC)])
PAUSE_IND_polIIH[is.finite(PAUSE_IND_polIIH)==F] <- max(PAUSE_IND_polIIH[is.finite(PAUSE_IND_polIIH)])
PAUSE_IND_RNAPolII_control[is.finite(PAUSE_IND_RNAPolII_control)==F] <- max(PAUSE_IND_RNAPolII_control[is.finite(PAUSE_IND_RNAPolII_control)])


DELTA_PAUSE_IND_pol2_hypbKD = PAUSE_IND_pol2_hypbKD-PAUSE_IND_pol2_ctrl
DELTA_PAUSE_IND_pol2_kdm4aKD = PAUSE_IND_pol2_kdm4aKD-PAUSE_IND_pol2_ctrl
DELTA_PAUSE_IND_pol2_mes4KD = PAUSE_IND_pol2_mes4KD-PAUSE_IND_pol2_ctrl
DELTA_PAUSE_IND_pol2_nelf_N = PAUSE_IND_pol2_nelf_N-PAUSE_IND_pol2_ctrl_N
DELTA_PAUSE_IND_pol2ser2P_hypbKD = PAUSE_IND_pol2ser2P_hypbKD-PAUSE_IND_pol2ser2P_ctrl
DELTA_PAUSE_IND_pol2ser2P_mes4KD = PAUSE_IND_pol2ser2P_mes4KD-PAUSE_IND_pol2ser2P_ctrl
DELTA_PAUSE_IND_pol2ser2P_nelf_N = PAUSE_IND_pol2ser2P_nelf_N-PAUSE_IND_pol2ser2P_ctrl_N


# ORDER decreasingly

DELTA_PAUSE_IND_pol2_hypbKD = DELTA_PAUSE_IND_pol2_hypbKD[order(DELTA_PAUSE_IND_pol2_hypbKD, decreasing=T)]
DELTA_PAUSE_IND_pol2_kdm4aKD = DELTA_PAUSE_IND_pol2_kdm4aKD[order(DELTA_PAUSE_IND_pol2_kdm4aKD, decreasing=T)]
DELTA_PAUSE_IND_pol2_mes4KD = DELTA_PAUSE_IND_pol2_mes4KD[order(DELTA_PAUSE_IND_pol2_mes4KD, decreasing=T)]
DELTA_PAUSE_IND_pol2_nelf_N = DELTA_PAUSE_IND_pol2_nelf_N[order(DELTA_PAUSE_IND_pol2_nelf_N, decreasing=T)]
DELTA_PAUSE_IND_pol2ser2P_hypbKD = DELTA_PAUSE_IND_pol2ser2P_hypbKD[order(DELTA_PAUSE_IND_pol2ser2P_hypbKD, decreasing=T)]
DELTA_PAUSE_IND_pol2ser2P_mes4KD = DELTA_PAUSE_IND_pol2ser2P_mes4KD[order(DELTA_PAUSE_IND_pol2ser2P_mes4KD, decreasing=T)]
DELTA_PAUSE_IND_pol2ser2P_nelf_N = DELTA_PAUSE_IND_pol2ser2P_nelf_N[order(DELTA_PAUSE_IND_pol2ser2P_nelf_N, decreasing=T)]
PAUSE_IND_pol2_ctrl_N = PAUSE_IND_pol2_ctrl_N[order(PAUSE_IND_pol2_ctrl_N, decreasing=T)]
PAUSE_IND_pol2_ctrl = PAUSE_IND_pol2_ctrl[order(PAUSE_IND_pol2_ctrl, decreasing=T)]
PAUSE_IND_pol2_hypbKD = PAUSE_IND_pol2_hypbKD[order(PAUSE_IND_pol2_hypbKD, decreasing=T)]
PAUSE_IND_pol2_kdm4aKD = PAUSE_IND_pol2_kdm4aKD[order(PAUSE_IND_pol2_kdm4aKD, decreasing=T)]
PAUSE_IND_pol2_mes4KD = PAUSE_IND_pol2_mes4KD[order(PAUSE_IND_pol2_mes4KD, decreasing=T)]
PAUSE_IND_pol2_nelf_N = PAUSE_IND_pol2_nelf_N[order(PAUSE_IND_pol2_nelf_N, decreasing=T)]
PAUSE_IND_pol2ser2P_ctrl_N = PAUSE_IND_pol2ser2P_ctrl_N[order(PAUSE_IND_pol2ser2P_ctrl_N, decreasing=T)]
PAUSE_IND_pol2ser2P_ctrl = PAUSE_IND_pol2ser2P_ctrl[order(PAUSE_IND_pol2ser2P_ctrl, decreasing=T)]
PAUSE_IND_pol2ser2P_hypbKD = PAUSE_IND_pol2ser2P_hypbKD[order(PAUSE_IND_pol2ser2P_hypbKD, decreasing=T)]
PAUSE_IND_pol2ser2P_mes4KD = PAUSE_IND_pol2ser2P_mes4KD[order(PAUSE_IND_pol2ser2P_mes4KD, decreasing=T)]
PAUSE_IND_pol2ser2P_nelf_N = PAUSE_IND_pol2ser2P_nelf_N[order(PAUSE_IND_pol2ser2P_nelf_N, decreasing=T)]
PAUSE_IND_polIIC = PAUSE_IND_polIIC[order(PAUSE_IND_polIIC, decreasing=T)]
PAUSE_IND_polIIH = PAUSE_IND_polIIH[order(PAUSE_IND_polIIH, decreasing=T)]
PAUSE_IND_RNAPolII_control = PAUSE_IND_RNAPolII_control[order(PAUSE_IND_RNAPolII_control, decreasing=T)]


# CREATE LIST

PAUSE_INDICE_VEC = list(
  DELTA_PAUSE_IND_pol2_hypbKD = DELTA_PAUSE_IND_pol2_hypbKD,
  DELTA_PAUSE_IND_pol2_kdm4aKD = DELTA_PAUSE_IND_pol2_kdm4aKD,
  DELTA_PAUSE_IND_pol2_mes4KD = DELTA_PAUSE_IND_pol2_mes4KD,
  DELTA_PAUSE_IND_pol2_nelf_N = DELTA_PAUSE_IND_pol2_nelf_N,
  DELTA_PAUSE_IND_pol2ser2P_hypbKD = DELTA_PAUSE_IND_pol2ser2P_hypbKD,
  DELTA_PAUSE_IND_pol2ser2P_mes4KD = DELTA_PAUSE_IND_pol2ser2P_mes4KD,
  DELTA_PAUSE_IND_pol2ser2P_nelf_N = DELTA_PAUSE_IND_pol2ser2P_nelf_N,
  PAUSE_IND_pol2_ctrl_N = PAUSE_IND_pol2_ctrl_N,
  PAUSE_IND_pol2_ctrl = PAUSE_IND_pol2_ctrl,
  PAUSE_IND_pol2_hypbKD = PAUSE_IND_pol2_hypbKD,
  PAUSE_IND_pol2_kdm4aKD = PAUSE_IND_pol2_kdm4aKD,
  PAUSE_IND_pol2_mes4KD = PAUSE_IND_pol2_mes4KD,
  PAUSE_IND_pol2_nelf_N = PAUSE_IND_pol2_nelf_N,
  PAUSE_IND_pol2ser2P_ctrl_N = PAUSE_IND_pol2ser2P_ctrl_N,
  PAUSE_IND_pol2ser2P_ctrl = PAUSE_IND_pol2ser2P_ctrl,
  PAUSE_IND_pol2ser2P_hypbKD = PAUSE_IND_pol2ser2P_hypbKD,
  PAUSE_IND_pol2ser2P_mes4KD = PAUSE_IND_pol2ser2P_mes4KD,
  PAUSE_IND_pol2ser2P_nelf_N = PAUSE_IND_pol2ser2P_nelf_N,
  PAUSE_IND_polIIC = PAUSE_IND_polIIC,
  PAUSE_IND_polIIH = PAUSE_IND_polIIH,
  PAUSE_IND_RNAPolII_control = PAUSE_IND_RNAPolII_control
)

saveRDS(PAUSE_INDICE_VEC, paste0(workdir, "PROJET_H2AV_2021/DATA/LIST_FEATURES/PAUSE_INDICE_VEC.RDS"))
PAUSE_INDICE_VEC = readRDS(paste0(workdir, "PROJET_H2AV_2021/DATA/LIST_FEATURES/PAUSE_INDICE_VEC.RDS"))



#########################################################################################################################
######################################  H2AV / POL2 QUANTIFICATION    ###################################################

# LOAD DATA

### Quantif TSS 

workdir="/home/cperrois/work/"
Q_H2AV_TSS_PW_R2=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/Q_H2AV_PW_2_L1_RPGC_upstr500_dnstr500.RDS"))
Q_H2AV_TSS_PW_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/Q_H2AV_PW_1_L1_RPGC_upstr500_dnstr500.RDS"))
Q_H2AV_TSS_PWH_inputnorm_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/PWH_inputnormRPGC_readsCounts_TSS.RDS"))
Q_H2AV_TSS_PWH_R1  =readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/PWH_RPGC_readsCounts_TSS.RDS"))
Q_H2AV_TSS_PWH_R2=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/Q_H2AV_PWH_A_L1_RPGC_upstr500_dnstr500.RDS"))
# nelf TSS_WT KD / KD + HU
Q_H2AV_TSS_PN_inputnorm_R1 =readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/PN_inputnormRPGC_readsCounts_TSS.RDS"))
Q_H2AV_TSS_PN_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/PN_RPGC_readsCounts_TSS.RDS"))
Q_H2AV_TSS_PNH_L1_inputnorm_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/PNH_L1_inputnormRPGC_readsCounts_TSS.RDS"))
Q_H2AV_TSS_PNH_L1_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/PNH_L1_RPGC_readsCounts_TSS.RDS"))
Q_H2AV_TSS_PN_REP_R2=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/Q_H2AV_PN_J_L1_RPGC_upstr500_dnstr500.RDS"))
Q_H2AV_TSS_PNH_REP_R2=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/Q_H2AV_PNH_J_L1_RPGC_upstr500_dnstr500.RDS"))
# Hypb TSS_WT/KD / KD + HU  
Q_H2AV_TSS_PHYPH_inputnorm_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/PHYPH_B_L1_inputnormRPGC_readsCounts_TSS.RDS"))
Q_H2AV_TSS_PHYPH_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/PHYPH_B_L1_RPGC_readsCounts_TSS.RDS"))
Q_H2AV_TSS_PHYP_inputnorm_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/PHYP_B_L1_inputnormRPGC_readsCounts_TSS.RDS"))
Q_H2AV_TSS_PHYP_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/PHYP_B_L1_RPGC_readsCounts_TSS.RDS"))
Q_H2AV_TSS_PHYPB_R2=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/Q_H2AV_PHYPB_A_L1_RPGC_upstr500_dnstr500.RDS"))
Q_H2AV_TSS_PHYPBH_R2=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/Q_H2AV_PHYPBH_A_L1_RPGC_upstr500_dnstr500.RDS"))

Q_POLII_TSS_ctrl_N=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/pol2_ctrl_N_RPGC_readsCounts_TSS.RDS"))
Q_POLII_TSS_ctrl=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/pol2_ctrl_RPGC_readsCounts_TSS.RDS"))
Q_POLII_TSS_hypbKD=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/pol2_hypbKD_RPGC_readsCounts_TSS.RDS"))
Q_POLII_TSS_nelf_N=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/pol2_nelf_N_RPGC_readsCounts_TSS.RDS"))
Q_POLII_TSS_ser2P_ctrl_N=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/pol2ser2P_ctrl_N_RPGC_readsCounts_TSS.RDS"))
Q_POLII_TSS_ser2P_ctrl=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/pol2ser2P_ctrl_RPGC_readsCounts_TSS.RDS"))
Q_POLII_TSS_ser2P_hypbKD=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/pol2ser2P_hypbKD_RPGC_readsCounts_TSS.RDS"))
Q_POLII_TSS_ser2P_nelf_N=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/TSS/pol2ser2P_nelf_N_RPGC_readsCounts_TSS.RDS"))

### Quantif GB

workdir="/home/cperrois/work/"
GNref = readRDS(paste0(workdir, "PROJET_H2AV/DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))
# WT + WT+HU 
Q_H2AV_GB_PW_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/Q_H2AV_PW_1_L1_RPGC_upstr500_dnstr500.RDS"))
Q_H2AV_GB_PW_R2=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/Q_H2AV_PW_2_L1_RPGC_upstr500_dnstr500.RDS"))

Q_H2AV_GB_PWH_inputnorm_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/PWH_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_H2AV_GB_PWH_R1  =readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/PWH_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H2AV_GB_PWH_REP2=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/Q_H2AV_PWH_A_L1_RPGC_upstr500_dnstr500.RDS"))
# nelf_ WT KD / KD + HU
Q_H2AV_GB_PN_inputnorm_R1 =readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/PN_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_H2AV_GB_PN_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/PN_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H2AV_GB_PNH_L1_inputnorm_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/PNH_L1_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_H2AV_GB_PNH_L1_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/PNH_L1_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H2AV_GB_PN_REP_R2=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/Q_H2AV_PN_J_L1_RPGC_upstr500_dnstr500.RDS"))
Q_H2AV_GB_PNH_REP_R2=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/Q_H2AV_PNH_J_L1_RPGC_upstr500_dnstr500.RDS"))
# Hypb_WT/KD / KD + HU  
Q_H2AV_GB_PHYPH_inputnorm_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/PHYPH_B_L1_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_H2AV_GB_PHYPH_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/PHYPH_B_L1_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H2AV_GB_PHYP_inputnorm_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/PHYP_B_L1_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_H2AV_GB_PHYP_R1=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/PHYP_B_L1_RPGC_readsCounts_GB_SCALED.RDS"))
Q_H2AV_GB_PHYPB_REP2=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/Q_H2AV_PHYPB_A_L1_RPGC_upstr500_dnstr500.RDS"))
Q_H2AV_GB_PHYPBH_REP2=readRDS(paste0(workdir,"PROJET_H2AV_2021/DATA/QUANTIF/GB/Q_H2AV_PHYPBH_A_L1_RPGC_upstr500_dnstr500.RDS"))

# ZSCORES

ZSCORE_H2AV_R3_HYPB_WT = computeZscore(Q_H2AV_hypb_R3,Q_H2AV_ctrl_R3)
ZSCORE_H2AV_NELF_WT = computeZscore(Q_H2AV_nelf,Q_H2AV_ctrl)
ZSCORE_H2AV_HYPB_WT = computeZscore(Q_H2AV_hypb,Q_H2AV_ctrl)
ZSCORE_RAD51_WT3_GB_N_WT = computeZscore(RAD51_N1_WT3,RAD51_L1_WT3)
ZSCORE_H2AV_WT3_GB_N_WT = computeZscore(Q_GB_WT3_NELF,Q_GB_WT3_WT)
ZSCORE_H2AV_SEA4_GB_N_WT = computeZscore(Q_GB_SEA4_NELF,Q_GB_SEA4_WT)
ZSCORE_RAD51_GB_WT_N = computeZscore(RAD51_ctrl,RAD51_nelf)
ZSCORE_H2AV_GB_WT_N = computeZscore(H2AV_ctrl,H2AV_nelf)
ZSCORE_PROFMAT_H2AV_WT_vs_NELF = computeZscore_PROFMAT(PROFMAT_H2AV_PW,PROFMAT_H2AV_PN)




### ADD '.1' to gene vector name for compatibility
for (i in 1: length(Q_H2AV_GB_PW_R1)){names(Q_H2AV_GB_PW_R1)[i]=paste0(names(Q_H2AV_GB_PW_R1)[i],".1")} 
for (i in 1: length(Q_H2AV_GB_PW_R2)){names(Q_H2AV_GB_PW_R2)[i]=paste0(names(Q_H2AV_GB_PW_R2)[i],".1")} 
for (i in 1: length(Q_H2AV_GB_PWH_REP2)){names(Q_H2AV_GB_PWH_REP2)[i]=paste0(names(Q_H2AV_GB_PWH_REP2)[i],".1")} 
for (i in 1: length(Q_H2AV_GB_PN_REP_R2)){names(Q_H2AV_GB_PN_REP_R2)[i]=paste0(names(Q_H2AV_GB_PN_REP_R2)[i],".1")} 
for (i in 1: length(Q_H2AV_GB_PNH_REP_R2)){names(Q_H2AV_GB_PNH_REP_R2)[i]=paste0(names(Q_H2AV_GB_PNH_REP_R2)[i],".1")} 
for (i in 1: length(Q_H2AV_GB_PHYPBH_REP2)){names(Q_H2AV_GB_PHYPBH_REP2)[i]=paste0(names(Q_H2AV_GB_PHYPBH_REP2)[i],".1")} 
for (i in 1: length(Q_H2AV_GB_PHYPB_REP2)){names(Q_H2AV_GB_PHYPB_REP2)[i]=paste0(names(Q_H2AV_GB_PHYPB_REP2)[i],".1")} 

# FILTER ACTIVE GENES
Q_H2AV_GB_PW_R1_f=Q_H2AV_GB_PW_R1[names(Q_H2AV_GB_PW_R1) %in% GNref]
Q_H2AV_GB_PW_R2_f=Q_H2AV_GB_PW_R2[names(Q_H2AV_GB_PW_R2) %in% GNref]
Q_H2AV_GB_PWH_inputnorm_R1_f=Q_H2AV_GB_PWH_inputnorm_R1[names(Q_H2AV_GB_PWH_inputnorm_R1) %in% GNref]
Q_H2AV_GB_PWH_R1_f=Q_H2AV_GB_PWH_R1[names(Q_H2AV_GB_PWH_R1) %in% GNref] 
Q_H2AV_GB_PWH_REP2_f=Q_H2AV_GB_PWH_REP2[names(Q_H2AV_GB_PWH_REP2) %in% GNref]
Q_H2AV_GB_PN_inputnorm_R1_f=Q_H2AV_GB_PN_inputnorm_R1[names(Q_H2AV_GB_PN_inputnorm_R1) %in% GNref]
Q_H2AV_GB_PN_R1_f=Q_H2AV_GB_PN_R1[names(Q_H2AV_GB_PN_R1) %in% GNref]
Q_H2AV_GB_PNH_L1_inputnorm_R1_f=Q_H2AV_GB_PNH_L1_inputnorm_R1[names(Q_H2AV_GB_PNH_L1_inputnorm_R1) %in% GNref]
Q_H2AV_GB_PNH_L1_R1_f=Q_H2AV_GB_PNH_L1_R1[names(Q_H2AV_GB_PNH_L1_R1) %in% GNref]
Q_H2AV_GB_PN_REP_R2_f=Q_H2AV_GB_PN_REP_R2[names(Q_H2AV_GB_PN_REP_R2) %in% GNref]
Q_H2AV_GB_PNH_REP_R2_f=Q_H2AV_GB_PNH_REP_R2[names(Q_H2AV_GB_PNH_REP_R2) %in% GNref]
Q_H2AV_GB_PHYPH_inputnorm_R1_f=Q_H2AV_GB_PHYPH_inputnorm_R1[names(Q_H2AV_GB_PHYPH_inputnorm_R1) %in% GNref]
Q_H2AV_GB_PHYPH_R1_f=Q_H2AV_GB_PHYPH_R1[names(Q_H2AV_GB_PHYPH_R1) %in% GNref]
Q_H2AV_GB_PHYP_inputnorm_R1_f=Q_H2AV_GB_PHYP_inputnorm_R1[names(Q_H2AV_GB_PHYP_inputnorm_R1) %in% GNref]
Q_H2AV_GB_PHYP_R1_f=Q_H2AV_GB_PHYP_R1[names(Q_H2AV_GB_PHYP_R1) %in% GNref]
Q_H2AV_GB_PHYPB_REP2_f=Q_H2AV_GB_PHYPB_REP2[names(Q_H2AV_GB_PHYPB_REP2) %in% GNref]
Q_H2AV_GB_PHYPBH_REP2_f=Q_H2AV_GB_PHYPBH_REP2[names(Q_H2AV_GB_PHYPBH_REP2) %in% GNref]

Q_H2AV_TSS_PWH_inputnorm_R1_f=Q_H2AV_TSS_PWH_inputnorm_R1[names(Q_H2AV_TSS_PWH_inputnorm_R1) %in% GNref]     
Q_H2AV_TSS_PWH_R1_f=Q_H2AV_TSS_PWH_R1[names(Q_H2AV_TSS_PWH_R1) %in% GNref]              
Q_H2AV_TSS_PWH_R2_f=Q_H2AV_TSS_PWH_R2[names(Q_H2AV_TSS_PWH_R2) %in% GNref]  
Q_H2AV_TSS_PN_inputnorm_R1_f=Q_H2AV_TSS_PN_inputnorm_R1[names(Q_H2AV_TSS_PN_inputnorm_R1) %in% GNref]     
Q_H2AV_TSS_PN_R1_f=Q_H2AV_TSS_PN_R1[names(Q_H2AV_TSS_PN_R1) %in% GNref]
Q_H2AV_TSS_PN_R2_f=Q_H2AV_TSS_PN_REP_R2[names(Q_H2AV_TSS_PN_REP_R2) %in% GNref]
Q_H2AV_TSS_PNH_L1_inputnorm_R1_f=Q_H2AV_TSS_PNH_L1_inputnorm_R1[names(Q_H2AV_TSS_PNH_L1_inputnorm_R1) %in% GNref] 
Q_H2AV_TSS_PNH_L1_R1_f=Q_H2AV_TSS_PNH_L1_R1[names(Q_H2AV_TSS_PNH_L1_R1) %in% GNref]            
Q_H2AV_TSS_PNH_R2_f=Q_H2AV_TSS_PNH_REP_R2[names(Q_H2AV_TSS_PNH_REP_R2) %in% GNref]
Q_H2AV_TSS_PHYPH_inputnorm_R1_f=Q_H2AV_TSS_PHYPH_inputnorm_R1[names(Q_H2AV_TSS_PHYPH_inputnorm_R1) %in% GNref]  
Q_H2AV_TSS_PHYPH_R1_f=Q_H2AV_TSS_PHYPH_R1[names(Q_H2AV_TSS_PHYPH_R1) %in% GNref]             
Q_H2AV_TSS_PHYPBH_R2_f=Q_H2AV_TSS_PHYPBH_R2[names(Q_H2AV_TSS_PHYPBH_R2) %in% GNref] 
Q_H2AV_TSS_PHYP_inputnorm_R1_f=Q_H2AV_TSS_PHYP_inputnorm_R1[names(Q_H2AV_TSS_PHYP_inputnorm_R1) %in% GNref]   
Q_H2AV_TSS_PHYP_R1_f=Q_H2AV_TSS_PHYP_R1[names(Q_H2AV_TSS_PHYP_R1) %in% GNref]              
Q_H2AV_TSS_PHYPB_R2_f=Q_H2AV_TSS_PHYPB_R2[names(Q_H2AV_TSS_PHYPB_R2) %in% GNref]               
Q_POLII_ctrl_N_f=Q_POLII_ctrl_N[names(Q_POLII_ctrl_N) %in% GNref]                 
Q_POLII_ctrl_f=Q_POLII_ctrl[names(Q_POLII_ctrl) %in% GNref]                    
Q_POLII_hypbKD_f=Q_POLII_hypbKD[names(Q_POLII_hypbKD) %in% GNref]                 
Q_POLII_nelf_N_f=Q_POLII_nelf_N[names(Q_POLII_nelf_N) %in% GNref]                  
Q_POLII_ser2P_ctrl_N_f=Q_POLII_ser2P_ctrl_N[names(Q_POLII_ser2P_ctrl_N) %in% GNref]           
Q_POLII_ser2P_ctrl_f=Q_POLII_ser2P_ctrl[names(Q_POLII_ser2P_ctrl) %in% GNref]              
Q_POLII_ser2P_hypbKD_f=Q_POLII_ser2P_hypbKD[names(Q_POLII_ser2P_hypbKD) %in% GNref]           
Q_POLII_ser2P_nelf_N_f=Q_POLII_ser2P_nelf_N[names(Q_POLII_ser2P_nelf_N) %in% GNref] 

ZSCORE_H2AV_R3_HYPB_WT = ZSCORE_H2AV_R3_HYPB_WT[names(ZSCORE_H2AV_R3_HYPB_WT) %in% GNref]
ZSCORE_H2AV_NELF_WT = ZSCORE_H2AV_NELF_WT[names(ZSCORE_H2AV_NELF_WT) %in% GNref]
ZSCORE_H2AV_HYPB_WT = ZSCORE_H2AV_HYPB_WT[names(ZSCORE_H2AV_HYPB_WT) %in% GNref]
ZSCORE_RAD51_WT3_GB_N_WT = ZSCORE_RAD51_WT3_GB_N_WT[names(ZSCORE_RAD51_WT3_GB_N_WT) %in% GNref]
ZSCORE_H2AV_WT3_GB_N_WT = ZSCORE_H2AV_WT3_GB_N_WT[names(ZSCORE_H2AV_WT3_GB_N_WT) %in% GNref]
ZSCORE_H2AV_SEA4_GB_N_WT = ZSCORE_H2AV_SEA4_GB_N_WT[names(ZSCORE_H2AV_SEA4_GB_N_WT) %in% GNref]
ZSCORE_RAD51_GB_WT_N = ZSCORE_RAD51_GB_WT_N[names(ZSCORE_RAD51_GB_WT_N) %in% GNref]
ZSCORE_H2AV_GB_WT_N = ZSCORE_H2AV_GB_WT_N[names(ZSCORE_H2AV_GB_WT_N) %in% GNref]


# CREATE LIST

LIST_QUANTIF = list(
  Q_H2AV_TSS_PW_R2=Q_H2AV_PW_R2,
  Q_H2AV_TSS_PW_R1=Q_H2AV_PW_R1,
  Q_H2AV_TSS_PWH_inputnorm_R1=Q_H2AV_PWH_inputnorm_R1,
  Q_H2AV_TSS_PWH_R1 = Q_H2AV_PWH_R1  ,
  Q_H2AV_TSS_PWH_R2=Q_H2AV_PWH_REP2,
  Q_H2AV_TSS_PN_inputnorm_R1= Q_H2AV_PN_inputnorm_R1 ,
  Q_H2AV_TSS_PN_R1=Q_H2AV_PN_R1,
  Q_H2AV_TSS_PNH_L1_inputnorm_R1=Q_H2AV_PNH_L1_inputnorm_R1,
  Q_H2AV_TSS_PNH_L1_R1=Q_H2AV_PNH_L1_R1,
  Q_H2AV_TSS_PN_REP_R2=Q_H2AV_PN_REP_R2,
  Q_H2AV_TSS_PNH_REP_R2=Q_H2AV_PNH_REP_R2,
  Q_H2AV_TSS_PHYPH_inputnorm_R1=Q_H2AV_PHYPH_inputnorm_R1,
  Q_H2AV_TSS_PHYPH_R1=Q_H2AV_PHYPH_R1,
  Q_H2AV_TSS_PHYP_inputnorm_R1=Q_H2AV_PHYP_inputnorm_R1,
  Q_H2AV_TSS_PHYP_R1=Q_H2AV_PHYP_R1,
  Q_H2AV_TSS_PHYPB_R2=Q_H2AV_PHYPB_REP2,
  Q_H2AV_TSS_PHYPBH_R2=Q_H2AV_PHYPBH_REP2,
  Q_POLII_ctrl_N=Q_POLII_ctrl_N,
  Q_POLII_ctrl=Q_POLII_ctrl,
  Q_POLII_hypbKD=Q_POLII_hypbKD,
  Q_POLII_nelf_N=Q_POLII_nelf_N,
  Q_POLII_ser2P_ctrl_N=Q_POLII_ser2P_ctrl_N,
  Q_POLII_ser2P_ctrl=Q_POLII_ser2P_ctrl,
  Q_POLII_ser2P_hypbKD=Q_POLII_ser2P_hypbKD,
  Q_POLII_ser2P_nelf_N=Q_POLII_ser2P_nelf_N,
  Q_H2AV_TSS_PWH_inputnorm_R1_f=Q_H2AV_TSS_PWH_inputnorm_R1_f,
  Q_H2AV_TSS_PWH_R1_f=Q_H2AV_TSS_PWH_R1_f,
  Q_H2AV_TSS_PWH_R2_f=Q_H2AV_TSS_PWH_R2_f,
  Q_H2AV_TSS_PN_inputnorm_R1_f=Q_H2AV_TSS_PN_inputnorm_R1_f,
  Q_H2AV_TSS_PN_R1_f=Q_H2AV_TSS_PN_R1_f,
  Q_H2AV_TSS_PN_REP_R2_f=Q_H2AV_TSS_PN_REP_R2_f,
  Q_H2AV_TSS_PNH_L1_inputnorm_R1_f=Q_H2AV_TSS_PNH_L1_inputnorm_R1_f,
  Q_H2AV_TSS_PNH_L1_R1_f=Q_H2AV_TSS_PNH_L1_R1_f,
  Q_H2AV_TSS_PNH_REP_R2_f=Q_H2AV_TSS_PNH_REP_R2_f,
  Q_H2AV_TSS_PHYPH_inputnorm_R1_f=Q_H2AV_TSS_PHYPH_inputnorm_R1_f,
  Q_H2AV_TSS_PHYPH_R1_f=Q_H2AV_TSS_PHYPH_R1_f,
  Q_H2AV_TSS_PHYPBH_R2_f=Q_H2AV_TSS_PHYPBH_R2_f,
  Q_H2AV_TSS_PHYP_inputnorm_R1_f=Q_H2AV_TSS_PHYP_inputnorm_R1_f,
  Q_H2AV_TSS_PHYP_R1_f=Q_H2AV_TSS_PHYP_R1_f,
  Q_H2AV_TSS_PHYPB_R2_f=Q_H2AV_TSS_PHYPB_R2_f,
  Q_POLII_ctrl_N_f=Q_POLII_ctrl_N_f,
  Q_POLII_ctrl_f=Q_POLII_ctrl_f,
  Q_POLII_hypbKD_f=Q_POLII_hypbKD_f,
  Q_POLII_nelf_N_f=Q_POLII_nelf_N_f,
  Q_POLII_ser2P_ctrl_N_f=Q_POLII_ser2P_ctrl_N_f,
  Q_POLII_ser2P_ctrl_f=Q_POLII_ser2P_ctrl_f,
  Q_POLII_ser2P_hypbKD_f=Q_POLII_ser2P_hypbKD_f,
  Q_POLII_ser2P_nelf_N_f=Q_POLII_ser2P_nelf_N_f,
  Q_H2AV_GB_PW_R2_f=Q_H2AV_GB_PW_R2_f,
  Q_H2AV_GB_PW_R1_f=Q_H2AV_GB_PW_R1_f,
  Q_H2AV_GB_PWH_inputnorm_R1_f=Q_H2AV_GB_PWH_inputnorm_R1_f,
  Q_H2AV_GB_PWH_R1_f=Q_H2AV_GB_PWH_R1_f, 
  Q_H2AV_GB_PWH_R2_f=Q_H2AV_GB_PWH_REP2_f,
  Q_H2AV_GB_PN_inputnorm_R1_f=Q_H2AV_GB_PN_inputnorm_R1_f,
  Q_H2AV_GB_PN_R1_f=Q_H2AV_GB_PN_R1_f,
  Q_H2AV_GB_PNH_L1_inputnorm_R1_f=Q_H2AV_GB_PNH_L1_inputnorm_R1_f,
  Q_H2AV_GB_PNH_L1_R1_f=Q_H2AV_GB_PNH_L1_R1_f,
  Q_H2AV_GB_PN_R2_f=Q_H2AV_GB_PN_REP_R2_f,
  Q_H2AV_GB_PNH_R2_f=Q_H2AV_GB_PNH_REP_R2_f,
  Q_H2AV_GB_PHYPH_inputnorm_R1_f=Q_H2AV_GB_PHYPH_inputnorm_R1_f,
  Q_H2AV_GB_PHYPH_R1_f=Q_H2AV_GB_PHYPH_R1_f,
  Q_H2AV_GB_PHYP_inputnorm_R1_f=Q_H2AV_GB_PHYP_inputnorm_R1_f,
  Q_H2AV_GB_PHYP_R1_f=Q_H2AV_GB_PHYP_R1_f,
  Q_H2AV_GB_PHYPB_R2_f=Q_H2AV_GB_PHYPB_REP2_f,
  Q_H2AV_GB_PHYPBH_R2_f=Q_H2AV_GB_PHYPBH_REP2_f,
  ZSCORE_H2AV_R3_HYPB_WT = ZSCORE_H2AV_R3_HYPB_WT,
  ZSCORE_H2AV_NELF_WT = ZSCORE_H2AV_NELF_WT,
  ZSCORE_H2AV_HYPB_WT = ZSCORE_H2AV_HYPB_WT,
  ZSCORE_RAD51_WT3_GB_N_WT = ZSCORE_RAD51_WT3_GB_N_WT,
  ZSCORE_H2AV_WT3_GB_N_WT = ZSCORE_H2AV_WT3_GB_N_WT,
  ZSCORE_H2AV_SEA4_GB_N_WT = ZSCORE_H2AV_SEA4_GB_N_WT,
  ZSCORE_RAD51_GB_WT_N = ZSCORE_RAD51_GB_WT_N,
  ZSCORE_H2AV_GB_WT_N = ZSCORE_H2AV_GB_WT_N,
  ZSCORE_PROFMAT_H2AV_WT_vs_NELF = ZSCORE_PROFMAT_H2AV_WT_vs_NELF,
)

saveRDS(LIST_QUANTIF,"/home/cperrois/work/PROJET_H2AV_2021/DATA/LIST_FEATURES/LIST_QUANTIF.RDS")
LIST_QUANTIF=readRDS("/home/cperrois/work/PROJET_H2AV_2021/DATA/LIST_FEATURES/LIST_QUANTIF.RDS")


#########################################################################################################################
###############################  LIST OF GENES CLOSE TO THE ORIGINES OF REPLICATION   ###################################

GSE65692_S2origins_dm6 = readRDS(paste0(workdir, "PROJET_H2AV/DATA/r6.13/GSE65692_S2origins_dm6.RDS"))

myfindOverlaps = function(genes_GR, ref_GR){
  ol1 = findOverlaps(genes_GR, ref_GR, minoverlap=200)
  ol1_genes = unique(ol1@from)
  return(genes_GR[ol1_genes])
}

GENES_dm6_ORI = myfindOverlaps(gene_dm6_gr_f, GSE65692_S2origins_dm6)
GENES_dm6_ORI = names(GENES_dm6_ORI)

saveRDS(GENES_dm6_ORI, paste0(workdir, "PROJET_H2AV/DATA/LIST_FEATURES/LIST_GN_ORI_dm6.RDS"))


# #########################################################################################################################
# ###########################################   CHIP H2AV T + P     ##########################################################
# LIST H2AV_VEC
# LOAD DATA

Q_PHYP_B_L1_RPGC_GB = readRDS(paste0(Quantifdir, "PHYP_B_L1_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PHYP_B_L2_RPGC_GB = readRDS(paste0(Quantifdir, "PHYP_B_L2_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PHYPH_B_L1_RPGC_GB = readRDS(paste0(Quantifdir, "PHYPH_B_L1_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PHYPH_B_L2_RPGC_GB = readRDS(paste0(Quantifdir, "PHYPH_B_L2_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PNCH_RPGC_GB = readRDS(paste0(Quantifdir, "PNCH_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PNH_L1_RPGC_GB = readRDS(paste0(Quantifdir, "PNH_L1_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PNH_L2_RPGC_GB = readRDS(paste0(Quantifdir, "PNH_L2_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PNHOA_L1_RPGC_GB = readRDS(paste0(Quantifdir, "PNHOA_L1_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PNHOA_L2_RPGC_GB = readRDS(paste0(Quantifdir, "PNHOA_L2_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PNH_RPGC_GB = readRDS(paste0(Quantifdir, "PNH_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PN_RPGC_GB = readRDS(paste0(Quantifdir, "PN_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PWCH_RPGC_GB = readRDS(paste0(Quantifdir, "PWCH_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PWH_RPGC_GB = readRDS(paste0(Quantifdir, "PWH_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PW_L1_RPGC_GB = readRDS(paste0(Quantifdir, "PW_L1_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PW_L2_RPGC_GB = readRDS(paste0(Quantifdir, "PW_L2_RPGC_readsCounts_GB_SCALED.RDS"))
Q_TNHOA_L1_RPGC_GB = readRDS(paste0(Quantifdir, "TNHOA_L1_RPGC_readsCounts_GB_SCALED.RDS"))
Q_TNHOA_L2_RPGC_GB = readRDS(paste0(Quantifdir, "TNHOA_L2_RPGC_readsCounts_GB_SCALED.RDS"))
Q_TWH_RPGC_GB = readRDS(paste0(Quantifdir, "TWH_RPGC_readsCounts_GB_SCALED.RDS"))
Q_PHYP_B_L1_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PHYP_B_L1_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PHYP_B_L2_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PHYP_B_L2_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PHYPH_B_L1_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PHYPH_B_L1_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PHYPH_B_L2_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PHYPH_B_L2_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PNCH_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PNCH_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PNH_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PNH_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PNH_L1_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PNH_L1_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PNH_L2_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PNH_L2_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PNHOA_L1_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PNHOA_L1_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PNHOA_L2_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PNHOA_L2_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PN_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PN_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PWCH_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PWCH_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PWH_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PWH_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PW_L1_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PW_L1_filt_sort_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PW_L2_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "PW_L2_filt_sort_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_TNHOA_L1_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "TNHOA_L1_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_TNHOA_L2_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "TNHOA_L2_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_TWH_inputnormRPGC_GB = readRDS(paste0(Quantifdir, "TWH_inputnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PNH_MNASEnormRPGC_GB = readRDS(paste0(Quantifdir, "PNH_MNASEnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PWCH_MNASEnormRPGC_GB = readRDS(paste0(Quantifdir, "PWCH_MNASEnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PNCH_MNASEnormRPGC_GB = readRDS(paste0(Quantifdir, "PNCH_MNASEnormRPGC_readsCounts_GB_SCALED.RDS"))
Q_PN_MNASEnormRPGC_GB = readRDS(paste0(Quantifdir, "PN_MNASEnormRPGC_readsCounts_GB_SCALED.RDS"))

Q_PHYP_B_L1_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PHYP_B_L1_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PHYP_B_L2_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PHYP_B_L2_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PHYPH_B_L1_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PHYPH_B_L1_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PHYPH_B_L2_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PHYPH_B_L2_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PNCH_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PNCH_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PNH_L1_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PNH_L1_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PNH_L2_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PNH_L2_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PNH_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PNH_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PNHOA_L1_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PNHOA_L1_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PNHOA_L2_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PNHOA_L2_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PN_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PN_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PWCH_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PWCH_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PWH_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PWH_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PW_L1_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PW_L1_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_PW_L2_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "PW_L2_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_TNHOA_L1_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "TNHOA_L1_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_TNHOA_L2_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "TNHOA_L2_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))
Q_TWH_MNASEnormRPKM_GB = readRDS(paste0(Quantifdir, "TWH_MNASEnormRPKM_readsCounts_GB_SCALED.RDS"))

Q_PNHOA_L1_normBy_TNHOA_L1_RPKM_GB = readRDS(paste0(Quantifdir, "PNHOA_L1_normBy_TNHOA_L1_RPKM_readsCounts_GB_SCALED.RDS"))
Q_PN_normBy_TNHOA_L1_RPKM_GB = readRDS(paste0(Quantifdir, "PN_normBy_TNHOA_L1_RPKM_readsCounts_GB_SCALED.RDS"))
Q_PWH_normBy_TWH_RPKM_GB = readRDS(paste0(Quantifdir, "PWH_normBy_TWH_RPKM_readsCounts_GB_SCALED.RDS"))


# FILTER ACTIVE GENES

Q_PHYP_B_L1_RPGC_GB_f = Q_PHYP_B_L1_RPGC_GB[names(Q_PHYP_B_L1_RPGC_GB) %in% GNref]
Q_PHYP_B_L2_RPGC_GB_f = Q_PHYP_B_L2_RPGC_GB[names(Q_PHYP_B_L2_RPGC_GB) %in% GNref]
Q_PHYPH_B_L1_RPGC_GB_f = Q_PHYPH_B_L1_RPGC_GB[names(Q_PHYPH_B_L1_RPGC_GB) %in% GNref]
Q_PHYPH_B_L2_RPGC_GB_f = Q_PHYPH_B_L2_RPGC_GB[names(Q_PHYPH_B_L2_RPGC_GB) %in% GNref]
Q_PNCH_RPGC_GB_f = Q_PNCH_RPGC_GB[names(Q_PNCH_RPGC_GB) %in% GNref]
Q_PNH_L1_RPGC_GB_f = Q_PNH_L1_RPGC_GB[names(Q_PNH_L1_RPGC_GB) %in% GNref]
Q_PNH_L2_RPGC_GB_f = Q_PNH_L2_RPGC_GB[names(Q_PNH_L2_RPGC_GB) %in% GNref]
Q_PNHOA_L1_RPGC_GB_f = Q_PNHOA_L1_RPGC_GB[names(Q_PNHOA_L1_RPGC_GB) %in% GNref]
Q_PNHOA_L2_RPGC_GB_f = Q_PNHOA_L2_RPGC_GB[names(Q_PNHOA_L2_RPGC_GB) %in% GNref]
Q_PNH_RPGC_GB_f = Q_PNH_RPGC_GB[names(Q_PNH_RPGC_GB) %in% GNref]
Q_PN_RPGC_GB_f = Q_PN_RPGC_GB[names(Q_PN_RPGC_GB) %in% GNref]
Q_PWCH_RPGC_GB_f = Q_PWCH_RPGC_GB[names(Q_PWCH_RPGC_GB) %in% GNref]
Q_PWH_RPGC_GB_f = Q_PWH_RPGC_GB[names(Q_PWH_RPGC_GB) %in% GNref]
Q_PW_L1_RPGC_GB_f = Q_PW_L1_RPGC_GB[names(Q_PW_L1_RPGC_GB) %in% GNref]
Q_PW_L2_RPGC_GB_f = Q_PW_L2_RPGC_GB[names(Q_PW_L2_RPGC_GB) %in% GNref]
Q_TNHOA_L1_RPGC_GB_f = Q_TNHOA_L1_RPGC_GB[names(Q_TNHOA_L1_RPGC_GB) %in% GNref]
Q_TNHOA_L2_RPGC_GB_f = Q_TNHOA_L2_RPGC_GB[names(Q_TNHOA_L2_RPGC_GB) %in% GNref]
Q_TWH_RPGC_GB_f = Q_TWH_RPGC_GB[names(Q_TWH_RPGC_GB) %in% GNref]
Q_PHYP_B_L1_inputnormRPGC_GB_f = Q_PHYP_B_L1_inputnormRPGC_GB[names(Q_PHYP_B_L1_inputnormRPGC_GB) %in% GNref]
Q_PHYP_B_L2_inputnormRPGC_GB_f = Q_PHYP_B_L2_inputnormRPGC_GB[names(Q_PHYP_B_L2_inputnormRPGC_GB) %in% GNref]
Q_PHYPH_B_L1_inputnormRPGC_GB_f = Q_PHYPH_B_L1_inputnormRPGC_GB[names(Q_PHYPH_B_L1_inputnormRPGC_GB) %in% GNref]
Q_PHYPH_B_L2_inputnormRPGC_GB_f = Q_PHYPH_B_L2_inputnormRPGC_GB[names(Q_PHYPH_B_L2_inputnormRPGC_GB) %in% GNref]
Q_PNCH_inputnormRPGC_GB_f = Q_PNCH_inputnormRPGC_GB[names(Q_PNCH_inputnormRPGC_GB) %in% GNref]
Q_PNH_inputnormRPGC_GB_f = Q_PNH_inputnormRPGC_GB[names(Q_PNH_inputnormRPGC_GB) %in% GNref]
Q_PNH_L1_inputnormRPGC_GB_f = Q_PNH_L1_inputnormRPGC_GB[names(Q_PNH_L1_inputnormRPGC_GB) %in% GNref]
Q_PNH_L2_inputnormRPGC_GB_f = Q_PNH_L2_inputnormRPGC_GB[names(Q_PNH_L2_inputnormRPGC_GB) %in% GNref]
Q_PNHOA_L1_inputnormRPGC_GB_f = Q_PNHOA_L1_inputnormRPGC_GB[names(Q_PNHOA_L1_inputnormRPGC_GB) %in% GNref]
Q_PNHOA_L2_inputnormRPGC_GB_f = Q_PNHOA_L2_inputnormRPGC_GB[names(Q_PNHOA_L2_inputnormRPGC_GB) %in% GNref]
Q_PN_inputnormRPGC_GB_f = Q_PN_inputnormRPGC_GB[names(Q_PN_inputnormRPGC_GB) %in% GNref]
Q_PWCH_inputnormRPGC_GB_f = Q_PWCH_inputnormRPGC_GB[names(Q_PWCH_inputnormRPGC_GB) %in% GNref]
Q_PWH_inputnormRPGC_GB_f = Q_PWH_inputnormRPGC_GB[names(Q_PWH_inputnormRPGC_GB) %in% GNref]
Q_PW_L1_inputnormRPGC_GB_f = Q_PW_L1_inputnormRPGC_GB[names(Q_PW_L1_inputnormRPGC_GB) %in% GNref]
Q_PW_L2_inputnormRPGC_GB_f = Q_PW_L2_inputnormRPGC_GB[names(Q_PW_L2_inputnormRPGC_GB) %in% GNref]
Q_TNHOA_L1_inputnormRPGC_GB_f = Q_TNHOA_L1_inputnormRPGC_GB[names(Q_TNHOA_L1_inputnormRPGC_GB) %in% GNref]
Q_TNHOA_L2_inputnormRPGC_GB_f = Q_TNHOA_L2_inputnormRPGC_GB[names(Q_TNHOA_L2_inputnormRPGC_GB) %in% GNref]
Q_TWH_inputnormRPGC_GB_f = Q_TWH_inputnormRPGC_GB[names(Q_TWH_inputnormRPGC_GB) %in% GNref]
Q_PNH_MNASEnormRPGC_GB_f = Q_PNH_MNASEnormRPGC_GB[names(Q_PNH_MNASEnormRPGC_GB) %in% GNref]
Q_PWCH_MNASEnormRPGC_GB_f = Q_PWCH_MNASEnormRPGC_GB[names(Q_PWCH_MNASEnormRPGC_GB) %in% GNref]
Q_PNCH_MNASEnormRPGC_GB_f = Q_PNCH_MNASEnormRPGC_GB[names(Q_PNCH_MNASEnormRPGC_GB) %in% GNref]
Q_PN_MNASEnormRPGC_GB_f = Q_PN_MNASEnormRPGC_GB[names(Q_PN_MNASEnormRPGC_GB) %in% GNref]

Q_PHYP_B_L1_MNASEnormRPKM_GB_f = Q_PHYP_B_L1_MNASEnormRPKM_GB[names(Q_PHYP_B_L1_MNASEnormRPKM_GB) %in% GNref]
Q_PHYP_B_L2_MNASEnormRPKM_GB_f = Q_PHYP_B_L2_MNASEnormRPKM_GB[names(Q_PHYP_B_L2_MNASEnormRPKM_GB) %in% GNref]
Q_PHYPH_B_L1_MNASEnormRPKM_GB_f = Q_PHYPH_B_L1_MNASEnormRPKM_GB[names(Q_PHYPH_B_L1_MNASEnormRPKM_GB) %in% GNref]
Q_PHYPH_B_L2_MNASEnormRPKM_GB_f = Q_PHYPH_B_L2_MNASEnormRPKM_GB[names(Q_PHYPH_B_L2_MNASEnormRPKM_GB) %in% GNref]
Q_PNCH_MNASEnormRPKM_GB_f = Q_PNCH_MNASEnormRPKM_GB[names(Q_PNCH_MNASEnormRPKM_GB) %in% GNref]
Q_PNH_L1_MNASEnormRPKM_GB_f = Q_PNH_L1_MNASEnormRPKM_GB[names(Q_PNH_L1_MNASEnormRPKM_GB) %in% GNref]
Q_PNH_L2_MNASEnormRPKM_GB_f = Q_PNH_L2_MNASEnormRPKM_GB[names(Q_PNH_L2_MNASEnormRPKM_GB) %in% GNref]
Q_PNH_MNASEnormRPKM_GB_f = Q_PNH_MNASEnormRPKM_GB[names(Q_PNH_MNASEnormRPKM_GB) %in% GNref]
Q_PNHOA_L1_MNASEnormRPKM_GB_f = Q_PNHOA_L1_MNASEnormRPKM_GB[names(Q_PNHOA_L1_MNASEnormRPKM_GB) %in% GNref]
Q_PNHOA_L2_MNASEnormRPKM_GB_f = Q_PNHOA_L2_MNASEnormRPKM_GB[names(Q_PNHOA_L2_MNASEnormRPKM_GB) %in% GNref]
Q_PN_MNASEnormRPKM_GB_f = Q_PN_MNASEnormRPKM_GB[names(Q_PN_MNASEnormRPKM_GB) %in% GNref]
Q_PWCH_MNASEnormRPKM_GB_f = Q_PWCH_MNASEnormRPKM_GB[names(Q_PWCH_MNASEnormRPKM_GB) %in% GNref]
Q_PWH_MNASEnormRPKM_GB_f = Q_PWH_MNASEnormRPKM_GB[names(Q_PWH_MNASEnormRPKM_GB) %in% GNref]
Q_PW_L1_MNASEnormRPKM_GB_f = Q_PW_L1_MNASEnormRPKM_GB[names(Q_PW_L1_MNASEnormRPKM_GB) %in% GNref]
Q_PW_L2_MNASEnormRPKM_GB_f = Q_PW_L2_MNASEnormRPKM_GB[names(Q_PW_L2_MNASEnormRPKM_GB) %in% GNref]
Q_TNHOA_L1_MNASEnormRPKM_GB_f = Q_TNHOA_L1_MNASEnormRPKM_GB[names(Q_TNHOA_L1_MNASEnormRPKM_GB) %in% GNref]
Q_TNHOA_L2_MNASEnormRPKM_GB_f = Q_TNHOA_L2_MNASEnormRPKM_GB[names(Q_TNHOA_L2_MNASEnormRPKM_GB) %in% GNref]
Q_TWH_MNASEnormRPKM_GB_f = Q_TWH_MNASEnormRPKM_GB[names(Q_TWH_MNASEnormRPKM_GB) %in% GNref]


ZSCORE_PN_PW_L1_f = computeZscore(Q_PN_RPGC_GB_f, Q_PW_L1_RPGC_GB_f)
ZSCORE_PN_PW_L2_f = computeZscore(Q_PN_RPGC_GB_f, Q_PW_L2_RPGC_GB_f)
ZSCORE_PNH_L1_PWH_f = computeZscore(Q_PNH_L1_RPGC_GB_f, Q_PWH_RPGC_GB_f)
ZSCORE_PNH_L2_PWH_f = computeZscore(Q_PNH_L2_RPGC_GB_f, Q_PWH_RPGC_GB_f)
ZSCORE_PNH_L1_PWCH_f = computeZscore(Q_PNH_L1_RPGC_GB_f, Q_PWCH_RPGC_GB_f)
ZSCORE_PNCH_PWCH_f = computeZscore(Q_PNCH_RPGC_GB_f, Q_PWCH_RPGC_GB_f)
ZSCORE_PNHOA_L1_PWH_f = computeZscore(Q_PNHOA_L1_RPGC_GB_f, Q_PWH_RPGC_GB_f)
ZSCORE_PNHOA_L2_PWH_f = computeZscore(Q_PNHOA_L2_RPGC_GB_f, Q_PWH_RPGC_GB_f)
ZSCORE_PNHOA_L1_TNHOA_L1_f = computeZscore(Q_PNHOA_L1_RPGC_GB_f, Q_TNHOA_L1_RPGC_GB_f)
ZSCORE_PNHOA_L2_TNHOA_L2_f = computeZscore(Q_PNHOA_L2_RPGC_GB_f, Q_TNHOA_L2_RPGC_GB_f)
ZSCORE_TNHOA_L1_TWH_f = computeZscore(Q_TNHOA_L1_RPGC_GB_f, Q_TWH_RPGC_GB_f)
ZSCORE_TNHOA_L2_TWH_f = computeZscore(Q_TNHOA_L2_RPGC_GB_f, Q_TWH_RPGC_GB_f)
ZSCORE_PHYP_B_L1_PW_L1_f = computeZscore(Q_PHYP_B_L1_RPGC_GB_f, Q_PW_L1_RPGC_GB_f)
ZSCORE_PHYP_B_L2_PW_L2_f = computeZscore(Q_PHYP_B_L2_RPGC_GB_f, Q_PW_L2_RPGC_GB_f)
ZSCORE_PHYPH_B_L1_PWH_f = computeZscore(Q_PHYPH_B_L1_RPGC_GB_f, Q_PWH_RPGC_GB_f)
ZSCORE_PHYPH_B_L2_PWH_f = computeZscore(Q_PHYPH_B_L2_RPGC_GB_f, Q_PWH_RPGC_GB_f)
ZSCORE_PHYPH_B_L1_PWCH_f = computeZscore(Q_PHYPH_B_L1_RPGC_GB_f, Q_PWCH_RPGC_GB_f)
ZSCORE_PHYPH_B_L2_PWCH_f = computeZscore(Q_PHYPH_B_L2_RPGC_GB_f, Q_PWCH_RPGC_GB_f)

ZSCORE_PN_PWH_f = computeZscore(Q_PN_RPGC_GB_f, Q_PWH_RPGC_GB_f)

## INPUT NORM
ZSCORE_PN_PW_L1_inputnorm_f = computeZscore(Q_PN_inputnormRPGC_GB_f, Q_PW_L1_inputnormRPGC_GB_f)
ZSCORE_PN_PW_L2_inputnorm_f = computeZscore(Q_PN_inputnormRPGC_GB_f, Q_PW_L2_inputnormRPGC_GB_f)
ZSCORE_PNH_L1_PWH_inputnorm_f = computeZscore(Q_PNH_L1_inputnormRPGC_GB_f, Q_PWH_inputnormRPGC_GB_f)
ZSCORE_PNH_L2_PWH_inputnorm_f = computeZscore(Q_PNH_L2_inputnormRPGC_GB_f, Q_PWH_inputnormRPGC_GB_f)
ZSCORE_PNH_L1_PWCH_inputnorm_f = computeZscore(Q_PNH_L1_inputnormRPGC_GB_f, Q_PWCH_inputnormRPGC_GB_f)
ZSCORE_PNCH_PWCH_inputnorm_f = computeZscore(Q_PNCH_inputnormRPGC_GB_f, Q_PWCH_inputnormRPGC_GB_f)
ZSCORE_PNHOA_L1_PWH_inputnorm_f = computeZscore(Q_PNHOA_L1_inputnormRPGC_GB_f, Q_PWH_inputnormRPGC_GB_f)
ZSCORE_PNHOA_L2_PWH_inputnorm_f = computeZscore(Q_PNHOA_L2_inputnormRPGC_GB_f, Q_PWH_inputnormRPGC_GB_f)
ZSCORE_PNHOA_L1_TNHOA_L1_inputnorm_f = computeZscore(Q_PNHOA_L1_inputnormRPGC_GB_f, Q_TNHOA_L1_inputnormRPGC_GB_f)
ZSCORE_PNHOA_L2_TNHOA_L2_inputnorm_f = computeZscore(Q_PNHOA_L2_inputnormRPGC_GB_f, Q_TNHOA_L2_inputnormRPGC_GB_f)
ZSCORE_TNHOA_L1_TWH_inputnorm_f = computeZscore(Q_TNHOA_L1_inputnormRPGC_GB_f, Q_TWH_inputnormRPGC_GB_f)
ZSCORE_TNHOA_L2_TWH_inputnorm_f = computeZscore(Q_TNHOA_L2_inputnormRPGC_GB_f, Q_TWH_inputnormRPGC_GB_f)
ZSCORE_PHYP_B_L1_PW_L1_inputnorm_f = computeZscore(Q_PHYP_B_L1_inputnormRPGC_GB_f, Q_PW_L1_inputnormRPGC_GB_f)
ZSCORE_PHYP_B_L2_PW_L2_inputnorm_f = computeZscore(Q_PHYP_B_L2_inputnormRPGC_GB_f, Q_PW_L2_inputnormRPGC_GB_f)
ZSCORE_PHYPH_B_L1_PWH_inputnorm_f = computeZscore(Q_PHYPH_B_L1_inputnormRPGC_GB_f, Q_PWH_inputnormRPGC_GB_f)
ZSCORE_PHYPH_B_L2_PWH_inputnorm_f = computeZscore(Q_PHYPH_B_L2_inputnormRPGC_GB_f, Q_PWH_inputnormRPGC_GB_f)
ZSCORE_PHYPH_B_L1_PWCH_inputnorm_f = computeZscore(Q_PHYPH_B_L1_inputnormRPGC_GB_f, Q_PWCH_inputnormRPGC_GB_f)
ZSCORE_PHYPH_B_L2_PWCH_inputnorm_f = computeZscore(Q_PHYPH_B_L2_inputnormRPGC_GB_f, Q_PWCH_inputnormRPGC_GB_f)
ZSCORE_PN_PWH_inputnorm_f = computeZscore(Q_PN_inputnormRPGC_GB_f, Q_PWH_inputnormRPGC_GB_f)



ZSCORE_PHYPH_B_L1_PW_L1_MNASEnormRPKM_f = computeZscore(Q_PHYP_B_L1_MNASEnormRPKM_GB_f, Q_PW_L1_MNASEnormRPKM_GB_f)
ZSCORE_PHYPH_B_L2_PW_L2_MNASEnormRPKM_f = computeZscore(Q_PHYP_B_L2_MNASEnormRPKM_GB_f, Q_PW_L2_MNASEnormRPKM_GB_f)
ZSCORE_PN_PW_L1_MNASEnormRPKM_f = computeZscore(Q_PN_MNASEnormRPKM_GB_f, Q_PW_L1_MNASEnormRPKM_GB_f)
ZSCORE_PN_PW_L2_MNASEnormRPKM_f = computeZscore(Q_PN_MNASEnormRPKM_GB_f, Q_PW_L2_MNASEnormRPKM_GB_f)
ZSCORE_PHYPH_B_L1_PWH_MNASEnormRPKM_f = computeZscore(Q_PHYPH_B_L1_MNASEnormRPKM_GB_f, Q_PWH_MNASEnormRPKM_GB_f)
ZSCORE_PHYPH_B_L2_PWH_MNASEnormRPKM_f = computeZscore(Q_PHYPH_B_L2_MNASEnormRPKM_GB_f, Q_PWH_MNASEnormRPKM_GB_f)
ZSCORE_PNH_L1_PWH_MNASEnormRPKM_f = computeZscore(Q_PNH_L1_MNASEnormRPKM_GB_f, Q_PWH_MNASEnormRPKM_GB_f)
ZSCORE_PNH_L2_PWH_MNASEnormRPKM_f = computeZscore(Q_PNH_L2_MNASEnormRPKM_GB_f, Q_PWH_MNASEnormRPKM_GB_f)
ZSCORE_PNCH_PWCH_MNASEnormRPKM_f = computeZscore(Q_PNCH_MNASEnormRPKM_GB_f, Q_PWCH_MNASEnormRPKM_GB_f)
ZSCORE_PNHOA_L1_PWH_MNASEnormRPKM_f = computeZscore(Q_PNHOA_L1_MNASEnormRPKM_GB_f, Q_PWH_MNASEnormRPKM_GB_f)
ZSCORE_PNHOA_L2_PWH_MNASEnormRPKM_f = computeZscore(Q_PNHOA_L2_MNASEnormRPKM_GB_f, Q_PWH_MNASEnormRPKM_GB_f)
ZSCORE_PNHOA_L1_TNHOA_L1_MNASEnormRPKM_f = computeZscore(Q_PNHOA_L1_MNASEnormRPKM_GB_f, Q_TNHOA_L1_MNASEnormRPKM_GB_f)
ZSCORE_PNHOA_L2_TNHOA_L2_MNASEnormRPKM_f = computeZscore(Q_PNHOA_L2_MNASEnormRPKM_GB_f, Q_TNHOA_L2_MNASEnormRPKM_GB_f)
ZSCORE_TNHOA_L1_TWH_MNASEnormRPKM_f = computeZscore(Q_TNHOA_L1_MNASEnormRPKM_GB_f, Q_TWH_MNASEnormRPKM_GB_f)
ZSCORE_TNHOA_L2_TWH_MNASEnormRPKM_f = computeZscore(Q_TNHOA_L2_MNASEnormRPKM_GB_f, Q_TWH_MNASEnormRPKM_GB_f)

# RATIO GAMMA / TOTAL
LOG2RATIO_PWH_TWH_f = log2(Q_PWH_RPGC_GB_f+1)-log2(Q_TWH_RPGC_GB_f[names(Q_PWH_RPGC_GB_f)]+1)
LOG2RATIO_PNHOA_TNHOA_f = log2(Q_PNHOA_L1_RPGC_GB_f+1)-log2(Q_TNHOA_L1_RPGC_GB_f[names(Q_PNHOA_L1_RPGC_GB_f)]+1)

LOG2RATIO_PWH_TWH_MNASEnormRPKM_f = log2(Q_PWH_MNASEnormRPKM_GB_f+1)-log2(Q_TWH_MNASEnormRPKM_GB_f[names(Q_PWH_MNASEnormRPKM_GB_f)]+1)
LOG2RATIO_PNHOA_TNHOA_MNASEnormRPKM_f = log2(Q_PNHOA_L1_MNASEnormRPKM_GB_f+1)-log2(Q_TNHOA_L1_MNASEnormRPKM_GB_f[names(Q_PNHOA_L1_MNASEnormRPKM_GB_f)]+1)


### ORDER deacreasingly
Q_PHYP_B_L1_RPGC_GB_f = Q_PHYP_B_L1_RPGC_GB_f[order(Q_PHYP_B_L1_RPGC_GB_f, decreasing=T)]
Q_PHYP_B_L2_RPGC_GB_f = Q_PHYP_B_L2_RPGC_GB_f[order(Q_PHYP_B_L2_RPGC_GB_f, decreasing=T)]
Q_PHYPH_B_L1_RPGC_GB_f = Q_PHYPH_B_L1_RPGC_GB_f[order(Q_PHYPH_B_L1_RPGC_GB_f, decreasing=T)]
Q_PHYPH_B_L2_RPGC_GB_f = Q_PHYPH_B_L2_RPGC_GB_f[order(Q_PHYPH_B_L2_RPGC_GB_f, decreasing=T)]
Q_PNCH_RPGC_GB_f = Q_PNCH_RPGC_GB_f[order(Q_PNCH_RPGC_GB_f, decreasing=T)]
Q_PNH_L1_RPGC_GB_f = Q_PNH_L1_RPGC_GB_f[order(Q_PNH_L1_RPGC_GB_f, decreasing=T)]
Q_PNH_L2_RPGC_GB_f = Q_PNH_L2_RPGC_GB_f[order(Q_PNH_L2_RPGC_GB_f, decreasing=T)]
Q_PNHOA_L1_RPGC_GB_f = Q_PNHOA_L1_RPGC_GB_f[order(Q_PNHOA_L1_RPGC_GB_f, decreasing=T)]
Q_PNHOA_L2_RPGC_GB_f = Q_PNHOA_L2_RPGC_GB_f[order(Q_PNHOA_L2_RPGC_GB_f, decreasing=T)]
Q_PNH_RPGC_GB_f = Q_PNH_RPGC_GB_f[order(Q_PNH_RPGC_GB_f, decreasing=T)]
Q_PN_RPGC_GB_f = Q_PN_RPGC_GB_f[order(Q_PN_RPGC_GB_f, decreasing=T)]
Q_PWCH_RPGC_GB_f = Q_PWCH_RPGC_GB_f[order(Q_PWCH_RPGC_GB_f, decreasing=T)]
Q_PWH_RPGC_GB_f = Q_PWH_RPGC_GB_f[order(Q_PWH_RPGC_GB_f, decreasing=T)]
Q_PW_L1_RPGC_GB_f = Q_PW_L1_RPGC_GB_f[order(Q_PW_L1_RPGC_GB_f, decreasing=T)]
Q_PW_L2_RPGC_GB_f = Q_PW_L2_RPGC_GB_f[order(Q_PW_L2_RPGC_GB_f, decreasing=T)]
Q_TNHOA_L1_RPGC_GB_f = Q_TNHOA_L1_RPGC_GB_f[order(Q_TNHOA_L1_RPGC_GB_f, decreasing=T)]
Q_TNHOA_L2_RPGC_GB_f = Q_TNHOA_L2_RPGC_GB_f[order(Q_TNHOA_L2_RPGC_GB_f, decreasing=T)]
Q_TWH_RPGC_GB_f = Q_TWH_RPGC_GB_f[order(Q_TWH_RPGC_GB_f, decreasing=T)]

LOG2RATIO_PWH_TWH_f = LOG2RATIO_PWH_TWH_f[order(LOG2RATIO_PWH_TWH_f, decreasing=T)]
LOG2RATIO_PNHOA_TNHOA_f = LOG2RATIO_PNHOA_TNHOA_f[order(LOG2RATIO_PNHOA_TNHOA_f, decreasing=T)]

Q_PHYP_B_L1_MNASEnormRPKM_GB_f = Q_PHYP_B_L1_MNASEnormRPKM_GB_f[order(Q_PHYP_B_L1_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PHYP_B_L2_MNASEnormRPKM_GB_f = Q_PHYP_B_L2_MNASEnormRPKM_GB_f[order(Q_PHYP_B_L2_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PHYPH_B_L1_MNASEnormRPKM_GB_f = Q_PHYPH_B_L1_MNASEnormRPKM_GB_f[order(Q_PHYPH_B_L1_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PHYPH_B_L2_MNASEnormRPKM_GB_f = Q_PHYPH_B_L2_MNASEnormRPKM_GB_f[order(Q_PHYPH_B_L2_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PNCH_MNASEnormRPKM_GB_f = Q_PNCH_MNASEnormRPKM_GB_f[order(Q_PNCH_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PNH_L1_MNASEnormRPKM_GB_f = Q_PNH_L1_MNASEnormRPKM_GB_f[order(Q_PNH_L1_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PNH_L2_MNASEnormRPKM_GB_f = Q_PNH_L2_MNASEnormRPKM_GB_f[order(Q_PNH_L2_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PNH_MNASEnormRPKM_GB_f = Q_PNH_MNASEnormRPKM_GB_f[order(Q_PNH_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PNHOA_L1_MNASEnormRPKM_GB_f = Q_PNHOA_L1_MNASEnormRPKM_GB_f[order(Q_PNHOA_L1_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PNHOA_L2_MNASEnormRPKM_GB_f = Q_PNHOA_L2_MNASEnormRPKM_GB_f[order(Q_PNHOA_L2_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PN_MNASEnormRPKM_GB_f = Q_PN_MNASEnormRPKM_GB_f[order(Q_PN_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PWCH_MNASEnormRPKM_GB_f = Q_PWCH_MNASEnormRPKM_GB_f[order(Q_PWCH_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PWH_MNASEnormRPKM_GB_f = Q_PWH_MNASEnormRPKM_GB_f[order(Q_PWH_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PW_L1_MNASEnormRPKM_GB_f = Q_PW_L1_MNASEnormRPKM_GB_f[order(Q_PW_L1_MNASEnormRPKM_GB_f, decreasing=T)]
Q_PW_L2_MNASEnormRPKM_GB_f = Q_PW_L2_MNASEnormRPKM_GB_f[order(Q_PW_L2_MNASEnormRPKM_GB_f, decreasing=T)]
Q_TNHOA_L1_MNASEnormRPKM_GB_f = Q_TNHOA_L1_MNASEnormRPKM_GB_f[order(Q_TNHOA_L1_MNASEnormRPKM_GB_f, decreasing=T)]
Q_TNHOA_L2_MNASEnormRPKM_GB_f = Q_TNHOA_L2_MNASEnormRPKM_GB_f[order(Q_TNHOA_L2_MNASEnormRPKM_GB_f, decreasing=T)]
Q_TWH_MNASEnormRPKM_GB_f = Q_TWH_MNASEnormRPKM_GB_f[order(Q_TWH_MNASEnormRPKM_GB_f, decreasing=T)]

LOG2RATIO_PWH_TWH_MNASEnormRPKM_f = LOG2RATIO_PWH_TWH_MNASEnormRPKM_f[order(LOG2RATIO_PWH_TWH_MNASEnormRPKM_f, decreasing=T)]
LOG2RATIO_PNHOA_TNHOA_MNASEnormRPKM_f = LOG2RATIO_PNHOA_TNHOA_MNASEnormRPKM_f[order(LOG2RATIO_PNHOA_TNHOA_MNASEnormRPKM_f, decreasing=T)]

Q_PHYP_B_L1_inputnormRPGC_GB_f = Q_PHYP_B_L1_inputnormRPGC_GB_f[order(Q_PHYP_B_L1_inputnormRPGC_GB_f, decreasing=T)]
Q_PHYP_B_L2_inputnormRPGC_GB_f = Q_PHYP_B_L2_inputnormRPGC_GB_f[order(Q_PHYP_B_L2_inputnormRPGC_GB_f, decreasing=T)]
Q_PHYPH_B_L1_inputnormRPGC_GB_f = Q_PHYPH_B_L1_inputnormRPGC_GB_f[order(Q_PHYPH_B_L1_inputnormRPGC_GB_f, decreasing=T)]
Q_PHYPH_B_L2_inputnormRPGC_GB_f = Q_PHYPH_B_L2_inputnormRPGC_GB_f[order(Q_PHYPH_B_L2_inputnormRPGC_GB_f, decreasing=T)]
Q_PNCH_inputnormRPGC_GB_f = Q_PNCH_inputnormRPGC_GB_f[order(Q_PNCH_inputnormRPGC_GB_f, decreasing=T)]
Q_PNH_inputnormRPGC_GB_f = Q_PNH_inputnormRPGC_GB_f[order(Q_PNH_inputnormRPGC_GB_f, decreasing=T)]
Q_PNH_L1_inputnormRPGC_GB_f = Q_PNH_L1_inputnormRPGC_GB_f[order(Q_PNH_L1_inputnormRPGC_GB_f, decreasing=T)]
Q_PNH_L2_inputnormRPGC_GB_f = Q_PNH_L2_inputnormRPGC_GB_f[order(Q_PNH_L2_inputnormRPGC_GB_f, decreasing=T)]
Q_PNHOA_L1_inputnormRPGC_GB_f = Q_PNHOA_L1_inputnormRPGC_GB_f[order(Q_PNHOA_L1_inputnormRPGC_GB_f, decreasing=T)]
Q_PNHOA_L2_inputnormRPGC_GB_f = Q_PNHOA_L2_inputnormRPGC_GB_f[order(Q_PNHOA_L2_inputnormRPGC_GB_f, decreasing=T)]
Q_PN_inputnormRPGC_GB_f = Q_PN_inputnormRPGC_GB_f[order(Q_PN_inputnormRPGC_GB_f, decreasing=T)]
Q_PWCH_inputnormRPGC_GB_f = Q_PWCH_inputnormRPGC_GB_f[order(Q_PWCH_inputnormRPGC_GB_f, decreasing=T)]
Q_PWH_inputnormRPGC_GB_f = Q_PWH_inputnormRPGC_GB_f[order(Q_PWH_inputnormRPGC_GB_f, decreasing=T)]
Q_PW_L1_inputnormRPGC_GB_f = Q_PW_L1_inputnormRPGC_GB_f[order(Q_PW_L1_inputnormRPGC_GB_f, decreasing=T)]
Q_PW_L2_inputnormRPGC_GB_f = Q_PW_L2_inputnormRPGC_GB_f[order(Q_PW_L2_inputnormRPGC_GB_f, decreasing=T)]
Q_TNHOA_L1_inputnormRPGC_GB_f = Q_TNHOA_L1_inputnormRPGC_GB_f[order(Q_TNHOA_L1_inputnormRPGC_GB_f, decreasing=T)]
Q_TNHOA_L2_inputnormRPGC_GB_f = Q_TNHOA_L2_inputnormRPGC_GB_f[order(Q_TNHOA_L2_inputnormRPGC_GB_f, decreasing=T)]
Q_TWH_inputnormRPGC_GB_f = Q_TWH_inputnormRPGC_GB_f[order(Q_TWH_inputnormRPGC_GB_f, decreasing=T)]



# CREATE LIST

LIST_H2AV_VEC = list(Q_PHYP_B_L1_RPGC_GB_f = Q_PHYP_B_L1_RPGC_GB_f,
                     Q_PHYP_B_L2_RPGC_GB_f = Q_PHYP_B_L2_RPGC_GB_f,
                     Q_PHYPH_B_L1_RPGC_GB_f = Q_PHYPH_B_L1_RPGC_GB_f,
                     Q_PHYPH_B_L2_RPGC_GB_f = Q_PHYPH_B_L2_RPGC_GB_f,
                     Q_PNCH_RPGC_GB_f = Q_PNCH_RPGC_GB_f,
                     Q_PNH_L1_RPGC_GB_f = Q_PNH_L1_RPGC_GB_f,
                     Q_PNH_L2_RPGC_GB_f = Q_PNH_L2_RPGC_GB_f,
                     Q_PNHOA_L1_RPGC_GB_f = Q_PNHOA_L1_RPGC_GB_f,
                     Q_PNHOA_L2_RPGC_GB_f = Q_PNHOA_L2_RPGC_GB_f,
                     Q_PNH_RPGC_GB_f = Q_PNH_RPGC_GB_f,
                     Q_PN_RPGC_GB_f = Q_PN_RPGC_GB_f,
                     Q_PWCH_RPGC_GB_f = Q_PWCH_RPGC_GB_f,
                     Q_PWH_RPGC_GB_f = Q_PWH_RPGC_GB_f,
                     Q_PW_L1_RPGC_GB_f = Q_PW_L1_RPGC_GB_f,
                     Q_PW_L2_RPGC_GB_f = Q_PW_L2_RPGC_GB_f,
                     Q_TNHOA_L1_RPGC_GB_f = Q_TNHOA_L1_RPGC_GB_f,
                     Q_TNHOA_L2_RPGC_GB_f = Q_TNHOA_L2_RPGC_GB_f,
                     Q_TWH_RPGC_GB_f = Q_TWH_RPGC_GB_f,
                     Q_PHYP_B_L1_inputnormRPGC_GB_f = Q_PHYP_B_L1_inputnormRPGC_GB_f,
                     Q_PHYP_B_L2_inputnormRPGC_GB_f = Q_PHYP_B_L2_inputnormRPGC_GB_f,
                     Q_PHYPH_B_L1_inputnormRPGC_GB_f = Q_PHYPH_B_L1_inputnormRPGC_GB_f,
                     Q_PHYPH_B_L2_inputnormRPGC_GB_f = Q_PHYPH_B_L2_inputnormRPGC_GB_f,
                     Q_PNCH_inputnormRPGC_GB_f = Q_PNCH_inputnormRPGC_GB_f,
                     Q_PNH_inputnormRPGC_GB_f = Q_PNH_inputnormRPGC_GB_f,
                     Q_PNH_L1_inputnormRPGC_GB_f = Q_PNH_L1_inputnormRPGC_GB_f,
                     Q_PNH_L2_inputnormRPGC_GB_f = Q_PNH_L2_inputnormRPGC_GB_f,
                     Q_PNHOA_L1_inputnormRPGC_GB_f = Q_PNHOA_L1_inputnormRPGC_GB_f,
                     Q_PNHOA_L2_inputnormRPGC_GB_f = Q_PNHOA_L2_inputnormRPGC_GB_f,
                     Q_PN_inputnormRPGC_GB_f = Q_PN_inputnormRPGC_GB_f,
                     Q_PWCH_inputnormRPGC_GB_f = Q_PWCH_inputnormRPGC_GB_f,
                     Q_PWH_inputnormRPGC_GB_f = Q_PWH_inputnormRPGC_GB_f,
                     Q_PW_L1_inputnormRPGC_GB_f = Q_PW_L1_inputnormRPGC_GB_f,
                     Q_PW_L2_inputnormRPGC_GB_f = Q_PW_L2_inputnormRPGC_GB_f,
                     Q_TNHOA_L1_inputnormRPGC_GB_f = Q_TNHOA_L1_inputnormRPGC_GB_f,
                     Q_TNHOA_L2_inputnormRPGC_GB_f = Q_TNHOA_L2_inputnormRPGC_GB_f,
                     Q_TWH_inputnormRPGC_GB_f = Q_TWH_inputnormRPGC_GB_f,
                     ZSCORE_PN_PW_L1_f = ZSCORE_PN_PW_L1_f,
                     ZSCORE_PN_PW_L2_f = ZSCORE_PN_PW_L2_f,
                     ZSCORE_PHYP_B_L1_PW_L1_f = ZSCORE_PHYP_B_L1_PW_L1_f,
                     ZSCORE_PHYP_B_L2_PW_L2_f = ZSCORE_PHYP_B_L2_PW_L2_f,
                     ZSCORE_PHYPH_B_L1_PWH_f = ZSCORE_PHYPH_B_L1_PWH_f,
                     ZSCORE_PHYPH_B_L2_PWH_f = ZSCORE_PHYPH_B_L2_PWH_f,
                     ZSCORE_PHYPH_B_L1_PWCH_f = ZSCORE_PHYPH_B_L1_PWCH_f,
                     ZSCORE_PHYPH_B_L2_PWCH_f = ZSCORE_PHYPH_B_L2_PWCH_f,
                     ZSCORE_PNH_L1_PWH_f = ZSCORE_PNH_L1_PWH_f,
                     ZSCORE_PNH_L2_PWH_f = ZSCORE_PNH_L2_PWH_f,
                     ZSCORE_PNH_L1_PWCH_f = ZSCORE_PNH_L1_PWCH_f,
                     ZSCORE_PNCH_PWCH_f = ZSCORE_PNCH_PWCH_f,
                     ZSCORE_PNHOA_L1_PWH_f = ZSCORE_PNHOA_L1_PWH_f,
                     ZSCORE_PNHOA_L2_PWH_f = ZSCORE_PNHOA_L2_PWH_f,
                     ZSCORE_PNHOA_L1_TNHOA_L1_f = ZSCORE_PNHOA_L1_TNHOA_L1_f,
                     ZSCORE_PNHOA_L2_TNHOA_L2_f = ZSCORE_PNHOA_L2_TNHOA_L2_f,
                     ZSCORE_TNHOA_L1_TWH_f = ZSCORE_TNHOA_L1_TWH_f,
                     ZSCORE_TNHOA_L2_TWH_f = ZSCORE_TNHOA_L2_TWH_f,
                     ZSCORE_PN_PWH_f = ZSCORE_PN_PWH_f,
                     ZSCORE_PN_PW_L1_inputnorm_f = ZSCORE_PN_PW_L1_inputnorm_f,
                     ZSCORE_PN_PW_L2_inputnorm_f = ZSCORE_PN_PW_L2_inputnorm_f,
                     ZSCORE_PNH_L1_PWH_inputnorm_f = ZSCORE_PNH_L1_PWH_inputnorm_f,
                     ZSCORE_PNH_L2_PWH_inputnorm_f = ZSCORE_PNH_L2_PWH_inputnorm_f,
                     ZSCORE_PNH_L1_PWCH_inputnorm_f = ZSCORE_PNH_L1_PWCH_inputnorm_f,
                     ZSCORE_PNCH_PWCH_inputnorm_f = ZSCORE_PNCH_PWCH_inputnorm_f,
                     ZSCORE_PNHOA_L1_PWH_inputnorm_f = ZSCORE_PNHOA_L1_PWH_inputnorm_f,
                     ZSCORE_PNHOA_L2_PWH_inputnorm_f = ZSCORE_PNHOA_L2_PWH_inputnorm_f,
                     ZSCORE_PNHOA_L1_TNHOA_L1_inputnorm_f = ZSCORE_PNHOA_L1_TNHOA_L1_inputnorm_f,
                     ZSCORE_PNHOA_L2_TNHOA_L2_inputnorm_f = ZSCORE_PNHOA_L2_TNHOA_L2_inputnorm_f,
                     ZSCORE_TNHOA_L1_TWH_inputnorm_f = ZSCORE_TNHOA_L1_TWH_inputnorm_f,
                     ZSCORE_TNHOA_L2_TWH_inputnorm_f = ZSCORE_TNHOA_L2_TWH_inputnorm_f,
                     ZSCORE_PHYP_B_L1_PW_L1_inputnorm_f = ZSCORE_PHYP_B_L1_PW_L1_inputnorm_f,
                     ZSCORE_PHYP_B_L2_PW_L2_inputnorm_f = ZSCORE_PHYP_B_L2_PW_L2_inputnorm_f,
                     ZSCORE_PHYPH_B_L1_PWH_inputnorm_f = ZSCORE_PHYPH_B_L1_PWH_inputnorm_f,
                     ZSCORE_PHYPH_B_L2_PWH_inputnorm_f = ZSCORE_PHYPH_B_L2_PWH_inputnorm_f,
                     ZSCORE_PHYPH_B_L1_PWCH_inputnorm_f = ZSCORE_PHYPH_B_L1_PWCH_inputnorm_f,
                     ZSCORE_PHYPH_B_L2_PWCH_inputnorm_f = ZSCORE_PHYPH_B_L2_PWCH_inputnorm_f,
                     ZSCORE_PN_PWH_inputnorm_f = ZSCORE_PN_PWH_inputnorm_f,
                     Q_PHYP_B_L1_MNASEnormRPKM_GB_f = Q_PHYP_B_L1_MNASEnormRPKM_GB_f,
                     Q_PHYP_B_L2_MNASEnormRPKM_GB_f = Q_PHYP_B_L2_MNASEnormRPKM_GB_f,
                     Q_PHYPH_B_L1_MNASEnormRPKM_GB_f = Q_PHYPH_B_L1_MNASEnormRPKM_GB_f,
                     Q_PHYPH_B_L2_MNASEnormRPKM_GB_f = Q_PHYPH_B_L2_MNASEnormRPKM_GB_f,
                     Q_PNCH_MNASEnormRPKM_GB_f = Q_PNCH_MNASEnormRPKM_GB_f,
                     Q_PNH_L1_MNASEnormRPKM_GB_f = Q_PNH_L1_MNASEnormRPKM_GB_f,
                     Q_PNH_L2_MNASEnormRPKM_GB_f = Q_PNH_L2_MNASEnormRPKM_GB_f,
                     Q_PNH_MNASEnormRPKM_GB_f = Q_PNH_MNASEnormRPKM_GB_f,
                     Q_PNHOA_L1_MNASEnormRPKM_GB_f = Q_PNHOA_L1_MNASEnormRPKM_GB_f,
                     Q_PNHOA_L2_MNASEnormRPKM_GB_f = Q_PNHOA_L2_MNASEnormRPKM_GB_f,
                     Q_PN_MNASEnormRPKM_GB_f = Q_PN_MNASEnormRPKM_GB_f,
                     Q_PWCH_MNASEnormRPKM_GB_f = Q_PWCH_MNASEnormRPKM_GB_f,
                     Q_PWH_MNASEnormRPKM_GB_f = Q_PWH_MNASEnormRPKM_GB_f,
                     Q_PW_L1_MNASEnormRPKM_GB_f = Q_PW_L1_MNASEnormRPKM_GB_f,
                     Q_PW_L2_MNASEnormRPKM_GB_f = Q_PW_L2_MNASEnormRPKM_GB_f,
                     Q_TNHOA_L1_MNASEnormRPKM_GB_f = Q_TNHOA_L1_MNASEnormRPKM_GB_f,
                     Q_TNHOA_L2_MNASEnormRPKM_GB_f = Q_TNHOA_L2_MNASEnormRPKM_GB_f,
                     Q_TWH_MNASEnormRPKM_GB_f = Q_TWH_MNASEnormRPKM_GB_f,
                     ZSCORE_PHYPH_B_L1_PW_L1_MNASEnormRPKM_f = ZSCORE_PHYPH_B_L1_PW_L1_MNASEnormRPKM_f,
                     ZSCORE_PHYPH_B_L2_PW_L2_MNASEnormRPKM_f = ZSCORE_PHYPH_B_L2_PW_L2_MNASEnormRPKM_f,
                     ZSCORE_PN_PW_L1_MNASEnormRPKM_f = ZSCORE_PN_PW_L1_MNASEnormRPKM_f,
                     ZSCORE_PN_PW_L2_MNASEnormRPKM_f = ZSCORE_PN_PW_L2_MNASEnormRPKM_f,
                     ZSCORE_PHYPH_B_L1_PWH_MNASEnormRPKM_f = ZSCORE_PHYPH_B_L1_PWH_MNASEnormRPKM_f,
                     ZSCORE_PHYPH_B_L2_PWH_MNASEnormRPKM_f = ZSCORE_PHYPH_B_L2_PWH_MNASEnormRPKM_f,
                     ZSCORE_PNH_L1_PWH_MNASEnormRPKM_f = ZSCORE_PNH_L1_PWH_MNASEnormRPKM_f,
                     ZSCORE_PNH_L2_PWH_MNASEnormRPKM_f = ZSCORE_PNH_L2_PWH_MNASEnormRPKM_f,
                     ZSCORE_PNCH_PWCH_MNASEnormRPKM_f = ZSCORE_PNCH_PWCH_MNASEnormRPKM_f,
                     ZSCORE_PNHOA_L1_PWH_MNASEnormRPKM_f = ZSCORE_PNHOA_L1_PWH_MNASEnormRPKM_f,
                     ZSCORE_PNHOA_L2_PWH_MNASEnormRPKM_f = ZSCORE_PNHOA_L2_PWH_MNASEnormRPKM_f,
                     ZSCORE_PNHOA_L1_TNHOA_L1_MNASEnormRPKM_f = ZSCORE_PNHOA_L1_TNHOA_L1_MNASEnormRPKM_f,
                     ZSCORE_PNHOA_L2_TNHOA_L2_MNASEnormRPKM_f = ZSCORE_PNHOA_L2_TNHOA_L2_MNASEnormRPKM_f,
                     ZSCORE_TNHOA_L1_TWH_MNASEnormRPKM_f = ZSCORE_TNHOA_L1_TWH_MNASEnormRPKM_f,
                     ZSCORE_TNHOA_L2_TWH_MNASEnormRPKM_f = ZSCORE_TNHOA_L2_TWH_MNASEnormRPKM_f,
                     LOG2RATIO_PWH_TWH_MNASEnormRPKM_f = LOG2RATIO_PWH_TWH_MNASEnormRPKM_f,
                     LOG2RATIO_PNHOA_TNHOA_MNASEnormRPKM_f = LOG2RATIO_PNHOA_TNHOA_MNASEnormRPKM_f,
                     LOG2RATIO_PWH_TWH_f = LOG2RATIO_PWH_TWH_f,
                     LOG2RATIO_PNHOA_TNHOA_f = LOG2RATIO_PNHOA_TNHOA_f)

saveRDS(LIST_H2AV_VEC, paste0(workdir, "PROJET_H2AV/DATA/LIST_FEATURES/LIST_H2AV_VEC.RDS"))
LIST_H2AV_VEC = readRDS(paste0(workdir, "PROJET_H2AV/DATA/LIST_FEATURES/LIST_H2AV_VEC.RDS"))
