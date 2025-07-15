#####################################################################################-
#          LOAD LIBRARIES  ----
####################################################################################-

library(gplots)
library(Vennerable)
library(SuperExactTest)
'%ni%' = Negate('%in%')



#####################################################################################-
#           FUNCTIONS  ----
#####################################################################################

workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"
outfig = paste0(workdir,"FIGURES/VENN/")
source(paste0(workdir, "functionR/FisherList_Hmap.R"))
source(paste0(workdir, "functionR/VENN_vennerable.R"))

#####################################################################################-
#           DATA  ----
#####################################################################################
## QUANTIF

LIST_H2AV_VEC = readRDS(paste0(workdir, "DATA/LIST_FEATURES/LIST_H2AV_VEC.RDS"))

GENES_dm6_ORI = readRDS(paste0(workdir, "DATA/LIST_FEATURES/LIST_GN_ORI_dm6.RDS"))

GNref = readRDS(paste0(workdir, "DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))


LIST_UP_DN_dereg_genes_nelfKD = readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_UP_DN_dereg_genes_nelfKD.RDS"))
UPdereg_genes_nelfKD = LIST_UP_DN_dereg_genes_nelfKD$UPdereg_genes_nelfKD
DNdereg_genes_nelfKD = LIST_UP_DN_dereg_genes_nelfKD$DNdereg_genes_nelfKD

UPdereg_genes_nelfKD = paste0(UPdereg_genes_nelfKD,".1")
DNdereg_genes_nelfKD = paste0(DNdereg_genes_nelfKD,".1")

###################################################################################

plotVENNerable = function(data_to_venn, Nb_REF, outdir){
	vennNAME = paste(names(data_to_venn), collapse="_")
	Venn_ovlp = Venn(data_to_venn)
	color_Venn = VennThemes(compute.Venn(Venn_ovlp))
  #set1
  color_Venn[["Set"]][["Set1"]]$col = "#016000"
  color_Venn[["SetText"]][["Set1"]]$col = "#016000"
  color_Venn[["Face"]][["100"]]$fill = "#63bf61"
  color_Venn[["Face"]][["100-1"]]$fill = "#63bf61"
  #set2
  color_Venn[["Set"]][["Set2"]]$col = "#ad6a05"
  color_Venn[["SetText"]][["Set2"]]$col = "#ad6a05"
  color_Venn[["Face"]][["010"]]$fill = "#fc9d44"
  color_Venn[["Face"]][["010-1"]]$fill = "#fc9d44"
  #set3
  color_Venn[["Set"]][["Set3"]]$col = "#232323"
  color_Venn[["SetText"]][["Set3"]]$col = "#232323"
  color_Venn[["Face"]][["001"]]$fill = "#a3a3a3"
  color_Venn[["Face"]][["001-1"]]$fill = "#a3a3a3"
  #set1 n set2
  color_Venn[["Face"]][["110"]]$fill = "#b5a029"
  color_Venn[["Face"]][["110-1"]]$fill = "#b5a029"
  #set1 n set3
  color_Venn[["Face"]][["101"]]$fill = "#276825"
  color_Venn[["Face"]][["101-1"]]$fill = "#276825"
  #set2 n set3
	color_Venn[["Face"]][["011"]]$fill = "#84603f"
  color_Venn[["Face"]][["011-1"]]$fill = "#84603f"
	#fisher exact test
	fishertest = fisher_namesList(data_to_venn,data_to_venn,Nb_REF)
	# Supertest for triple intersection
	ResSuperTest = supertest(data_to_venn, n = Nb_REF)
	SupetTestToPrint = as.matrix(ResSuperTest$P.value)
	colnames(SupetTestToPrint) = "intersect_pval"
	########### PLOT #########################
	pdf(paste0(outdir, "PLOTvenn_",vennNAME,".pdf"))
	plot(Venn_ovlp, gp = color_Venn, doWeights = TRUE)
	par(mfrow=c(2,1))
	textplot(fishertest$p.value,  valign="top")
	textplot(fishertest$odds.ratio)
	textplot(SupetTestToPrint)
	plot(Venn_ovlp, gp = color_Venn, doWeights = TRUE, show = list(SetLabels = FALSE, FaceText = ""))
	dev.off()
}


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


# UP reg / TOP 5% H2AV
List1 = list(TOP5_H2AV = getNameList(LIST_H2AV_VEC$ZSCORE_PN_PW_L1_f, topdown = "top", prct = 5), UP_reg_NELF_KD = UPdereg_genes_nelfKD)
plotVENNerable(data_to_venn = List1, Nb_REF = length(GNref), outdir = outfig)

# DOWN reg / TOP 5% H2AV
List1 = list(TOP5_H2AV = getNameList(LIST_H2AV_VEC$ZSCORE_PN_PW_L1_f, topdown = "top", prct = 5), DN_reg_NELF_KD = DNdereg_genes_nelfKD)
plotVENNerable(data_to_venn = List1, Nb_REF = length(GNref), outdir = outfig)

# UP/DOWN reg / TOP 5% H2AV
List1 = list(TOP5_H2AV = getNameList(LIST_H2AV_VEC$ZSCORE_PN_PW_L1_f, topdown = "top", prct = 5), UP_reg_NELF_KD = UPdereg_genes_nelfKD, DN_reg_NELF_KD = DNdereg_genes_nelfKD)
plotVENNerable(data_to_venn = List1, Nb_REF = length(GNref), outdir = outfig)






# TOP 10 / BOT 10% ZSCORE_PNH_L1_PWH_f / ORI and NO ORI

List1 = list(UP_ZSCORE_PNH_L1_PWH_f = getNameList(LIST_H2AV_VEC$ZSCORE_PNH_L1_PWH_f, topdown = "top", prct = 10), DN_ZSCORE_PNH_L1_PWH_f = getNameList(LIST_H2AV_VEC$ZSCORE_PNH_L1_PWH_f, topdown = "down", prct = 10), GENES_dm6_ORI = GENES_dm6_ORI)
plotVENNerable(data_to_venn = List1, Nb_REF = length(GNref), outdir = outfig)

List1 = list(UP_ZSCORE_PNH_L1_PWH_f = getNameList(LIST_H2AV_VEC$ZSCORE_PNH_L1_PWH_f, topdown = "top", prct = 10), DN_ZSCORE_PNH_L1_PWH_f = getNameList(LIST_H2AV_VEC$ZSCORE_PNH_L1_PWH_f, topdown = "down", prct = 10), GENES_dm6_NO_ORI = GNref[GNref %ni% GENES_dm6_ORI])
plotVENNerable(data_to_venn = List1, Nb_REF = length(GNref), outdir = outfig)




# TOP 10% ZSCORE_PNH_L1_PWH_f / GN ACTIF / ORI and NO ORI

List1 = list(UP_ZSCORE_PNH_L1_PWH_f = getNameList(LIST_H2AV_VEC$ZSCORE_PNH_L1_PWH_f, topdown = "top", prct = 10), GN_ACTIF = GNref[0:6966], GENES_dm6_ORI = GENES_dm6_ORI)
plotVENNerable(data_to_venn = List1, Nb_REF = length(GNref), outdir = outfig)

List1 = list(UP_ZSCORE_PNH_L1_PWH_f = getNameList(LIST_H2AV_VEC$ZSCORE_PNH_L1_PWH_f, topdown = "top", prct = 10), GN_ACTIF = GNref[0:6954], GENES_dm6_NO_ORI = GNref[GNref %ni% GENES_dm6_ORI])
plotVENNerable(data_to_venn = List1, Nb_REF = length(GNref), outdir = outfig)



# TOP 10% ZSCORE_PNH_L1_PWH_f / ORI and NO ORI

List1 = list(UP_ZSCORE_PNH_L1_PWH_f = getNameList(LIST_H2AV_VEC$ZSCORE_PNH_L1_PWH_f, topdown = "top", prct = 10), GENES_dm6_ORI = GENES_dm6_ORI)
plotVENNerable(data_to_venn = List1, Nb_REF = length(GNref), outdir = outfig)

List1 = list(UP_ZSCORE_PNH_L1_PWH_f = getNameList(LIST_H2AV_VEC$ZSCORE_PNH_L1_PWH_f, topdown = "top", prct = 10), GENES_dm6_NO_ORI = GNref[GNref %ni% GENES_dm6_ORI])
plotVENNerable(data_to_venn = List1, Nb_REF = length(GNref), outdir = outfig)



# TOP 10% ZSCORE_PN_L1_PW_f / ORI and NO ORI

List1 = list(UP_ZSCORE_PN_PW_L1_f = getNameList(LIST_H2AV_VEC$ZSCORE_PN_PW_L1_f, topdown = "top", prct = 10), GENES_dm6_ORI = GENES_dm6_ORI)
plotVENNerable(data_to_venn = List1, Nb_REF = length(GNref), outdir = outfig)

List1 = list(UP_ZSCORE_PN_PW_L1_f = getNameList(LIST_H2AV_VEC$ZSCORE_PN_PW_L1_f, topdown = "top", prct = 10), GENES_dm6_NO_ORI = GNref[GNref %ni% GENES_dm6_ORI])
plotVENNerable(data_to_venn = List1, Nb_REF = length(GNref), outdir = outfig)

