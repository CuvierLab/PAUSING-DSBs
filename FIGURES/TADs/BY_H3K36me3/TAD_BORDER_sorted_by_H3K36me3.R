
#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

library(GenomicFeatures)

.libPaths(c(.libPaths(),"/home/cperrois/work/Rpackages/R-3.4.3/"))
library(dplyr) # used for recover names(coverages)
library(gsubfn) # used in bam2coverage
library(Rsamtools)
library(GenomicRanges)
library(GenomeInfoDb)
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library(normr)
library("htmlwidgets",lib="/home/cperrois/work/Rpackages/R-3.4.3/")
library(seqplots)
library("BSgenome.Dmelanogaster.UCSC.dm6")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library(rtracklayer)
library(devtools)
'%ni%' = Negate('%in%')
#library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome")
"%ni%" = Negate("%in%")
# Load FUNCTION
workdir = "/home/cperrois/work/PROJET_H2AV_2021/"
fonctiondir="/work/cperrois/functionR/"
source(paste0(fonctiondir, "AVG_PROFILE.R"))


#######################################################

LIST_QUANTIF_K36=readRDS("/home/cperrois/work/PROJET_H2AV/DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS")
Q_H3K36me3_2C4_GB_f=LIST_QUANTIF_K36$Q_H3K36me3_2C4_GB_f
ZSCORE_H3K36me3_2N4_H3K36me3_2C4_f=LIST_QUANTIF_K36$ZSCORE_H3K36me3_2N4_H3K36me3_2C4_f
workdir2 = "/home/cperrois/work/"

PAUSE_INDICE_VEC = readRDS(paste0(workdir2, "PROJET_H2AV/DATA/LIST_FEATURES/PAUSE_INDICE_VEC.RDS"))
PAUSE_IND_pol2_ctrl_N = PAUSE_INDICE_VEC$PAUSE_IND_pol2_ctrl_N
PAUSE_IND_pol2_ctrl_N_f=PAUSE_IND_pol2_ctrl_N[names(PAUSE_IND_pol2_ctrl_N) %in% GNref]
DELTA_PAUSE_IND_pol2_nelf_N = PAUSE_INDICE_VEC$DELTA_PAUSE_IND_pol2_nelf_N


GNref = readRDS(paste0(workdir, "../PROJET_H2AV_2021/DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))
gene_dm6_gr = readRDS(paste0(workdir, "DATA/r6.13/TxDb.Dmelanogaster.Ensembl.dm6.RDS"))
seqlevelsStyle(gene_dm6_gr) <- "UCSC"
names(gene_dm6_gr)=paste0(names(gene_dm6_gr),".1")
gene_dm6_gr_TSS = resize(gene_dm6_gr, 1, fix="start")
gene_dm6_gr_TSS_act=gene_dm6_gr_TSS[names(gene_dm6_gr_TSS) %in% GNref]


### TD actif 
TAD_ramirez_2015_dm6_active = rtracklayer::import.bed( paste0(workdir,"../PROJET_H2AV/DATA/HIC/TAD_RAMIREZ/TAD_ramirez_2015_dm6_active.bed"))


GR_TAD_ramirez_2015_dm6_active_5_3=TAD_ramirez_2015_dm6_active
GR_TAD_ramirez_2015_dm6_active_5_3$name = paste0("TAD_",seqnames(GR_TAD_ramirez_2015_dm6_active_5_3) ,"_",start(GR_TAD_ramirez_2015_dm6_active_5_3))
strand(GR_TAD_ramirez_2015_dm6_active_5_3)="+"

###

GR_TAD_ramirez_2015_dm6_active_3_5=GRanges( 
                  seqnames=seqnames(GR_TAD_ramirez_2015_dm6_active_5_3),
                  name=paste0("TAD_",seqnames(GR_TAD_ramirez_2015_dm6_active_5_3),"_",end(GR_TAD_ramirez_2015_dm6_active_5_3)),
                  ranges = IRanges(start=start(GR_TAD_ramirez_2015_dm6_active_5_3),end=end(GR_TAD_ramirez_2015_dm6_active_5_3)),
                  strand="-",
                  score=score(GR_TAD_ramirez_2015_dm6_active_5_3))


### TAD  sens : 5_3 
GR_TAD_active_5_3=resize(GR_TAD_ramirez_2015_dm6_active_5_3, width=1, fix="start")
dist_GR_TAD_active_5_3_to_gene=distanceToNearest(GR_TAD_active_5_3, gene_dm6_gr_TSS_act, ignore.strand=FALSE)
df_dist_GR_TAD_active_5_3_to_gene = as.data.frame(distanceToNearest(GR_TAD_active_5_3, gene_dm6_gr_TSS_act))
colnames=c("TAD" ,"Gene", "distance")
colnames(df_dist_GR_TAD_active_5_3_to_gene)=colnames
df_dist_GR_TAD_active_5_3_to_gene$Gene=names(gene_dm6_gr_TSS_act[dist_GR_TAD_active_5_3_to_gene@to])
df_dist_GR_TAD_active_5_3_to_gene$start=start(GR_TAD_ramirez_2015_dm6_active_5_3)
df_dist_GR_TAD_active_5_3_to_gene$end=end(GR_TAD_ramirez_2015_dm6_active_5_3)
df_dist_GR_TAD_active_5_3_to_gene$strand="+"
df_dist_GR_TAD_active_5_3_to_gene$name=GR_TAD_active_5_3$name
df_dist_GR_TAD_active_5_3_to_gene$TAD=paste0(df_dist_GR_TAD_active_5_3_to_gene$TAD,"_R")

### TAD  sens : 3_5 
GR_TAD_active_3_5=resize(GR_TAD_ramirez_2015_dm6_active_3_5, width=1, fix="start")
dist_GR_TAD_active_3_5_to_gene=distanceToNearest(GR_TAD_active_3_5, gene_dm6_gr_TSS_act, ignore.strand=FALSE)
df_dist_GR_TAD_active_3_5_to_gene = as.data.frame(distanceToNearest(GR_TAD_active_3_5, gene_dm6_gr_TSS_act))
colnames=c("TAD" ,"Gene", "distance")
colnames(df_dist_GR_TAD_active_3_5_to_gene)=colnames
df_dist_GR_TAD_active_3_5_to_gene$Gene=names(gene_dm6_gr_TSS_act[dist_GR_TAD_active_3_5_to_gene@to])
df_dist_GR_TAD_active_3_5_to_gene$start=start(GR_TAD_ramirez_2015_dm6_active_5_3)
df_dist_GR_TAD_active_3_5_to_gene$end=end(GR_TAD_ramirez_2015_dm6_active_5_3)
df_dist_GR_TAD_active_3_5_to_gene$strand="-"
df_dist_GR_TAD_active_3_5_to_gene$name=GR_TAD_active_3_5$name
df_dist_GR_TAD_active_3_5_to_gene$TAD=paste0(df_dist_GR_TAD_active_3_5_to_gene$TAD,"_L")


#####################
df_dist_GR_TAD_active_to_gene_dup_flip=rbind(df_dist_GR_TAD_active_5_3_to_gene,df_dist_GR_TAD_active_3_5_to_gene)
#####################
for (elm in (1:dim(df_dist_GR_TAD_active_to_gene_dup_flip)[1])){
df_dist_GR_TAD_active_to_gene_dup_flip$Q_H3K36me3_2C4_GB_f[elm]=as.numeric(Q_H3K36me3_2C4_GB_f[df_dist_GR_TAD_active_to_gene_dup_flip$Gene[elm]])
df_dist_GR_TAD_active_to_gene_dup_flip$ZSCORE_H3K36me3_2N4_H3K36me3_2C4_f[elm]=as.numeric(ZSCORE_H3K36me3_2N4_H3K36me3_2C4_f[df_dist_GR_TAD_active_to_gene_dup_flip$Gene[elm]])
df_dist_GR_TAD_active_to_gene_dup_flip$PAUSE_IND_pol2_ctrl_N_f[elm]=as.numeric(PAUSE_IND_pol2_ctrl_N_f[df_dist_GR_TAD_active_to_gene_dup_flip$Gene[elm]])
df_dist_GR_TAD_active_to_gene_dup_flip$DELTA_PAUSE_IND_pol2_nelf_N[elm]=as.numeric(DELTA_PAUSE_IND_pol2_nelf_N[df_dist_GR_TAD_active_to_gene_dup_flip$Gene[elm]])
}


saveRDS(df_dist_GR_TAD_active_to_gene_dup_flip,"/home/cperrois/work/PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_TAD_BORDER/TAD_ACTIF_BORDER/TAD_BORDER_nearest_dis_GN_H3K36me3_oriented_posit_start_order/df_dist_GR_TAD_active_to_gene_dup_flip.RDS")

row_to_keep = which(df_dist_GR_TAD_active_to_gene_dup_flip$distance < 5000)
df_dist_GR_TAD_active_to_gene_dup_flip = df_dist_GR_TAD_active_to_gene_dup_flip[row_to_keep,]
summary(as.numeric(df_dist_GR_TAD_active_to_gene_dup_flip$distance))

# Q_H3K36me3_2C4_GB_f
df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3<-df_dist_GR_TAD_active_to_gene_dup_flip[order(df_dist_GR_TAD_active_to_gene_dup_flip[,8],decreasing=T),]
head(df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3)


#### Quartile construction
Q1_TAD_B_odr_Q_H3K36me3_start_posit_name=head(df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3$name,round(dim(df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3)[1]/4))
D50UP_TAD_B_odr_Q_H3K36me3_start_posit=head(df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3,round(dim(df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3)[1]/2))
Q2_TAD_B_odr_Q_H3K36me3_start_posit_name=tail(D50UP_TAD_B_odr_Q_H3K36me3_start_posit$name,round(dim(D50UP_TAD_B_odr_Q_H3K36me3_start_posit)[1]/2))
D50DN_TAD_B_odr_Q_H3K36me3_start_posit=tail(df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3,round(dim(df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3)[1]/2))
Q3_TAD_B_odr_Q_H3K36me3_start_posit_name=head(D50DN_TAD_B_odr_Q_H3K36me3_start_posit$name,round(dim(D50DN_TAD_B_odr_Q_H3K36me3_start_posit)[1]/2))
Q4_TAD_B_odr_Q_H3K36me3_start_posit_name=tail(df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3$name,round(dim(df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3)[1]/4))

# quartile Data frame 
row_to_keep = which(df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3$name %in% Q1_TAD_B_odr_Q_H3K36me3_start_posit_name )
df_Q1_TAD_B_odr_Q_H3K36me3_start_posit = df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3[row_to_keep,]
row_to_keep = which(df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3$name %in% Q2_TAD_B_odr_Q_H3K36me3_start_posit_name )
df_Q2_TAD_B_odr_Q_H3K36me3_start_posit = df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3[row_to_keep,]

row_to_keep = which(df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3$name %in% Q3_TAD_B_odr_Q_H3K36me3_start_posit_name )
df_Q3_TAD_B_odr_Q_H3K36me3_start_posit = df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3[row_to_keep,]
row_to_keep = which(df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3$name %in% Q4_TAD_B_odr_Q_H3K36me3_start_posit_name )
df_Q4_TAD_B_odr_Q_H3K36me3_start_posit = df_dist_GR_TAD_active_to_gene_dup_flip_odr_Q_H3K36me3[row_to_keep,]


##############################################################"""

vec_ch_Q1=c()
for(elm in 1:dim(df_Q1_TAD_B_odr_Q_H3K36me3_start_posit)[1]){
ch_spli=strsplit(df_Q1_TAD_B_odr_Q_H3K36me3_start_posit$name,"_")
ch=ch_spli[[elm]][2]
vec_ch_Q1=c(vec_ch_Q1,ch)
}
GR_Q1_TAD_B_odr_Q_H3K36me_start_posit= GRanges(      
                  seqnames=vec_ch_Q1,
                  names=df_Q1_TAD_B_odr_Q_H3K36me3_start_posit$name,
                  ranges = IRanges(start=df_Q1_TAD_B_odr_Q_H3K36me3_start_posit$start,end=df_Q1_TAD_B_odr_Q_H3K36me3_start_posit$end),
                  nearest_GN_id= df_Q1_TAD_B_odr_Q_H3K36me3_start_posit$Gene,
                  TAD_B_ID = df_Q1_TAD_B_odr_Q_H3K36me3_start_posit$TAD,
                  distance = df_Q1_TAD_B_odr_Q_H3K36me3_start_posit$distance,
                  strand=df_Q1_TAD_B_odr_Q_H3K36me3_start_posit$strand
                  )


#####################""                
vec_ch_Q2=c()
for(elm in 1:dim(df_Q2_TAD_B_odr_Q_H3K36me3_start_posit)[1]){
ch_spli=strsplit(df_Q2_TAD_B_odr_Q_H3K36me3_start_posit$name,"_")
ch=ch_spli[[elm]][2]
vec_ch_Q2=c(vec_ch_Q2,ch)
}
GR_Q2_TAD_B_odr_Q_H3K36me_start_posit= GRanges(      
                  seqnames=vec_ch_Q2,
                  names=df_Q2_TAD_B_odr_Q_H3K36me3_start_posit$name,
                  ranges = IRanges(start=df_Q2_TAD_B_odr_Q_H3K36me3_start_posit$start,end=df_Q2_TAD_B_odr_Q_H3K36me3_start_posit$end),
                  nearest_GN_id= df_Q2_TAD_B_odr_Q_H3K36me3_start_posit$Gene,
                  TAD_B_ID = df_Q2_TAD_B_odr_Q_H3K36me3_start_posit$TAD,
                  distance = df_Q2_TAD_B_odr_Q_H3K36me3_start_posit$distance,
                  strand=df_Q2_TAD_B_odr_Q_H3K36me3_start_posit$strand
                  )

###########################

vec_ch_Q3=c()
for(elm in 1:dim(df_Q3_TAD_B_odr_Q_H3K36me3_start_posit)[1]){
ch_spli=strsplit(df_Q3_TAD_B_odr_Q_H3K36me3_start_posit$name,"_")
ch=ch_spli[[elm]][2]
vec_ch_Q3=c(vec_ch_Q3,ch)
}
GR_Q3_TAD_B_odr_Q_H3K36me_start_posit= GRanges(      
                  seqnames=vec_ch_Q3,
                  names=df_Q3_TAD_B_odr_Q_H3K36me3_start_posit$name,
                  ranges = IRanges(start=df_Q3_TAD_B_odr_Q_H3K36me3_start_posit$start,end=df_Q3_TAD_B_odr_Q_H3K36me3_start_posit$end),
                  nearest_GN_id= df_Q3_TAD_B_odr_Q_H3K36me3_start_posit$Gene,
                  TAD_B_ID = df_Q3_TAD_B_odr_Q_H3K36me3_start_posit$TAD,
                  distance = df_Q3_TAD_B_odr_Q_H3K36me3_start_posit$distance,
                  strand=df_Q3_TAD_B_odr_Q_H3K36me3_start_posit$strand
                  )
#############################################""
vec_ch_Q4=c()
for(elm in 1:dim(df_Q4_TAD_B_odr_Q_H3K36me3_start_posit)[1]){
ch_spli=strsplit(df_Q4_TAD_B_odr_Q_H3K36me3_start_posit$name,"_")
ch=ch_spli[[elm]][2]
vec_ch_Q4=c(vec_ch_Q4,ch)
}
GR_Q4_TAD_B_odr_Q_H3K36me_start_posit= GRanges(      
                  seqnames=vec_ch_Q4,
                  names=df_Q4_TAD_B_odr_Q_H3K36me3_start_posit$name,
                  ranges = IRanges(start=df_Q4_TAD_B_odr_Q_H3K36me3_start_posit$start,end=df_Q4_TAD_B_odr_Q_H3K36me3_start_posit$end),
                  nearest_GN_id= df_Q4_TAD_B_odr_Q_H3K36me3_start_posit$Gene,
                  TAD_B_ID = df_Q4_TAD_B_odr_Q_H3K36me3_start_posit$TAD,
                  distance = df_Q4_TAD_B_odr_Q_H3K36me3_start_posit$distance,
                  strand=df_Q4_TAD_B_odr_Q_H3K36me3_start_posit$strand
                  )
########################################

vec_bord_Q1=c()
for(elm in 1:length(GR_Q1_TAD_B_odr_Q_H3K36me_start_posit)){
ch_spli=strsplit(GR_Q1_TAD_B_odr_Q_H3K36me_start_posit$TAD_B_ID,"_")
#print(ch_spli[[elm]][3])
ch=ch_spli[[elm]][2]
vec_bord_Q1=c(vec_bord_Q1,ch)}
index_TAD_R=as.vector(which(vec_bord_Q1=="R"))
index_TAD_L=as.vector(which(vec_bord_Q1=="L"))
GR_Q1_TAD_B_R_odr_Q_H3K36me_start_posit=GR_Q1_TAD_B_odr_Q_H3K36me_start_posit[index_TAD_R]
GR_Q1_TAD_B_L_odr_Q_H3K36me_start_posit=GR_Q1_TAD_B_odr_Q_H3K36me_start_posit[index_TAD_L]



vec_bord_Q2=c()
for(elm in 1:length(GR_Q2_TAD_B_odr_Q_H3K36me_start_posit)){
ch_spli=strsplit(GR_Q2_TAD_B_odr_Q_H3K36me_start_posit$TAD_B_ID,"_")
ch=ch_spli[[elm]][2]
vec_bord_Q2=c(vec_bord_Q2,ch)}
index_TAD_R=as.vector(which(vec_bord_Q2=="R"))
index_TAD_L=as.vector(which(vec_bord_Q2=="L"))
GR_Q2_TAD_B_R_odr_Q_H3K36me_start_posit=GR_Q2_TAD_B_odr_Q_H3K36me_start_posit[index_TAD_R]
GR_Q2_TAD_B_L_odr_Q_H3K36me_start_posit=GR_Q2_TAD_B_odr_Q_H3K36me_start_posit[index_TAD_L]


vec_bord_Q3=c()
for(elm in 1:length(GR_Q3_TAD_B_odr_Q_H3K36me_start_posit)){
ch_spli=strsplit(GR_Q3_TAD_B_odr_Q_H3K36me_start_posit$TAD_B_ID,"_")
ch=ch_spli[[elm]][2]
vec_bord_Q3=c(vec_bord_Q3,ch)}
index_TAD_R=as.vector(which(vec_bord_Q3=="R"))
index_TAD_L=as.vector(which(vec_bord_Q3=="L"))
GR_Q3_TAD_B_R_odr_Q_H3K36me_start_posit=GR_Q3_TAD_B_odr_Q_H3K36me_start_posit[index_TAD_R]
GR_Q3_TAD_B_L_odr_Q_H3K36me_start_posit=GR_Q3_TAD_B_odr_Q_H3K36me_start_posit[index_TAD_L]




vec_bord_Q4=c()
for(elm in 1:length(GR_Q4_TAD_B_odr_Q_H3K36me_start_posit)){
ch_spli=strsplit(GR_Q4_TAD_B_odr_Q_H3K36me_start_posit$TAD_B_ID,"_")
ch=ch_spli[[elm]][2]
vec_bord_Q4=c(vec_bord_Q4,ch)}
index_TAD_R=as.vector(which(vec_bord_Q4=="R"))
index_TAD_L=as.vector(which(vec_bord_Q4=="L"))
GR_Q4_TAD_B_R_odr_Q_H3K36me_start_posit=GR_Q4_TAD_B_odr_Q_H3K36me_start_posit[index_TAD_R]
GR_Q4_TAD_B_L_odr_Q_H3K36me_start_posit=GR_Q4_TAD_B_odr_Q_H3K36me_start_posit[index_TAD_L]



##### AVG 
  workdir2="/home/cperrois/work/"
  workdir="/home/cperrois/work/"
# BIGWIG
tmp <- paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_TAD_BORDER/TAD_ACTIF_BORDER/TMPgetPlotSetArray")
 bwdir = paste0(workdir2, "PROJET_H2AV_2021/DATA/BIGWIG/")
  BW_PW_R1=paste0(bwdir,"H2AV_PW_1_L1_RPGC.bw")
  BW_PW_R2=paste0(bwdir,"H2AV_PW_2_L1_RPGC.bw")
  BW_PHYP_R1=paste0(bwdir,"PHYP_B_L1_trimmed_filt_sort_RPGC.bw")
  BW_PHYPB_R2=paste0(bwdir,"H2AV_PHYPB_A_L1_RPGC.bw")
  BW_PN_R1=paste0(bwdir,"PN_trimmed_filt_sort_RPGC.bw")
  BW_PN_R2=paste0(bwdir,"H2AV_PN_J_L1_RPGC.bw")
  BW_PWH_R2=paste0(bwdir,"H2AV_PWH_A_L1_RPGC.bw")
  BW_PWHOA_L1="/home/cperrois/work/ALIGN_CHIPSEQ_H2AV_BGI26062021/05_BIGWIG/H2AV_PWHOA_L1_RPGC.bw"
  BW_PHYPBH_R1=paste0(bwdir,"PHYPH_B_L2_trimmed_filt_sort_RPGC.bw")
  BW_PHYPBH_R2=paste0(bwdir,"H2AV_PHYPBH_A_L1_RPGC.bw")
  BW_PNH_R2_L1="/home/cperrois/work/ALIGN_CHIPSEQ_H2AV_BGI26062021/05_BIGWIG/H2AV_PNH_J_L1_RPGC.bw"
  BW_PNHOA_L1="/home/cperrois/work/PROJET_H2AV/DATA/BIGWIG_H2AV/PNHOA_L1_trimmed_filt_sort_RPGC.bw"
  BW_PWH_R1=paste0(bwdir,"PWH_trimmed_filt_sort_RPGC.bw")

  BW_PHYPH_R1=paste0(bwdir,"PHYPH_B_L2_trimmed_filt_sort_RPGC.bw")
  BW_PHYPBH_R2=paste0(bwdir,"H2AV_PHYPBH_A_L1_RPGC.bw")

  BW_PHYPB_R1=paste0(bwdir,"PHYP_B_L1_trimmed_filt_sort_RPGC.bw")
  BW_PHYPB_R2=paste0(bwdir,"H2AV_PHYPB_A_L1_RPGC.bw")

  bwdir = paste0(workdir2, "PROJET_H2AV_2021/DATA/snakeMake_ALIGN_CHIPSEQ_Cuvier_dec2021/05_BIGWIG/")
  BW_Rad51_N_bis=paste0(bwdir,"Rad51_N_bis_RPGC.bw")
  BW_Rad51_N_HU_bis=paste0(bwdir,"Rad51_N_HU_bis_RPGC.bw")
  BW_Rad51_WT_bis=paste0(bwdir,"Rad51_WT_bis_RPGC.bw")
  BW_Rad51_WT_HU_bis=paste0(bwdir,"Rad51_WT_HU_bis_RPGC.bw")

##########  ALL BORDURE    ############################
  GR_list_toPlot=c(
  "GR_Q1_TAD_B_odr_Q_H3K36me_start_posit",
  "GR_Q2_TAD_B_odr_Q_H3K36me_start_posit",
  "GR_Q3_TAD_B_odr_Q_H3K36me_start_posit",
  "GR_Q4_TAD_B_odr_Q_H3K36me_start_posit"
  )


   for(GR in GR_list_toPlot)
    {
      pdf(paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_TAD_BORDER/TAD_ACTIF_BORDER/TAD_BORDER_CONTROLE/AVG_PROF_K36_" ,GR,".pdf"))
      par(lwd=2)
      seqPlotSDoutliers_scaleFact(c(BW_H3K36me3_2C4_RPGC_R1,BW_H3K36me3_2N4_RPGC_R1),tmp,GR,c(0,15),c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#1a8856", "#77ebb6")) 

      dev.off()
    }


    for(GR in GR_list_toPlot)
    {
      pdf(paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_TAD_BORDER/TAD_ACTIF_BORDER/TAD_BORDER_nearest_dis_GN_H3K36me3_oriented_posit_start_order/AVG_PROF_H2AV_N_" ,GR,".pdf"))
      par(lwd=2)
      # H2AV 
        seqPlotSDoutliers_scaleFact(c(BW_PW_R1,BW_PN_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,3),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20,scalingF = c(1,1), sd=c(T,3),gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PW_R2,BW_PN_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,3),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PWH_R2,BW_PNH_R2_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PWHOA_L1,BW_PNHOA_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PW_R2,BW_PWH_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PW_R1,BW_PWH_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PN_R2,BW_PNHOA_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PN_R1,BW_PNH_R2_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
      dev.off()
    }

 for(GR in GR_list_toPlot)
    {
      pdf(paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_TAD_BORDER/TAD_ACTIF_BORDER/TAD_BORDER_nearest_dis_GN_H3K36me3_oriented_posit_start_order/AVG_PROF_H2AV_HYPB_" ,GR,".pdf"))
      par(lwd=2)
      # H2AV 
        seqPlotSDoutliers_scaleFact(c(BW_PW_R1, BW_PHYPB_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,3),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PW_R2, BW_PHYPB_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,4),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PWH_R1, BW_PHYPH_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PWH_R2, BW_PHYPBH_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,3),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PHYPB_R1, BW_PHYPH_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,3),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PHYPB_R2, BW_PHYPBH_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,3),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
      dev.off()
    }


 for(GR in GR_list_toPlot)
    {
      pdf(paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_TAD_BORDER/TAD_ACTIF_BORDER/TAD_BORDER_nearest_dis_GN_H3K36me3_oriented_posit_start_order/AVG_PROF_RAD51_" ,GR,".pdf"))
      par(lwd=2)
        seqPlotSDoutliers_scaleFact(c(BW_Rad51_WT_bis,BW_Rad51_WT_HU_bis),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),c(1,3),c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
        seqPlotSDoutliers_scaleFact(c(BW_Rad51_N_bis,BW_Rad51_N_HU_bis),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),c(1,3),c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
        seqPlotSDoutliers_scaleFact(c(BW_Rad51_WT_bis,BW_Rad51_N_bis),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),c(1,3),c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
        seqPlotSDoutliers_scaleFact(c(BW_Rad51_WT_HU_bis,BW_Rad51_N_HU_bis),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),c(1,3),c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      dev.off()
    }

################################################"""""

 GR_list_toPlot=c(
"GR_Q1_TAD_B_R_odr_Q_H3K36me_start_posit",
"GR_Q1_TAD_B_L_odr_Q_H3K36me_start_posit",
"GR_Q2_TAD_B_R_odr_Q_H3K36me_start_posit",
"GR_Q2_TAD_B_L_odr_Q_H3K36me_start_posit",
"GR_Q3_TAD_B_R_odr_Q_H3K36me_start_posit",
"GR_Q3_TAD_B_L_odr_Q_H3K36me_start_posit",
"GR_Q4_TAD_B_R_odr_Q_H3K36me_start_posit",
"GR_Q4_TAD_B_L_odr_Q_H3K36me_start_posit")
for(GR in GR_list_toPlot)
    {
      pdf(paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_TAD_BORDER/TAD_ACTIF_BORDER/TAD_BORDER_CONTROLE/AVG_PROF_K36_" ,GR,".pdf"))
      par(lwd=2)
  		#seqPlotSDoutliers_scaleFact(c(BW_H3K36me3_2C4_RPGC_R1,BW_me3_L_bis_RPGC_R2),tmp,GR,c(0,15),c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#1a8856", "#77ebb6")) 
 		  #seqPlotSDoutliers_scaleFact(c(BW_H3K36me3_2N4_RPGC_R1,BW_me3_N_bis_RPGC_R2),tmp,GR,c(0,15),c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#1a8856", "#77ebb6")) 
      seqPlotSDoutliers_scaleFact(c(BW_H3K36me3_2C4_RPGC_R1,BW_H3K36me3_2N4_RPGC_R1),tmp,GR,c(0,15),c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#1a8856", "#77ebb6")) 

      dev.off()
    }


    for(GR in GR_list_toPlot)
    {
      pdf(paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_TAD_BORDER/TAD_ACTIF_BORDER/TAD_BORDER_nearest_dis_GN_H3K36me3_oriented_posit_start_order/AVG_PROF_H2AV_N_" ,GR,".pdf"))
      par(lwd=2)
      # H2AV 
        seqPlotSDoutliers_scaleFact(c(BW_PW_R1,BW_PN_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,3),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20,scalingF = c(1,1), sd=c(T,3),gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PW_R2,BW_PN_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,3),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PWH_R2,BW_PNH_R2_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PWHOA_L1,BW_PNHOA_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PW_R2,BW_PWH_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PW_R1,BW_PWH_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PN_R2,BW_PNHOA_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PN_R1,BW_PNH_R2_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
      dev.off()
    }

 for(GR in GR_list_toPlot)
    {
      pdf(paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_TAD_BORDER/TAD_ACTIF_BORDER/TAD_BORDER_nearest_dis_GN_H3K36me3_oriented_posit_start_order/AVG_PROF_H2AV_HYPB_" ,GR,".pdf"))
      par(lwd=2)
      # H2AV 
        seqPlotSDoutliers_scaleFact(c(BW_PW_R1, BW_PHYPB_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,3),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PW_R2, BW_PHYPB_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,4),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PWH_R1, BW_PHYPH_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,5),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PWH_R2, BW_PHYPBH_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,3),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PHYPB_R1, BW_PHYPH_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,3),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
        seqPlotSDoutliers_scaleFact(c(BW_PHYPB_R2, BW_PHYPBH_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,3),xlim=c(5000,5000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434")) 
      dev.off()
    }


#end 