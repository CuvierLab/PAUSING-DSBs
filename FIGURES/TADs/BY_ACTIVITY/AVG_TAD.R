#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-

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
#####################################################################################-
#        GENE LEVEL 
#####################################################################################-
# Load FUNCTION
  workdir="/home/cperrois/work/"
  fonctiondir="/work/cperrois/functionR/"
  source(paste0(fonctiondir, "AVG_PROFILE.R"))

# BIGWIG
  #tmp <- create(paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_TAD/TMPgetPlotSetArray"))
  tmp <- paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_TAD/TMPgetPlotSetArray")
  bwdir = paste0(workdir, "PROJET_H2AV_2021/DATA/BIGWIG/")
  BW_PW_R1=paste0(bwdir,"H2AV_PW_1_L1_RPGC.bw")
  BW_PW_R2=paste0(bwdir,"H2AV_PW_2_L1_RPGC.bw")
  BW_PHYPB_R1=paste0(bwdir,"PHYP_B_L1_trimmed_filt_sort_RPGC.bw")
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



#TAD 
  TAD_ramirez_2015_dm6_active = rtracklayer::import.bed( paste0(workdir,"PROJET_H2AV/DATA/HIC/TAD_RAMIREZ/TAD_ramirez_2015_dm6_active.bed"))
  TAD_ramirez_2015_dm6_inactive = rtracklayer::import.bed( paste0(workdir,"PROJET_H2AV/DATA/HIC/TAD_RAMIREZ/TAD_ramirez_2015_dm6_inactive.bed"))

  TAD_ramirez_2015_dm6_active = reduce(TAD_ramirez_2015_dm6_active)
  TAD_ramirez_2015_dm6_inactive = reduce(TAD_ramirez_2015_dm6_inactive)
  
# Refrence Genome 
  r6_ref_genes=readRDS('/work/cperrois/PROJET_H2AV_2021/DATA/r6.13/TxDb.GR.dm6.RDS')
  GNref = readRDS(paste0(workdir, "PROJET_H2AV/DATA/LIST_FEATURES/LIST_GN_ACTIFS.RDS"))


# PLOT AVG 
  GR_list_toPlot=c(
  "TAD_ramirez_2015_dm6_active"
  )
  for(GR in GR_list_toPlot){
    pdf(paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_TAD/AVG_TAD_ALL/AVG_PLOT_" ,GR,".pdf"))
    par(lwd=2)
    # H2AV 
      seqPlotSDoutliers_scaleFact(c(BW_PW_R1,BW_PN_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PW_R2,BW_PN_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PWH_R2,BW_PNH_R2_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,4),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PWHOA_L1,BW_PNHOA_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,3),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PW_R2,BW_PWH_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,4),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PW_R1,BW_PWH_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PN_R2,BW_PNHOA_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PN_R1,BW_PNH_R2_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,4),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PWH_R1,BW_PHYPBH_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PWH_R2,BW_PHYPBH_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,4),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PW_R1,BW_PHYPB_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PW_R2,BW_PHYPB_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,4),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PHYPB_R1,BW_PHYPBH_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PHYPB_R2,BW_PHYPBH_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_active)),ylim=c(1,4),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1"))
   
    dev.off()
    }


  GR_list_toPlot=c(
  "TAD_ramirez_2015_dm6_inactive"
  )
  for(GR in GR_list_toPlot){
    pdf(paste0(workdir,"PROJET_H2AV_2021/FIGURE/AVG_PLOT/AVG_TAD/AVG_TAD_ALL/AVG_PLOT_H2AV_" ,GR,".pdf"))
    par(lwd=2)
      seqPlotSDoutliers_scaleFact(c(BW_PW_R1,BW_PN_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PW_R2,BW_PN_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PWH_R2,BW_PNH_R2_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,4),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PWHOA_L1,BW_PNHOA_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,3),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PW_R2,BW_PWH_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,4),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PW_R1,BW_PWH_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PN_R2,BW_PNHOA_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PN_R1,BW_PNH_R2_L1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,4),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PWH_R1, BW_PHYPBH_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PWH_R2, BW_PHYPBH_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,4),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PW_R1, BW_PHYPB_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PW_R2, BW_PHYPB_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,4),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PHYPB_R1, BW_PHYPBH_R1),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,2),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
      seqPlotSDoutliers_scaleFact(c(BW_PHYPB_R2, BW_PHYPBH_R2),tmp,GR,anchor=median(width(TAD_ramirez_2015_dm6_inactive)),ylim=c(1,4),xlim=c(25000,25000),type="af",bin=10,smooth=T,spar=0.20, scalingF = c(1,1), sd=c(T,3), gnme="dm6", colvec = c("#285bad", "#eb3434","#2bce58","#ad07c1")) 
    dev.off()
}

