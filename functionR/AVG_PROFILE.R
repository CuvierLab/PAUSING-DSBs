


#####################################################################################-
#          LOAD LIBRARIES  ----
#####################################################################################-
#library(dplyr) # used for recover names(coverages)
##library(gsubfn) # used in bam2coverage
#library(Rsamtools)
#library(GenomicRanges)
#library(GenomeInfoDb)
#library("GenomicFeatures")
#library("GenomicAlignments")
#library("BiocParallel")
#library(normr)
#library(seqplots)
#library(rtracklayer)
#library(devtools)
#'%ni%' = Negate('%in%')
#library("BSgenome.Hsapiens.UCSC.hg38")
#library("BSgenome")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#
#BiocManager::install("seqplots")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)


#####################################################################################-
#          INSIDE FUNCTIONS  ----
#####################################################################################-
 myplotfun <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin){
   for(mygr in gr.v){
     sze <- length(get(mygr))
     # print(mygr)
     o.tmp <- toString(rtracklayer::export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed")))
     # o.tmp <- toString(export.bed(mygr,paste0(tmp,"/",mygr,"_#",sze,"peaks.bed")))
     #plotTMP <- getPlotSetArray(bw.l,o.tmp,"hg38",bin = bin,ignore_strand = F,xmin = xlim[1],xmax=xlim[2],rm0 = F, type = "af")
      plotTMP <- getPlotSetArray(bw.l,o.tmp,"dm6",bin = bin,ignore_strand = F,xmin = xlim[1],xmax=xlim[2],rm0 = F, type = "af")
file.remove(paste0(o.tmp))
     print(plotAverage(plotTMP,xlab='Relative position [bp]', ylim=ylim, ylab='number of reads',keepratio = F,error.estimates = T,colvec =c("black","firebrick2")))
     print(plotAverage(plotTMP,xlab='Relative position [bp]', ylim=ylim, ylab='number of reads',keepratio = F,error.estimates = T,colvec =c("black","firebrick2"), legend=F))
   }
 }

#####################################################################################-

seqPlotSDoutliers_scaleFact <- function(bw.l,tmp,gr.v,ylim,xlim=xlim,bin=bin,anchor=10000,sd=c(F,10),err=T,type="pf",smooth=T,spar=0.7, KR = F,gnme=NA,ignore.strand=F, scalingF = c(1,1), colvec = c("black","firebrick2"))
{
  bw.n <- NULL
  o.tmp <- NULL
  for(n in 1:length(bw.l)){
    bw.c <- bw.l[n]
    bw.n[n] <- gsub("(.*).bw","\\1",basename(bw.c))
  }
  for(mygr in gr.v){
    sze <- length(get(mygr))
    print(mygr)
    o.tmp <- c(o.tmp,toString(rtracklayer::export.bed(get(mygr),paste0(tmp,"/",mygr,"_#",sze,"peaks.bed"))))
  }

  gpsa <- getPlotSetArray(bw.l,o.tmp,gnme,bin = bin,ignore_strand = ignore.strand,xmin = xlim[1],xmax=xlim[2],rm0 = F,type=type, xanchored=anchor)
  gpsa.data <- gpsa$data
  for(mygr in gr.v){
    for(my.bw in bw.n){
      sze <- length(get(mygr))
      gpsa.mtx <- data.frame(gpsa.data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["heatmap"]])
      Value_per_bin = c(gpsa.mtx)
      means <- colMeans(gpsa.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
      means[which(is.nan(means))] <- 0 # change NA in 0
      if(sd[1] %in% T){
        gpsa.scl.mtx <- gpsa.mtx %>% mutate_all(scale) # scale the data (center reduce)
        gpsa.scl.mtx[abs(gpsa.scl.mtx) > sd[2]] <- NA # Remove value X SD away (sd = 3 by default ~ 98% of the data)
        means <- colMeans(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,na.rm=T) # Now you can do the mean on original data without 3 SD away outliers
        stderror <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx,2,function(n){
          sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
        })
        conint <- apply(gpsa.mtx + gpsa.scl.mtx - gpsa.scl.mtx, 2, function(n) {
          qt(0.95, sum(!is.na(n))) * sd(n, na.rm = TRUE)/sqrt(sum(!is.na(n)))
        })
        stderror[is.na(stderror)] <- 0
        conint[is.na(conint)] <- 0
        if(smooth){
          means = smooth.spline(1:(length(means)), means, spar=spar)$y
       
        }
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["stderror"]] <- stderror # change the means vector from getPlotSetArray object
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["conint"]] <- conint  # change the means vector from getPlotSetArray object
      }
      if(smooth){
        means = smooth.spline(1:(length(means)), means, spar=spar)$y
      }
      if(my.bw == bw.n[[scalingF[1]]]){
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means*scalingF[2]      #16.33333*1.761243 # change the means vector from getPlotSetArray object + scale fact for comparison

      }else{
        gpsa$data[[paste0(mygr,"_#",sze,"peaks")]][[my.bw]][["means"]] <- means # change the means vector from getPlotSetArray object
      }
    }

  }
  file.remove(paste0(o.tmp))
  par(mfrow=c(1,1))
  plotAverage(gpsa,xlab='Relative position [bp]', ylim=ylim, ylab='Signal', main = paste0("Plot profile \n SD Removed_",sd[1]," ",sd[2]), keepratio = KR,error.estimates = err, colvec = colvec, pointsize = 20, legend_ext = T)
  plotAverage(gpsa,xlab='Relative position [bp]', ylim=ylim, ylab='Signal', main = paste0("Plot profile \n SD Removed_",sd[1]," ",sd[2]), keepratio = KR,error.estimates = err, colvec = colvec, pointsize = 20, legend=F)
  # par(mfrow=c(2,2))
  # plot(density(Value_per_bin))

}



#####################################################################################-

split_GRanges_inList <- function(GR, NnamesSplit, Nsplit = NULL){
  # namesSplit is either a character vector (then  decreasinglyordered splitted with Nsplit)
  # either a list of character vectors (then a list of granges is created according to list of names )
  namesSplit = get(NnamesSplit)
  if(is.numeric(namesSplit)){
    namesSplit = names(namesSplit[order(namesSplit, decreasing=T)])
    GR = GR[namesSplit]
    GRList = split(GR, ceiling(seq_along(GR)/ceiling(length(GR)/Nsplit)))
    names(GRList) =  paste0("GR_",NnamesSplit, "_Q" ,seq(1,Nsplit,1))
  }
  if(is.character(namesSplit)){
    GR = GR[namesSplit]
    GRList = split(GR, ceiling(seq_along(GR)/ceiling(length(GR)/Nsplit)))
    names(GRList) =  paste0("GR_",NnamesSplit, "_Q" ,seq(1,Nsplit,1))
  }
  if(is.list(namesSplit)){
    GRList = list()
    GRList = unlist(lapply(namesSplit, function(subnames){GRList = c(GRList, GR[names(GR) %in% subnames])}))
    names(GRList) =  paste0("GR_",names(namesSplit))
  }
  return(GRList)
}

## Sort each list of genes in POSDOM according to a quantif
sortListGenes = function(GNList, Quantif){
## EXEMPLE ::  List_genes_DOM_pause_indice_ctrl = sortListGenes(GNList = List_genes_DOM, Quantif = pause_indice_ctrl)
  lapply(GNList, function(GN){
    QuantifGN = Quantif[names(Quantif) %in% GN]
    GN = names(QuantifGN[order(QuantifGN, decreasing=T)])
  })
}



#### GET GN LIST of a vector decile or else
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
