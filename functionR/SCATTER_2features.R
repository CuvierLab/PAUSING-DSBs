
# Single Cell RNA seq project
# Cuvier''s Team
# Schaak - Heurteau - Depierre*
# 2017

#function to plot two vector and compute peasron R

'%ni%' = Negate('%in%')

library(MASS)
library(ggplot2)
library(viridis)

print("Use example : compareFunction(feat1 = , feat2 = ,  namefeat1 = , namefeat2 = , info = NULL, outpath = /path/to/out, yl=c(-1,1), xl=c(-1,1))")

theme_set(theme_bw(base_size = 16))

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


# compareFunction_listCol(feat1 = log(Q_Notch1_2+1), feat2 = log(Q_HDAC1_DMSO+1), namefeat1 = "Q_Notch1_2", namefeat2 = "Q_HDAC1_DMSO", listCol = colorList, info = "LOG_FilterHighLow", outpath = paste0(workdir,"FIGURE/SCATTER/"),  yl=c(0,15),  xl=c(0,15))
#
# colorList = list(NO_NOTCH = genespm500_ovlpPEAKS_ALLfeature$ALL_TSS[genespm500_ovlpPEAKS_ALLfeature$ALL_TSS %ni% unique(c(genespm500_ovlpPEAKS_ALLfeature$Notch1_1, genespm500_ovlpPEAKS_ALLfeature$Notch1_2))],
#               LOW_NOTCH = genespm500_ovlpPEAKS_ALLfeature$LOW_NOTCH, HIGH_NOTCH = genespm500_ovlpPEAKS_ALLfeature$HIGH_NOTCH)

print("Use example : compareFunction_listCol2(feat1, feat2,  namefeat1, namefeat2, listCol = NULL, vecCol = c(\"#839159\", \"#451c00\"), logT = T, info = NULL, outpath, yl=c(0,15), xl=c(0,15))")


compareFunction_listCol3 = function(feat1, feat2,  namefeat1, namefeat2, listCol = NULL, info = NULL, outpath, yl=c(0,15), xl=c(0,15)){
  ########################  NO FILTER  ##########################################################
  feat1 = as.matrix(feat1)
  feat2 = as.matrix(feat2)
	comRownames = Reduce(intersect, list(rownames(feat1),rownames(feat2)))
	feat1 = feat1[comRownames,1]
	feat2 = feat2[comRownames,1]
  data=data.frame(feat1=feat1, feat2=feat2)
  data$grpCol = NA
  data$NCol = NA
  data$grpCol[rownames(data) %in% listCol[[1]]] = "grey50"
  data$grpCol[rownames(data) %in% listCol[[2]]] = "turquoise"
  data$grpCol[rownames(data) %in% listCol[[3]]] = "blue"
  data$NCol[rownames(data) %in% listCol[[1]]] = names(listCol)[1]
  data$NCol[rownames(data) %in% listCol[[2]]] = names(listCol)[2]
  data$NCol[rownames(data) %in% listCol[[3]]] = names(listCol)[3]
	pearson = cor.test(data$feat1, data$feat2)
	pearsonT = round(pearson$estimate, 2)
  pearsonP = pearson$p.value
  data$density <- get_density(data$feat1,data$feat2, n=500)
  res.lm = lm(feat2~feat1)
  print(summary(res.lm))
  summary(res.lm)$r.squared
  # GGPLOT = ggplot(data = data) + geom_point(aes(x=feat1, y=feat2, color = density)) +
  #          scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
  #          #geom_smooth(method = "lm", formula=data$feat2~data$feat1, se = FALSE)
  GGPLOT = ggplot(data = data) + geom_point(aes(x=feat1, y=feat2, color = density), shape = 4) +
           scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)

           #geom_smooth(method = "lm", formula=data$feat2~data$feat1, se = FALSE)
  GGPLOT_colList = ggplot(data = data) + geom_point(aes(x=feat1, y=feat2, color = NCol), shape = 4) +
           scale_color_manual(values=c( "blue",  "turquoise", "grey50")) + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data$feat2~data$feat1, se = FALSE)
  ################### FILTER 1-2 ##################################################################
  feat1_f1 = feat1[names(feat1) %in% unique(c(listCol[[1]], listCol[[2]]))]
  feat2_f1 = feat2[names(feat2) %in% unique(c(listCol[[1]], listCol[[2]]))]
  feat2_f1 = as.matrix(feat2_f1)
  feat1_f1 = as.matrix(feat1_f1)
  feat2_f1 = as.matrix(feat2_f1)
	comRownames = Reduce(intersect, list(rownames(feat1_f1),rownames(feat2_f1)))
	feat1_f1 = feat1_f1[comRownames,1]
	feat2_f1 = feat2_f1[comRownames,1]
  data_f1=data.frame(feat1_f1=feat1_f1, feat2_f1=feat2_f1)
  data_f1$grpCol = NA
  data_f1$NCol = NA
  data_f1$grpCol[rownames(data_f1) %in% listCol[[1]]] = "grey50"
  data_f1$grpCol[rownames(data_f1) %in% listCol[[2]]] = "turquoise"
  data_f1$grpCol[rownames(data_f1) %in% listCol[[3]]] = "blue"
  data_f1$NCol[rownames(data_f1) %in% listCol[[1]]] = names(listCol)[1]
  data_f1$NCol[rownames(data_f1) %in% listCol[[2]]] = names(listCol)[2]
  data_f1$NCol[rownames(data_f1) %in% listCol[[3]]] = names(listCol)[3]
	pearson_f1 = cor.test(data_f1$feat1_f1, data_f1$feat2_f1)
	pearsonT_f1 = round(pearson_f1$estimate, 2)
  pearsonP_f1 = pearson_f1$p.value
  data_f1$density <- get_density(data_f1$feat1_f1,data_f1$feat2_f1, n=500)
  res.lm_f1 = lm(feat2_f1~feat1_f1)
  print(summary(res.lm_f1))
  summary(res.lm_f1)$r.squared
  GGPLOT_f1 = ggplot(data = data_f1) + geom_point(aes(x=feat1_f1, y=feat2_f1, color = density), shape = 4) +
           scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data_f1$feat2_f1~data$feat1_f1, se = FALSE)
  GGPLOT_colList_f1 = ggplot(data = data_f1) + geom_point(aes(x=feat1_f1, y=feat2_f1, color = NCol), shape = 4) +
           scale_color_manual(values=c("turquoise", "grey50")) + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data_f1$feat2_f1~data$feat1_f1, se = FALSE)
  ################### _f2 FILTER 2-3 ##################################################################
  feat1_f2 = feat1[names(feat1) %in% unique(c(listCol[[2]], listCol[[3]]))]
  feat2_f2 = feat2[names(feat2) %in% unique(c(listCol[[2]], listCol[[3]]))]
  feat2_f2 = as.matrix(feat2_f2)
  feat1_f2 = as.matrix(feat1_f2)
  feat2_f2 = as.matrix(feat2_f2)
	comRownames = Reduce(intersect, list(rownames(feat1_f2),rownames(feat2_f2)))
	feat1_f2 = feat1_f2[comRownames,1]
	feat2_f2 = feat2_f2[comRownames,1]
  data_f2=data.frame(feat1_f2=feat1_f2, feat2_f2=feat2_f2)
  data_f2$grpCol = NA
  data_f2$NCol = NA
  data_f2$grpCol[rownames(data_f2) %in% listCol[[1]]] = "grey50"
  data_f2$grpCol[rownames(data_f2) %in% listCol[[2]]] = "turquoise"
  data_f2$grpCol[rownames(data_f2) %in% listCol[[3]]] = "blue"
  data_f2$NCol[rownames(data_f2) %in% listCol[[1]]] = names(listCol)[1]
  data_f2$NCol[rownames(data_f2) %in% listCol[[2]]] = names(listCol)[2]
  data_f2$NCol[rownames(data_f2) %in% listCol[[3]]] = names(listCol)[3]
	pearson_f2 = cor.test(data_f2$feat1_f2, data_f2$feat2_f2)
	pearsonT_f2 = round(pearson_f2$estimate, 2)
  pearsonP_f2 = pearson_f2$p.value
  data_f2$density <- get_density(data_f2$feat1_f2,data_f2$feat2_f2, n=500)
  res.lm_f2 = lm(feat2_f2~feat1_f2)
  print(summary(res.lm_f2))
  summary(res.lm_f2)$r.squared
  GGPLOT_f2 = ggplot(data = data_f2) + geom_point(aes(x=feat1_f2, y=feat2_f2, color = density), shape = 4) +
           scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data_f2$feat2_f2~data$feat1_f2, se = FALSE)
  GGPLOT_colList_f2 = ggplot(data = data_f2) + geom_point(aes(x=feat1_f2, y=feat2_f2, color = NCol), shape = 4) +
           scale_color_manual(values=c( "blue",  "turquoise")) + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data_f1$feat2_f1~data$feat1_f1, se = FALSE)

	pdf(paste0(outpath,"scatter_",namefeat1,"_VS_", namefeat2, "_",info,".pdf"))
  # without x & y lim
    ########################  NO FILTER  ##########################################################
    print(GGPLOT)
    print(GGPLOT_colList)
    plot(feat1, feat2, pch=20, col=data$grpCol, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : ", pearsonT," - pv: ",pearsonP))
    abline(lm(feat2~feat1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
  	plot(feat1, feat2, pch=20, col=data$grpCol,ylim=yl, xlim=xl, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : R2: ",summary(res.lm)$r.squared))
    abline(lm(feat2~feat1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
    ################### _f1 _FILTER 1-2 ##################################################################
    print(GGPLOT_f1)
    print(GGPLOT_colList_f1)
    plot(feat1_f1, feat2_f1, pch=20, col=data_f1$grpCol, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : ", pearsonT_f1," - pv: ",pearsonP_f1))
    abline(lm(feat2_f1~feat1_f1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
    plot(feat1_f1, feat2_f1, pch=20, col=data_f1$grpCol, ylim=yl, xlim=xl, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : R2: ",summary(res.lm_f1)$r.squared))
    abline(lm(feat2_f1~feat1_f1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
    ################### _f2 FILTER 2-3 ##################################################################
    print(GGPLOT_f2)
    print(GGPLOT_colList_f2)
    plot(feat1_f2, feat2_f2, pch=20, col=data_f2$grpCol, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : ", pearsonT_f2," - pv: ",pearsonP_f2))
    abline(lm(feat2_f2~feat1_f2), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
    plot(feat1_f2, feat2_f2, pch=20, col=data_f2$grpCol, ylim=yl, xlim=xl, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : R2: ",summary(res.lm_f2)$r.squared))
    abline(lm(feat2_f2~feat1_f2), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))


#	add_loess(feat1, feat2, col="red")
	dev.off()

}



compareFunction_listCol2 = function(feat1, feat2,  namefeat1, namefeat2, listCol = NULL, namelistCol = NULL, vecCol = c("#839159", "#451c00"), logT = T, info = NULL, outpath, yl=c(0,15), xl=c(0,15)){
  ## LOG TRANS
  if(logT %in% T){
    feat1 = log2(feat1+1)
    feat2 = log2(feat2+1)
  }
  ## listCol
  # i.e. si listCol est un element a 1 vecteur, prendre le complement (La partie FILTER 1-2 sera alors inutile)
  listUnique = F
  if(length(listCol) > 1){
    listUnique = TRUE
  }
  if(length(listCol) %in% 1){
    listCol_temp = listCol[[1]]
    listCol[[1]] = names(feat1)[names(feat1) %ni% listCol_temp]
    listCol[[2]] = listCol_temp
    names(listCol) = c(paste0("1_NO_",names(listCol)[1]), paste0("2_",names(listCol)[1]))
  }
  ########################  NO FILTER  ##########################################################
  feat1 = as.matrix(feat1)
  feat2 = as.matrix(feat2)
	comRownames = Reduce(intersect, list(rownames(feat1),rownames(feat2)))
	feat1 = feat1[comRownames,1]
	feat2 = feat2[comRownames,1]
  data=data.frame(feat1=feat1, feat2=feat2)
  data$grpCol = NA
  data$NCol = NA
  data$grpCol[rownames(data) %in% listCol[[1]]] = vecCol[1]
  data$grpCol[rownames(data) %in% listCol[[2]]] = vecCol[2]
  data$NCol[rownames(data) %in% listCol[[1]]] = names(listCol)[1]
  data$NCol[rownames(data) %in% listCol[[2]]] = names(listCol)[2]
	pearson = cor.test(data$feat1, data$feat2)
	pearsonT = round(pearson$estimate, 2)
  pearsonP = pearson$p.value
  # if(summary(data$feat1)[2] == summary(data$feat1)[3]){
  #
  # }
  data$density <- get_density(data$feat1,data$feat2, n=500)
  res.lm = lm(feat2~feat1)
  print(summary(res.lm))
  summary(res.lm)$r.squared
  GGPLOT = ggplot(data = data) + geom_point(aes(x=feat1, y=feat2, color = density), shape = 1) +
           scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data$feat2~data$feat1, se = FALSE)
  GGPLOT_colList = ggplot(data = data) + geom_point(aes(x=feat1, y=feat2, color = NCol), shape = 1) +
           scale_color_manual(values=c(vecCol[1], vecCol[2])) + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data$feat2~data$feat1, se = FALSE)
  ################### FILTER 1-2 ##################################################################
  # i.e. pour avoir uniquement ceux dans union(listCol 1 et 2) si total est plus grand
  if(listUnique == T){
    feat1_f1 = feat1[names(feat1) %in% unique(c(listCol[[1]], listCol[[2]]))]
    feat2_f1 = feat2[names(feat2) %in% unique(c(listCol[[1]], listCol[[2]]))]
    feat2_f1 = as.matrix(feat2_f1)
    feat1_f1 = as.matrix(feat1_f1)
    feat2_f1 = as.matrix(feat2_f1)
    comRownames = Reduce(intersect, list(rownames(feat1_f1),rownames(feat2_f1)))
    feat1_f1 = feat1_f1[comRownames,1]
    feat2_f1 = feat2_f1[comRownames,1]
    data_f1=data.frame(feat1_f1=feat1_f1, feat2_f1=feat2_f1)
    data_f1$grpCol = NA
    data_f1$NCol = NA
    data_f1$grpCol[rownames(data_f1) %in% listCol[[1]]] = vecCol[1]
    data_f1$grpCol[rownames(data_f1) %in% listCol[[2]]] = vecCol[2]
    data_f1$NCol[rownames(data_f1) %in% listCol[[1]]] = names(listCol)[1]
    data_f1$NCol[rownames(data_f1) %in% listCol[[2]]] = names(listCol)[2]
    pearson_f1 = cor.test(data_f1$feat1_f1, data_f1$feat2_f1)
    pearsonT_f1 = round(pearson_f1$estimate, 2)
    pearsonP_f1 = pearson_f1$p.value
    data_f1$density <- get_density(data_f1$feat1_f1,data_f1$feat2_f1, n=500)
    res.lm_f1 = lm(feat2_f1~feat1_f1)
    print(summary(res.lm_f1))
    summary(res.lm_f1)$r.squared
    GGPLOT_f1 = ggplot(data = data_f1) + geom_point(aes(x=feat1_f1, y=feat2_f1, color = density), shape = 1) +
    scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
    #geom_smooth(method = "lm", formula=data_f1$feat2_f1~data$feat1_f1, se = FALSE)
    GGPLOT_colList_f1 = ggplot(data = data_f1) + geom_point(aes(x=feat1_f1, y=feat2_f1, color = NCol), shape = 1) +
    scale_color_manual(values=c(vecCol[1], vecCol[2])) + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
    #geom_smooth(method = "lm", formula=data_f1$feat2_f1~data$feat1_f1, se = FALSE)
  }
################# PLOT in PDF
  pdf(paste0(outpath,"scatter_",namefeat1,"_VS_", namefeat2,"_by_",namelistCol, "_",info,".pdf"))
  # without x & y lim
    ########################  NO FILTER  ##########################################################
    print(GGPLOT)
    print(GGPLOT_colList)
    plot(feat1, feat2, pch=0, col=data$grpCol, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : ", pearsonT," - pv: ",pearsonP))
    abline(lm(feat2~feat1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
  	plot(feat1, feat2, pch=0, col=data$grpCol,ylim=yl, xlim=xl, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : R2: ",summary(res.lm)$r.squared))
    abline(lm(feat2~feat1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
    ################### _f1 _FILTER 1-2 ##################################################################
  if(listUnique == T){
    print(GGPLOT_f1)
    print(GGPLOT_colList_f1)
    plot(feat1_f1, feat2_f1, pch=0, col=data_f1$grpCol, xlab = namefeat1, ylab = namefeat2,
      main=paste0(namefeat1 , " & ", namefeat2, " : ", pearsonT_f1," - pv: ",pearsonP_f1))
      abline(lm(feat2_f1~feat1_f1), col="red")
      abline(h=0, v=0, col="grey")
      lines(x = c(0,100), y = c(0,100))
      plot(feat1_f1, feat2_f1, pch=0, col=data_f1$grpCol, ylim=yl, xlim=xl, xlab = namefeat1, ylab = namefeat2,
        main=paste0(namefeat1 , " & ", namefeat2, " : R2: ",summary(res.lm_f1)$r.squared))
        abline(lm(feat2_f1~feat1_f1), col="red")
        abline(h=0, v=0, col="grey")
        lines(x = c(0,100), y = c(0,100))

  }

#	add_loess(feat1, feat2, col="red")
	dev.off()

}


compareFunction_listColRV = function(feat1, feat2,  namefeat1, namefeat2, listCol = NULL, info = NULL, outpath, yl=c(0,15), xl=c(0,15)){
  ########################  NO FILTER  ##########################################################
  feat1 = as.matrix(feat1)
  feat2 = as.matrix(feat2)
	comRownames = Reduce(intersect, list(rownames(feat1),rownames(feat2)))
	feat1 = feat1[comRownames,1]
	feat2 = feat2[comRownames,1]
  data=data.frame(feat1=feat1, feat2=feat2)
  data$grpCol = NA
  data$NCol = NA
  data$grpCol[rownames(data) %in% listCol[[1]]] = "grey50"
  data$grpCol[rownames(data) %in% listCol[[2]]] = "green"
  data$grpCol[rownames(data) %in% listCol[[3]]] = "red"
  data$NCol[rownames(data) %in% listCol[[1]]] = names(listCol)[1]
  data$NCol[rownames(data) %in% listCol[[2]]] = names(listCol)[2]
  data$NCol[rownames(data) %in% listCol[[3]]] = names(listCol)[3]
  data$NCol = as.character(data$NCol)
  data$NCol = factor(data$NCol, levels=unique(data$NCol))
	pearson = cor.test(data$feat1, data$feat2)
	pearsonT = round(pearson$estimate, 2)
  pearsonP = pearson$p.value
  data$density <- get_density(data$feat1,data$feat2, n=500)
  res.lm = lm(feat2~feat1)
  print(summary(res.lm))
  summary(res.lm)$r.squared
  GGPLOT = ggplot(data = data) + geom_point(aes(x=feat1, y=feat2, color = density), shape = 4) +
           scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data$feat2~data$feat1, se = FALSE)
  GGPLOT_colList = ggplot(data = data) + geom_point(aes(x=feat1, y=feat2, color = NCol), shape = 4) +
           scale_color_manual(values=c("grey50",  "green", "red" )) + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data$feat2~data$feat1, se = FALSE)
  ################### FILTER 1-2 ##################################################################
  feat1_f1 = feat1[names(feat1) %in% unique(c(listCol[[1]], listCol[[2]]))]
  feat2_f1 = feat2[names(feat2) %in% unique(c(listCol[[1]], listCol[[2]]))]
  feat2_f1 = as.matrix(feat2_f1)
  feat1_f1 = as.matrix(feat1_f1)
  feat2_f1 = as.matrix(feat2_f1)
	comRownames = Reduce(intersect, list(rownames(feat1_f1),rownames(feat2_f1)))
	feat1_f1 = feat1_f1[comRownames,1]
	feat2_f1 = feat2_f1[comRownames,1]
  data_f1=data.frame(feat1_f1=feat1_f1, feat2_f1=feat2_f1)
  data_f1$grpCol = NA
  data_f1$NCol = NA
  data_f1$grpCol[rownames(data_f1) %in% listCol[[1]]] = "grey50"
  data_f1$grpCol[rownames(data_f1) %in% listCol[[2]]] = "green"
  data_f1$grpCol[rownames(data_f1) %in% listCol[[3]]] = "red"
  data_f1$NCol[rownames(data_f1) %in% listCol[[1]]] = names(listCol)[1]
  data_f1$NCol[rownames(data_f1) %in% listCol[[2]]] = names(listCol)[2]
  data_f1$NCol[rownames(data_f1) %in% listCol[[3]]] = names(listCol)[3]
  data_f1$NCol = as.character(data_f1$NCol)
  data_f1$NCol = factor(data_f1$NCol, levels=unique(data_f1$NCol))
	pearson_f1 = cor.test(data_f1$feat1_f1, data_f1$feat2_f1)
	pearsonT_f1 = round(pearson_f1$estimate, 2)
  pearsonP_f1 = pearson_f1$p.value
  data_f1$density <- get_density(data_f1$feat1_f1,data_f1$feat2_f1, n=500)
  res.lm_f1 = lm(feat2_f1~feat1_f1)
  print(summary(res.lm_f1))
  summary(res.lm_f1)$r.squared
  GGPLOT_f1 = ggplot(data = data_f1) + geom_point(aes(x=feat1_f1, y=feat2_f1, color = density), shape = 4) +
           scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data_f1$feat2_f1~data$feat1_f1, se = FALSE)
  GGPLOT_colList_f1 = ggplot(data = data_f1) + geom_point(aes(x=feat1_f1, y=feat2_f1, color = NCol), shape = 4) +
           scale_color_manual(values=c("grey50", "green")) + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data_f1$feat2_f1~data$feat1_f1, se = FALSE)
  ################### _f2 FILTER 2-3 ##################################################################
  feat1_f2 = feat1[names(feat1) %in% unique(c(listCol[[2]], listCol[[3]]))]
  feat2_f2 = feat2[names(feat2) %in% unique(c(listCol[[2]], listCol[[3]]))]
  feat2_f2 = as.matrix(feat2_f2)
  feat1_f2 = as.matrix(feat1_f2)
  feat2_f2 = as.matrix(feat2_f2)
	comRownames = Reduce(intersect, list(rownames(feat1_f2),rownames(feat2_f2)))
	feat1_f2 = feat1_f2[comRownames,1]
	feat2_f2 = feat2_f2[comRownames,1]
  data_f2=data.frame(feat1_f2=feat1_f2, feat2_f2=feat2_f2)
  data_f2$grpCol = NA
  data_f2$NCol = NA
  data_f2$grpCol[rownames(data_f2) %in% listCol[[1]]] = "grey50"
  data_f2$grpCol[rownames(data_f2) %in% listCol[[2]]] = "green"
  data_f2$grpCol[rownames(data_f2) %in% listCol[[3]]] = "red"
  data_f2$NCol[rownames(data_f2) %in% listCol[[1]]] = names(listCol)[1]
  data_f2$NCol[rownames(data_f2) %in% listCol[[2]]] = names(listCol)[2]
  data_f2$NCol[rownames(data_f2) %in% listCol[[3]]] = names(listCol)[3]
  data_f2$NCol = as.character(data_f2$NCol)
  data_f2$NCol = factor(data_f2$NCol, levels=unique(data_f2$NCol))
	pearson_f2 = cor.test(data_f2$feat1_f2, data_f2$feat2_f2)
	pearsonT_f2 = round(pearson_f2$estimate, 2)
  pearsonP_f2 = pearson_f2$p.value
  data_f2$density <- get_density(data_f2$feat1_f2,data_f2$feat2_f2, n=500)
  res.lm_f2 = lm(feat2_f2~feat1_f2)
  print(summary(res.lm_f2))
  summary(res.lm_f2)$r.squared
  GGPLOT_f2 = ggplot(data = data_f2) + geom_point(aes(x=feat1_f2, y=feat2_f2, color = density), shape = 4) +
           scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data_f2$feat2_f2~data$feat1_f2, se = FALSE)
  GGPLOT_colList_f2 = ggplot(data = data_f2) + geom_point(aes(x=feat1_f2, y=feat2_f2, color = NCol), shape = 4) +
           scale_color_manual(values=c( "green", "red")) + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data_f1$feat2_f1~data$feat1_f1, se = FALSE)

	pdf(paste0(outpath,"scatter_",namefeat1,"_VS_", namefeat2, "_",info,".pdf"))
  # without x & y lim
    ########################  NO FILTER  ##########################################################
    print(GGPLOT)
    print(GGPLOT_colList)
    plot(feat1, feat2, pch=4, col=data$grpCol, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : ", pearsonT," - pv: ",pearsonP))
    abline(lm(feat2~feat1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
  	plot(feat1, feat2, pch=4, col=data$grpCol,ylim=yl, xlim=xl, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : R2: ",summary(res.lm)$r.squared))
    abline(lm(feat2~feat1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
    ################### _f1 _FILTER 1-2 ##################################################################
    print(GGPLOT_f1)
    print(GGPLOT_colList_f1)
    plot(feat1_f1, feat2_f1, pch=4, col=data_f1$grpCol, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : ", pearsonT_f1," - pv: ",pearsonP_f1))
    abline(lm(feat2_f1~feat1_f1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
    plot(feat1_f1, feat2_f1, pch=4, col=data_f1$grpCol, ylim=yl, xlim=xl, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : R2: ",summary(res.lm_f1)$r.squared))
    abline(lm(feat2_f1~feat1_f1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
    ################### _f2 FILTER 2-3 ##################################################################
    print(GGPLOT_f2)
    print(GGPLOT_colList_f2)
    plot(feat1_f2, feat2_f2, pch=4, col=data_f2$grpCol, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : ", pearsonT_f2," - pv: ",pearsonP_f2))
    abline(lm(feat2_f2~feat1_f2), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
    plot(feat1_f2, feat2_f2, pch=4, col=data_f2$grpCol, ylim=yl, xlim=xl, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : R2: ",summary(res.lm_f2)$r.squared))
    abline(lm(feat2_f2~feat1_f2), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))


#	add_loess(feat1, feat2, col="red")
	dev.off()

}


compareFunction = function(feat1, feat2,  namefeat1, namefeat2, logT=F, info = NULL, outpath, yl=c(-1,1), xl=c(-1,1)){
  if(logT %in% T){
    feat1 = log2(feat1+1)
    feat2 = log2(feat2+1)
  }
  feat1 = as.matrix(feat1)
  feat2 = as.matrix(feat2)
	comRownames = Reduce(intersect, list(rownames(feat1),rownames(feat2)))
	feat1 = feat1[comRownames,1]
	feat2 = feat2[comRownames,1]
  data=data.frame(feat1=feat1, feat2=feat2)
	pearson = cor.test(data$feat1, data$feat2)
	pearsonT = round(pearson$estimate, 2)
  pearsonP = pearson$p.value
  data$density <- get_density(data$feat1,data$feat2, n=500)
  res.lm = lm(feat2~feat1)
  print(summary(res.lm))
  summary(res.lm)$r.squared
  GGPLOT = ggplot(data = data) + geom_point(aes(x=feat1, y=feat2, color = density)) + xlab(namefeat1) + ylab(namefeat2) +
           scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data$feat2~data$feat1, se = FALSE)
	pdf(paste0(outpath,"scatter_",namefeat1,"_VS_", namefeat2, "_",info,".pdf"))
  # without x & y lim
    print(GGPLOT)
    textplot("R_squared")
    textplot(summary(res.lm)$r.squared)
    textplot("P-val Intercept")
    textplot(summary(res.lm)$coefficients[7])
    textplot("P-val feat1 ")
    textplot(summary(res.lm)$coefficients[8])
  	plot(feat1, feat2, pch=20, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : ", pearsonT," - pv: ",pearsonP))
    abline(lm(feat2~feat1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
  	plot(feat1, feat2, pch=20,ylim=yl, xlim=xl, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : R2: ",summary(res.lm)$r.squared))
    abline(lm(feat2~feat1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
#	add_loess(feat1, feat2, col="red")
	dev.off()

}



##################

# ADD by REFKA 2021
#need : library("ggforce")
#1)cette fonction permet ploter / scatter de deux vecteurs
# exemple de vecteur : 

# str(Q_TSS_PW_R2_dwn_2K_f)
        # Named num [1:6967] 301 577 608 249 505 ...
        # - attr(*, "names")= chr [1:6967] "FBgn0000017.1" "FBgn0000018.1" "FBgn0000032.1" "FBgn0000042.1" ...

#2)de colorer des groupes de gènes sités dans la liste "LIST_clusters" 
#chaque groupe de gènes est mis sous forme de vecteur de GN 
#> head(names(Q_TSS_PW_R2_dwn_2K_f))
#[1] "FBgn0000017.1" "FBgn0000018.1" "FBgn0000032.1" "FBgn0000042.1"
#[5] "FBgn0000043.1" "FBgn0000052.1"

# exemple : LIST_clusters
#LIST_clusters=list(UP_10_Q_H3K36me3_2C4=UP_10_Q_H3K36me3_2C4,
#                    DN_10_Q_H3K36me3_2C4=DN_10_Q_H3K36me3_2C4)
#> str(LIST_clusters)
#List of 2
# $ UP_10_Q_H3K36me3_2C4: chr [1:697] "FBgn0020235.1" "FBgn0260990.1" "FBgn0039697.1" "FBgn0027518.1" ...
# $ DN_10_Q_H3K36me3_2C4: chr [1:697] "FBgn0262035.1" "FBgn0002543.1" "FBgn0025111.1" "FBgn0259168.1" ...
#> 


 
compare_cluster_Function = function(feat1, feat2,LIST_clusters,namefeat1, namefeat2,pdf_height=NULL, pdf_width=NULL, logT=F, info = NULL, outpath, yl=c(-1,1), xl=c(-1,1)){
  if(logT %in% T){
    feat1 = log2(feat1+1)
    feat2 = log2(feat2+1)
  }
  feat1 = as.matrix(feat1)
  feat2 = as.matrix(feat2)
	comRownames = Reduce(intersect, list(rownames(feat1),rownames(feat2)))
	feat1 = feat1[comRownames,1]
	feat2 = feat2[comRownames,1]
  data=data.frame(feat1=feat1, feat2=feat2)
	pearson = cor.test(data$feat1, data$feat2)
	pearsonT = round(pearson$estimate, 2)
  pearsonP = pearson$p.value
  data$density <- get_density(data$feat1,data$feat2, n=500)
  res.lm = lm(feat2~feat1)
  print(summary(res.lm))
  summary(res.lm)$r.squared
  data$clusters=0
  l=1
  while (l <= length(LIST_clusters))
  {
      print(names(LIST_clusters[l]))
      for(e in (1:dim(data)[1])){
          if (rownames(data)[e] %in% LIST_clusters[[l]])
              {data$clusters[e] <-names(LIST_clusters[l])}
  }
  l=l+1
  }
  GGPLOT = ggplot(data = data) + geom_point(aes(x=feat1, y=feat2, color = density)) + xlab(namefeat1) + ylab(namefeat2) +
           scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
  GGPLOT_2 = ggplot(data = data) + geom_point(aes(x=feat1, y=feat2, color = clusters)) + xlab(namefeat1) + ylab(namefeat2)+ scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)

  #GGPLOT_3= ggplot(data = data) + geom_mark_ellipse(aes(fill = clusters,label = clusters),expand = unit(0.5,"mm"),label.buffer = unit(-5, 'mm'))+ theme(legend.position = "none") +geom_point(aes(x=feat1, y=feat2, color = clusters)) + xlab(namefeat1) + ylab(namefeat2)+ scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
  GGPLOT_3= ggplot(data = data) + geom_mark_ellipse(aes(x=feat1, y=feat2,fill = clusters,label = clusters),expand = unit(1,"mm"),label.buffer = unit(-1, 'mm'))+geom_point(aes(x=feat1, y=feat2,color = clusters)) + xlab(namefeat1) + ylab(namefeat2)  


           #geom_smooth(method = "lm", formula=data$feat2~data$feat1, se = FALSE)
	pdf(paste0(outpath,"scatter_",namefeat1,"_VS_", namefeat2, "_",info,".pdf"), height=pdf_height, width=pdf_height)
  # without x & y lim
    print(GGPLOT)
    print(GGPLOT_2)
    print(GGPLOT_3)
    textplot("R_squared")
    textplot(summary(res.lm)$r.squared)
    textplot("P-val Intercept")
    textplot(summary(res.lm)$coefficients[7])
    textplot("P-val feat1 ")
    textplot(summary(res.lm)$coefficients[8])
  	plot(feat1, feat2, pch=20, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : ", pearsonT," - pv: ",pearsonP))
    abline(lm(feat2~feat1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
  	plot(feat1, feat2, pch=20,ylim=yl, xlim=xl, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : R2: ",summary(res.lm)$r.squared))
    abline(lm(feat2~feat1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
#	add_loess(feat1, feat2, col="red")
	dev.off()

}




























compareFunctionNOPDF = function(feat1, feat2,  yl=c(-1,1), xl=c(-1,1)){
  namefeat1 = "coco"
  namefeat = "cucu"
  feat1 = as.matrix(feat1)
  feat2 = as.matrix(feat2)
	comRownames = Reduce(intersect, list(rownames(feat1),rownames(feat2)))
	feat1 = feat1[comRownames,1]
	feat2 = feat2[comRownames,1]
  data=data.frame(feat1=feat1, feat2=feat2)
	pearson = cor.test(data$feat1, data$feat2)
	pearsonT = round(pearson$estimate, 2)
  pearsonP = pearson$p.value
  data$density <- get_density(data$feat1,data$feat2, n=500)
  GGPLOT = ggplot(data = data) + geom_point(aes(x=feat1, y=feat2, color = density)) +
           scale_color_viridis() + scale_x_continuous(limits = xl) + scale_y_continuous(limits = yl)
           #geom_smooth(method = "lm", formula=data$feat2~data$feat1, se = FALSE)
	# pdf(paste0(outpath,"scatter_",namefeat1,"_VS_", namefeat2,".pdf"))
  # without x & y lim
    print(GGPLOT)
  	plot(feat1, feat2, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : ", pearsonT," - pv: ",pearsonP))
    abline(lm(feat2~feat1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
  	plot(feat1, feat2,ylim=yl, xlim=xl, xlab = namefeat1, ylab = namefeat2,
    main=paste0(namefeat1 , " & ", namefeat2, " : ", pearsonT," - pv: ",pearsonP))
    abline(lm(feat2~feat1), col="red")
    abline(h=0, v=0, col="grey")
    lines(x = c(0,100), y = c(0,100))
#	add_loess(feat1, feat2, col="red")
	# dev.off()

}
