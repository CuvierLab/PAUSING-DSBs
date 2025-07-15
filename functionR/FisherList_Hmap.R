
# Single Cell RNA seq project
# Cuvier''s Team
# Schaak - Heurteau - Depierre*
# 2017

# modified function from pascal/gael
# fisher_namesList :
# fisher function adapted to test intersection genes list against genes list
# list 1 and 2 in input are list of char vector, total correspond to the total number of genes, statistical "univers"

# Hmap_pval :
# plot pval from fisher_namesList
# input :
# mat = output return from fisher_namesList
# graph = "/path/to/output/mygraph.pdf"
# bool = T to exprt pdf grap according to graph input, if False, plot directed in X11
# threshold=1e-8 correspond to pval threshold for the beginning of the blue colorRamp

library("pheatmap")


# List_genesDOM_LIST_UP_K27_SHARED_GN = fisher_namesList(List_genes_DOM_f, LIST_UP_K27_SHARED_GN, 6831)
# # Sel2_DESeq2__ExprDiff_CUTLL1_GN_splitted = fisher_namesList(Sel2_DESeq2,ExprDiff_CUTLL1_GN_splitted,length(rownames(ExprDiff_CUTLL1)))
#

# list1 = List_genes_DOM_f
# list2 = LIST_UP_K27_SHARED_GN
# total = 6831


fisher_namesList=function(list1, list2, total, test.side = "greater")
 {
 require(GenomicRanges)
 ## create matrices to fill
 m_pval=matrix(nrow=length(list1),ncol=length(list2))
 m_or=matrix(nrow=length(list1),ncol=length(list2))
 m_eff=matrix(nrow=length(list1),ncol=length(list2))
 m_pct_list1=matrix(nrow=length(list1),ncol=length(list2))
 m_pct_list2=matrix(nrow=length(list1),ncol=length(list2))

 colnames(m_pval)=names(list2)
 rownames(m_pval)=names(list1)
 colnames(m_or)=names(list2)
 rownames(m_or)=names(list1)
 colnames(m_eff)=names(list2)
 rownames(m_eff)=names(list1)
 colnames(m_pct_list1)=names(list2)
 rownames(m_pct_list1)=names(list1)
 colnames(m_pct_list2)=names(list2)
 rownames(m_pct_list2)=names(list1)

 ## fill matrices with either fisher pval, fisher odds, effectives or percentage
 for (i in 1:length(list1)){
   liste1=list1[[i]]
    for (j in 1:length(list2)){
      liste2=list2[[j]]
      inter=length(Reduce(intersect, list(liste1,liste2)))
      m_pval[i,j]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(liste2)-inter,max(total-length(liste1)-length(liste2)+inter,0)),ncol=2),alternative=test.side)$p.value,2)
      m_or[i,j]=signif(fisher.test(matrix(c(inter,length(liste1)-inter,length(liste2)-inter,max(total-length(liste1)-length(liste2)+inter,0)),ncol=2),alternative=test.side)$estimate,2)
      m_eff[i,j]=inter
      m_pct_list1[i,j] = inter/length(liste1)
      m_pct_list2[i,j] = inter/length(liste2)

    }
 }

## get list of matrices to return
 l=list()
 l[[1]]=m_pval
 l[[2]]=m_or
 l[[3]]=m_eff
 l[[4]]=m_pct_list1
 l[[5]]=m_pct_list2
 names(l)=c("p.value", "odds.ratio", "effectives", "pct_list1", "pct_list2")
 return(l)
 }


########################################################################################################################################################################
###    FROM MATRICES RETURN FROM fisher_namesList(), PLOT HEATMAP
########################################################################################################################################################################
 Hmap_pval=function(mat, values_ordonnée, values_abscisse, graph,bool, threshold=1e-8, odds_scale = c(0,2)){
   if (bool==T){pdf(graph)} # should we plot in pdf -> yes if bool = T
   # Prepare pval matrix
   mat_pv=mat$p.value
   mat_pv=mat_pv+10^-300 # to avoid pval = 0
   # ADDED on 23/10/2019
   mat_pv[mat_pv > 0.05] = 1 # set pval to 1 if >0.05 (basis sup threshold)
   if (threshold>1 | threshold<=0) (stop("threshold must be in ]0,1]"))
   mat_pv[mat_pv < threshold] = threshold # set pval to threshold if threshold (basis inferior threshold)
   mat_pv=log10(mat_pv) # transform pval to log10

# ####### OLD WAY TO PLOT PVALUES
#    if (bool==T){postscript(graph)}
#    # mat$p.value
#    ramp <- colorRamp(c("blue", "white"))
#    palette(rgb(ramp((200:0)/200),max=255))
#    plot.new()
#    plot.window(c(-0.4,1.2),c(-0.2,1.2))
#
#    long=ncol(mat_pv)
#    long2=nrow(mat_pv)
#    for (j in 1:ncol(mat_pv)){
#      for (k in 1:nrow(mat_pv)){
#        rect((j-1)/long,(k-1)/long2,j/long,k/long2,col=((mat_pv[k,j])*(mat_pv[k,j]>logt)*(200/logt))+200*(mat_pv[k,j]<logt))
#        text(-0.13,(k-1/2)/long2,rownames(mat_pv)[k],col='black',cex=0.7)
#      }
#      text((j-1/2)/long,1.13,colnames(mat_pv)[j],col='black',srt=45,cex=0.7)
#    }
#    palette(rgb(ramp((200:0)/200),max=255))
#    plot.new()
#    plot.window(c(-0.4,1.2),c(-0.2,1.2))
#    for (i in 0:200){
#      rect((200-i)/200,0,(199-i)/200,1,col=(i),border=i)
#    }
##########################
####   PLOT PVALUE   #####
##########################
   par(mfrow=c(2,1))
   logt=log10(threshold)
   pval_info = matrix(c(1, threshold))
   rownames(pval_info) = c("white", "blue")
   textplot(pval_info)
   textplot(mat$p.value)
   title("PVALUES FISHER EXACT TEST")

   par(mfrow=c(1,1))

  v_breaks = c(seq(min(mat_pv), 0.1, length.out=25), seq(0.09, max(mat_pv), length.out=25))
  v_cols = c(colorRampPalette(c("blue", "white"))(24), "white", colorRampPalette(c("white", "white"))(24))
    pheatmap(mat_pv, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols,  display_numbers = mat$p.value, fontsize_number = 18, main="PVALUES FISHER EXACT TEST")
  v_breaks = c(seq(log10(threshold), 0.1, length.out=25), seq(0.09, max(mat_pv), length.out=25))
  v_cols = c(colorRampPalette(c("blue", "white"))(24), "white", colorRampPalette(c("white", "white"))(24))
    pheatmap(mat_pv, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, main="PVALUES FISHER EXACT TEST")
    pheatmap(mat_pv, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols,  display_numbers = mat$p.value, fontsize_number = 18, main="PVALUES FISHER EXACT TEST")
    pheatmap(mat_pv, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols,  display_numbers = mat$effectives, fontsize_number = 18, main="PVALUES FISHER EXACT TEST - EFFECTIVES")
    pheatmap(mat_pv, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols,  display_numbers = round(mat$pct_list1,2), fontsize_number = 18, main="% of list2 in list1 / i.e. column feat. in row feat.")
    pheatmap(mat_pv, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols,  display_numbers = round(mat$pct_list2,2), fontsize_number = 18, main="% of list1 in list2 / i.e. row feat. in comumn feat.")

####### OLD WAY TO PLOT ODDS RATIO, REMOVED ON 18/06/2020
   #  mat_or = mat$odds.ratio
   # mat_or = (log(mat_or)+1)*100
   # mat_or[mat_or > 200] = 200
   # mat_or[mat_or < 0] = 0
   # ramp <- colorRamp(c("darkblue", "white", "yellow"))
   # palette(rgb(ramp((200:0)/200),max=255))
   # plot.new()
   # plot.window(c(-0.4,1.2),c(-0.2,1.2))
   # long=ncol(mat_or)
   # long2=nrow(mat_or)
   # for (j in 1:ncol(mat_or)){
   #   for (k in 1:nrow(mat_or)){
   #     rect((j-1)/long,(k-1)/long2,j/long,k/long2,col=((mat_or[k,j])))
   #     # rect((j-1)/long,(k-1)/long2,j/long,k/long2,col=((mat_or[k,j])*(mat_or[k,j]>ThOdds)*(200/ThOdds))+200*(mat_or[k,j]<ThOdds))
   #     text(-0.13,(k-1/2)/long2,rownames(mat_or)[k],col='black',cex=0.7)
   #   }
   #   text((j-1/2)/long,1.13,colnames(mat_or)[j],col='black',srt=45,cex=0.7)
   # }
   # palette(rgb(ramp((200:0)/200),max=255))
   # plot.new()
   # plot.window(c(-0.4,1.2),c(-0.2,1.2))
   # for (i in 0:200){
   #   rect((200-i)/200,0,(199-i)/200,1,col=(i),border=i)
   #  }
#############################
####   PLOT ODDS RATIO   ####
#############################
    textplot(mat$odds.ratio) # PLOT ODDS RATIO TABLE
    title("ODDS RATIO FISHER EXACT TEST")

    plot.new()
    ## RAW ODDS RATIO
    M = mat$odds.ratio
    v_breaks = c(seq(min(M), 0.9, length.out=25), seq(1.1, max(M), length.out=25))
    v_cols = c(colorRampPalette(c("yellow", "white"))(24), "white", colorRampPalette(c("white", "darkblue"))(24))
    pheatmap(M, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, main="ODDS RATIO FISHER EXACT TEST")
    pheatmap(M, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, display_numbers = mat$odds.ratio, fontsize_number = 18, main="ODDS RATIO FISHER EXACT TEST")
    pheatmap(M, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, display_numbers = mat$p.value, fontsize_number = 18, main="ODDS RATIO FISHER EXACT TEST - pval numbers displayed")
    pheatmap(M, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols,  display_numbers = mat$effectives, fontsize_number = 18, main="ODDS RATIO FISHER EXACT TEST - EFFECTIVES")


    ## LIMITED ODDS RATIO to borders defined with "odds_scale" param
    M_v2 = M
    M_v2[M_v2 > odds_scale[2]] = odds_scale[2]
    v_breaks = c(seq(odds_scale[1], 0.8, length.out=25), seq(1.2, odds_scale[2], length.out=25))
    v_cols = c(colorRampPalette(c("yellow", "white"))(24), "white", colorRampPalette(c("white", "darkblue"))(24))
    pheatmap(M_v2, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, main="ODDS RATIO FISHER EXACT TEST")
    pheatmap(M_v2, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, display_numbers = mat$odds.ratio, fontsize_number = 18, main="ODDS RATIO FISHER EXACT TEST")
    pheatmap(M_v2, cluster_rows=F, cluster_cols=F, breaks= v_breaks, color=v_cols, display_numbers = mat$p.value, fontsize_number = 18, main="ODDS RATIO FISHER EXACT TEST - pval numbers displayed")


    # textplot(lapply(values_ordonnée, length))
    # textplot(lapply(values_abscisse, length))
    if(is.character(values_ordonnée[[1]])){ #if the split of genes nmaes is not done on quantitatif value, plot effectif of each feature
      barplot(unlist(lapply(values_ordonnée, length)), cex.names=0.5)
    }else{
      boxplot(values_ordonnée, main=deparse(substitute(values_ordonnée)), outline=F)
      barplot(unlist(values_ordonnée), main=deparse(substitute(values_ordonnée)))
      plot(density(unlist(values_ordonnée)), main=deparse(substitute(values_ordonnée)))
      hist(unlist(values_ordonnée), 50, main=deparse(substitute(values_ordonnée)))
    }
    if(is.character(values_abscisse[[1]])){ #if the split of genes nmaes is not done on quantitatif value, plot effectif of each feature
      barplot(unlist(lapply(values_abscisse, length)), cex.names=0.5)
    }else{
      boxplot(values_abscisse, main=deparse(substitute(values_abscisse)), outline=F)
      barplot(unlist(values_abscisse), main=deparse(substitute(values_abscisse)))
      plot(density(unlist(values_abscisse)), main=deparse(substitute(values_abscisse)))
      hist(unlist(values_abscisse), 50, main=deparse(substitute(values_abscisse)))
    }
   if (bool==T){
     dev.off()
   }
 }


# end
