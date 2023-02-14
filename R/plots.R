#' Plot of cell occurrences
#'
#' Plot of occurrences of each cell(cell type) in the whole spatial transcriptome.
#'
#' @param nums Results of \code{deconvoluting}.
#'
#' @return This function will return a plot showing the occurrence of cells as well as a global variable named "occurred.cells"
#' showing the names of the cells occurred.
#'
#' @export
cell.occur = function(nums){
  nums = as.matrix(nums)
  counts = nums>0
  occurs = apply(counts,1,sum)
  occurred.cells <<- names(occurs[which(occurs>0)])
  occurs = occurs[order(occurs,decreasing = T)]
  ncells = length(occurs)
  options(warn = -1)
  totoccurs = min(which(occurs==0))-1
  p <- ggplot2::ggplot(data=NULL,aes(x=1:ncells,y=occurs))+
    ggplot2::geom_point(size=0.2)+
    ggplot2::geom_line()+
    ggplot2::geom_hline(aes(yintercept=0),color="red")+
    ggplot2::labs(x=" ",y="occurrence")+
    ggplot2::theme_classic()
  if(totoccurs != Inf){
    p <- p +
      ggplot2::geom_vline(aes(xintercept=totoccurs),color="red",linetype="dashed")+
      ggplot2::scale_x_discrete(limits=c(1,totoccurs))
  }
  else{totoccurs = ncells}
  options(warn = 1)
  print(paste0(totoccurs," in ",ncells," cells are used."))
  return(p)
}

#' #' Heatmap of any correlation
#' #'
#' #' This function
#' #'
#' #' @param x The data whose correlation will be calculated and visualized. One column represents one variable.
#' #' @param group
#' #' @param method Method of calculating correlation coefficient.
#' #'
#' #' @return A global variable named "corr", represents for correlation of cell type
#' #' proportion across spatial locations and its visualization.
#' #'
#' #' @export
#' corr.heatmap = function(x, group, method = c("pearson","spearman")){
#'   x = as.matrix(x)
#'   corr <<- Hmisc::rcorr(x,type=method)[["r"]]
#'   corrtp = reshape2::melt(corr)
#'   ggplot(corrtp,aes(x=Var1,y=Var2,fill=value))+
#'     geom_tile()+
#'     scale_fill_gradient2(low="White", high="Blue")+
#'     labs(x=NULL,y=NULL,title = "correlation")+
#'     theme_bw(base_size = 15)+
#'     theme(axis.text.x = element_blank(),
#'           axis.text.y = element_blank())
#' }

#' Converting single cell results to type results
#'
#' This function converts single-cell abundance to cell-type abundance.
#'
#' @param nums The results of single-cell resolution
#' @param annotations The annotation of the single cells.
#'
#' @export
sc2type = function(nums,annotations){
  ord = order(annotations[,2])
  annotations = annotations[ord,]
  tab = table(annotations[,2])
  nums = nums[ord,]
  ntypes = dim(tab)
  typenums = matrix(nrow = ntypes, ncol = ncol(nums))
  rownames(typenums) = dimnames(tab)[[1]]
  colnames(typenums) = colnames(nums)
  begin=1
  for(i in 1:ntypes){
    if(tab[i]==1)typenums[i,]=nums[begin,]
    else{
      end = begin+tab[i]-1
      names(end) = NULL
      typenums[i,] = apply(nums[begin:end,],2,sum)
    }
    begin = end+1
  }
  return(typenums)
}

#' Plot of weights for each cell type
#'
#' Plot of weights for each cell type
#'
#' @param nums Results of \code{deconvoluting}.
#' @param cell.type Whether to merge \code{nums} by cell type.
#' @param annotations If \code{cell.type = F}, this should be the annotation of the single cells.
#' @param cell.names Name of cells (or cell types) to be plotted. Default is all.
#' @param coords Coordinates of spatial spots. \code{rownames} should be the barcode(name) of each spot
#' and \code{colnames} should be "x" and "y".
#' @param name The name of this pdf file.
#'
#' @export
cell.type.weight = function(nums,cell.type,annotations,cell.names=NULL,coords,name){
  if(cell.type == T) tnums = sc2type(nums,annotations)
  else tnums = nums
  ntnums = apply(tnums,2,function(x){if(sum(x)!=0)x=x/sum(x);return(x)})
  if(length(cell.names)!=0) ntnums = ntnums[cell.names,]
  ntypes = nrow(ntnums)
  plots = vector(mode = "list",length = ntypes)
  celltypes = rownames(ntnums)
  for(i in 1:ntypes){
    weight = data.frame(ntnums[i,])
    toplot = cbind.data.frame(coords,weight)
    colnames(toplot) = c("x","y","weight")
    plots[[i]] = ggplot(toplot)+
      geom_point(aes(x=x,y=y,color=weight),size=1)+
      scale_color_gradient(limits=c(0,1),low="#F5F5F5",high="blue")+
      labs(title = celltypes[i])+
      theme_classic()
  }
  pdf(name)
  invisible(lapply(plots, print))
  dev.off()
}

#' Spatial expression
#'
#' Visualization of spatial expression of some genes.
#'
#' @param st Spatial transcriptomics.
#' @param coords Coordinates of spatial spots.
#' @param gene.list List of genes.
#'
#' @return A pdf file named "spatial_genes".
#'
#' @export
spatial_gene = function(st,coords,gene.list){
  if(sum(! gene.list %in% rownames(st))>0){
    warning("Some genes are not in st. Such genes are ignored.")
  }
  genelist = intersect(gene.list,rownames(st))
  ngenes = length(gene.list)
  plots = vector(mode = "list",length = ngenes)
  for(i in 1:ngenes){
    toplot = cbind.data.frame(coords,st[genelist[i],])
    colnames(toplot) = c("x","y","expression")
    plots[[i]] = ggplot(toplot)+
      geom_point(aes(x=x,y=y,color=expression),size=1)+
      scale_color_gradient(low="#F5F5F5",high="blue")+
      labs(title = genelist[i])+
      theme_classic()
  }
  pdf("spatial_genes")
  invisible(lapply(plots, print))
  dev.off()
}

#' Visualization of colocalization
#'
#' This function draws a network to show colocalization of cells (cell types).
#'
#' @param corr Correlation matrix returned by \code{coloc.corr}
#' @param thresh Threshold of correlation coefficient. Value large than this value will be regarded as an edge.
#' @param cell.type Whether the result is by cell-type (\code{TRUE}) or by single-cell (\code{FALSE}).
#' @param annotations If \code{cell.type = T}, this should be the names of the cell types.
#' If \code{cell.type = F}, this should be the annotation of the single cells.
#' @param ntypes Number of cell types.
#' @param color The colors of cell types. Default is \code{rainbow(ntypes)}.
#' @param name The name of this pdf file.

coloc.network=function(corr,thresh,cell.type,annotations,ntypes,color=NULL,name){
  diag(corr)=0
  g = graph_from_adjacency_matrix(corr>thresh,mode="undirected",weighted=T)
  V(g)$annotations = annotations[,2]
  if(length(color)==0)color = grDevices::rainbow(ntypes)
  V(g)$color = color[factor(V(g)$annotations)]
  deg = degree(g,mode = "all")
  ncells = nrow(corr)
  if(!dir.exists(dir))dir.create(dir)
  pdf(file.path(dir, name))
  plot(g,vertex.frame.color="white",vertex.size=100/ncells*deg,vertex.label.cex=0.2,layout=layout_with_fr,edge.width=0.2)
  legend(x=-1.5,y=1.5,levels(factor(V(g)$annotations)),pch=21,col="#777777",pt.bg=color,cex=0.5)
  dev.off()
}

#' This is the CARD function
#' @export
spatial.piechart=function(nums,coords,colors=NULL,title=NULL){
  ntnums = apply(nums,2,function(x){if(sum(x)!=0)x=x/sum(x);return(x)})
  data = cbind.data.frame(t(ntnums),coords)
  ntypes = nrow(ntnums)
  if(length(colors)==0) colors = grDevices::rainbow(ntypes)
  ct.select = rownames(ntnums)
  occurred = which(apply(ntnums,1,sum)!=0)
  ggplot() +
    geom_scatterpie(aes(x=x, y=y,r = 0.52),data=data,
                                              cols=ct.select,color=NA) + coord_fixed(ratio = 1) +
                     scale_fill_manual(values = colors[occurred])+
                     labs(title=title)+
                     theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                           panel.background = element_blank(),
                           plot.background = element_blank(),
                           panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                           axis.text =element_blank(),
                           axis.ticks =element_blank(),
                           axis.title =element_blank(),
                           legend.title=element_text(size = 12,face="bold"),
                           legend.text=element_text(size = 10),
                           legend.key = element_rect(colour = "transparent", fill = "white"),
                           legend.key.size = unit(0.45, 'cm'),
                           strip.text = element_text(size = 12,face="bold"),
                           legend.position="bottom")+
                     guides(fill=guide_legend(title="Cell Type"))
}

#' Corrected profile v.s. original profile
#'
#' Comparison of the corrected and original profile.
#'
#' @param ref Original reference.
#' @param ests Corrected expression profile.
#' @param barcodes (Optional) the barcodes of the cells whose profile is to compared. Default is all.
#'
#' @return A pdf file named "profile_comparison".
#'
#' @export
profile_comparison = function(ref,ests,barcodes=NULL){
  ref = ref+0.5
  ref = apply(ref,2,function(x){x/sum(x)*1e06})
  ref = ref[rownames(ests),colnames(ests)]
  if(is.null(barcodes)){
    barcodes=colnames(ests)
  }
  ncells = length(barcodes)
  plots = vector(mode = "list",length = ncells)
  for(i in 1:ncells){
    geneslabel = correctedgenes(ref[,i],ests[,i])
    toplot = cbind.data.frame(log10(ref[,i]),log10(ests[,i]),rownames(ests))
    colnames(toplot) = c("x","y","genes")
    plots[[i]] = ggplot()+
      geom_point(data=toplot,aes(x=x,y=y,alpha=0.8))+
      ggrepel::geom_text_repel(data=toplot[geneslabel,],aes(x=x,y=y,
                                                   label=genes,
                          fontface="italic"),size=3)+
      labs(x="origninal",y="corrected",title=barcodes[i])+
      geom_abline(slope = 1,intercept=0,linetype="dashed",color="red")+
      theme_classic()
  }
  pdf("profile_comparison.pdf")
  invisible(lapply(plots, print))
  dev.off()
}

#
correctedgenes = function(ref1,est1,top=20,thresh=50){
  diff = est1-ref1
  ord = order(diff,decreasing = T)
  tops = ord[1:top]
  totake = tops[which(tops>thresh)]
  if(length(totake)==0)totake =tops[which(tops>0)]
  return(totake)
}


#' Pie chart of a single spot
#'
#' Pie chart of a single spot.
#'
#' @param num One column of the results.
#' @param colors The colors of cell types. Default is \code{rainbow(ntypes)}.
#' @param title Title of this plot.
#'
#' @export
spot.pie = function(num,colors=NULL,title){
  ntypes = length(num)
  if(length(colors)==0)colors = grDevices::rainbow(ntypes)
  ntnum = num/sum(num)
  toplot = data.frame(types = names(ntnum), prop = ntnum)
  colors = colors[toplot$prop!=0]
  toplot = toplot[toplot$prop!=0,]
  noccurtypes = nrow(toplot)
  toplot$label = paste(round(toplot$prop*100,2),"%",sep="")

  sum = 0
  for(i in 1:noccurtypes){
    toplot$ypos[i] = 180*sum + 90*toplot$prop[i]
    sum = sum + toplot$prop[i]
  }
  toplot$xpos = 1.4
  toplot$z = 1

  ggplot(toplot)+
    geom_segment(aes(x = xpos+0.05,y = ypos,xend = xpos+0.1,yend = ypos,alpha = 0.7),color = colors,size = 1.5)+
    geom_col(aes(z,rev(prop)*180,fill = types,alpha = 0.8),size = 0.5)+
    #geom_point(aes(x = xpos+0.1,y = ypos),size = 1,color = colors,alpha = 0.5)+
    coord_polar(theta = "y")+
    scale_fill_manual(values = rev(colors),guide = NULL)+
    #geom_text(aes(x = xpos+0.15,y = ypos,label = label),size = 5,adj = 0.5)+
    scale_alpha(guide = NULL)+
    labs(title=title)+
    theme(panel.background = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "top")

}

#' Showing cells of interest
#'
#' Summary of the result, for users to identify cells of interest.
#'
#' @param nums The result of \code{deconvoluting}
#' @param occurred.cells Barcodes of cells whose abundance is not zero. This value can be generated by \code{cell.occur}.
#'
#' @return a plot showing the mean abundance and variation of all occurred cells. Cells with high mean abundance and/or high variance will be labeled.
#'
#' @export
show.cellsofinterest = function(nums,occurred.cells){
  occurrednums = nums[occurred.cells,]
  aver = apply(occurrednums,1,mean)
  varia = apply(occurrednums,1,var)
  logratio = log(varia/aver)

  aver20 = order(aver,decreasing = T)[1:20]
  varia20 = order(varia,decreasing = T)[1:20]
  ratio20 = order(abs(logratio),decreasing = T)[1:20]
  tolabel = union(union(aver20,varia20),ratio20)

  p = ggplot()+
    geom_point(aes(x=aver,y=varia))+
    geom_point(aes(x=aver[tolabel],y=varia[tolabel]),color="red")+
    geom_text_repel(aes(x=aver[tolabel],y=varia[tolabel],label=names(aver)[tolabel]), size=3,
                    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"))+
    labs(x="mean abundance",y="variance")+
    theme_classic()
  return(p)
}
