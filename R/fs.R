#' Processing deconvolution
#'
#' Processing deconvolution, getting the number of each reference cell(cell type)
#' in each spot of spatial transcriptomics.
#'
#' @param ref Reference, the scRNA-seq data. Pseudo-count will be added and counts will be normalized to TPM.
#' @param st Spatial transcriptomics to be deconvoluted
#' @param cellnames (optional) Cells(Cell types) to be used for deconvolution in reference.
#' If missing, all cells will be used.
#' @param genemode Mode of genes to be included in deconvolution. Three modes are available: "\code{default}", using the intersect of \code{ref} and \code{st};
#' "\code{customized}", indicating the genes yourself; "\code{filtered}", use some parameters to filter the genes (see \code{var_thresh} and \code{exp_thresh}
#' for details).
#' @param gene.list If \code{genemode = "customized"}, this will be the genes you want to use in deconvolution.
#' @param var_thresh If \code{genemode = "filtered"}, this will be the threshold of variance of gene expression in reference.
#' @param exp_thresh If \code{genemode = "filtered"}, this will be the lowest average gene expression per spot in spatial transcriptomics.
#' @param hpmode Mode of choosing hyperparameter adjusting the proportion of
#' LS term and penalty term. Three modes are available: "\code{default}",  using the default
#' hyperparameter we set according to the number of cells and genes; "\code{customized}",
#' setting the hyperparameter yourself; "\code{autoselection}", automatically calculating
#' and selecting the optimal hyperparameter(may take a long time).
#' @param hp If \code{hpmode = "customized"}, this will be the hyperparameter you want to use.
#' @param normalize Whether to normalize the reference or not. `TRUE` is the default value and is recommended when deconvoluting spatial transcriptomics.
#' @param thre Threshold, numbers less than this value in the result will be regarded as zero.
#' Default is \code{1e-10}.
#' @param dopar Whether to use parallel computing. \code{TRUE} is recommended when
#' \code{hpmode = "autoselection"}.
#' @param ncores Number of cores to be used when doing parallel computing.
#' @param realtime Whether to write results into disk in real time.
#'
#' @return A matrix, showing the number of estimated absolute abundance of each cell(cell type) in each spot.
#'
#' @export
deconvoluting = function(ref, st, cellnames=NULL, genemode, gene.list, var_thresh=0.025, exp_thresh=0.03, hpmode, hp, normalize=T, thre=1e-10, dopar=T, ncores, realtime=F, dir=NULL){
  if((dopar==T) & (missing(ncores))){stop("Parameter \"ncores\" is required to avoid latent errors.")}
  if(missing(aver_cell)){stop("Average number of cells not indicated")}
  ref = as.matrix(ref); st = as.matrix(st)
  if(!is.null(cellnames)){
    ncells = length(cellnames)
    ref = ref[,cellnames]
  }
  else ncells = ncol(ref)

  genemode = match.arg(genemode,c("default","customized","filtered"))
  if(genemode=="default"){
    gene.list = intersect(rownames(ref),rownames(st))
  }
  else if(genemode=="customized"){
    intersection = intersect(rownames(ref),rownames(st))
    if(sum(! gene.list %in% intersection)>0){
      warning("Some genes are not in ref or st. Such genes are ignored.")
    }
    gene.list = intersect(intersection,gene.list)
  }
  else if(genemode=="filtered"){
    gene.list = gene.filter(ref,st,var_thresh,exp_thresh)
  }

  ngenes = length(gene.list)
  print(paste0("Number of selected genes: ",ngenes))
  ref = ref[gene.list,];st = st[gene.list,]

  nspots = ncol(st)

  if(normalize==T){
    ref = ref + 0.5
    ref = apply(ref,2,function(x){x/sum(x)*1e06})
  }
  print("Calculating correlation matrix ...")
  r = Hmisc::rcorr(ref)[["r"]]
  r[r<0] = 0

  hpmode = match.arg(hpmode, c("default","customized","autoselection"))
  if(hpmode=="default"){
    hp = round(1e05/ncells^2*ngenes)
    print(paste0("Set hyperparameter as: ",hp))
    nums = rundeconv(ref,st,ncells,nspots,r,hp,dopar,ncores,realtime,dir)
  }
  else if(hpmode=="customized"){
    if(missing(hp)) stop("missing hyperparameter!")
    nums = rundeconv(ref,st,ncells,nspots,r,hp,dopar,ncores,realtime,dir)
  }
  else if(hpmode=="autoselection"){
    hp = 1e05/ncells^2*ngenes
    hplist = c(round(hp/100),round(hp/10),round(hp),round(hp*10),round(hp*100),round(hp*1000))
    numlist = list();residlist = list()
    for(i in 1:6){
      print(paste0("iteration ",i,": hp = ",hplist[i]))
      numlist[[i]] = rundeconv(ref,st,ncells,nspots,r,hplist[i],dopar,ncores,realtime,dir)
      residlist[[i]] = resids(numlist[[i]],ref,st,mode = "spot",inside = "in")
    }
    diff = matrix(nrow = 5, ncol = nspots)
    for(i in 1:5)diff[i,]=(residlist[[i+1]]-residlist[[i]])/hplist[i]
    ord = apply(diff,2,order,decreasing=T)
    totake = as.numeric(names(which.max(table(ord[1,]))))
    print(paste0("Select the hyperparameter as: ",hplist[totake]))
    nums = numlist[[totake]]
  }
  nums[nums<thre] = 0

  nums = as(nums,"dgCMatrix")
  return(nums)
}

rundeconv = function(ref,st,ncells,nspots,r,hp,dopar,ncores,realtime,dir){
  G = Hessian(ref,r,hp)
  sc = norm(G,"2")
  print("Running deconvolution ...")
  if(dopar==F){
    nums = matrix(nrow = ncells,ncol = nspots)
    pb <- txtProgressBar(style=3)
    if(realtime==T){
      if(length(dir)==0)dir = "./real_time_results"
      if(!dir.exists(dir))dir.create(dir)
      barcodes = colnames(st)
      for(i in 1:nspots){
        nums[,i] = solveqp(ncells,r = r,x = ref,y = st[,i],G = G,sc = sc)
        write.csv(nums[,i],file = paste0(dir,"/",barcodes[i],".csv"),row.names = F)
        setTxtProgressBar(pb, i/nspots)
      }
    }
    else{
      for(i in 1:nspots){
        nums[,i] = solveqp(ncells,r = r,x = ref,y = st[,i],G = G,sc = sc)
        setTxtProgressBar(pb, i/nspots)
      }
    }
    close(pb)
  }
  else{
    cl = snow::makeCluster(ncores)
    doSNOW::registerDoSNOW(cl)
    pb = txtProgressBar(max = nspots, style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    opts = list(progress = progress)
    if(realtime==T){
      if(length(dir)==0)dir = "./real_time_results"
      if(!dir.exists(dir))dir.create(dir)
      nums = foreach::foreach(y=iterators::iter(st,by="col"),.combine=cbind,.inorder=T,.packages = "quadprog",.export = "solveqp",.options.snow = opts) %dopar% {
        num = solveqp(ncells,r = r,x = ref,y = y,G = G,sc = sc)
        write.csv(num,file = paste0(dir,"/",colnames(y),".csv"),row.names = F)
        return(num)
      }
    }
    else{
      nums = foreach::foreach(y=iterators::iter(st,by="col"),.combine=cbind,.inorder=T,.packages = "quadprog",.export = "solveqp",.options.snow = opts) %dopar% {
        num = solveqp(ncells,r = r,x = ref,y = y,G = G,sc = sc)
        return(num)
      }
    }
    close(pb)
    snow::stopCluster(cl)
  }
  rownames(nums) = colnames(ref)
  colnames(nums) = colnames(st)
  return(nums)
}

#' @export
Hessian = function(ref,r,hp){
  # Calculating the Hessian matrix for this quadratic programming problem
  print("Calculating Hessian matrix ...")
  # ngenes = nrow(ref)
  # G = .Fortran("Hessian_f",ref,ngenes,ncells,r,hp)

  # G=matrix(data=0,nrow=ncells,ncol=ncells)
  # for(i in 1:ncells) G[i,i]=sum(ref[,i]^2)+hp*sum(r[i,])
  # for(i in 1:(ncells-1)){
  #   for(j in (i+1):ncells){
  #     G[i,j]=2*(crossprod(ref[,i],ref[,j])-hp*r[i,j])
  #   }
  # }
  # G=G+t(G)

  diag(r) = 0
  p1 = 2* t(ref) %*% ref
  p2 = 2*hp* (diag(apply(r,1,sum)) - r)
  G = p1+p2
  return(as.matrix(G))
}

solveqp = function(ncells,r,x,y,G,sc){
  # Solving the QP problem
  d=2*t(x) %*% y #t(d)=2*t(y) %*% x
  a=diag(rep(1,ncells))
  b=rep(0,ncells)
  result=quadprog::solve.QP(G/sc,d/sc,a,b)
  return(result[["solution"]])
}

#' Converting single cell results to cell type results
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
  names(tab)=NULL
  begin=1
  for(i in 1:ntypes){
    end = begin+tab[i]-1
    if(tab[i]==1){
      typenums[i,]=as.matrix(nums[begin,])
    }
    else{
      typenums[i,] = apply(nums[begin:end,],2,sum)
    }
    begin = end+1
  }
  return(typenums)
}

#' Giving results interpretability
#'
#' One way of gaining interpretability, estimating absolute abundance
#'
#' @param res Result of `deconvoluting`.
#' @param aver.cell A priori value indicating the average cell number of one spot. This varies according to the platform.
#'
#' @return A cell-by-spot matrix, each place is the estimated absolute abundance of the cell in the spot.
#'
#' @export
to.absolute.abundance = function(res,aver.cell){
  nspots = ncol(res)
  totalcell = nspots*aver.cell
  sumnum = sum(res)
  res = res/sumnum*totalcell
  return(res)
}

#' Giving results interpretability
#'
#' One way of gaining interpretability, normalize results to proportion
#'
#' @param res Result of `deconvoluting`.
#'
#' @return A cell-by-spot matrix, with the sum of each column equaling to 1.
#'
#' @export
to.proportion = function(res){
  res = apply(res,2,function(x){x/sum(x)})
  return(res)
}

#' Getting the gene expression profile of each cell type.
#'
#' If you want to use the gene expression profile of each cell type as reference,
#' this function can calculate the average expression of genes within all cells
#' of every cell type.
#'
#' @param sc Raw scRNA-seq data.
#' @param annotations A \code{data.frame} with 2 columns: the first column should be barcodes,
#' the second column should be the corresponding cell type.
#' @param gene.list (optional) Genes to be calculated. If missing, all genes
#' will be calculated.
#'
#' @return A matrix showing the gene expression profile of each cell type,
#' one row represents one cell type and one column represents one gene.
#' If there are cells missing annotations, this function will also write a global variable named `new.annotation`,
#' where cells missing annotation are labeled as "_Unknown".
#'
#' @export
get.ref = function(sc,annotations,gene.list = NULL){
  sc = as(sc,"dgCMatrix")
  if(!is.null(gene.list))sc = sc[gene.list,]
  ngenes=nrow(sc)
  barcodes=colnames(sc)
  ncells=length(barcodes)
  shared=intersect(barcodes,annotations[,1])
  if(length(shared)<ncells){
    warning("There are cells missing annotation. Labeled as '_Unknown'.")
  }
  new.annotation=as.matrix(matrix(nrow=ncells,ncol=1))
  for(i in 1:ncells){
    idx=which(annotations[,1]==barcodes[i])
    if(length(idx)==0){
      new.annotation[i,1]="_Unknown"
    }
    else{new.annotation[i,1]=annotations[i,2]}
  }
  ords = order(new.annotation)
  new.annotation = new.annotation[ords,]
  tab = table(new.annotation)
  sc = sc[,ords]
  ntypes=length(tab)
  ref=matrix(nrow=ngenes,ncol=ntypes)
  rownames(ref)=rownames(sc)
  colnames(ref)=dimnames(tab)[[1]]
  begin=1
  for(i in 1:ntypes){
    if(tab[i]==1)ref[,i]=sc[,begin]
    else{
      end = begin+tab[i]-1
      names(end) = NULL
      ref[,i] = apply(sc[,begin:end],1,mean)
    }
    begin = end+1
  }
  if(length(shared)<ncells){
    new.annotation <<- cbind.data.frame(colnames(sc),new.annotation)
  }
  return(ref)
}

#' Getting differentially expressed genes
#'
#' Getting differentially expressed genes across cells(cell types).
#' The expression of DE genes shouldn't be too low in both scRNA-seq data
#' and spatial transcriptomics data, and should have some difference across cells(cell types).
#'
#' @param ref Reference.
#' @param st Spatial transcriptomics data.
#' @param var_thresh The threshold of variance of gene expression in reference.
#' @param exp_thresh Lowest average gene expression per spot in spatial transcriptomics.
#'
#' @export
gene.filter = function(ref,st,var_thresh = 0.025,exp_thresh = 0.003){
  ref = as.matrix(ref); st = as.matrix(st)
  nspots = ncol(st)
  MIN.OBS = exp_thresh*nspots
  bulk.vec = rowSums(st)
  gene.list = rownames(ref)
  if(length(grep("[Mm][Tt]-",gene.list))>0)gene.list = gene.list[-grep("mt-", gene.list)]
  gene.list = intersect(gene.list, names(bulk.vec))
  gene.list = gene.list[bulk.vec[gene.list] >= MIN.OBS]
  vars=apply(ref[gene.list,],1,var)
  gene.list = names(vars[which(vars>=var_thresh)])
  return(gene.list)
}

#' Calculating residuals
#'
#' Calculating residuals for the whole spatial transciptome or some genes.
#'
#' @param nums Result generated by \code{deconvoluting}.
#' @param ref Reference.
#' @param st Spatial transcriptomics data.
#' @param mode Mode of calculating residuals. "\code{spot}" for spot-wise and "\code{gene}" for gene-wise.
#' @param genes (optional) If \code{mode = "gene"}, this will be a gene or a list of genes whose
#' residuals you want to calculate. The default is all genes.
#' @param inside Indicates whether normalization is required. You don't need to adjust it when using.
#'
#' @export
resids = function(nums,ref,st,mode,genes,inside="out"){
  ref = as.matrix(ref); nums = as.matrix(nums); st = as.matrix(st)
  mode = match.arg(mode, c("spot", "gene"))
  if(mode=="gene"){
    if(!missing(genes))ref = ref[genes,];st = st[genes,]
  }
  inside = match.arg(inside, c("in","out"))
  if(inside == "out"){
    ref = ref+0.5
    ref = apply(ref,2,function(x){x/sum(x)*1e06})
  }
  ress = ref %*% nums - st
  if(mode=="gene")res = apply(ress,1,function(x){norm(x,"2")^2})
  else if(mode=="spot")res = apply(ress,2,function(x){norm(x,"2")^2})
  return(res)
}

#' Calculating Entropy and Perplexity
#'
#' Calculating Entropy and Perplexity for each spot.
#'
#' @param nums Results.
#'
#' @return A matrix with 2 rows: the first is Entropy, the second is Perplexity.
#'
#' @export
hNpp = function(nums, thre = 1e-10){
  nums = as.matrix(nums)
  ncells = nrow(nums); nspots = ncol(ncells)
  nums = apply(nums,2,function(x){x/sum(x)})
  s = ifelse(nums==0,0,-nums*log(nums,2))
  entropy = apply(s,2,sum)
  perplexity = 2^entropy
  return(rbind(entropy,perplexity))
}

#' Sample single cells.
#'
#' Sample some single cells from original scRNA-seq data. If cell annotations are provided,
#' this function will run stratified sampling; otherwise it will run random sampling.
#'
#' @param ncells Number of cells in raw scRNA-seq data.
#' @param annotations A \code{data.frame} with 2 columns: the first column should be barcodes,
#' the second column should be the corresponding cell type. If missing, complete random sampling will be run.
#' @param size Total number of cells you want to take out.
#' @param prot Protection of cell types with low occurrence in raw data. If \code{prot = T},
#' it will be ensured that there is at least one cell of each cell type in the sampling results;
#' otherwise due to \code{size}, some cell types may not appear. Default is \code{T}.
#'
#' @return If annotations are provided, this function will return a \code{data.frame} containing
#' barcodes of selected cells and their annotations, and its \code{rownames} are the original indices.
#' Otherwise only indices.
#'
#' @export
cell.sampling=function(ncells,annotations = NULL,size,prot = T){
  if(ncells<size){
    warning("ncells < size. All single cells will be returned.")
    return(annotations)
  }
  if(length(annotations)==0){
    print("Missing annotations. Random sampling will be carried out.")
    refcells = sample(1:ncells,size)
    return(refcells)
  }
  else{
    if(ncells!=dim(annotations)[1]) stop("The dimension of 'annotations' do not match 'ncells'.")
    annotations = annotations[order(annotations[,2]),]
    counts = as.data.frame(table(annotations[,2]))
    cell_type_names = counts[,1]
    counts = as.matrix(counts[,2])
    totake = round(counts*size/sum(counts))
    if(prot==T)totake[totake==0]=totake[totake==0]+1
    print(paste0("cells actually taken: ",sum(totake)))
    labs = cbind.data.frame(annotations[,2],1:ncells)
    colnames(labs) = c("annotations","index")
    refcells = sampling::strata(labs,stratanames ="annotations",size = totake,method = "srswor")
    return(annotations[refcells[,2],])
  }
}

#' Correction of single-cell expression profile.
#'
#' Use results and spatial transcriptomics to reversely estimate the spatial-specific expression profile.
#'
#' @param num One column of the results.
#' @param st1 One column of Spatial transcriptomics data. The spot in which the expression profile is to correct.
#' @param ref Reference.
#' @param gene.list Genes to be estimated. Default is all.
#' @param dopar Whether to use parallel computing.
#'
#' @export
profile_correction = function(num, st1, ref, gene.list=NULL, dopar=T, ncores){
  num = as.matrix(num); st1 = as.matrix(st1); ref = as.matrix(ref)
  ref = ref + 0.5
  ref = apply(ref,2,function(x){x/sum(x)*1e06})
  cells = which(num!=0)
  if(length(gene.list)!=0){
    st1 = as.matrix(st1[gene.list,]);ref = as.matrix(ref[gene.list,])
  }
  ref = as.matrix(ref[,cells]); num = as.matrix(num[cells])
  ngenes = dim(st1)[1]; ncells = dim(num)[1]
  if(dopar==T){
    cl = snow::makeCluster(ncores)
    doSNOW::registerDoSNOW(cl)
    pb = txtProgressBar(max = ngenes, style = 3)
    progress = function(n) setTxtProgressBar(pb, n)
    opts = list(progress = progress)
    ests = foreach::foreach(i=1:ngenes,.combine=rbind,.inorder=T,.packages = "quadprog",.export = "solveqp2",.options.snow = opts) %dopar% {
      est = solveqp2(st1[i,],num,ref[i,],ncells)
      return(est)
    }
    close(pb)
    snow::stopCluster(cl)
    rownames(ests) = rownames(ref)
    colnames(ests) = colnames(ref)
  }
  else{
    ests = matrix(nrow=ngenes,ncol=ncells)
    rownames(ests) = rownames(ref)
    colnames(ests) = colnames(ref)
    for(gene in gene.list){
      ests[gene,] = solveqp2(st1[gene,],num,ref[gene,],ncells)
    }
  }
  return(ests)
}

#' Solving QP in correcting expression profile
solveqp2 = function(y,num,ref,ncells){
  G = diag(2,nrow=ncells,ncol=ncells)
  d = 2*t(ref) # d^T = 2*ref
  a = num # a^T = num^T
  b = y
  result=quadprog::solve.QP(G,d,a,b)
  return(t(result[["solution"]]))
}

#' Seurat interface 1
#'
#' Seurat interface for single-cell transcriptomics
#'
#' @param sc_seu Seurat object of single-cell transcriptomics
#' @param anno_slot Character, name of slot storing annotations information, under "meta.data" slot.
#'
#' @return A list, containing single-cell expression matrix and annotations.
#'
#' @export
extract.seurat.sc = function(sc_seu, anno_slot){
  sc = list()
  sc[["expr"]] = sc_seu@assays$RNA@counts
  sc[["annotations"]] = eval(parse(text=paste0(substitute(sc_seu),"$",anno_slot)))
  return(sc)
}

#' Seurat interface 2
#'
#' Seurat interface for spatial transcriptomics
#'
#' @param sc_seu Seurat object of spatial transcriptomics
#' @param coord_slot Character, name of slot storing annotations information, under "image" slot. Default is "image".
#'
#' @return A list, containing single-cell expression matrix and annotations.
#'
#' @export
extract.seurat.st = function(st_seu, coord_slot="image"){
  st = list()
  st[["expr"]] = st_seu@assays$Spatial@counts
  coords = eval(parse(text=paste0(substitute(st_seu),"@images$",coord_slot,"@coordinates")))
  st[["coords"]] = as.matrix(coords[,2:3])
  return(st)
}
