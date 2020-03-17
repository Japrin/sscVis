#' dispaly message with time stamp
#' @param msg characters; message to display
#' @export
loginfo <- function(msg) {
  timestamp <- sprintf("%s", Sys.time())
  msg <- paste0("[",timestamp, "] ", msg,"\n")
  cat(msg)
}

#' Correlation calculation. use BLAS and data.table to speed up.
#'
#' @importFrom data.table frank
#' @importFrom RhpcBLASctl omp_get_num_procs omp_set_num_threads
#' @param x matrix; input data, rows for variable (genes), columns for observations (cells).
#' @param y matrix; input data, rows for variable (genes), columns for observations (cells) (default: NULL)
#' @param method character; method used. (default: "pearson")
#' @param nthreads integer; number of threads to use. if NULL, automatically detect the number. (default: NULL)
#' @param na.rm logical; remove missing values. (default: T)
#' @details calcualte the correlation among variables(rows)
#' @return correlation coefficient matrix among rows
#' @export
cor.BLAS <- function(x,y=NULL,method="pearson",nthreads=NULL,na.rm=T)
{
  if(is.null(nthreads))
  {
    nprocs <- RhpcBLASctl::omp_get_num_procs()
    RhpcBLASctl::omp_set_num_threads(max(nprocs-1,1))
  }else{
    RhpcBLASctl::omp_set_num_threads(nthreads)
  }
  cor.pearson <- function(x,y=NULL)
  {
    if(is.null(y)){
	  ### x = x - rowMeans(t(na.omit(t(x))))
      x = x - rowMeans(x,na.rm=na.rm)
      x = x / sqrt(rowSums(x^2,na.rm=na.rm))
      ### cause 'memory not mapped' :( ; and slower in my evaluation: 38 sec .vs. 12 sec.
      #x.cor = tcrossprod(x)
	  if(na.rm){ x[is.na(x)] <- 0 }
      x.cor = x %*% t(x)
      return(x.cor)
    }else{
      x = x - rowMeans(x,na.rm=na.rm)
      x = x / sqrt(rowSums(x^2,na.rm=na.rm))
      y = y - rowMeans(y,na.rm=na.rm)
      y = y / sqrt(rowSums(y^2,na.rm=na.rm))
      #xy.cor <- tcrossprod(x,y)
	  if(na.rm){ x[is.na(x)] <- 0 }
	  if(na.rm){ y[is.na(y)] <- 0 }
      xy.cor <- x %*% t(y)
      return(xy.cor)
    }
  }
  x <- as.matrix(x)
  if(!is.matrix(x)){
    warning("x is not like a matrix")
    return(NULL)
  }
  if(!is.null(y)){
    y <- as.matrix(y)
    if(!is.matrix(y)){
      warning("y is not like a matrix")
      return(NULL)
    }
  }
  if(method=="pearson"){
    return(cor.pearson(x,y))
  }else if(method=="spearman"){
    if(is.null(y)){
      return(cor.pearson(t(apply(x, 1, data.table::frank, na.last="keep"))))
    }else{
      return(cor.pearson(t(apply(x, 1, data.table::frank, na.last="keep")),
                         t(apply(y, 1, data.table::frank, na.last="keep"))))
    }
  }else{
    warning("method must be pearson or spearman")
    return(NULL)
  }
}


#' Modified from limma::removeBatchEffect, a little different design matrix
#' @param x matrix; samples in columns and variables in rows
#' @param batch character; batch vector (default: NULL)
#' @param covariates double; other covariates to adjust (default: NULL)
#' @param ... parameters passed to lmFit
#' @details Modified from limma::removeBatchEffect, a little different design matrix
#' @return a matrix with dimention as input ( samples in rows and variables in columns)
#' @importFrom limma lmFit
#' @importFrom stats model.matrix
#' @export
simple.removeBatchEffect <- function (x, batch = NULL, covariates = NULL, ...)
{
    if (is.null(batch) && is.null(covariates))
        return(as.matrix(x))
	if(!is.null(batch) && length(unique(batch))==1){
		return(t(scale(t(as.matrix(x)),scale=F)))
	}
    if (!is.null(batch)) {
        batch <- as.factor(batch)
        batch <- model.matrix(~batch)
    }
    if (!is.null(covariates))
        covariates <- as.matrix(covariates)
    X.batch <- cbind(batch, covariates)
    fit <- limma::lmFit(x, X.batch, ...)
    beta <- fit$coefficients
    beta[is.na(beta)] <- 0
    ret.V <- as.matrix(x) - beta %*% t(X.batch)
    return(ret.V)
}

#' run dynamicTreeCut::cutreeDynamic on rows of the given matrix
#' @param dat data frame or matrix;
#' @param method.hclust character; clustering method for hclust [default: "ward.D2"]
#' @param method.distance character; distance method for hclust [default: "spearman"]
#' @param ... parameters passed to dynamicTreeCut::cutreeDynamic
#' @details dynamicTreeCut::cutreeDynamic on rows of the given matrix
#' @return a matrix with dimention as input ( samples in rows and variables in columns)
#' @importFrom stats dist hclust as.dist as.dendrogram
#' @importFrom dynamicTreeCut cutreeDynamic
#' @export
run.cutreeDynamic <- function(dat,method.hclust="ward.D2",method.distance="spearman",
							  #deepSplit=4, minClusterSize=2,
							  ...)
{
	#requireNamespace("dendextend")
    obj.hclust <- NULL
    if(method.distance=="spearman" || method.distance=="pearson"){
		tryCatch({
			###obj.distM <- as.dist(1-sscClust:::cor.BLAS((dat),method=method.distance,nthreads=1))
			obj.distM <- as.dist(1-cor.BLAS((dat),method=method.distance,nthreads=1))
			obj.hclust <- stats::hclust(obj.distM, method=method.hclust)
		},error = function(e){
			cat("using spearman/pearson as distance failed;try to fall back to use euler distance ... \n");
		})
    }
    if(is.null(obj.hclust)){
		obj.distM <- stats::dist(dat)
		obj.hclust <- stats::hclust(obj.distM,method.hclust)
    }
    ##### if method.hclust=="complete", some clusters from cutreeDynamic are not consistent with the dendrogram
    cluster.label <- dynamicTreeCut::cutreeDynamic(obj.hclust,distM=as.matrix(obj.distM),
                                                #method = "hybrid",
                                                #deepSplit=deepSplit,minClusterSize=minClusterSize,
                                                ...)
    obj.dend <- as.dendrogram(obj.hclust)
    ncls <- length(unique(cluster.label))
    colSet.cls <- auto.colSet(ncls,"Paired")
    ###colSet.cls <- sscClust:::auto.colSet(ncls)
    names(colSet.cls) <- unique(cluster.label[order.dendrogram(obj.dend)])
    col.cls <- data.frame("k0"=sapply(cluster.label,function(x){ colSet.cls[as.character(x)] }))

#    print(colSet.cls)
#    print(table(cluster.label))
#	pdf("test.pdf")
#	opar <- par(mar=c(15,4,4,2))
#	plot(obj.dend)
#	dev.off()
#	par(opar)

    obj.branch <- dendextend::color_branches(obj.dend,
                                 clusters=cluster.label[order.dendrogram(obj.dend)],
                                 col=colSet.cls)
	obj.branch <- dendextend::set(obj.branch,"branches_lwd", 1.5)

#	pdf("test.01.pdf")
#	opar <- par(mar=c(15,4,4,2))
#	plot(obj.branch)
#	dev.off()
#	par(opar)
    ###aa <- plot.branch(obj.hclust,"test.01",cluster=cluster.label)

    return(list("hclust"=obj.hclust,"dist"=obj.distM,"cluster"=cluster.label,"branch"=obj.branch))
}

#' run cutree on rows of the given matrix
#' @param dat data frame or matrix;
#' @param method.hclust character; clustering method for hclust [default: "ward.D2"]
#' @param method.distance character; distance method for hclust [default: "spearman"]
#' @param k integer; number of clusters [default: 1]
#' @param ... parameters passed to cutree
#' @details cutree on rows of the given matrix
#' @return a matrix with dimention as input ( samples in rows and variables in columns)
#' @importFrom stats dist hclust as.dendrogram as.dist
#' @importFrom dendextend color_branches cutree
#' @export
run.cutree <- function(dat,method.hclust="ward.D2",method.distance="spearman",k=1,
							  ...)
{
	#requireNamespace("dendextend")
	ret <- list()
    {
		branch <- FALSE
		obj.hclust <- NULL
		if(method.distance=="spearman" || method.distance=="pearson"){
			tryCatch({
				obj.distM <- as.dist(1-cor.BLAS((dat),method=method.distance,nthreads=1))
				#obj.distM <- as.dist(1-cor.BLAS((dat),method=method.distance,nthreads=1))
				obj.hclust <- stats::hclust(obj.distM, method=method.hclust)
			},error = function(e){
				cat("using spearman/pearson as distance failed;try to fall back to use euler distance ... \n");
			})
		}
		if(is.null(obj.hclust)){
			obj.distM <- stats::dist(dat)
			obj.hclust <- stats::hclust(obj.distM,method.hclust)
		}
		obj.dend <- as.dendrogram(obj.hclust)
		cluster <- cutree(obj.hclust,k=k,...)
        colSet.cls <- auto.colSet(length(unique(cluster)),"Paired")
		branch <- dendextend::color_branches(obj.dend,clusters=cluster[order.dendrogram(obj.dend)],col=colSet.cls)
		branch <- dendextend::set(branch,"branches_lwd", 1.5)

		ret[["hclust"]] <- obj.hclust
		ret[["dist"]] <- obj.distM
		ret[["cluster"]] <- cluster
		ret[["branch"]] <- branch
    }
    return(ret)
}

####### operations on sce object ######

#' Build an SingleCellExperiment object
#'
#' Build an SingleCellExperiment object from a matrix or data frame
#'
#' @importFrom SingleCellExperiment SingleCellExperiment rowData
#' @importFrom S4Vectors metadata `metadata<-`
#' @importFrom stats setNames
#' @param x matrix/data.frame or SingleCellExperiment; input expression data
#' @param display.name a vector, should be human readable gene name
#' @param assay.name assay name (default "exprs")
#' @details if x is an object of SingleCellExperiment, just clear the metadata;
#' if x is matrix/data.frame, convert it to an object of SingleCellExperiment.
#' Also a vector `display.name` can be provided, which would be used in some plots,
#' such as geneOnTSNE, heatmap. The row names of SingleCellExperiment object usually be gene id
#' (e.g. entrez ID, Ensemble ID), the `display.name` should be human readable gene name (
#' e.g. HGNC gene symbol). If `display.name` is NULL (default), the row names of SingleCellExperiment object
#' would be used.
#' @return an object of \code{SingleCellExperiment} class
#' @export
ssc.build <- function(x,display.name=NULL,assay.name="exprs")
{
  obj <- NULL
  if(class(x)=="SingleCellExperiment")
  {
    obj <- x
    metadata(obj)$ssc <- list()
  }else if(class(x) %in% c("matrix","data.frame")){
    obj <- SingleCellExperiment(assays = setNames(list(as.matrix(x)),assay.name))
  }else if(class(x) %in% c("dgCMatrix","dgTMatrix")){
    obj <- SingleCellExperiment(assays = setNames(list(x),assay.name))
  }
  metadata(obj)$ssc$colSet <- list()
  if(is.null(rowData(obj)[["display.name"]])){
    if(!is.null(display.name)){
      f.na <- is.na(display.name)
      display.name[f.na] <- row.names(obj)[f.na]
      rowData(obj)[,"display.name"] <- display.name
    }else{
      rowData(obj)[,"display.name"] <- row.names(obj)
    }
    if(is.null(names(rowData(obj)$display.name))){ names(rowData(obj)$display.name) <- row.names(obj) }
  }
  if(is.null(obj)){
    warning("SingleCellExperiment object building failed!")
  }
  return(obj)
}


#' sort by hierarchical clustering
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param order.col logical; wheter order columns (default: FALSE)
#' @param order.row logical; wheter order row (default: FALSE)
#' @param clustering.distance character; one of spearmn, pearson and euclidean (default: "spearman")
#' @param clustering.method character; method for hclust (default: "complete")
#' @param k.row integer; number of clusters in the rows (default: 1)
#' @param k.col integer; number of clusters in the columns (default: 1)
#' @importFrom stats hclust
#' @importFrom dendextend color_branches
#' @details order genes or cells according the assay, using hclust.
#' @export
ssc.assay.hclust <- function(obj,assay.name="exprs",
                             order.col=FALSE,order.row=FALSE,
                             clustering.distance="spearman",clustering.method="complete",
                             k.row=1,k.col=1)
{
    if(order.col && ncol(obj)>2)
    {
		ret.col <- run.cutree(t(assay(obj,assay.name)),method.hclust=clustering.method,
							   method.distance=clustering.distance,k=k.col)

        obj <- obj[,ret.col$hclust$order]
        metadata(obj)$assay.hclust$col <- ret.col$hclust
        metadata(obj)$assay.hclust$branch.col <- ret.col$branch
#        branch.col <- FALSE
#        obj.hclust.col <- NULL
#        if(clustering.distance=="spearman" || clustering.distance=="pearson"){
#            tryCatch({
#                dist.out <- cor.BLAS(t(assay(obj,assay.name)),method=clustering.distance,nthreads=1)
#                obj.hclust.col <- hclust(as.dist(1-dist.out), method=clustering.method)
#            },error = function(e){
#                cat("using spearman/pearson as distance failed;try to fall back to use euler distance ... \n");
#            })
#        }
#        if(is.null(obj.hclust.col)){
#            obj.hclust.col <- hclust(dist(t(assay(obj,assay.name))),method=clustering.method)
#        }
#        branch.col <- dendextend::color_branches(as.dendrogram(obj.hclust.col),k=k.col)
#        obj <- obj[,obj.hclust.col$order]
#        metadata(obj)$assay.hclust$col <- obj.hclust.col
#        metadata(obj)$assay.hclust$branch.col <- branch.col
    }
    if(order.row && nrow(obj)>2)
	{
		ret.row <- run.cutree((assay(obj,assay.name)),method.hclust=clustering.method,
							   method.distance=clustering.distance,k=k.row)

        obj <- obj[ret.row$hclust$order,]
        metadata(obj)$assay.hclust$row <- ret.row$hclust
        metadata(obj)$assay.hclust$branch.row <- ret.row$branch
#        branch.row <- FALSE
#        obj.hclust.row <- NULL
#        if(clustering.distance=="spearman" || clustering.distance=="pearson"){
#            tryCatch({
#                dist.out <- cor.BLAS(assay(obj,assay.name),method=clustering.distance,nthreads=1)
#                obj.hclust.row <- hclust(as.dist(1-dist.out), method=clustering.method)
#            },error = function(e){
#                cat("using spearman/pearson as distance failed;try to fall back to use euler distance ... \n");
#            })
#        }
#        if(is.null(obj.hclust.row)){
#            obj.hclust.row <- hclust(dist(assay(obj,assay.name)),method=clustering.method)
#        }
#        branch.row <- dendextend::color_branches(as.dendrogram(obj.hclust.row),k=k.row)
#        obj <- obj[obj.hclust.row$order,]
#        metadata(obj)$assay.hclust$row <- obj.hclust.row
#        metadata(obj)$assay.hclust$branch.row <- branch.row
    }
    return(obj)
}

#' order genes and cells
#' @param obj object of \code{singleCellExperiment} class
#' @param columns.order character; columns of colData(obj) used for ordering (default: NULL)
#' @param gene.desc data.frame; it must contain columns geneID and Group (default: NULL)
#' @importFrom data.table as.data.table
#' @export
ssc.order <- function(obj,columns.order=NULL,gene.desc=NULL)
{
    if(!is.null(gene.desc) && ("Group" %in% colnames(gene.desc)) && ("geneID" %in% colnames(gene.desc))){
		gene.desc <- as.data.table(gene.desc)
		gene.desc <- gene.desc[order(Group),]
        obj <- obj[gene.desc$geneID,]
    }
    if(!is.null(columns.order)){
        annDF <- as.data.frame(colData(obj)[columns.order])
        annDF <- annDF[eval(parse(text=sprintf("with(annDF,order(%s))", columns.order))),,drop=F]
        obj <- obj[,rownames(annDF)]
    }
    return(obj)
}

#' calcualte the average expression of cells of each group
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param gene character; only consider the specified gnees (default: NULL)
#' @param column character; columns in colData(obj) to be averaged. (default: "majorCluster")
#' @param ncell.downsample integer; for each group, number of cells downsample to. (default: NULL)
#' @param avg character; average method. can be one of "mean", "diff", "zscore" . (default: "mean")
#' @param ret.type character; return type. can be one of "data.melt", "data.cast", "data.mtx". (default: "data.melt")
#' @importFrom plyr ldply
#' @importFrom Matrix rowMeans
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats rowSds
#' @importFrom data.table dcast setkey
#' @importFrom S4Vectors DataFrame
#' @details multiple average methods are implemented
#' @export
ssc.average.cell <- function(obj,assay.name="exprs",gene=NULL,column="majorCluster",ncell.downsample=NULL,
                             avg="mean",ret.type="data.melt")
{
  if(!all(column %in% colnames(colData(obj)))){
    warning(sprintf("some column(s) not in the obj: %s \n",column))
    return(NULL)
  }
  if(!is.null(gene)){
    obj <- obj[gene,]
  }
  
  #### downsample cells
  if(!is.null(ncell.downsample)){
    clust <- colData(obj)[,column]
    names(clust) <- colnames(obj)
    grp.list <- unique(clust)
    f.cell <- unlist(sapply(grp.list,function(x){
      x <- names(clust[clust==x])
      sample(x,min(length(x),ncell.downsample)) }))
    obj <- obj[,f.cell]
  }
  
  cls <- unique(colData(obj)[,column,drop=F])
  data.melt.df <- as.data.table(ldply(seq_len(nrow(cls)),function(i){
    r.filter <- data.table(as.data.frame(cls[i,,drop=F]))
    f.in <- colSums(laply(names(r.filter),function(x){ obj[[x]]==r.filter[[x]][1] },.drop = F))==length(r.filter)
    f.out <- !f.in
    obj.in <- obj[,f.in]
    avg.in <- NULL
    avg.in <- Matrix::rowMeans(assay(obj.in,assay.name),na.rm = T)
    if(avg=="mean"){
      dat.ret <- data.table(geneID=names(avg.in))
      for(x in names(r.filter)){
        dat.ret[[x]] <- r.filter[[x]][1]
      }
      dat.ret[["avg"]] <- avg.in
      return(dat.ret)
    }else if (avg=="diff"){
      obj.out <- obj[,f.out]
      avg.out <- Matrix::rowMeans(assay(obj.out,assay.name),na.rm = T)
      dat.ret <- data.table(geneID=names(avg.in))
      for(x in names(r.filter)){
        dat.ret[[x]] <- r.filter[[x]][1]
      }
      dat.ret[["avg"]] <- avg.in - avg.out
      return(dat.ret)
    }else if (avg=="zscore"){
      obj.out <- obj[,f.out]
      avg.out <- Matrix::rowMeans(assay(obj.out,assay.name),na.rm = T)
      ##sd.r <- matrixStats::rowSds(assay(obj,assay.name))
      sd.r <- DelayedMatrixStats::rowSds(DelayedArray(assay(obj,assay.name)),na.rm = T)
      dat.ret <- data.table(geneID=names(avg.in))
      for(x in names(r.filter)){
        dat.ret[[x]] <- r.filter[[x]][1]
      }
      dat.ret[["avg"]] <- (avg.in-avg.out)/sd.r
      dat.ret$avg[is.na(dat.ret$avg)] <- 0
      return(dat.ret)
    }
  }))
  
  if(ret.type=="data.melt"){
    return(data.melt.df)
  }else if(ret.type=="data.dcast"){
    dat.df <- dcast(data.melt.df,sprintf("geneID~%s",paste0(column,collapse = "+")),value.var="avg")
    ##rownames(dat.df) <- dat.df[,1]
    setkey(dat.df,"geneID")
    dat.df <- dat.df[rownames(obj),]
    return(dat.df)
  }else if(ret.type=="data.mtx"){
    dat.df <- dcast(data.melt.df,sprintf("geneID~%s",paste0(column,collapse = "+")),value.var="avg")
    dat.mtx <- as.matrix(dat.df[,-1])
    rownames(dat.mtx) <- dat.df[[1]]
    dat.mtx <- dat.mtx[rownames(obj),]
    return(dat.mtx)
  }else if(ret.type=="sce"){
    dat.df <- dcast(data.melt.df,sprintf("geneID~%s",paste0(column,collapse = "+")),value.var="avg")
    dat.mtx <- as.matrix(dat.df[,-1])
    rownames(dat.mtx) <- dat.df[[1]]
    dat.mtx <- dat.mtx[rownames(obj),]
    obj.ret <- ssc.build(dat.mtx,assay.name=assay.name,display.name=rowData(obj)$display.name)
    obj.ret.colDat <- unique(data.melt.df[,column,with=F])
    obj.ret.colDat$cid <- do.call(paste, c(obj.ret.colDat, list(sep = '_')))
    setkey(obj.ret.colDat,"cid")
    obj.ret.colDat.df <- as.data.frame(obj.ret.colDat[colnames(obj.ret),])
    rownames(obj.ret.colDat.df) <- obj.ret.colDat.df$cid
    colData(obj.ret) <- DataFrame(obj.ret.colDat.df)
    return(obj.ret)
  }
}


#' calcualte the average expression of genes of each cell
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param gene list; named list; NULL for all genes (default: NULL)
#' @param ncell.downsample integer; for each group, number of cells downsample to. (default: NULL)
#' @param avg character; average method. can be one of "mean", "diff", "zscore" . (default: "mean")
#' @param ret.type character; return type. can be one of "data.melt", "data.cast", "data.mtx". (default: "data.melt")
#' @importFrom plyr laply llply
#' @importFrom Matrix colMeans
#' @importFrom SummarizedExperiment rowData
#' @importFrom SingleCellExperiment `reducedDims<-` reducedDims
#' @details multiple average methods are implemented
#' @export
ssc.average.gene <- function(obj,assay.name="exprs",gene=NULL,ncell.downsample=NULL,
                             avg="mean",ret.type="sce")
{
  if(is.null(gene)){
    gene <- list("Score"=seq_len(nrow(obj)))
  }else{
    gene <- llply(gene,function(x){ match(x,rowData(obj)$display.name)  })
  }
  if(avg=="mean"){
    mean.mtx <- laply(gene,function(x){ colMeans(assay(obj,assay.name)[x,,drop=F],na.rm = T) },.drop = F)
    rownames(mean.mtx) <- names(gene)
    obj.mean <- ssc.build(mean.mtx)
    colData(obj.mean) <- colData(obj)
    reducedDims(obj.mean) <- reducedDims(obj)
    return(obj.mean)
  }else{
    return(NULL)
  }
}

#' Convert display name (gene symbol) to gene id
#' @importFrom SingleCellExperiment rowData
#' @param obj object of \code{SingleCellExperiment} class
#' @param display.name character; disaply name
#' @return If successfull vector contains the gene id; otherwise NULL
#' @export
ssc.displayName2id <- function(obj,display.name)
{
  ret <- NULL
  if("display.name" %in% colnames(rowData(obj)))
  {
    lookup.table <- structure(rowData(obj)[,"display.name"],
                              names=rownames(obj))
    #ret <- lookup.table[ids]
    ret <- lookup.table[lookup.table %in% display.name]
    ret <- structure(names(ret),names=lookup.table[names(ret)])
    ret <- ret[display.name]
  }
  return(ret)
}

#' Convert gene id to display name (gene symbol)
#' @importFrom SingleCellExperiment rowData
#' @param obj object of \code{SingleCellExperiment} class
#' @param ids character; gene ids
#' @return If successfull vector contains the disaply name; otherwise NULL
#' @export
ssc.id2displayName <- function(obj,ids)
{
  ret <- NULL
  if("display.name" %in% colnames(rowData(obj)))
  {
    lookup.table <- structure(rowData(obj)[,"display.name"],
                              names=rownames(obj))
    ret <- lookup.table[ids]
  }
  return(ret)
}

#' downsample
#' @param obj object of \code{singleCellExperiment} class
#' @param ncell.downsample integer; for each group, number of cells downsample to. (default: NULL)
#' @param group.var character; column in the colData(obj) used for grouping. (default: "majorCluster")
#' @param rn.seed integer; random number seed (default: 9999)
#' @param priority character; priority. (default: "")
#' @param rd character; data reducedDim(obj, rd) will be used for distance calculation (default: "pca")
#' @details return a downsampled object of \code{singleCellExperiment} class.
#' @importFrom data.table `:=` as.data.table
#' @importFrom matrixStats colMedians
#' @export
ssc.downsample <- function(obj, ncell.downsample=NULL, group.var="majorCluster",rn.seed=9999,
						   priority="",rd="pca")
{
    #### downsample cells
    set.seed(rn.seed)
	if(is.null(colnames(obj))){
		stop(sprintf("no colnames for obj\n"))
	}
    clust <- structure(obj[[group.var]],names=colnames(obj))
    grp.list <- unique(clust)
	if(priority=="silhouette"){
		sil <- ssc.plot.silhouette(obj,group.var,reducedDim.name=rd,do.plot=F)
		p.tb <- as.data.table(sil[,c(1:3)])
		p.tb$cellID <- colnames(obj)
		p.tb$group <- obj[[group.var]]
		p.tb[,silhouette.rank:=rank(-sil_width),by=c("group")]
		obj[[priority]] <- p.tb[["sil_width"]]
		obj[[sprintf("%s.rank",priority)]] <- p.tb[[sprintf("%s.rank",priority)]]
	}else if(priority=="distCenter"){
		dat.map <- reducedDim(obj,rd)
		p.tb <- as.data.table(ldply(grp.list,function(x){
						   xx <- names(clust[clust==x])
						   dat.block <- dat.map[xx,]
						   cc <- colMedians(dat.block)
						   out.tb <- data.table(cellID=xx, group=x,
												distCenter=sqrt(rowSums(sweep(dat.block,2,cc,"-")^2)))
						   out.tb[,distCenter.rank:=rank(distCenter)]
						   return(out.tb)
						   }))
		obj[[priority]] <- p.tb[[priority]][match(colnames(obj),p.tb$cellID)]
		obj[[sprintf("%s.rank",priority)]] <- p.tb[[sprintf("%s.rank",priority)]][match(colnames(obj),p.tb$cellID)]
	}
    if(!is.null(ncell.downsample)){
		if(priority==""){
			f.cell <- unlist(sapply(grp.list,function(x){
								 x <- names(clust[clust==x])
								 sample(x,min(length(x),ncell.downsample)) }))
		}else{
			f.cell <- p.tb$cellID[ p.tb[[sprintf("%s.rank",priority)]] <= ncell.downsample ]
		}
        obj <- obj[,f.cell]
    }
    return(obj)
}


#' convert assay data to long format
#' @param obj object of \code{singleCellExperiment} class
#' @param gene.id should be in rownames(obj); genes to plot
#' @param gene.symbol should be in rowData(obj)[,"display.name"]; genes to plot
#' @param assay.name character; which assay. NULL for all assays. (default: NULL)
#' @param col.idx character; output extra columns of the colData(obj). NULL for donnot output extra columns. (default: NULL)
#' @importFrom data.table as.data.table melt
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment rowData
#' @importFrom data.table `:=` as.data.table
#' @return Returns a data.table
#' @details convert assay data to table of long format. It can be specified which genes and assays will be contained in the long table.
#' @export
ssc.toLongTable <- function(obj,gene.id,gene.symbol,assay.name=NULL,col.idx=NULL)
{
	#requireNamespace("reshape2")
	if(missing(gene.id) && missing(gene.symbol)){
		warning("No gene.id or gene.symbol provided!")
		return(NULL)
	}
	if(is.null(names(rowData(obj)[,"display.name"]))){
		names(rowData(obj)[,"display.name"]) <- rownames(obj)
	}
	if(missing(gene.id) && !is.null(gene.symbol)){
		
		gene.list <- rowData(obj)[,"display.name"][which(rowData(obj)[,"display.name"] %in% gene.symbol)]
	}else{
		if(is.null(gene.id)){
			gene.id <- rownames(obj)
		}else{
			gene.id <- intersect(gene.id,rownames(obj))
		}
		gene.list <- rowData(obj)[,"display.name"][gene.id]
	}
	if(length(gene.list)==0){ 
		warning("No data found for the provided genes!")
		return(NULL)
	}
	obj <- obj[names(gene.list),]
	if(is.null(assay.name)){
		assay.name <- assayNames(obj)
	}

	dat.long <- NULL
	for(aname in intersect(assayNames(obj),assay.name)){
		dat.i <- as.data.table(reshape2::melt(assay(obj,aname)))
		colnames(dat.i) <- c("geneID","aid",aname)
		dat.i[,geneID:=as.character(geneID)]
		dat.i[,aid:=as.character(aid)]
		if(is.null(dat.long)){
			dat.long <- dat.i
		}else{
			dat.long <- merge(dat.long,dat.i,by=c("geneID","aid"))
		}
	}

	if(!is.null(dat.long) && !is.null(col.idx)){
		dat.extra.tb <- cbind(data.table(aid=colnames(obj)),
							  as.data.frame(colData(obj)[,col.idx,drop=F]))
		dat.long <- merge(dat.long,dat.extra.tb,by="aid")
	}

	return(dat.long)
}

#' Add module scores for gene expression programs in single cells. Mofidy from Seurat::AddModuleScore
#'
#' Calculate the average expression levels of each program (cluster) on single cell level,
#' subtracted by the aggregated expression of control feature sets.
#' All analyzed features are binned based on averaged expression, and the control features are
#' randomly selected from each bin.
#'
#' @param obj object of SingleCellExperiment
#' @param features Feature expression programs in named list
#' @param pool List of features to check expression levels agains, defaults to \code{rownames(x = object)}
#' @param nbin Number of bins of aggregate expression levels for all analyzed features
#' @param ctrl Number of control features selected from the same bin per analyzed feature
#' @param assay.name Name of assay to use
#' @param adjB character; batch column of the colData(obj). (default: NULL)
#' @param do.scale logical; scale the data (default: TRUE)
#' @param seed Set a random seed
#' @return Returns a SingleCellExperiment object with module scores added to object meta data
#' @importFrom ggplot2 cut_number
#' @importFrom Matrix rowMeans colMeans
#' @importFrom plyr llply laply
#' @importFrom stats rnorm
#' @references Tirosh et al, Science (2016); Seurat's source code
#' @export
ssc.moduleScore <- function(obj, features, pool = NULL,
							nbin = 24, ctrl = 100, assay.name = "exprs",
							adjB=NULL,do.scale=T,seed = 1)
{
	set.seed(seed = seed)
	if(is.null(names(features))){
		names(features) <- sprintf("M%02d",seq_along(features))
	}
	features <- llply(features,function(x){ intersect(x,rowData(obj)$display.name) })
	if(any(sapply(features,length)<1)){
		warning(sprintf("no expression data of the provided genes!\n"))
		return(obj)
	}
	features <- llply(features,function(x){ rownames(obj)[match(x,rowData(obj)$display.name)] })
	if(is.null(pool)){
		pool <- rownames(obj)
	}else{
		pool <- rownames(obj)[match(pool,rowData(obj)$display.name)]
	}
	assay.data <- assay(obj,assay.name)
	data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
	data.avg <- data.avg[order(data.avg)]
	data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30,
						   n = nbin, labels = FALSE, right = FALSE)
	###data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg) / (nbin + 1))))
	names(x = data.cut) <- names(x = data.avg)

	ctrl.use <- llply(features,function(x){
							  feat.rnd <- c()
							  for(i in seq_along(x)){
									feat.rnd <- c(feat.rnd,
												  names(sample(data.cut[which(data.cut==data.cut[x[i]])],size=ctrl,replace=F))
												  )
							  }
							  return(feat.rnd)
						   })
    .calScore <- function(x){
		 #### adjB and  scale
		 if(!is.null(adjB)){
			 dat.use <- simple.removeBatchEffect(assay.data[x,,drop=F],batch = obj[[adjB]])
		 }else{
			 dat.use <- t(scale(t(assay.data[x,,drop=F]),scale=F))
		 }
		 if(do.scale){
			 dat.use <- t(scale(t(dat.use)))
		 }
		 ###
		 Matrix::colMeans(dat.use)
	}
	ctrl.scores <- laply(ctrl.use,.calScore,.drop=F)
	features.scores <- laply(features,.calScore,.drop=F)
	features.scores.use <- features.scores - ctrl.scores
	rownames(features.scores.use) <- names(features)
	features.scores.use <- as.data.frame(x = t(x = features.scores.use))
	obj[[colnames(x = features.scores.use)]] <- features.scores.use[[1]]
	return(obj)
}

#' scale the data
#' @param obj object of \code{singleCellExperiment} class
#' @param gene.id should be in rownames(obj); genes to plot
#' @param gene.symbol should be in rowData(obj)[,"display.name"]; genes to plot
#' @param assay.name character; which assay (default: "exprs")
#' @param adjB character; batch column of the colData(obj). (default: NULL)
#' @param do.scale logical; whether scale the data. (default: F)
#' @details scale the data specified in assay.name, then store the scaled data to ${assay.name}.scale. One of gene.id and gene.symbol must be provided.
#' @importFrom SummarizedExperiment `assay<-` assay
#' @export
ssc.scale <- function(obj,gene.id,gene.symbol,assay.name="norm_exprs",adjB=NULL,do.scale=F)
{
	if(missing(gene.id) && missing(gene.symbol)){
		warning("No gene.id or gene.symbol provided!")
		return(NULL)
	}
	if(is.null(names(rowData(obj)[,"display.name"]))){
		names(rowData(obj)[,"display.name"]) <- rownames(obj)
	}
	if(missing(gene.id) && !is.null(gene.symbol)){
		
		gene.list <- rowData(obj)[,"display.name"][which(rowData(obj)[,"display.name"] %in% gene.symbol)]
	}else{
		gene.id <- intersect(gene.id,rownames(obj))
		gene.list <- rowData(obj)[,"display.name"][gene.id]
	}
	if(length(gene.list)==0){ 
		warning("No data found for the provided genes!")
		return(NULL)
	}
	obj <- obj[names(gene.list),]
	dat.block <- assay(obj,assay.name)
	if(!is.null(adjB)){
		dat.block <- simple.removeBatchEffect(dat.block,batch=obj[[adjB]])
	}
	if(do.scale){
		dat.block <- t(scale(t(dat.block)))
	}
	assay(obj,sprintf("%s.scale",assay.name)) <- dat.block
	return(obj)
}


#' scale the assay per gene
#' @param obj object of \code{singleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param assay.new character; assay name to store the scaled- expression;if NULL, will be (assay.name).z (default: NULL)
#' @param covar character; perform scale in cells grouping by covar (default: "patient")
#' @param n.cores integer; number of cores used, if NULL it will be determined automatically (default: NULL)
#' @param z.lo double; z-score lower boundary; if set, z-score lower than this will be set to this (default: NULL)
#' @param z.hi double; z-score higher boundary; if set, z-score higher than this will be set to this (default: NULL)
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom doParallel registerDoParallel
#' @importFrom data.table as.data.table
#' @importFrom plyr ldply
#' @importFrom stats sd
#' @importFrom SummarizedExperiment `assay<-` assay
#' @export
ssc.assay.zscore <- function(obj,assay.name="exprs",assay.new=NULL,covar="patient",n.cores=NULL,
                             z.lo=NULL,z.hi=NULL)
{
    #requireNamespace("data.table")

    dat.plot <- as.matrix(assay(obj,assay.name))

    cell.info <- as.data.table(colData(obj)[covar])
    cell.info$cellID <- colnames(obj)
    cell.info.split <- split(cell.info,by=covar,sorted=T)

    RhpcBLASctl::omp_set_num_threads(1)
    doParallel::registerDoParallel(cores = n.cores)
    dat.plot.z.df <- data.table(ldply(names(cell.info.split),function(x){
                dat.block <- dat.plot[,cell.info.split[[x]][["cellID"]],drop=F]
                #dat.block <- scale(t(dat.block))
                rowM <- rowMeans(dat.block, na.rm = T)
                rowSD <- apply(dat.block, 1, sd, na.rm = T)
                dat.block <- sweep(dat.block, 1, rowM)
                dat.block <- sweep(dat.block, 1, rowSD, "/")
                dat.block <- t(dat.block)
                dat.block.df <- cbind(cellID=rownames(dat.block),as.data.frame(dat.block),stringsAsFactors=F)
                return(dat.block.df)
            },.progress = "none",.parallel=T))
    dat.plot.z <- as.matrix(dat.plot.z.df[,-c("cellID")])
    rownames(dat.plot.z) <- dat.plot.z.df$cellID
    dat.plot.z <- t(dat.plot.z)

    if(!is.null(z.lo) && !is.null(z.hi)){
        dat.plot.z[dat.plot.z < z.lo] <- z.lo
        dat.plot.z[dat.plot.z > z.hi] <- z.hi
    }
    assay(obj,if(is.null(assay.new)) sprintf("%s.z",assay.name) else assay.new) <- dat.plot.z[,colnames(obj)]
    return(obj)
}

