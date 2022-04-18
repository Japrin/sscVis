#' Wraper for silhouette()
#' @importFrom cluster silhouette
#' @importFrom graphics plot
#' @importFrom stats dist
#' @param obj object of \code{SingleCellExperiment}
#' @param cluster.label character; which column of colData of obj to used as cluster label.
#' @param reducedDim.name character; which reducedDim to use. (default: "iCor.tsne")
#' @param do.plot logical; whether plot
#' @param ... Arguments to be passed to plot()
#' @return an object, sil, of class \code{silhouette}
#' @export
ssc.plot.silhouette <- function(obj,cluster.label,reducedDim.name="iCor.tsne",do.plot=T, ...){
  dist.obj <- dist(reducedDim(obj,reducedDim.name))
  sil <- cluster::silhouette(as.numeric(as.factor(colData(obj)[,cluster.label])),dist.obj)
  if(do.plot){
    plot(sil, ...)
  }
  return(sil)
}


#' Plot on tSNE map
#' @param obj object of \code{SingleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param gene, character; genes to be showed. (default: NULL)
#' @param columns character; columns in colData(obj) to be showd. (default: NULL)
#' @param splitBy character; columns in colData(obj). Split the dataset to mupltiple subset then plot them one by one (default: NULL)
#' @param plotDensity logical; whether plot 2D density. (default F)
#' @param colSet list; mapping iterms in the names to colors in the values. (default: list())
#' @param reduced.name character; names in the reducedDimNames. (default: "iCor.tsne")
#' @param reduced.dim integer; which dimensions of the reduced data to be used. (default: c(1,2))
#' @param out.prefix character; output prefix. (default: NULL)
#' @param p.ncol integer; number of columns in the figure layout. (default: 3)
#' @param width numeric; width of the plot, used for geneOnTSNE. (default: NA)
#' @param height numeric; height of the plot, used for geneOnTSNE. (default: NA)
#' @param base_aspect_ratio numeric; base_aspect_ratio, used for plotting metadata. (default 1.1)
#' @param peaks integer or character; index or names of the peaks. (default: NULL)
#' @param xlim integer or NULL; only draw points lie in the ragne specified by xlim and ylim (default NULL)
#' @param ylim integer or NULL; only draw points lie in the ragne specified by xlim and ylim (default NULL)
#' @param size double; points' size. If NULL, infer from number of points (default NULL)
#' @param palette.name character; which palette to use. (default: "YlOrRd")
#' @param adjB character; batch column of the colData(obj). (default: NULL)
#' @param clamp integer vector; expression values will be clamped to the range defined by this parameter, such as c(0,15). (default: "none" )
#' @param do.scale logical; whether scale the expression value. (default: FALSE)
#' @param label double; label size. if NULL, no label showed. (default: NULL )
#' @param par.repel list; passed to geom_text_repel
#' @param par.geneOnTSNE character; other parameters of geneOnTSNE
#' @param my.ggPoint function; used to plot scatter plot, such as geom_point, geom_point_rast, geom_scattermore, geom_scattermost (default: geom_point)
#' @param par.geom_point list; extra parameters for geom_point/geom_point_rast; (default: list())
#' @param par.legend list; lengend parameters, used to overwrite the default setting; (default: list())
#' @param show.legend logical; if NULL, determined automatically; (default: NULL)
#' @param theme.use function; which theme to use (default: theme_bw)
#' @param legend.w numeric; adjust legend width (default: 1)
#' @param fun.extra function;  (default: NULL)
#' @param verbose logical;  (default: FALSE)
#' @importFrom SingleCellExperiment colData reducedDim
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual theme_bw aes_string guides guide_legend coord_cartesian
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggrastr geom_point_rast
#' @importFrom scattermore geom_scattermore geom_scattermost
#' @importFrom cowplot save_plot plot_grid
#' @importFrom utils read.table
#' @importFrom RColorBrewer brewer.pal
#' @importFrom data.table as.data.table
#' @importFrom SummarizedExperiment assay
#' @importFrom stats median
#' @importFrom grDevices pdf png dev.off
#' @importFrom graphics par plot.new title
#' @details If `gene` is not NULL, expression of the specified genes will be plot on the tSNE map; if columns in not
#' NULL, colData of obj with names in `columns` will be plot on the tSNE map. The tSNE map used is specified by option
#' `reduced.name` and `reduced.dim`. Both `gene` and `columns` can be non-NULL. For list `colSet`, each element define
#' a color mapping for the responding iterm in the `column`; if not specifed, automatically generated color mapping will
#' be used.
#' @export
ssc.plot.tsne <- function(obj, assay.name="exprs", gene=NULL, columns=NULL,splitBy=NULL,
                             plotDensity=F, colSet=list(),
                             reduced.name="iCor.tsne",reduced.dim=c(1,2),xlim=NULL,ylim=NULL,size=NULL,
                             palette.name="YlOrRd",adjB=NULL,clamp="none",do.scale=FALSE,
                             label=NULL,par.repel=list(force=1),
                             #vector.friendly=F,par.geom_point=list(),
                             my.ggPoint=geom_point,par.geom_point=list(),
                             par.legend=list(),show.legend=NULL,
                             theme.use=theme_bw,legend.w=1,verbose=F,fun.extra=NULL,
                             par.geneOnTSNE=list(scales="free",pt.order="value",pt.alpha=0.1),
                             out.prefix=NULL,p.ncol=3,width=NA,height=NA,base_aspect_ratio=1.1,peaks=NULL)
{
  #requireNamespace("ggplot2")
  #requireNamespace("cowplot")
  if(length(reduced.dim)!=2){ warning(sprintf("Wrong parameter, reduced.dim!!")); return(); }
  if(is.null(reducedDim(obj,reduced.name))){
    warning(sprintf("No reducedDim: %s\n",reduced.name))
    return()
  }
  dat.map <- reducedDim(obj,reduced.name)[,reduced.dim]
  if(!is.null(columns))
  {
    if(all(columns %in% colnames(colData(obj))))
    {
      if(is.list(colSet)){
        multi.p <- lapply(columns,function(cc){
          if(is.null(colSet[[cc]])){
            if(is.null(metadata(obj)$ssc$colSet[[cc]])){
                cc.values <- sort(unique(colData(obj)[,cc]))
                colSet[[cc]] <- structure(auto.colSet(length(cc.values),name = "Paired"),
                                          names=as.character(cc.values))
            }else{
                colSet[[cc]] <- metadata(obj)$ssc$colSet[[cc]]
            }
          }
          dat.plot <- data.frame(sample=rownames(dat.map),stringsAsFactors = F)
          if(!is.null(splitBy)){
            dat.plot <- as.data.frame(cbind(dat.plot,dat.map,colData(obj)[,c(cc,splitBy),drop=F]))
            colnames(dat.plot) <- c("sample","Dim1","Dim2",cc,"splitBy")
          }else{
            dat.plot <- as.data.frame(cbind(dat.plot,dat.map,colData(obj)[,cc,drop=F]))
            colnames(dat.plot) <- c("sample","Dim1","Dim2",cc)
          }
	  if(par.geneOnTSNE$pt.order=="value"){
	      dat.plot <- dat.plot[order(dat.plot[,cc]),]
	  }else if(par.geneOnTSNE$pt.order=="random"){
	      dat.plot <- dat.plot[sample(nrow(dat.plot),nrow(dat.plot)),]
	  }
          npts <- nrow(dat.plot)
          if(is.numeric(dat.plot[,cc])){
            nvalues <- Inf
            if(clamp!="none"){
                dat.plot[[cc]][ dat.plot[[cc]] < clamp[1] ] <- clamp[1]
                dat.plot[[cc]][ dat.plot[[cc]] > clamp[2] ] <- clamp[2]
            }
          }else{
            nvalues <- length(unique(dat.plot[,cc]))
          }
	  #if(vector.friendly){
	  #    my.ggPoint <- geom_point_rast
	  #    ##my.ggPoint <- scattermore::geom_scattermore
	  #}else{
	  #    my.ggPoint <- geom_point
	  #}
          p <- ggplot2::ggplot(dat.plot,aes(Dim1,Dim2)) +
			  do.call(my.ggPoint,c(list(mapping=aes_string(colour=cc),
									  show.legend=if(is.null(show.legend)) { if(!is.numeric(dat.plot[,cc]) && nvalues>40) F else NA } else show.legend,
									  size=if(is.null(size)) auto.point.size(npts)*1.1 else size),
								   par.geom_point)) +
            #my.ggPoint(aes_string(colour=cc),
            #           show.legend=if(!is.numeric(dat.plot[,cc]) && nvalues>40) F else NA,
            #           size=if(is.null(size)) auto.point.size(npts)*1.1 else size) +
            labs(x=sprintf("Dim%d",reduced.dim[1]),y=sprintf("Dim%d",reduced.dim[2]))
          if(!is.null(label)){
              dat.plot.label <- as.data.table(dat.plot)[,.(Dim1=median(.SD$Dim1),
                                                           Dim2=median(.SD$Dim2)),
                                                        by=cc]
              p <- p + do.call(ggrepel::geom_text_repel,c(list(aes_string("Dim1","Dim2",label = cc),
                                                               size=if(length(label) > 1) label[cc] else label,data=dat.plot.label),
                                                          par.repel))
          }
          if(!is.null(splitBy)){
            p <- p + ggplot2::facet_wrap(~splitBy,ncol = if(length(p.ncol)>1) p.ncol[2] else NULL)
          }
          if(is.numeric(dat.plot[,cc])){
            p <- p + do.call(scale_colour_gradientn,
                             c(list(colours = getColorPaletteFromNameContinuous(palette.name)),
			       par.legend))
          }else{
            p <- p + do.call(scale_colour_manual,c(list(values = colSet[[cc]]),par.legend))
          }
	  p <- p + theme.use() + labs(title=cc) + theme(plot.title = element_text(hjust = 0.5))
	  legend.ncol <- if(nvalues>10 && !is.infinite(nvalues)) ceiling(nvalues/10) else NULL
          p <- p + coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE)
	  if(!is.numeric(dat.plot[,cc])){
	      p <- p + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=if(nvalues<=26) 4 else 2.0),
								      label.theme = element_text(size=8)), ncol=legend.ncol)
	  }
	  if(legend.w==0){ p <- p + theme(legend.position="none") }

	  if(!is.null(fun.extra)){
	      p <- fun.extra(p)
	      #p <- do.call(`+`,list(p,fun.extra()))
	  }

          return(p)
        })
        pp <- cowplot::plot_grid(plotlist=multi.p,
                                 ncol = if(length(columns)>1) p.ncol[1] else 1,align = "hv")
        if(!is.null(out.prefix)){
          cowplot::save_plot(sprintf("%s.columnsOntSNE.pdf",out.prefix),pp,
                             ncol = if(length(columns)>1) p.ncol[1] else 1,
                             base_aspect_ratio=base_aspect_ratio,
							 base_height=if(!is.na(height)) height else 3.71)
        }else{
	    if(verbose){
		    return(list("plot"=pp,"list"=multi.p))
	    }else{
		    return(pp)
	    }
        }
      }else{
        warning(sprintf("invalidate parameter: colSet. Please check that!"))
      }
    }else{
      warning(sprintf("some columns not in the data. Not plot be produced!"))
    }
  }
  if(!is.null(gene)){
    if(length(gene)==1 && file.exists(gene)){ gene <- read.table(gene,header = T)[,1] }
    #if(all(!(gene %in% rownames(obj)))){
    #  ### gene symbol?
    #  gene <- ssc.displayName2id(obj,display.name = gene)
    #}
    #if(is.null(names(gene))){
    #  names(gene) <- gene
    #}
    gene <- ssc.displayName2id(obj,display.name = gene)

    dat.onTSNE <- assay(obj,assay.name)[gene,,drop=F]
    if(!is.null(adjB)){
      dat.onTSNE <- simple.removeBatchEffect(dat.onTSNE,batch=colData(obj)[[adjB]])
    }
    if(do.scale){
      dat.onTSNE <- t(scale(t(dat.onTSNE)))
    }
    p <- do.call(ggGeneOnTSNE,c(list(Y=dat.onTSNE, dat.map=dat.map, gene.to.show=gene,
                                     p.ncol=p.ncol,xlim=xlim,ylim=ylim,
                                     size=size,width=width,height=height,palette.name=palette.name,
                                     clamp=clamp,my.ggPoint=my.ggPoint,theme.use=theme.use,
                                     par.geom_point=par.geom_point,par.legend=par.legend,
                                     splitBy=if(is.null(splitBy)) NULL else obj[[splitBy]],
				     fun.extra=fun.extra,
				     verbose=verbose,
                                     out.prefix=out.prefix),
                                par.geneOnTSNE))
    if(is.null(out.prefix)){
      #print(p)
      return(p)
    }
  }
  if(plotDensity){
    if(is.null(out.prefix)){
      plotDensity2D(dat.map,peaks = peaks)
    }else{
      pdf(sprintf("%s.density.pdf",out.prefix),width = 5,height = 5)
      plotDensity2D(dat.map,peaks = peaks)
      dev.off()
    }
  }
}


#' Plot violin
#' @param obj object of \code{SingleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param gene character; genes to be showed. (default: NULL)
#' @param columns character; columns in colData(obj) to be showd. (default: NULL)
#' @param group.var character; column in the colData(obj) used for grouping. (default: "majorCluster")
#' @param group.in character; only thoes groups to be shown. NULL for all groups. (default: NULL)
#' @param splitBy character; columns in colData(obj). Split the dataset to mupltiple subset then plot them one by one (default: NULL)
#' @param clamp integer vector; expression values will be clamped to the range defined by this parameter. (default: c(0,12))
#' @param out.prefix character; output prefix. (default: NULL)
#' @param p.ncol integer; number of columns in the figure layout. (default: 3)
#' @param base_aspect_ratio numeric; base_aspect_ratio, used for plotting metadata. (default 1.1)
#' @param adjB character; batch column of the colData(obj). (default: NULL)
#' @param do.scale logical; whether scale the expression value. (default: FALSE)
#' @param par.legend list; lengend parameters, used to overwrite the default setting; (default: list())
#' @param par.boxplot list; geom_boxplot parameters. (default: list(outlier.shape = NA,width=0.25,alpha=0.8))
#' @param par.text list; geom_text parameters. (default: list(vjust=0))
#' @param palette.name character; which palette to use. (default: "YlOrRd")
#' @param angle.axis.x numeric; rotation angle. (default 60)
#' @param add character; other plots to add. (default NULL)
#' @param ... parameter passed to cowplot::save_plot
#' @importFrom SingleCellExperiment colData reducedDim
#' @importFrom ggplot2 ggplot aes geom_violin scale_fill_gradientn theme_bw theme aes_string facet_grid element_text geom_boxplot scale_colour_brewer
#' @importFrom cowplot save_plot plot_grid
#' @importFrom data.table melt data.table
#' @importFrom SummarizedExperiment assay
#' @importFrom grDevices pdf png dev.off
#' @importFrom graphics par plot.new title
#' @importFrom S4Vectors `metadata<-` metadata
#' @importFrom stats quantile
#' @importFrom Matrix t
#' @details If `gene` is not NULL, violin of the genes' expression will be plot; if columns in not
#' NULL, colData of obj with names in `columns` will be plot in violin.
#' @export
ssc.plot.violin <- function(obj, assay.name="exprs", gene=NULL, columns=NULL,par.legend=list(),
                            group.var="majorCluster",group.in=NULL,splitBy=NULL,
                            clamp=c(0,12),adjB=NULL,do.scale=F,
                            angle.axis.x=60,add=NULL,
                            par.boxplot=list(outlier.shape = NA,width=0.25,alpha=0.8),
                            par.text=list(vjust = 0),
                            palette.name="YlOrRd",
                            out.prefix=NULL,p.ncol=1,base_aspect_ratio=1.1,...)
{
  #requireNamespace("ggplot2")
  #requireNamespace("data.table")
  if(!is.null(gene))
  {
	  gene <- ssc.displayName2id(obj,display.name = gene)

	  dat.violin <- assay(obj,assay.name)[gene,,drop=F]
	  if(!is.null(adjB)){
		dat.violin <- simple.removeBatchEffect(dat.violin,batch=colData(obj)[[adjB]])
	  }
      if(do.scale){
        dat.violin <- t(scale(Matrix::t(dat.violin)))
      }

	  dat.plot <- as.matrix(Matrix::t(dat.violin))
	  colnames(dat.plot) <- ssc.id2displayName(obj,colnames(dat.plot))
	  dat.plot.df <- data.table::data.table(sample=rownames(dat.plot),stringsAsFactors = F)
	  dat.plot.df <- cbind(dat.plot.df,as.data.frame(colData(obj)[,group.var,drop=F]))
	  dat.plot.df <- cbind(dat.plot.df,dat.plot)
	  if(!is.null(splitBy)){
	    dat.plot.df <- cbind(dat.plot.df,data.table(splitBy=obj[[splitBy]]))
	    dat.plot.df <- data.table::melt(dat.plot.df,id.vars=c("sample",group.var,"splitBy"),
	                                    variable.name="gene",value.name=assay.name)
	    dat.plot.df.grpMean <- dat.plot.df[,lapply(.SD,mean,na.rm=T),by=c("gene",group.var,"splitBy"),.SDcols=assay.name]
	    colnames(dat.plot.df.grpMean) <- c("gene",group.var,"splitBy","meanExp")
	    dat.plot.df <- dat.plot.df.grpMean[dat.plot.df,,on=c("gene",group.var,"splitBy")]
	  }else{
	    dat.plot.df <- data.table::melt(dat.plot.df,id.vars=c("sample",group.var),
	                                    variable.name="gene",value.name=assay.name)  
	    dat.plot.df.grpMean <- dat.plot.df[,lapply(.SD,mean,na.rm=T),by=c("gene",group.var),.SDcols=assay.name]
	    colnames(dat.plot.df.grpMean) <- c("gene",group.var,"meanExp")
	    dat.plot.df <- dat.plot.df.grpMean[dat.plot.df,,on=c("gene",group.var)]
	  }
	  dat.plot.df.grpMean$meanExp.label <- sprintf("%0.2f",dat.plot.df.grpMean$meanExp)
      dat.plot.df[,gene:=factor(gene,levels=colnames(dat.plot),ordered=T)]

      if(is.null(clamp)){
          clamp <- quantile(dat.plot.df[[assay.name]],c(0.05,0.95))
      }
      dat.plot.df[meanExp<clamp[1],meanExp:=clamp[1],]
      dat.plot.df[meanExp>clamp[2],meanExp:=clamp[2],]
      dat.plot.df[[assay.name]][ dat.plot.df[[assay.name]] < clamp[1] ] <- clamp[1]
      dat.plot.df[[assay.name]][ dat.plot.df[[assay.name]] > clamp[2] ] <- clamp[2]

	  if(!is.null(group.in)){
		  dat.plot.df <- dat.plot.df[dat.plot.df[[group.var[1]]] %in% group.in,]
	  }
	  p <- ggplot(dat.plot.df, aes_string(group.var[1], assay.name))
	  if(length(group.var)==1){
  		p <- p +
  		  geom_violin(scale = "width",aes(fill=meanExp),color=NA,show.legend = T) +
  		  do.call(scale_fill_gradientn,
  		                 c(list(colours = getColorPaletteFromNameContinuous(palette.name),
  		                        limits=clamp),
  		                   par.legend))
  		#   do.call(scale_fill_gradient2,
  		#           c(list(low = "yellow",mid = "red",high = "black",midpoint = mean(clamp), limits=clamp),
  		# 								 par.legend))
	  }else if(length(group.var)==2)
	  {
  		p <- p +
  			geom_boxplot(aes_string(colour = group.var[2])) +
  			scale_colour_brewer(palette = "Set1")
  			#geom_violin(scale = "width",aes_string(fill="meanExp",linetype=group.var[2],color=group.var[2]),
  			#            show.legend = T) +
  			#scale_fill_gradient2(low = "yellow",mid = "red",high = "black",midpoint = mean(clamp),limits=clamp)
	  }
	  if(!is.null(add)){
	    if("boxplot" %in% add){
	      p <- p + do.call(geom_boxplot,par.boxplot)
	    }
	    if("text" %in% add){
	      p <- p + do.call(geom_text,c(list(data=dat.plot.df.grpMean,
	                                        mapping=aes_string(y="meanExp",label="meanExp.label")),
	                                   par.text))
	    }
	  }
	  if(!is.null(splitBy)){
	    p <- p + facet_grid(splitBy~gene,scales="free_y")
	  }else{
	    p <- p + facet_wrap(gene~.,strip.position = "left",scales="free_y",dir="v",ncol=p.ncol)
	      #		facet_grid(gene ~ .,switch = "y",scales = "free_y")
	  }
	  p <- p + theme_bw(base_size = 12) + 
	    theme(axis.text.x = element_text(angle = angle.axis.x, hjust = 1),
	          #strip.background = element_blank(),
	          strip.placement = "inside")
  } else if(!is.null(columns)){
	  dat.plot.df <- as.data.table(cbind(data.frame(cellID=colnames(obj),stringsAsFactors=F),
						  as.data.frame(colData(obj)[,c(group.var,columns),drop=F])))
	  dat.plot.df <- melt(dat.plot.df,id.vars=c("cellID",group.var))
	  if(!is.null(group.in)){
		  dat.plot.df <- dat.plot.df[dat.plot.df[[group.var[1]]] %in% group.in,]
	  }
	  p <- ggplot(dat.plot.df, aes_string(group.var[1], "value"))
      if(length(group.var)==1){
          p <- p + geom_boxplot()
      }else if(length(group.var)==2){
          p <- p + geom_boxplot(aes_string(colour=group.var[2]))
          if(is.null(metadata(obj)$ssc$colSet)){
              p <- p + scale_colour_brewer(palette = "Set1")
          }else{
              p <- p + scale_colour_manual(values = metadata(obj)$ssc$colSet[[group.var[2]]])
          }
      }
      p <- p + theme_bw(base_size = 12) +
			facet_grid(variable ~ .,switch = "y",scales = "free_y") +
			theme(axis.text.x = element_text(angle = 60, hjust = 1),strip.placement = "inside")
  }
  if(!is.null(out.prefix)){
    cowplot::save_plot(sprintf("%s.violin.%s.pdf",out.prefix,if(!is.null(gene)) "gene" else "columns"),p,
                       ncol = p.ncol,
                       base_aspect_ratio=base_aspect_ratio,...)
  }else{
    return(p)
  }
}


#' Plot pca result, such as scree plot.
#' @param obj object of \code{SingleCellExperiment} class
#' @param out.prefix character; output prefix. (default: NULL)
#' @param p.ncol integer; number of columns in the figure layout. (default: 2)
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual theme_bw ylab
#' @export
ssc.plot.pca <- function(obj, out.prefix=NULL,p.ncol=2)
{
  #requireNamespace("ggplot2")
  eigenv <- metadata(obj)$ssc$pca.res$eigenv.prop * 100
  dat.plot.eigenv <- data.frame(PC=seq_along(eigenv),
                                eigenv=eigenv,
                                isKneePts=as.character(seq_along(eigenv)==metadata(obj)$ssc$pca.res$kneePts),
                                stringsAsFactors = F)
  p <- ggplot2::ggplot(head(dat.plot.eigenv,n=30),mapping = aes(PC,eigenv)) +
    geom_point(aes(colour=isKneePts),show.legend=F) + ylab("Variation explained (%)") +
    scale_colour_manual(values = c("TRUE"="#E41A1C","FALSE"="#377EB8")) +
    theme_bw()
  print(p)
}

#' Plot correlation.
#' @param obj object of \code{SingleCellExperiment} class
#' @param feat1 character; feature on x axis
#' @param feat2 character; feature on y axis
#' @param type1 character; feature type on x axis. (default: "gene")
#' @param type2 character; feature type on y axis. (default: "gene")
#' @param assay.name character; which assay (default: "exprs")
#' @param legend.w numeric; adjust legend width (default: 1)
#' @param legend.ncol integer;  (default: NULL)
#' @param add.legend logical;  (default: F)
#' @param out.prefix character; output prefix. (default: NULL)
#' @param vector.friendly logical; output vector friendly figure (default: FALSE)
#' @param ... parameter passed to geom_point
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab geom_smooth
#' @importFrom ggpubr stat_cor
#' @importFrom ggrastr geom_point_rast
#' @importFrom SummarizedExperiment assay
#' @export
ssc.plot.cor <- function(obj,feat1,feat2,type1="gene",type2="gene", assay.name="exprs",
						 legend.w=1,legend.ncol=NULL,add.legend=F, out.prefix=NULL,
						 vector.friendly=F,...)
{
  #requireNamespace("ggplot2")
  dat.plot <- data.table(cellID=colnames(obj),
						  feat1=if(type1=="gene") assay(obj,assay.name)[feat1,] else obj[[feat1]],
						  feat2=if(type2=="gene") assay(obj,assay.name)[feat2,] else obj[[feat2]])
  dat.ext <- colData(obj)
  dat.ext <- dat.ext[,setdiff(colnames(dat.ext),c(feat1,feat2)),drop=F]
  dat.plot <- cbind(dat.plot,as.data.frame(dat.ext))
  dat.plot <- dat.plot[!is.na(feat1) & !is.na(feat2),]

  ##p <- ggscatter(dat.plot,x="feat1",y="feat2",...) +
  if(vector.friendly){
	  my.ggPoint <- geom_point
  }else{
	  my.ggPoint <- geom_point_rast
  }
  p <- ggplot(dat.plot,aes_string(x="feat1",y="feat2")) +
	  my.ggPoint(...) +
    ylab(feat1) + xlab(feat2)
  if(add.legend==T){
	  p <- p + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=4.0),
														   ncol=legend.ncol,
                                                           label.theme = element_text(size=8)))
  }
  p <- p + geom_smooth(method='lm') + stat_cor()
###  if(!vector.friendly){
###	p <- p + geom_smooth(method='lm') + stat_cor()
###  }else{
###	  tmpfilename <- sprintf("%s.tmp.cor.%s.%s.png",if(!is.null(out.prefix)) out.prefix else "",feat1,feat2)
###	  ggsave(filename=tmpfilename,
###			 plot = p + theme_void() + theme(legend.position = "none",
###											 axis.line.x = element_blank(), axis.line.y = element_blank(),
###											 axis.title.x = element_blank(), axis.title.y = element_blank(),
###											 axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
###											 plot.title = element_blank()),
###			 width=8,height=6)
###	  pbuild.params <- ggplot_build(plot = p)$layout$panel_params[[1]]
###	  range.values <- c( pbuild.params$x.range, pbuild.params$y.range)
###	  img <- png::readPNG(source = tmpfilename)
###	  blank <- ggplot(data = p$data,mapping = aes(feat1,feat2)) + geom_blank()
###	  blank <- blank + p$theme +
###		  ylab(feat1) + xlab(feat2) +
###		  coord_cartesian(xlim = range.values[1:2], ylim = range.values[3:4], expand = F)
###	  blank <- blank + annotation_raster(raster = img,
###										 xmin = range.values[1], xmax = range.values[2],
###										 ymin = range.values[3], ymax = range.values[4])
###	  blank <- blank + geom_smooth(method='lm') + stat_cor()
###	  legend.blank <- cowplot::get_legend(p)
###	  if(legend.w==0){
###		  p <- blank + theme(legend.position="none")
###	  }else{
###		  #p <- cowplot::plot_grid(blank, legend.blank, rel_widths = c(4, 1*legend.w))
###		  p <- cowplot::plot_grid(blank, legend.blank, rel_widths = c(4, if(is.null(legend.ncol)) 1*legend.w else legend.ncol*legend.w))
###	  }
####	  p <- blank
###	  file.remove(tmpfilename)
###  }
  return(p)
}


#' plot heatmap
#' @param obj object of \code{SingleCellExperiment} class
#' @param assay.name character; which assay (default: "exprs")
#' @param out.prefix character; output prefix. (default: NULL)
#' @param ncell.downsample integer; number of cells downsample to. (default: NULL)
#' @param ave.by character; average the expression profile grouping by this. (default: NULL)
#' @param columns character; columns in colData(obj) to be showd. must be subset of columns of colData(obj) and ave.by (if it's not NULL) (default: NULL)
#' @param columns.order character; columns of colData(obj) used for ordering (default: NULL)
#' @param gene.desc data.frame; it must contain columns geneID and Group (default: NULL)
#' @param colSet list; mapping iterms in the names to colors in the values. (default: list())
#' @param pdf.width double; width of the pdf file. (default:16)
#' @param pdf.height double; height of the pdf file. (default:15)
#' @param do.scale logical; wheter scale the rows, just for visualization. (default: TRUE)
#' @param z.lo double; z-score lower boundary; z-score lower than this will be set to this (default: -2.5)
#' @param z.hi double; z-score higher boundary; z-score higher than this will be set to this (default: 2.5)
#' @param z.step double; z-score step, used for coloring the expression value (default: 1)
#' @param exp.title character; title for the expression legend (default: "Exp")
#' @param do.clustering.row logical; wheter order row (default: TRUE)
#' @param do.clustering.col logical; wheter order columns (default: TRUE)
#' @param dend.col dendrogram of the columns, 'cluster_columns' of ComplexHeatmap::Heatmap  (default: FALSE)
#' @param dend.row dendrogram of the rows, 'cluster_rows' of ComplexHeatmap::Heatmap (default: FALSE)
#' @param clustering.distance character; one of spearmn, pearson, cosine and euclidean (default: "spearman")
#' @param clustering.method character; method for hclust (default: "complete")
#' @param k.row integer; number of clusters in the rows (default: 1)
#' @param k.col integer; number of clusters in the columns (default: 1)
#' @param returnHT logical; whether return HT; (default: FALSE)
#' @param palette.name character; which palette to use, such as "RdBu","RdYlBu" (default: NULL)
#' @param palette.ann.numeric character; which palette to use, such as "RdBu","RdYlBu" (default: "RdYlBu")
#' @param Y.level.ann.numeric vector; value range for numeric annotation (default: NULL)
#' @param row.split vector; used for row; must be named or is corresponding to the rows of obj (default: NULL)
#' @param column.split vector; used for column; (default: NULL)
#' @param annotation_legend_param list; (default: list())
#' @param ann.bar.height double; height of the top annotation (default: 1.5)
#' @param mytitle character; title of the figure (default: "")
#' @param ... parameters pass to plotMatrix.simple()
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap decorate_annotation
#' @importFrom circlize colorRamp2
#' @importFrom gridBase baseViewports
#' @importFrom grid pushViewport grid.text gpar unit
#' @importFrom RColorBrewer brewer.pal
#' @importFrom SingleCellExperiment rowData colData
#' @importFrom SummarizedExperiment `rowData<-` `colData<-`
#' @importFrom grDevices pdf png dev.off
#' @importFrom graphics par plot.new title
#' @importFrom stats sd
#' @importFrom utils packageVersion
#' @details identify marker genes based on aov and AUC.
#' @export
ssc.plot.heatmap <- function(obj, assay.name="exprs",out.prefix=NULL,
                             ncell.downsample=NULL,ave.by=NULL,
                             columns=NULL,columns.order=NULL,gene.desc=NULL,
                             colSet=list(), pdf.width=16,pdf.height=15,
                             do.scale=TRUE,z.lo=-2.5,z.hi=2.5,z.step=1,exp.title="Exp",
                             do.clustering.row=T,do.clustering.col=T,
                             dend.col=FALSE,dend.row=FALSE,
                             clustering.distance="spearman",clustering.method="complete",
							               k.row=1,k.col=1,
							               returnHT=FALSE,
                             palette.name=NULL,
			     palette.ann.numeric="RdYlBu",Y.level.ann.numeric=NULL,
			     row.split=NULL,column.split=NULL,
                             annotation_legend_param=list(),ann.bar.height=1.5, mytitle="",...)
{
    #requireNamespace("ComplexHeatmap")
    #requireNamespace("circlize")
    #requireNamespace("gridBase")
    #requireNamespace("grid")
    #requireNamespace("RColorBrewer")

    if(!is.null(gene.desc) && ("Group" %in% colnames(gene.desc)) && ("geneID" %in% colnames(gene.desc))){
        obj <- obj[gene.desc$geneID,]
    }
    if(!is.null(ncell.downsample) && ncell.downsample < ncol(obj) ){
        obj <- obj[,sample(seq_len(ncol(obj)),ncell.downsample)]
    }

    n <- nrow(obj)
    m <- ncol(obj)
    if(n<2) { loginfo(sprintf("Too few genes: n=%s",n)); return(NULL) }
    if(m<2) { loginfo(sprintf("Too few samples: m=%s",m)); return(NULL) }

    if (!is.null(row.split) && is.null(names(row.split))) {
        names(row.split) <- unname(rowData(obj)$display.name)
    }
    if (!is.null(column.split) && is.null(names(column.split))) {
        names(column.split) <- colnames(obj)
    }

    #### sort
    if(is.null(ave.by)){
        obj <- ssc.assay.hclust(obj,assay.name=assay.name,
                                # order.col=if(is.logical(dend.col) && FALSE==dend.col) do.clustering.col else FALSE,
                                # order.row=if(is.logical(dend.row) && FALSE==dend.row) do.clustering.row else FALSE,
                                order.col = do.clustering.col,
                                order.row = do.clustering.row,
                                clustering.distance=clustering.distance,clustering.method=clustering.method,
                                k.row=1,k.col=1)
    }else{
        avg.colDat <- unique((colData(obj)[,unique(c(ave.by,columns,columns.order)),drop=F]))
        obj <- ssc.average.cell(obj,assay.name=assay.name,column=ave.by,ret.type="sce")
        ####columns <- intersect(ave.by,columns)
        ####columns.order <- intersect(ave.by,columns.order)
        # if(nrow(avg.colDat)==ncol(obj) && all(avg.colDat[[ave.by]] %in% colnames(obj)) ){
        #     rownames(avg.colDat) <- avg.colDat[[ave.by]]
        #     colData(obj) <- avg.colDat[colnames(obj),,drop=F]
        # }
        obj <- ssc.assay.hclust(obj,assay.name,order.col=do.clustering.col,order.row=do.clustering.row,
                                clustering.distance=clustering.distance,clustering.method=clustering.method)
    }

    #### visualization of annotation on top of heatmap
    ha.col <- NULL
    annDF <- data.frame()
    if(!is.null(columns))
    {
        if(!is.null(columns.order)){
            obj <- ssc.order(obj,columns.order=columns.order)
        }
        annDF <- as.data.frame(colData(obj)[columns])
        ###if(length(colSet)==0)
	for(i in seq_along(columns))
	{
	    x <- columns[i]
	    if(!x %in% names(colSet))
	    {
		if(class(colData(obj)[,x])=="numeric"){
		    if(!is.null(Y.level.ann.numeric)){
			Y.level <- Y.level.ann.numeric
		    }else{
			if(all(colData(obj)[,x]<=1) && all(colData(obj)[,x]>=0)){
			    Y.level <- c(0,1)
			}else{
			    Y.level <- pretty(colData(obj)[,x],n=8)
			}
		    }
		    # continious version
		    palette.ann.numeric.values <- getColorPaletteFromNameContinuous(palette.ann.numeric)
		    colSet[[x]] <- colorRamp2(seq(Y.level[1],Y.level[length(Y.level)],
						  length=length(palette.ann.numeric.values)),
					      ##rev(brewer.pal(n = 7, name = "RdYlBu")),
					      palette.ann.numeric.values,
					      space="LAB")
		    if(is.null(annotation_legend_param[[x]])){
			annotation_legend_param[[x]] <- list(color_bar="continuous",
							     legend_direction="horizontal",
							     legend_width=unit(4, "cm"),
							     legend_height=unit(2, "cm"))
		    }else{
			annotation_legend_param[[x]] <- c(annotation_legend_param[[x]],
							  list(color_bar="continuous",
							     legend_direction="horizontal",
							     legend_width=unit(4, "cm"),
							     legend_height=unit(2, "cm")))
		    }
		}else{
		    group.value <- sort(unique(colData(obj)[,x]))
		    colSet[[x]] <- structure(auto.colSet(length(group.value),name="Accent"),
					     names=as.character(group.value))
		}
	    }
	}

        g.show.legend <- T
	#print(str(annDF))
	#print(str(colSet))
        ha.col <- ComplexHeatmap::HeatmapAnnotation(df = annDF, col = colSet,
                                    show_legend = g.show.legend,
				    #simple_anno_size = unit(ann.bar.height, "cm"),
				    annotation_height = unit(rep(ann.bar.height,ncol(annDF)), "cm"),
                                    annotation_legend_param = annotation_legend_param)
        ###top_annotation_height <- unit(ann.bar.height * ncol(annDF), "cm")
    }

    obj <- ssc.order(obj,columns.order=NULL,gene.desc=gene.desc)

    #### scale data for visualization
    dat.plot <- as.matrix(assay(obj,assay.name))
    rownames(dat.plot) <- unname(rowData(obj)$display.name)
    #### scale by row
    if(do.scale)
    {
        rowM <- rowMeans(dat.plot, na.rm = T)
        rowSD <- apply(dat.plot, 1, sd, na.rm = T)
        dat.plot <- sweep(dat.plot, 1, rowM)
        dat.plot <- sweep(dat.plot, 1, rowSD, "/")
        if(!is.null(z.lo)){ dat.plot[dat.plot < z.lo] <- z.lo }
        if(!is.null(z.hi)){ dat.plot[dat.plot > z.hi] <- z.hi }
    }else{
        ###tmp.var <- pretty(abs(dat.plot),n=8)
        tmp.var <- pretty((dat.plot),n=8)
        if(is.null(z.lo)){ z.lo <- tmp.var[1] }
        if(is.null(z.hi)) { z.hi <- tmp.var[length(tmp.var)] }
        if(is.null(z.step)) { z.step <- tmp.var[2]-tmp.var[1] }
    }

    ##### plot
	if(!is.null(out.prefix))
	{
		pdf(sprintf("%s.pdf",out.prefix),width=pdf.width,height=pdf.height)
		par(mar=c(4,12,4,4))
		plot.new()
		title(main = mytitle,cex.main=2)
		##legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5,inset=c(-0.03,0),xpd=T)
		### Integrating Grid Graphics Output with Base Graphics Output
		vps <- gridBase::baseViewports()
		grid::pushViewport(vps$inner, vps$figure, vps$plot)
	}

    if(is.null(palette.name)){
        exp.palette <- rev(brewer.pal(n = 7, name = ifelse(do.scale,"RdBu","RdYlBu")))
    }else{
        ###exp.palette <- rev(brewer.pal(n = 7, name = palette.name))
        exp.palette <- getColorPaletteFromNameContinuous(palette.name)
    }

    if(!is.null(row.split)){
        row.split <- row.split[rownames(dat.plot)]
    }
    if(!is.null(column.split)){
        column.split <- column.split[colnames(dat.plot)]
    }
    ht <- do.call(plotMatrix.simple,c(list(dat=dat.plot,out.prefix=NULL,exp.name=exp.title,show.number=F,
                                           do.clust=NULL,z.lo=z.lo,z.hi=z.hi,palatte=exp.palette,
                                           ###clust.row=FALSE,clust.column=FALSE,
                                           returnHT=TRUE,
                                           par.legend=list(at = seq(z.lo,z.hi,z.step)),
                                           mytitle=mytitle, top_annotation = ha.col),
                                      if(as.character(packageVersion("ComplexHeatmap")) %in% c("1.17.1")) list() else list(column.split=column.split),
                                      list(...)))

	if(!is.null(out.prefix)){
		ComplexHeatmap::draw(ht, newpage= FALSE,merge_legends = TRUE,split=row.split)
#		if(!is.null(ha.col)){
#			for(i in seq_along(names(ha.col@anno_list))){
#			  ComplexHeatmap::decorate_annotation(names(ha.col@anno_list)[i],
#									{grid.text(names(ha.col@anno_list)[i], unit(-4, "mm"),
#											   gp=grid::gpar(fontsize=14),just = "right")})
#			}
#		}
		dev.off()
	}
	if(returnHT){ return(ht) }
}

#' plot gene expression density
#' @param obj object of \code{SingleCellExperiment} class
#' @param out.prefix character; output prefix. required
#' @param gene.id should be in rownames(obj); genes to plot
#' @param gene.symbol should be in rowData(obj)[,"display.name"]; genes to plot
#' @param assay.name character; which assay (default: "exprs")
#' @param pallete.name character; pallete to use. (default: "heat")
#' @param expT double; expression threshold of the genes. (default: c(0.3,0.3))
#' @param ann.txt.dis double; adjust the position of the annotation text (default: 0.3)
#' @param ann.txt.cex double; cex for annotation text (default: 1.2)
#' @param my.title character; title of the figure (default: "")
#' @param par.title list; parameters for drawing title by mtext (default: list(side=3,cex=1.8,line=-3,adj=0.5))
#' @importFrom ks kde
#' @importFrom fields image.plot
#' @importFrom scales viridis_pal
#' @importFrom stats density
#' @importFrom grDevices pdf png dev.off heat.colors rainbow
#' @importFrom graphics par plot.new title lines axis abline text box mtext layout
#' @details make density plot of genes. Note, density estimation from ggplot2 is different (and not so 'effective' as that from ks::kde). One of gene.id and gene.symbol must be provided.
#' @export
ssc.plot.GeneDensity <- function(obj,out.prefix,gene.id,gene.symbol,assay.name="norm_exprs",
								expT=c(0.3,0.3),pallete.name="heat",
								#adjB=NULL,do.scale=F,
								ann.txt.dis=0.3,ann.txt.cex=1.2,
                                my.title="",par.title=list(side=3,cex=1.8,line=-3,adj=0.5))
        
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
		gene.id <- intersect(rownames(obj),gene.id)
		gene.list <- rowData(obj)[,"display.name"][gene.id]
	}
	if(length(gene.list)==0){ 
		warning("No data found for the provided genes!")
		return(NULL)
	}
	obj <- obj[names(gene.list),]
	dat.block <- as.matrix(assay(obj,sprintf("%s",assay.name)))
	dat.plot <- cbind(data.table(cellID=colnames(obj)),(t(dat.block)))
	gene.list <- rowData(obj)[,"display.name"]

	if(length(gene.list)==2){
		colnames(dat.plot)[c(2,3)] <- c("x","y")
		par.old = par(no.readonly = T)

		.density <- ks::kde(dat.plot[,c("x","y")])
		.density.x <- density(dat.plot$x)
		.density.y <- density(dat.plot$y)
		zz <- c(5,10,20,30,40,50,60,70,80,90,95)
		if(pallete.name=="heat"){
			col.den2d <- c("transparent", rev(heat.colors(length(zz))))
		}else if(pallete.name=="rainbow"){
			col.den2d <- c("transparent", rev(rainbow(length(zz), end = 4/6)))
		}else if(pallete.name=="viridis"){
			col.den2d <- c("transparent", (scales::viridis_pal()(length(zz))))
		}

		dat.plot.x.range <- pretty(range(dat.plot$x))
		dat.plot.y.range <- pretty(range(dat.plot$y))
		dat.plot.x.range <- c(dat.plot.x.range[1],dat.plot.x.range[length(dat.plot.x.range)])
		dat.plot.y.range <- c(dat.plot.y.range[1],dat.plot.y.range[length(dat.plot.y.range)])

		mar.left <- 7
		pdf(sprintf("%s.dens2D.pdf",out.prefix),width = 8,height = 7)
		par(mar = c(0.8,mar.left,1,0),cex.lab=2,cex.axis=1.5)
		#layout(matrix(1:6, nrow = 2, byrow = T), widths = c(10,3,2.5), heights = c(3,10))
		nf <- layout(matrix(c(4,4,4,1,0,0,2,3,5),3,3,byrow = TRUE), c(10,3,2.5), c(1.5,3,10))
		#layout.show(nf)
		#dev.off()

		### upper
		plot(NULL, type = "n",  ylab = "",xlab="", xlim = dat.plot.x.range, ylim = c(0, max(.density.x$y)),
			 main = NA, axes = F, xaxs="i")
		lines(.density.x, col = "darkblue", lwd = 2)
		abline(v=expT[1],col="red4",lty=2,lwd=1.5)
		title(ylab = "density", line = 4.5)
		axis(2, las = 1)

		### main
        ### difference between ks_1.11.7 and ks_1.13.1
		par(mar = c(6,mar.left,0,0))
	#	plot(NULL, type = "n",  ylab = "", xlim = dat.plot.x.range, ylim = dat.plot.y.range,
	#		 main = NA, axes = F, xaxs="i",yaxs="i",xlab="")
		plot(.density,display="filled.contour2", cont=zz,xlab="", ylab="",col=col.den2d,
			 xlim=dat.plot.x.range,ylim=dat.plot.y.range)
		title(ylab = gene.list[2], line = 4.5)
		title(xlab = gene.list[1])
		.addAnn <- function()
		{
		    abline(v=expT[1],col="red4",lty=2,lwd=1.5)
		    abline(h=expT[2],col="red4",lty=2,lwd=1.5)
		    ## contingency table
		    cont.tb <- c()
		    ### DP
		    nn <- sum(dat.plot$x >= expT[1] & dat.plot$y >= expT[2])
		    cont.tb <- c(cont.tb,nn)
		    text(par('usr')[2],par('usr')[4]-ann.txt.dis,
			     labels = sprintf("%4.2f %%\n(%d/%d)", nn*100/nrow(dat.plot),nn,nrow(dat.plot)),
			     adj = 1.1,xpd=T,cex=ann.txt.cex)
		    ### SP.1
		    nn <- sum(dat.plot$x >= expT[1] & dat.plot$y < expT[2])
		    cont.tb <- c(cont.tb,nn)
		    text(par('usr')[2],expT[2]-ann.txt.dis,
			     labels = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),
			     adj=1.1,xpd=T,cex=ann.txt.cex)
		    ### SP.2
		    nn <- sum(dat.plot$x < expT[1] & dat.plot$y >= expT[2])
		    cont.tb <- c(cont.tb,nn)
		    text(expT[1],par('usr')[4]-ann.txt.dis,
			     labels = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),
			     adj = 1.1,xpd=T,cex=ann.txt.cex)
		    ### DN
		    nn <- sum(dat.plot$x < expT[1] & dat.plot$y < expT[2])
		    cont.tb <- c(cont.tb,nn)
		    text(expT[1],expT[2]-ann.txt.dis,
			     labels = sprintf("%4.2f %%\n(%d/%d)",nn*100/nrow(dat.plot),nn,nrow(dat.plot)),
			     adj=1.1,xpd=T,cex=ann.txt.cex)
		    return(matrix(cont.tb,ncol=2))
		}
		cont.mtx <- .addAnn()
		res.fisher <- fisher.test(cont.mtx)
		if(is.null(my.title)){
		    my.title <- sprintf("OR:%4.4f, p: %4.4f", res.fisher$estimate,res.fisher$p.value)
		}
		box()

		# right density plot
		par(mar = c(6,0.8,0,1))
		plot(NULL, type = "n", xlab = "density",ylab="",
			 ylim = dat.plot.y.range, xlim = c(0, max(.density.y$y)),
			 main = NA, axes = F, yaxs="i")
		lines(x=.density.y$y,y = .density.y$x, col = "darkblue", lwd = 2)
		abline(h=expT[2],col="red4",lty=2,lwd=1.5)
		axis(1,hadj = 0)

		### title
		par(mar=c(0,mar.left,0,4))
		plot.new()
        do.call(mtext,c(list(text=my.title),par.title))
		
		### legend
		par(mar = c(6,0.25,0,1))
		plot.new()
		fields::image.plot(zlim=c(0,zz[length(zz)]),legend.only=TRUE,add=T,
				   col = col.den2d,
				   axis.args=list( at=zz, labels=sprintf("%s%%",100-zz)),
				   legend.width=12,legend.mar=28)
		par(par.old)
		dev.off()

		###### density estimation from ggplot2 is different (and not so 'effective' as that from ks::kde)
		#	p <- ggplot(dat.plot, aes(x, y)) +
		#		geom_point(alpha=0.5,size=0,color="white") +
		#		#stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE,bins=150) +
		#		#geom_density2d(bins=150) +
		#		#ggalt::stat_bkde2d(bandwidth=c(0.1,0.1),aes(fill = ..nlevel..), geom = "polygon")+
		#        scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9,"YlOrRd"))+
		#		#geom_contour() +
		#		#guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10)) +
		#		theme_bw()
		#	p2 <- ggExtra::ggMarginal(p, type = 'density')
		#	ggsave(sprintf("%s.test.00.pdf",out.prefix),p2,width=6,height=4)
	}
}

