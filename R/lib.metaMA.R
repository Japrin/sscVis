##### code modified from metaMA

#' effect size from moderate t test
#' @param tstat vector; t value of limma's moderate t test
#' @param ntilde doule; .
#' @param m doule; . 
#' @details calculate the effect size
#' @return a matrix 
#' @export
effectsize <- function(tstat,ntilde,m)
{
	# cm=gamma(m/2)/(sqrt(m/2)*gamma((m-1)/2))
	ln.cm<-lgamma(m/2)-log(sqrt(m/2))-lgamma((m-1)/2)  
	cm<-exp(ln.cm)
	d=tstat/sqrt(ntilde)
	dprime=cm*d
	terme1=m/((m-2)*ntilde)
	vard=terme1+d^2*(terme1*ntilde-1/cm^2)
	vardprime=cm^2*(terme1+dprime^2*(terme1*ntilde-1/cm^2))
	result=cbind(d,vard,dprime,vardprime)
	colnames(result)=c("d","vard","dprime","vardprime")
	result
}

#' directional Effect Size combination
#' @param ES data frame or matrix; row for genes, column for studies
#' @param varES data frame or matrix; row for genes, column for studies
#' @param useREM logical; . [default TRUE]
#' @importFrom plyr ldply
#' @importFrom stats pnorm
#' @details calculate the weight averaged effect size from multiple studies
#' @return a data.frame
#' @export
directEScombi <- function(ES,varES,useREM=TRUE)
{
	listres <- vector("list",2)
	f.Q.NA <- function (dadj, varadj)
	{
		#w <- 1/varadj
		#tmp1 <- w * dadj
		#mu <- rowSums(tmp1,na.rm=TRUE)/rowSums(w,na.rm=TRUE)
		#Q <- rowSums(w * (dadj - mu)^2,na.rm=TRUE)
		Q <- sapply(seq_len(nrow(dadj)),function(i){
					   f.use <- varadj[i,]!=0
					   w <- 1/varadj[i,f.use]
					   tmp1 <- w * dadj[i,f.use]
					   mu <- sum(tmp1,na.rm=TRUE)/sum(w,na.rm=TRUE)
					   sum(w * (dadj[i,f.use] - mu)^2,na.rm=TRUE)
					})
	}
	tau2.NA <- function(Q, num.studies, my.weights)
	{
			vwts <- rowSums(my.weights,na.rm=TRUE)
			tmp2 <- rowSums(my.weights^2,na.rm=TRUE)
			tau2 <- pmax(0, (Q - (num.studies - 1))/(vwts - tmp2/vwts))
			return(tau2)
	}
	tau2.NA.oneRow <- function(Q, num.studies, my.weights)
	{
		f.use <- is.finite(my.weights)
		vwts <- sum(my.weights[f.use],na.rm=TRUE)
		tmp2 <- sum(my.weights[f.use]^2,na.rm=TRUE)
		tau2 <- pmax(0, (Q - (num.studies - 1))/(vwts - tmp2/vwts))
		return(tau2)
	}
	Qvals <- f.Q.NA(ES, varES)

#	num.studies <- dim(ES)[2]
#	if (useREM) { varES <- varES + tau2.NA(Qvals, num.studies, my.weights = 1/varES) }
#	wt <- 1/varES
#	MUvals <- rowSums(ES * wt,na.rm=TRUE)/rowSums(wt,na.rm=TRUE)
#	MUsES <- sqrt(1/rowSums(wt,na.rm=TRUE))
#	zSco <- MUvals/MUsES
#	###rpvalESc <- 2*(1-pnorm(abs(zSco)))
#	rpvalESc <- 2*(pnorm(-abs(zSco)))
#	rpadjEsc <- p.adjust(rpvalESc,method="BH")
#	###res=which(p.adjust(rpvalESc,method="BH")<=BHth)
#	result=cbind(MUvals,MUsES,zSco,rpvalESc,rpadjEsc)

	result <- ldply(seq_len(nrow(varES)),function(i){
						ES.i <- ES[i,]
						varES.i <- varES[i,]
                        hasNA.ES.i <- which(is.na(ES.i))
                        hasNA.varES.i <- which(is.na(varES.i))
                        if(length(hasNA.ES.i) > 0 || length(hasNA.varES.i) > 0){
                            if(!all(hasNA.ES.i==hasNA.varES.i)){
							    return(cbind(MUvals=0,MUsES=0,zSco=0,rpvalESc=1))
                            }
                        }
						if(all(varES.i==0,na.rm=T) && all(ES.i==0,na.rm=T)){
							return(cbind(MUvals=0,MUsES=0,zSco=0,rpvalESc=1))
						}
						f.use <- varES.i!=0
						if(useREM && sum(f.use,na.rm=T)>1){
							varES.i <- varES.i + tau2.NA.oneRow(Qvals[i],sum(f.use,na.rm=T),1/varES.i[f.use])
						}
						wt.i <- 1/varES.i[f.use]
						MUvals <- sum(ES.i[f.use] * wt.i,na.rm=TRUE)/sum(wt.i,na.rm=TRUE)
						MUsES <- sqrt(1/sum(wt.i,na.rm=TRUE))
						zSco <- MUvals/MUsES
						rpvalESc <- 2*(pnorm(-abs(zSco)))
						return(cbind(MUvals,MUsES,zSco,rpvalESc))
					})
	colnames(result) <- c("comb.ES","comb.ES.sd","comb.Z","comb.p")
	rownames(result) <- rownames(varES)
	result[,"comb.padj"] <- p.adjust(result[,"comb.p"],method="BH")

	result
}

#' directional Effect Size combination
#' @param dat.long data.table; input data
#' @param idx.geneID character; column for geneID (default: "geneID")
#' @param idx.aid character; column for aid (default: "aid")
#' @param idx.ES character; column for ES (default: "dprime")
#' @param idx.varES character; column for varES (default: "vardprime")
#' @param ... parameters passed to directEScombi()
#' @importFrom data.table dcast
#' @details wrapper function to run directEScombi()
#' @return a data.table
#' @export
directEScombiFromLongTable <- function(dat.long,idx.geneID="geneID",idx.aid="aid",
									   idx.ES="dprime",idx.varES="vardprime",...)
{
	### ES
	dat.x.dprime.tb <- dcast(dat.long, sprintf("%s ~ %s",idx.geneID,idx.aid),value.var=idx.ES)
	dat.x.dprime.mtx <- as.matrix(dat.x.dprime.tb[,-c("geneID"),with=F])
	rownames(dat.x.dprime.mtx) <- dat.x.dprime.tb[["geneID"]]
	### varES
	dat.x.vardprime.tb <- dcast(dat.long, sprintf("%s ~ %s",idx.geneID,idx.aid),value.var=idx.varES)
	dat.x.vardprime.mtx <- as.matrix(dat.x.vardprime.tb[,-c("geneID"),with=F])
	rownames(dat.x.vardprime.mtx) <- dat.x.vardprime.tb[["geneID"]]

	if(!all(rownames(dat.x.dprime.mtx)==rownames(dat.x.vardprime.mtx))){
		stop(sprintf("strange thing: !all(rownames(dat.x.dprime.mtx)==rownames(dat.x.vardprime.mtx)) (%s)\n",x))
	}
	dat.x.es.combi <- directEScombi(dat.x.dprime.mtx, dat.x.vardprime.mtx,...)
	dat.x.es.combi.tb <- cbind(data.table(geneID=rownames(dat.x.es.combi)),
							   dat.x.es.combi)
	return(dat.x.es.combi.tb)
}


#' collapse effect size, give input table and group variate
#' @param dat.long data.table; these columns are required: "geneID","aid","P.Value","adj.P.Val","sig","dprime","vardprime","group.2nd"
#' @param group.2nd character; column of obj, used for collapse [default: NULL]
#' @param mode.collapse character; collapse method. one of "min", "comb" [default: "comb"]
#' @param th.adj.P double; threshold [default: 0.01]
#' @param th.dprime double; threshold [default: 0.15]
#' @importFrom data.table as.data.table data.table setorderv
#' @importFrom plyr ldply
#' @details collapse effect size
#' @return a gene table
#' @export
collapseEffectSizeLong <- function(dat.long,mode.collapse="comb",group.2nd="group.2nd",
                                   th.adj.P=0.01,th.dprime=0.15)
{
  #dat.long[["group.2nd"]] <- dat.long[[group.2nd]]
  dat.collapse <- NULL
  if(mode.collapse=="min"){
    dat.collapse <- dat.long[,.(dprime=dprime[which.min(abs(dprime))],
                                vardprime=vardprime[which.min(abs(dprime))],
                                P.Value=P.Value[which.min(abs(dprime))],
                                adj.P.Val=adj.P.Val[which.min(abs(dprime))],
                                sig=all(adj.P.Val<th.adj.P) & all(dprime>th.dprime) ),
                             by=c("geneID",group.2nd)]
  }else if(mode.collapse=="comb"){
    aid.mapping.tb <- unique(dat.long[,c("aid",group.2nd),with=F])
    n.group.2nd <- aid.mapping.tb[,.N,by=group.2nd]
    ########
    ##
    colidx.freq <- c()
    if(all(c("freq._case","cluster.size") %in% colnames(dat.long))){
        colidx.freq <- c("freq._case","cluster.size")
    }
    dat.collapse.nOne <- dat.long[dat.long[[group.2nd]] %in% n.group.2nd[N==1,][[group.2nd]],
                                  c("geneID",group.2nd,"dprime","vardprime","P.Value","adj.P.Val","sig",colidx.freq),with=F]
    dat.collapse.nMul <- as.data.table(ldply(n.group.2nd[N>1,][[group.2nd]], function(x) {
				   dat.x <- dat.long[dat.long[[group.2nd]]==x,]
				   dat.x.sig.min <- dat.x[,.(gene.sig.min=(all(adj.P.Val<th.adj.P) &
								       all(dprime>th.dprime))),
							    by="geneID"]
				   dat.x.es.combi.tb <- directEScombiFromLongTable(dat.x)
				   dat.x.es.combi.tb <- merge(dat.x.es.combi.tb,dat.x.sig.min,by="geneID")
				   ret.tb <-data.table(geneID=dat.x.es.combi.tb$geneID,
						       group.2nd=x,
						       dprime=dat.x.es.combi.tb[["comb.ES"]],
						       vardprime=dat.x.es.combi.tb[["comb.ES.sd"]]^2,
						       P.Value=dat.x.es.combi.tb[["comb.p"]],
						       adj.P.Val=dat.x.es.combi.tb[["comb.padj"]],
						       #sig.all=dat.x.es.combi.tb[["gene.sig.min"]],
						       sig=(dat.x.es.combi.tb[["gene.sig.min"]]) |
							   (dat.x.es.combi.tb[["comb.ES"]]>th.dprime &
							      dat.x.es.combi.tb[["comb.padj"]]<th.adj.P)
						       )
				   colnames(ret.tb)[2] <- group.2nd
				   if(all(c("freq._case","cluster.size") %in% colnames(dat.long)))
				   {
				       freq.collapsed.tb <- dat.x[,.(freq._case=(sum(.SD$freq._case*.SD$cluster.size))/sum(.SD$cluster.size),
									  cluster.size=sum(.SD$cluster.size)),
								    by="geneID"]
				       ret.tb <- merge(ret.tb,freq.collapsed.tb,by="geneID")
				   }
				   return(ret.tb)
				}))
    dat.collapse <- rbind(dat.collapse.nOne,dat.collapse.nMul)
    setorderv(dat.collapse,c("geneID",group.2nd))
  }
  return(dat.collapse)
}
