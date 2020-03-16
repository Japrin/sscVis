.onLoad <- function(libname , pkgname) {
	if(getRversion() >= "2.15.1"){
		utils::globalVariables(c(".",".N",".SD","variable","value","Group","x",
								 "PC","isKneePts",
								 "meanExp",
								 "geneID","aid",
								 "Dim1","Dim2",
								 "freq","freq.med","N","NTotal","NOther",
								 "p.value","p.adj","p.signif","y_pos",
								 "dprime","vardprime","P.Value","adj.P.Val",
								 "silhouette.rank","sil_width",
								 "distCenter.rank","distCenter"))
	}
}
