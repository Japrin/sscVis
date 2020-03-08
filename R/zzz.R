.onLoad <- function(libname , pkgname) {
	if(getRversion() >= "2.15.1"){
		utils::globalVariables(c(".",".N",".SD","variable","value","Group",
								 "PC","isKneePts",
								 "meanExp",
								 "geneID","aid",
								 "Dim1","Dim2",
								 "freq","freq.med","N","NTotal","NOther",
								 "p.value","p.adj","p.signif","y_pos",
								 "silhouette.rank","sil_width",
								 "distCenter.rank","distCenter"))
	}
}
