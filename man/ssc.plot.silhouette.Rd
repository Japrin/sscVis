% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sscVis.R
\name{ssc.plot.silhouette}
\alias{ssc.plot.silhouette}
\title{Wraper for silhouette()}
\usage{
ssc.plot.silhouette(
  obj,
  cluster.label,
  reducedDim.name = "iCor.tsne",
  do.plot = T,
  ...
)
}
\arguments{
\item{obj}{object of \code{SingleCellExperiment}}

\item{cluster.label}{character; which column of colData of obj to used as cluster label.}

\item{reducedDim.name}{character; which reducedDim to use. (default: "iCor.tsne")}

\item{do.plot}{logical; whether plot}

\item{...}{Arguments to be passed to plot()}
}
\value{
an object, sil, of class \code{silhouette}
}
\description{
Wraper for silhouette()
}
