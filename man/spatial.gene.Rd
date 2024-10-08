% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{spatial.gene}
\alias{spatial.gene}
\title{Spatial expression profile}
\usage{
spatial.gene(
  st,
  coords,
  gene.list,
  size = 1,
  img.backgound = F,
  img,
  pdfout = T
)
}
\arguments{
\item{st}{Spatial transcriptomics.}

\item{coords}{Coordinates of spatial spots.}

\item{gene.list}{List of genes.}

\item{size}{Size of spot. Default is 1.}

\item{img}{Image to be used as background.}

\item{pdfout}{Whether to draw the plot in a pdf file named "spatial_gene". When more than one gene is plotted, this is recommended to be \code{TRUE}.}

\item{img.background}{Whether to use an image as background. Default is \code{FALSE}.}
}
\value{
A pdf file named "spatial_genes".
}
\description{
Visualization of spatial expression of some genes.
}
