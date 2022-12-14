% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs.R
\name{profile_correction}
\alias{profile_correction}
\title{Correction of single-cell expression profile.}
\usage{
profile_correction(num, st1, ref, gene.list = NULL, dopar = T, ncores)
}
\arguments{
\item{num}{One column of the results.}

\item{st1}{One column of Spatial transcriptomics data. The spot in which the expression profile is to correct.}

\item{ref}{Reference.}

\item{gene.list}{Genes to be estimated. Default is all.}

\item{dopar}{Whether to use parallel computing.}
}
\description{
Use results and spatial transcriptomics to reversely estimate the spatial-specific expression profile.
}
