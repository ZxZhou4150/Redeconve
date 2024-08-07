% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fs.R
\name{to.absolute.abundance}
\alias{to.absolute.abundance}
\title{Giving results interpretability}
\usage{
to.absolute.abundance(res, aver.cell)
}
\arguments{
\item{res}{Result of \code{deconvoluting}.}

\item{aver.cell}{A priori value indicating the average cell number of one spot. This varies according to the platform.}
}
\value{
A cell-by-spot matrix, each place is the estimated absolute abundance of the cell in the spot.
}
\description{
One way of gaining interpretability, estimating absolute abundance
}
\seealso{
\code{\link[=to.proportion]{to.proportion()}}
}
