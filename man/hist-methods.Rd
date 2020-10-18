\name{hist-methods}
\docType{methods}
\alias{hist-methods}
\alias{hist}
\title{Plot Density Estimate}
\description{
 Plot the density estimates for each sample.

\emph{Usage}

\code{
hist(x,
     which      = "",
     size       = 0,
     transfo    = log2,
     xlab       = "log intensity",
     ylab       = "density",
     names      = "namepart",
     type       = "l",
     col        = 1:6,
     lty        = 1:5,
     add.legend = FALSE,
     verbose    = TRUE,
     ...)
}
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}} or \code{\link{ExprTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{size}{length of sequence to be generated as subset.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{xlab}{a title for the x axis.}
  \item{ylab}{a title for the y axis.}
  \item{names}{optional vector of sample names.}
  \item{type}{type for the plot.}
  \item{col}{colors to use for the different arrays.}
  \item{lty}{line types to use for the different arrays.}
  \item{add.legend}{logical, if TRUE then a legend will be drawn.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{\dots}{optional arguments to be passed to \code{plot}.}
}
\details{
 Plots the non-parametric density estimates for each sample.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as callplot.
}
\author{Christian Stratowa}
\note{
 For objects of class \code{DataTreeSet} it is no longer necessary to \code{attachInten} since
 each data tree will be imported separately.
}
\seealso{\code{\link{plotDensity}}}
\examples{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

if (interactive()) {
hist(data.test3)
}
}
\keyword{methods}
