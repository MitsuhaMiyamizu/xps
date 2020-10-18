\name{mboxplot-methods}
\docType{methods}
\alias{mboxplot-methods}
\alias{mboxplot}
\title{Box Plots of Relative M Values}
\description{
  Produce boxplots of relative M values for the set of arrays.

\emph{Usage}

\code{
mboxplot(x, which = "", size = 0, transfo = log2, method  = "mean", range = 0, ylim = c(-1,1), outline = FALSE, names = "namepart", ...)
}
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}} or \code{\link{ExprTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{size}{length of sequence to be generated as subset.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{method}{method to create the reference data, \dQuote{mean} or \dQuote{median}.}
  \item{range}{determines how far the plot whiskers extend out from the box.}
  \item{ylim}{range for the plotted y values.}
  \item{outline}{if \code{outline} is not true, the outliers are not drawn.}
  \item{names}{optional vector of sample names.}
  \item{\dots}{optional arguments to be passed to \code{boxplot}.}
}
\details{
 Create boxplots of M plots, where M is determined relative to a pseudo-mean reference chip.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as boxplot.
}
\author{Christian Stratowa}
\note{
 For a \code{\link{DataTreeSet}} object, data must first be attached using method 
 \code{\link{attachInten}}.
}
\seealso{\code{\link{mvaplot}}, \code{\link{boxplot}}}
\examples{
# load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

# need to attach scheme mask and probe intensities
data.test3 <- attachMask(data.test3)
data.test3 <- attachInten(data.test3)

if (interactive()) {
mboxplot(data.test3, ylim=c(-6,6))
}

# optionally remove mask and data to free memory
data.test3 <- removeInten(data.test3)
data.test3 <- removeMask(data.test3)
}
\keyword{methods}
