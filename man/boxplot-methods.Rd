\name{boxplot-methods}
\docType{methods}
\alias{boxplot-methods}
\alias{boxplot}
\title{Box Plots}
\description{
  Produce box-and-whisker plot(s) of the samples.

\emph{Usage}

\code{
boxplot(x, which = "", size = 0, transfo = log2, range = 0, names = "namepart", bmar = NULL, ...)
}
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}}, \code{\link{ExprTreeSet}} or \code{\link{QualTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{size}{length of sequence to be generated as subset.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{range}{determines how far the plot whiskers extend out from the box.}
  \item{names}{optional vector of sample names.}
  \item{bmar}{optional \code{list} for bottom margin and axis label magnification \code{cex.axis}.}
  \item{\dots}{optional arguments to be passed to \code{boxplot}.}
}
\details{
 Creates a boxplot for slot \code{data} for an object of class \code{\link{DataTreeSet}},
 \code{\link{ExprTreeSet}} or \code{\link{QualTreeSet}}.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as boxplot.

 For \code{bmar=NULL} the default list \code{bmar = list(b=6, cex.axis=1.0)} will be used initially.
 However, both bottom margin and axis label magnification will be adjusted depending on the number
 of label characters and the number of smaples.  
}
\author{Christian Stratowa}
\note{
 For a \code{\link{DataTreeSet}} object, data must first be attached using method \code{\link{attachInten}}.

 Alternatively it is possible to use the pre-calculated quantiles stored in the \code{userinfo}
 of the data trees by calling \code{which="userinfo:varlist"}, where the \code{varlist} to call
 is described in method \code{\link{treeInfo}}.
}
\seealso{\code{\link{plotBoxplot}}, \code{\link[graphics]{boxplot}}}
\examples{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## need to attach scheme mask and probe intensities only if "userinfo" is not used
data.test3 <- attachMask(data.test3)
data.test3 <- attachInten(data.test3)

if (interactive()) {
boxplot(data.test3)
}

## optionally remove mask and data to free memory
data.test3 <- removeInten(data.test3)
data.test3 <- removeMask(data.test3)

## alternatively use the quantiles stored in userinfo of trees
if (interactive()) {
boxplot(data.test3, which="userinfo:fIntenQuant")
}

rm(scheme.test3, data.test3)
gc()
}
\keyword{methods}
