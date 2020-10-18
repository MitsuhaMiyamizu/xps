\name{pmplot-methods}
\docType{methods}
\alias{pmplot-methods}
\alias{pmplot}
\title{Barplot of PM and MM Intensities.}
\description{
Creates a barplot of mean perfect match and mismatch intensities.

\emph{Usage}

\code{
pmplot(x,
       which   = "",
       size    = 0,
       transfo = NULL,
       method  = "mean",
       names   = "namepart",
       beside  = TRUE,
       col     = c("red", "blue"),
       legend  = c("PM","MM"),
       las     = 2,
       ylab    = "mean intensities",
       ...)
}
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{size}{length of sequence to be generated as subset.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{method}{method to compute average intensities, \dQuote{mean} or \dQuote{median}.}
  \item{names}{optional vector of sample names.}
  \item{beside}{logical. If \code{FALSE}, mean intensities are portrayed as stacked bars, 
   and if \code{TRUE} the columns are portrayed as juxtaposed bars.}
  \item{col}{color of PM, MM bars.}
  \item{legend}{a vector of text used to construct a legend for the plot, 
   or a logical indicating whether a legend should be included.}
  \item{las}{the style of axis labels.}
  \item{ylab}{a title for the y axis.}
  \item{\dots}{optional arguments to be passed to \code{barplot}.}
}
\details{
 Produces barplots of mean perfect match and mismatch intensities for slot \code{data}
 for an object of class \code{\link{DataTreeSet}}.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as pmplot.
}
\author{Christian Stratowa}
\note{
 Data must first be attached to class \code{\link{DataTreeSet}} using method \code{\link{attachInten}}.
}
\seealso{\code{\link{plotPM}}, \code{\link{boxplot}}, \code{\link{barplot}}}
\examples{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## need to attach scheme mask and probe intensities
data.test3 <- attachMask(data.test3)
data.test3 <- attachInten(data.test3)

if (interactive()) {
pmplot(data.test3)
}

## optionally remove mask and data to free memory
data.test3 <- removeInten(data.test3)
data.test3 <- removeMask(data.test3)
}
\keyword{methods}
