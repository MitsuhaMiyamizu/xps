\name{RLE-methods}
\docType{methods}
\alias{RLE-methods}
\alias{RLE}
\title{Relative Log Expression (RLE)}
\description{
  Produce boxplots of Relative Log Expression (RLE) values for the set of arrays.
  Alternatively, summary statistics or RLE values can be extracted as data.frame.

\emph{Usage}

\code{
RLE(x,
    treename = "*",
    type     = c("plot", "stats", "values"),
    qualopt  = NULL,
         ...) 
}
}
\arguments{
  \item{x}{object of class \code{\link{QualTreeSet}}.}
  \item{treename}{vector of tree names to export.}
  \item{type}{type of output, \code{plot}, \code{stats} or \code{values}.}
  \item{qualopt}{quality control option, i.e. \sQuote{raw}, \sQuote{adjusted}, \sQuote{normalized} or \sQuote{all}.}
  \item{\dots}{optional arguments to be passed to \code{rleplot}.}
}
\details{
 Create boxplots of Relative Log Expression (RLE) values for the set of arrays, i.e. of M plots,
 where M is determined relative to a pseudo-median reference chip.

 Alternatively it is possible to extract either the summary statistics as \code{data.frame} 
 (\code{type="stats"}) or all RLE values as \code{data.frame} (\code{type="values"}).

 If an object of class \code{\link{QualTreeSet}} was created by fitting a probe level model 
 with \code{qualopt="all"} then \code{RLE} will plot or extract RLE for "all" quality options. If 
 you want to plot or extract RLE for a certain quality option only, e.g. "normalized" data only, then 
 you can  use parameter \code{qualopt} with \code{qualopt="<qualopt>"}. 
}
\author{Christian Stratowa}
\seealso{\code{\link{plotRLE}}, \code{\link{rleplot}}}
\examples{
\dontrun{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## qualification - rlm
rlm.all <- rmaPLM(data.test3, "tmp_Test3RLMall", filedir=getwd(), tmpdir="", qualopt="all", option="transcript", add.data=FALSE)

## plot expression levels
if (interactive()) {
RLE(rlm.all)
RLE(rlm.all, qualopt="normalized")
qcRLE <- RLE(rlm.all, type="stats")
qcRLE <- RLE(rlm.all, type="values")
qcRLE <- RLE(rlm.all, treename="TestA1_normalized.rlm", type="stats")
qcRLE <- RLE(rlm.all, treename="TestA1_normalized.rlm", type="values")
}
}
}
\keyword{methods}
