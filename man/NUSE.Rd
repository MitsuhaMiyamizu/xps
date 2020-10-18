\name{NUSE-methods}
\docType{methods}
\alias{NUSE-methods}
\alias{NUSE}
\title{Normalized Unscaled Standard Errors (NUSE)}
\description{
  Produce boxplot of Normalized Unscaled Standard Errors (NUSE) for the set of arrays. 
  Alternatively, summary statistics or NUSE values can be extracted as data.frame.

\emph{Usage}

\code{
NUSE(x,
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
  \item{\dots}{optional arguments to be passed to \code{nuseplot}.}
}
\details{
 Create boxplots of Normalized Unscaled Standard Errors (NUSE) for the set of arrays.

 Alternatively it is possible to extract either the summary statistics as \code{data.frame} 
 (\code{type="stats"}) or all NUSE values as \code{data.frame} (\code{type="values"}).

 If an object of class \code{\link{QualTreeSet}} was created by fitting a probe level model 
 with \code{qualopt="all"} then \code{NUSE} will plot or extract NUSE for "all" quality options. If 
 you want to plot or extract NUSE for a certain quality option only, e.g. "normalized" data only, then 
 you can  use parameter \code{qualopt} with \code{qualopt="<qualopt>"}. 
}
\author{Christian Stratowa}
\seealso{\code{\link{plotNUSE}}, \code{\link{nuseplot}}}
\examples{
\dontrun{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## qualification - rlm
rlm.all <- rmaPLM(data.test3, "tmp_Test3RLMall", filedir=getwd(), tmpdir="", qualopt="all", option="transcript", add.data=FALSE)

## plot expression levels
if (interactive()) {
NUSE(rlm.all)
NUSE(rlm.all, qualopt="normalized")
qcNUSE <- NUSE(rlm.all, type="stats")
qcNUSE <- NUSE(rlm.all, type="values")
qcNUSE <- NUSE(rlm.all, treename="TestA1_normalized.rlm", type="stats")
qcNUSE <- NUSE(rlm.all, treename="TestA1_normalized.rlm", type="values")
}
}
}
\keyword{methods}
