\name{normalize}
\alias{normalize.constant}
\alias{normalize.lowess}
\alias{normalize.quantiles}
\alias{normalize.supsmu}
\alias{normalize}
\alias{xpsNormalize-methods}
\alias{xpsNormalize}
\title{Normalization on Affymetrix Probe Level Data or on Expression Levels}
\description{
 Functions that allow to normalize Affymetrix arrays both at the probe level
 (\dQuote{low-level normalization}) and/or at the expression level
 (\dQuote{high-level normalization}).
}
\usage{
normalize(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, select = "all", method = "mean", option = "transcript:all", logbase = "0", exonlevel = "", refindex = 0, refmethod = "mean", params = list(0.02, 0), add.data = TRUE, verbose = TRUE)

normalize.constant(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, method = "mean", logbase = "0", exonlevel = "", refindex = 0, refmethod = "mean", params = list(0.02, 0), add.data = TRUE, verbose = TRUE)

normalize.lowess(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, logbase = "log2", exonlevel = "", refindex = 0, refmethod = "mean", params = list(0.67, 3, 0.0, 0.0), add.data = TRUE, verbose = TRUE)

normalize.quantiles(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, exonlevel = "", add.data = TRUE, verbose = TRUE)

normalize.supsmu(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, logbase = "log2", exonlevel = "", refindex = 0, refmethod = "mean", params = list(0.0, 0.0, 0.0, 0.0), add.data = TRUE, verbose = TRUE)

xpsNormalize(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{DataTreeSet} or \code{\link{ExprTreeSet}}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{update}{logical. If \code{TRUE} the existing ROOT data file \code{filename} will be updated.}
  \item{select}{type of probes to select for normalization.}
  \item{method}{normalization method to use.}
  \item{option}{option determining the grouping of probes for normalization, and the selection
   of the probes.}
  \item{logbase}{logarithm base as character, one of \sQuote{0}, \sQuote{log},
   \sQuote{log2}, \sQuote{log10}.}
  \item{exonlevel}{exon annotation level determining which probes should be used
   for summarization; exon/genome arrays only.}
  \item{refindex}{index of reference tree to use, or 0.}
  \item{refmethod}{for \code{refindex=0}, either trimmed mean or median of trees.}
  \item{params}{vector of parameters for normalization method.}
  \item{add.data}{logical. If \code{TRUE} expression data will be included as slot \code{data}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{DataTreeSet} or \code{\link{ExprTreeSet}}.}
  \item{\dots}{the arguments described above.}
}
\details{
 Functions that allow to normalize Affymetrix arrays both at the probe level
 (\dQuote{low-level normalization}) and/or at the expression level
 (\dQuote{high-level normalization}).

 Please have a look at vignette \dQuote{xpsPreprocess.pdf} for details on how to use
 function \code{normalize}. 

 \code{xpsNormalize} are the \code{DataTreeSet} or \code{\link{ExprTreeSet}} methods,
 respectively, called by function \code{normalize}, containing the same parameters.
}
\value{
 An object of type \code{\link{DataTreeSet}} or \code{\link{ExprTreeSet}}.
}
\author{Christian Stratowa}
\section{Warning }{
Functions \code{normalize.lowess} and \code{normalize.supsmu} have only be tested
for \code{object}s of type \code{\link{ExprTreeSet}} but not for objects of type 
\code{\link{DataTreeSet}}, i.e. for probe level intensities.
}
\seealso{\code{\link{express}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## RMA background
data.bg.rma <- bgcorrect.rma(data.test3,"tmp_Test3NormRMA",filedir=getwd(),tmpdir="",verbose=FALSE)
## normalize quantiles
data.qu.rma <- normalize.quantiles(data.bg.rma,"tmp_Test3NormRMA",filedir=getwd(),tmpdir="",update=TRUE,verbose=FALSE)
## summarize medianpolish
data.mp.rma <- summarize.rma(data.qu.rma,"tmp_Test3NormRMA",filedir=getwd(),tmpdir="",update=TRUE,verbose=FALSE)
}
\keyword{manip}
