\name{summarize}
\alias{summarize}
\alias{summarize.mas4}
\alias{summarize.mas5}
\alias{summarize.rma}
\alias{xpsSummarize-methods}
\alias{xpsSummarize}
\title{Probe Set Summarizing Functions}
\description{
 Converts Affymetrix probe level data to expression levels by summarizing the probe set
 values into one expression measure and a standard error for this summary.
}
\usage{
summarize(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, select = "none", method = character(), option = "transcript", logbase = "0", exonlevel = "", params = list(), xps.scheme = NULL, add.data = TRUE, verbose = TRUE)

summarize.mas4(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, option = "transcript", exonlevel = "", xps.scheme = NULL, add.data = TRUE, verbose = TRUE)

summarize.mas5(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, option = "transcript", exonlevel = "", xps.scheme = NULL, add.data = TRUE, verbose = TRUE)

summarize.rma(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, option = "transcript", exonlevel = "", xps.scheme = NULL, add.data = TRUE, verbose = TRUE)

xpsSummarize(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{DataTreeSet}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{update}{logical. If \code{TRUE} the existing ROOT data file \code{filename} will be updated.}
  \item{select}{type of probes to select for summarization.}
  \item{method}{summarization method to use.}
  \item{option}{option determining the grouping of probes for summarization, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon arrays only.}
  \item{logbase}{logarithm base as character, one of \sQuote{0}, \sQuote{log},
   \sQuote{log2}, \sQuote{log10}.}
  \item{exonlevel}{exon annotation level determining which probes should be used
   for summarization; exon/genome arrays only.}
  \item{params}{vector of parameters for summarization method.}
  \item{xps.scheme}{optional alternative \code{SchemeTreeSet}.}
  \item{add.data}{logical. If \code{TRUE} expression data will be included as slot \code{data}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{DataTreeSet}.}
  \item{\dots}{the arguments described above.}
}
\details{
 Converts Affymetrix probe level data to expression levels by summarizing the probe set
 values into one expression measure and a standard error for this summary.

 \code{xpsSummarize} is the \code{DataTreeSet} method called by function \code{summarize},
 containing the same parameters.
}
\value{
 An \code{\link{ExprTreeSet}}.
}
\author{Christian Stratowa}
\seealso{\code{\link{express}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## RMA background
data.bg.rma <- bgcorrect.rma(data.test3,"tmp_Test3RMA",filedir=getwd(),tmpdir="",verbose=FALSE)
## normalize quantiles
data.qu.rma <- normalize.quantiles(data.bg.rma,"tmp_Test3RMA",filedir=getwd(),tmpdir="",update=TRUE,verbose=FALSE)
## summarize medianpolish
data.mp.rma <- summarize.rma(data.qu.rma,"tmp_Test3RMA",filedir=getwd(),tmpdir="",update=TRUE,verbose=FALSE)

## get expression data.frame
expr.rma <- exprs(data.mp.rma)
head(expr.rma)

## plot expression levels
if (interactive()) {
boxplot(data.mp.rma)
boxplot(log2(expr.rma[,3:6]))
}
}
\keyword{manip}
