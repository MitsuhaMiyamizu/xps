\name{express}
\alias{express}
\alias{xpsPreprocess-methods}
\alias{xpsPreprocess}
\title{Compute expression levels from raw data}
\description{
 This function allows to combine different algorithms to compute expression levels,
 or to return the result for different algorithms only.
}
\usage{
express(xps.data,
        filename = character(),
        filedir = getwd(),
        tmpdir = "",
        update = FALSE,
      # background correction
        bgcorrect.method = NULL,
        bgcorrect.select = character(),
        bgcorrect.option = character(),
        bgcorrect.params = list(),
      # normalization
        normalize.method = NULL,
        normalize.select = character(),
        normalize.option = character(),
        normalize.logbase = character(),
        normalize.params = list(),
      # expression values
        summarize.method = NULL,
        summarize.select = character(),
        summarize.option = character(),
        summarize.logbase = character(),
        summarize.params = list(),
      # reference values
        reference.index  = 0,
        reference.method = "mean",
        reference.params = list(0),
      # misc.
        exonlevel  = "",
        xps.scheme = NULL,
        add.data   = TRUE,
        bufsize    = 32000,
        verbose    = TRUE)

xpsPreprocess(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{DataTreeSet}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{update}{logical. If \code{TRUE} the existing ROOT data file \code{filename} will be updated.}
  \item{bgcorrect.method}{background method to use.}
  \item{bgcorrect.select}{type of probes to select for background correction.}
  \item{bgcorrect.option}{type of background correction to use.}
  \item{bgcorrect.params}{vector of parameters for background method.}
  \item{normalize.method}{normalization method to use.}
  \item{normalize.select}{type of probes to select for normalization.}
  \item{normalize.option}{normalization option.}
  \item{normalize.logbase}{logarithm base as character, one of \sQuote{0}, \sQuote{log},
   \sQuote{log2}, \sQuote{log10}.}
  \item{normalize.params}{vector of parameters for normalization method.}
  \item{summarize.method}{summarization method to use.}
  \item{summarize.select}{type of probes to select for summarization.}
  \item{summarize.option}{option determining the grouping of probes for summarization, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon arrays only.}
  \item{summarize.logbase}{logarithm base as character, one of \sQuote{0}, \sQuote{log},
   \sQuote{log2}, \sQuote{log10}.}
  \item{summarize.params}{vector of parameters for summarization method.}
  \item{reference.index}{index of reference tree to use, or 0.}
  \item{reference.method}{for \code{refindex=0}, either trimmed mean or median of trees.}
  \item{reference.params}{vector of parameters for reference method.}
  \item{exonlevel}{exon annotation level determining which probes should be used for summarization; exon/genome arrays only.}
  \item{xps.scheme}{optional alternative \code{SchemeSet}.}
  \item{add.data}{logical. If \code{TRUE} expression data will be included as slot \code{data}.}
  \item{bufsize}{integer which sets the buffer size of the tree branch baskets (default is 32000).}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{DataTreeSet}.}
  \item{\dots}{the arguments described above.}
}
\details{
 This function allows to combine different algorithms to compute expression levels,
 or to return the result for different algorithms only.

 Please have a look at vignette \dQuote{xpsPreprocess.pdf} for details on how to use
 function \code{express}. 

 \code{xpsPreprocess} is the \code{DataTreeSet} method called by function \code{express},
 containing the same parameters.
}
\value{
 An object of type \code{\link{DataTreeSet}} or \code{\link{ExprTreeSet}}.
}
\author{Christian Stratowa}
\seealso{\code{\link{bgcorrect}}, \code{\link{normalize}}, \code{\link{summarize}}}
\examples{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## compute rma with a single call to express()
expr.rma <- express(data.test3,"tmp_Test3Exprs",filedir=getwd(),tmpdir="",update=FALSE,
            bgcorrect.method="rma",bgcorrect.select="none",bgcorrect.option="pmonly:epanechnikov",bgcorrect.params=c(16384),
            normalize.method="quantile",normalize.select="pmonly",normalize.option="transcript:together:none",normalize.logbase="0",normalize.params=c(0.0),
            summarize.method="medianpolish",summarize.select="pmonly",summarize.option="transcript",summarize.logbase="log2",summarize.params=c(10, 0.01, 1.0),
            verbose=FALSE)

## get expression data.frame
expr <- exprs(expr.rma)
head(expr)

## plot expression levels
if (interactive()) {
boxplot(expr.rma)
boxplot(log2(expr[,3:6]))
}

\dontrun{
## examples using Affymetrix human tissue dataset (see also xps/examples/script4exon.R)

## example - exon array, e.g. HuEx-1_0-st-v2:
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
scheme.exon <- root.scheme(paste(scmdir,"Scheme_HuEx10stv2r2_na25.root",sep="/"))
data.exon   <- root.data(scheme.exon, paste(datdir,"HuTissuesExon_cel.root",sep="/"))

workdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/exon"
expr.rma <- express(data.exon,"HuExonExprs",filedir=workdir,tmpdir="",update=F,
            bgcorrect.method="rma",bgcorrect.select="antigenomic",bgcorrect.option="pmonly:epanechnikov",bgcorrect.params=c(16384),
            normalize.method="quantile",normalize.select="pmonly",normalize.option="transcript:together:none",normalize.logbase="0",normalize.params=c(0.0),
            summarize.method="medianpolish",summarize.select="pmonly",summarize.option="transcript",summarize.logbase="log2",summarize.params=c(10, 0.01, 1.0),
            exonlevel="metacore+affx")
}
}
\keyword{manip}
