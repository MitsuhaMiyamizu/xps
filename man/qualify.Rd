\name{qualify}
\alias{qualify}
\alias{qualify.rlm}
\alias{xpsQualify-methods}
\alias{xpsQualify}
\title{Probe Set Quality Control Functions}
\description{
 Converts Affymetrix probe level data to expression levels by fitting a multichip model.
}
\usage{
qualify(xps.data,
        filename   = character(0),
        filedir    = getwd(),
        tmpdir     = "",
        update     = FALSE,
        select     = "none",
        method     = character(),
        option     = "transcript",
        logbase    = "log2",
        exonlevel  = "",
        params     = list(),
        xps.scheme = NULL,
        add.data   = TRUE,
        verbose    = TRUE)

qualify.rlm(xps.data,
            filename   = character(0),
            filedir    = getwd(),
            tmpdir     = "",
            update     = FALSE,
            option     = "transcript",
            exonlevel  = "",
            xps.scheme = NULL,
            add.data   = TRUE,
            verbose    = TRUE)

xpsQualify(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{DataTreeSet}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{update}{logical. If \code{TRUE} the existing ROOT data file \code{filename} will be updated.}
  \item{select}{type of probes to select for summarization.}
  \item{method}{qualification method to use, currently \code{rlm}.}
  \item{option}{option determining the grouping of probes for summarization, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon/genome arrays only.}
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
 Converts Affymetrix probe level data to expression levels by fitting a multichip model.

 This function stores three types of ROOT trees in \code{filename}: \cr
 - quality trees containing expression levels, normalized unscaled standard errors (NUSE),
   relative log expressions (RLE) \cr
 - residual trees containing the residual SE and the model fit weights \cr
 - border trees containing the border intensities, mean border intensities and center of intensities (COI)

 \code{xpsQualify} is the \code{DataTreeSet} method called by function \code{qualify},
 containing the same parameters.
}
\value{
 An \code{\link{QualTreeSet}}.
}
\author{Christian Stratowa}
\note{
 This function takes any \code{DataTreeSet} and computes expression levels by summarizing the probe set
 values into one expression measure. It does NOT do any further preprocessing such as background 
 correction or (quantile) normalization. If you want to do background correction and/or normalization
 first then you need to use function \code{\link{fitQC}}.
}
\seealso{\code{\link{fitQC}}}
\examples{
\dontrun{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## compute RMA stepwise

## background correction
data.bg.rma <- bgcorrect.rma(data.test3,"tmp_Test3RMABgrd",filedir=getwd())

## normalize quantiles
data.qu.rma <- normalize.quantiles(data.bg.rma,"tmp_Test3RMANorm",filedir=getwd())

## summarize medianpolish
data.mp.rma <- summarize.rma(data.qu.rma,"tmp_Test3RMAExpr",filedir=getwd(),tmpdir="")

## qualification - rlm

## fit model on raw data
data.raw.rlm <- qualify.rlm(data.test3, "tmp_Test3RawQual", filedir=getwd(), tmpdir="", option="transcript", add.data=TRUE)

## fit model on background adjusted data
data.adj.rlm <- qualify.rlm(data.bg.rma, "tmp_Test3AdjQual", filedir=getwd(), tmpdir="", option="transcript", add.data=TRUE)

## fit model on normalized data
data.nrm.rlm <- qualify.rlm(data.qu.rma, "tmp_Test3NormQual", filedir=getwd(), tmpdir="", option="transcript", add.data=TTRUE)

## get expression levels
expr.raw.rlm <- validData(data.raw.rlm)
expr.adj.rlm <- validData(data.adj.rlm)
expr.nrm.rlm <- validData(data.nrm.rlm)

## get borders
brd.raw <- borders(data.raw.rlm)
brd.adj <- borders(data.adj.rlm)

## get residuals
res.raw <- residuals(data.raw.rlm)
res.adj <- residuals(data.adj.rlm)

## get weights
w.raw <- weights(data.raw.rlm)
w.adj <- weights(data.adj.rlm)
}
}
\keyword{manip}
