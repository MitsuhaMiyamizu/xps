\name{fitQC}
\alias{fitQC}
\alias{xpsQualityControl-methods}
\alias{xpsQualityControl}
\title{Functions for fitting probe-level models}
\description{
 This function allows to combine different algorithms to compute background correction, normalization
 and fit a multichip model for summarization.
}
\usage{
fitQC(xps.data,
      filename = character(),
      filedir = getwd(),
      tmpdir = "",
      update = FALSE,
    # background correction
      bgcorrect.method  = "rma",
      bgcorrect.select  = "none",
      bgcorrect.option  = "pmonly:epanechnikov",
      bgcorrect.params  = c(16384),
    # normalization
      normalize.method  = "quantile",
      normalize.select  = "pmonly",
      normalize.option  = "transcript:together:none",
      normalize.logbase = "0",
      normalize.params  = c(0.0),
    # quality control
      qualify.method    = "rlm",
      qualify.select    = "pmonly",
      qualify.qualopt   = "all",
      qualify.option    = "transcript",
      qualify.estimator = "huber",
      qualify.logbase   = "log2",
      qualify.params    = list(10, 0.01, 1.0),
    # reference values
      reference.index   = 0,
      reference.method  = "mean",
      reference.params  = list(0.0),
    # misc.
      exonlevel  = "",
      xps.scheme = NULL,
      add.data   = FALSE,
      bufsize    = 32000,
      verbose    = TRUE)

xpsQualityControl(object, ...)
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
  \item{qualify.method}{qualification method to use, currently \code{rlm}.}
  \item{qualify.select}{type of probes to select for qualification.}
  \item{qualify.qualopt}{option determining the data to which to apply qualification, one of 
   \sQuote{raw}, \sQuote{adjusted}, \sQuote{normalized}, \sQuote{all}.}
  \item{qualify.option}{option determining the grouping of probes for qualification, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon arrays only.}
  \item{qualify.estimator}{option determining the M-estimator to use, one of \sQuote{huber},
   \sQuote{fair}, \sQuote{cauchy}, \sQuote{gemanmcclure}, \sQuote{welsch}, \sQuote{tukey}, \sQuote{andrew}.}
  \item{qualify.logbase}{logarithm base as character, one of \sQuote{0}, \sQuote{log},
   \sQuote{log2}, \sQuote{log10}.}
  \item{qualify.params}{vector of parameters for qualification method.}
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
 This function allows to combine different algorithms to compute background correction, normalization
 and fit a multichip model for summarization.

 \code{xpsQualityControl} is the \code{DataTreeSet} method called by function \code{fitQC},
 containing the same parameters.
}
\value{
 An object of type \code{\link{QualTreeSet}}.
}
\author{Christian Stratowa}
\seealso{\code{\link{fitRLM}}, \code{\link{qualify}}, \code{\link{express}}}
\examples{
\dontrun{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## qualification - rlm
rlm.all <- fitQC(data.test3, "tmp_Test3RLMall", filedir=getwd(), tmpdir="",
                 qualify.method="rlm", qualify.qualopt="all", qualify.option="transcript", add.data=FALSE)

## get expression data.frame
expr.rlm.all <- validData(rlm.all)

## get borders
brd.rlm.all <- borders(rlm.all)

## get residuals
res.rlm.all <- residuals(rlm.all)

## get weights
w.rlm.all <- weights(rlm.all)

## plot expression levels
if (interactive()) {
coiplot(rlm.all)
borderplot(rlm.all)
nuseplot(rlm.all)
rleplot(rlm.all)
image(rlm.all, type="resids")
}
}
}
\keyword{manip}
