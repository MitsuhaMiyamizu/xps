\name{fitRLM}
\alias{fitRLM}
\alias{rmaPLM}
\title{Functions for fitting RMA as probe-level model}
\description{
 Convert Affymetrix probe level data to expression levels by fitting RMA as multichip model.
}
\usage{
fitRLM(xps.data,
       filename   = character(),
       filedir    = getwd(),
       tmpdir     = "",
       background = "pmonly",
       normalize  = TRUE,
       qualopt    = "all",
       option     = "transcript",
       exonlevel  = "",
       params     = list(16384, 0.0, 1.0, 10, 0.01, 1),
       xps.scheme = NULL,
       add.data   = FALSE,
       bufsize    = 32000,
       verbose    = TRUE)

rmaPLM(xps.data,
       filename   = character(),
       filedir    = getwd(),
       tmpdir     = "",
       background = "pmonly",
       normalize  = TRUE,
       qualopt    = "all",
       option     = "transcript",
       exonlevel  = "",
       params     = list(16384, 0.0, 1.0, 10, 0.01, 1),
       xps.scheme = NULL,
       add.data   = FALSE,
       bufsize    = 32000,
       verbose    = TRUE)
}
\arguments{
  \item{xps.data}{object of class \code{DataTreeSet}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{background}{probes used to compute background, one of \sQuote{pmonly}, \sQuote{mmonly},
   \sQuote{both}; for genome/exon arrays one of \sQuote{genomic}, \sQuote{antigenomic}}
  \item{normalize}{logical. If \code{TRUE} normalize data using quantile normalization.}
  \item{qualopt}{option determining the data to which to apply qualification, one of 
   \sQuote{raw}, \sQuote{adjusted}, \sQuote{normalized}, \sQuote{all}.}
  \item{option}{option determining the grouping of probes for qualification, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon arrays only.}
  \item{exonlevel}{exon annotation level determining which probes should be used for summarization; exon/genome arrays only.}
  \item{params}{list of (default) parameters for \code{rma}.}
  \item{xps.scheme}{optional alternative \code{SchemeSet}.}
  \item{add.data}{logical. If \code{TRUE} expression data will be included as slot \code{data}.}
  \item{bufsize}{integer which sets the buffer size of the tree branch baskets (default is 32000).}
  \item{verbose}{logical, if \code{TRUE} print status information.}
}
\details{
 Convert Affymetrix probe level data to expression levels by fitting RMA as multichip model.
}
\value{
 An object of type \code{\link{QualTreeSet}}.
}
\author{Christian Stratowa}
\seealso{\code{\link{fitQC}}, \code{\link{qualify}}, \code{\link{express}}}
\examples{
\dontrun{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## qualification - rlm
rlm.all <- rmaPLM(data.test3, "tmp_Test3RLMall", filedir=getwd(), tmpdir="", qualopt="all", option="transcript", add.data=FALSE)

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
