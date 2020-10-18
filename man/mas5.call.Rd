\name{mas5.call}
\alias{mas5.call}
\alias{xpsMAS5Call-methods}
\alias{xpsMAS5Call}
\title{MAS 5.0 Absolute Detection Call}
\description{
 Performs the Wilcoxon signed rank-based gene expression presence/absence 
 detection algorithm first implemented in the Affymetrix Microarray Suite
 version 5.
}
\usage{
mas5.call(xps.data,
          filename = character(0), filedir = getwd(), tmpdir = "",
          tau = 0.015, alpha1 = 0.04, alpha2 = 0.06, ignore.saturated = TRUE, bgcorrect.option = "none", 
          option = "transcript", exonlevel = "", xps.scheme = NULL, add.data = TRUE, verbose = TRUE)

xpsMAS5Call(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{DataTreeSet}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{tau}{a small positive constant.}
  \item{alpha1}{a significance threshold in (0,alpha2).}
  \item{alpha2}{a significance threshold in (alpha1,0.5).}
  \item{ignore.saturated}{logical. If \code{TRUE} do the saturation correction described in the paper, with a saturation level of 46000.}
  \item{bgcorrect.option}{bgcorrect option determining wether to subtract background first, one of 
   \sQuote{none} or \sQuote{correctbg}.}
  \item{option}{option determining the grouping of probes for summarization, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon arrays only.}
  \item{exonlevel}{exon annotation level determining which probes should be used for summarization;
   exon/genome arrays only.}
  \item{xps.scheme}{optional alternative \code{SchemeTreeSet}.}
  \item{add.data}{logical. If \code{TRUE} call data will be added to slots \code{data} and \code{detcall}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{DataTreeSet}.}
  \item{\dots}{the arguments described above.}
}
\details{
 This function performs the hypothesis test:

  H0: median(Ri) = tau, corresponding to absence of transcript
  H1: median(Ri) > tau, corresponding to presence of transcript

 where Ri = (PMi - MMi) / (PMi + MMi) for each i a probe-pair in the probe-set
 represented by data.

 The p-value that is returned estimates the usual quantity:

   Pr(observing a more "present looking" probe-set than data | data is absent)

 Small p-values imply presence while large ones imply absence of transcript. 
 The detection call is computed by thresholding the p-value as in:

 call \dQuote{P} if p-value < alpha1 \cr
 call \dQuote{M} if alpha1 <= p-value < alpha2 \cr
 call \dQuote{A} if alpha2 <= p-value

 The defaults for \code{tau}, \code{alpha1} and \code{alpha2} correspond to those in MAS5.0
 for expression arrays. However, when using this function for exon or whole genome arrays,
 new values for \code{alpha1} and \code{alpha2} must be determined.Furthermore, in these
 cases it may be better to use \code{bgcorrect.option = "correctbg"} to get reasonable present
 calls. Note that the recommended function for exon/genome arrays is \code{\link{dabg.call}}.

 In order to use an alternative \code{\link{SchemeTreeSet}} set the corresponding
 SchemeTreeSet \code{xps.scheme}.

 \code{xpsMAS5Call} is the \code{DataTreeSet} method called by function \code{mas5.call}, 
 containing the same parameters.
}
\value{
 A \code{\link{CallTreeSet}}
}
\author{Christian Stratowa}
\references{
Liu, W. M. and Mei, R. and Di, X. and Ryder, T. B. and Hubbell, E. and Dee,
S. and Webster, T. A. and Harrington, C. A. and Ho, M. H. and Baid, J. and
Smeekens, S. P. (2002) Analysis of high density expression microarrays with
signed-rank call algorithms, Bioinformatics, 18(12), pp. 1593-1599.

Liu, W. and Mei, R. and Bartell, D. M. and Di, X. and Webster, T. A. and
Ryder, T. (2001) Rank-based algorithms for analysis of microarrays,
Proceedings of SPIE, Microarrays: Optical Technologies and Informatics, 4266.

Affymetrix (2002) Statistical Algorithms Description Document, Affymetrix Inc., 
Santa Clara, CA, whitepaper.
\url{http://www.affymetrix.com/support/technical/whitepapers/sadd_whitepaper.pdf}
}
\seealso{\code{\link{dabg.call}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## MAS5 detection call
call.mas5 <- mas5.call(data.test3,"tmp_Test3Call",tmpdir="",verbose=FALSE)

## get data.frames
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)
head(pval.mas5)
head(pres.mas5)

## plot results
if (interactive()) {
callplot(call.mas5, beside=FALSE, ylim=c(0,125))
}

rm(scheme.test3, data.test3)
gc()
}
\keyword{manip}
