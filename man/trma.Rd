\name{trma}
\alias{trma}
\title{transposed Robust Multi-Array Average Expression Measure}
\description{
 This function converts a \code{\link{DataTreeSet}} into an \code{\link{ExprTreeSet}}
 using the transposed robust multi-array average (RMA) expression measure.
}
\usage{
trma(xps.data,
     filename   = character(0),
     filedir    = getwd(),
     tmpdir     = "",
     background = "pmonly",
     normalize  = TRUE,
     option     = "transcript",
     exonlevel  = "",
     params     = list(16384, 0.0, 1.0, 10, 0.01, 2),
     xps.scheme = NULL,
     add.data   = TRUE,
     verbose    = TRUE)
}
\arguments{
  \item{xps.data}{object of class \code{\link{DataTreeSet}}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{background}{probes used to compute background, one of \sQuote{pmonly}, \sQuote{mmonly},
   \sQuote{both}; for genome/exon arrays one of \sQuote{genomic}, \sQuote{antigenomic}}
  \item{normalize}{logical. If \code{TRUE} normalize data using quantile normalization.}
  \item{option}{option determining the grouping of probes for summarization, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon arrays only.}
  \item{exonlevel}{exon annotation level determining which probes should be used for summarization;
   exon/genome arrays only.}
  \item{params}{list of (default) parameters for \code{rma}.}
  \item{xps.scheme}{optional alternative \code{SchemeTreeSet}.}
  \item{add.data}{logical. If \code{TRUE} expression data will be included as slot \code{data}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
}
\details{
This function computes the tRMA (transposed Robust Multichip Average) expression measure 
described in Giorgi et al. for both expression arrays and exon arrays. 

To use method \code{xpsRMA} or function \code{express} to compute \code{trma} you need 
to set \code{params = list(16384, 0.0, 1.0, 10, 0.01, 2)}.

For further details please see \code{\link{rma}}
}
\value{
 An \code{\link{ExprTreeSet}}
}
\author{Christian Stratowa}
\references{
 Federico M. Giorgi, Anthony M. Bolger, Marc Lohse and Bjoern Usadel (2010), Algorithm-driven
 Artifacts in median polish summarization of Microarray data. BMC Bioinformatics 11:553
}
\seealso{\code{\link{rma}}, \code{\link{xpsRMA}}, \code{\link{express}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

data.trma <- trma(data.test3,"tmp_Test3tRMA",tmpdir="",background="pmonly",normalize=TRUE,verbose=FALSE)

## get data.frame
expr.trma <- validData(data.trma)
head(expr.trma)

rm(scheme.test3, data.test3)
gc()
}
\keyword{manip}
