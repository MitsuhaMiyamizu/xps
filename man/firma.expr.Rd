\name{firma.expr}
\alias{firma.expr}
\title{Get Expression Levels from FIRMA}
\description{
 Extracts FIRMA expression levels from data.frame \code{data}.
}
\usage{
firma.expr(xps.data,
           probeset = NULL,
           option   = "probeset")
}
\arguments{
  \item{xps.data}{object of class \code{\link{ExprTreeSet}}.}
  \item{probeset}{transcriptID or probesetID or NULL.}
  \item{option}{option determining the \code{probeset} type for which to extract expression levels, 
   one of \sQuote{transcript}, \sQuote{probeset}, \sQuote{exon}.}
}
\details{
Function \code{firma.expr} returns the expression levels from slot \code{data} for a given \code{probeset},
or for all probesets or transcripts in case of \code{probeset=NULL}. Row names will be the Affymetrix 
transcriptIDs, probesetIDs or exonIDs, respectively, dependent on the selected \code{option}.
}
\value{
  A \code{\link{data.frame}}.
}
\author{Christian Stratowa}
\note{
For \code{option="probeset"} parameter \code{probeset} should usually be the transcriptID in order to get 
the expression levels for all probesetIDs of the corresponding transcriptID.
}
\seealso{\code{\link{firma}}}
\examples{
\dontrun{
## get transcript expression levels for all transcripts or for transcript=2429277
expr.firma <- firma.expr(sub.firma.ps, probeset=NULL, option="transcript")
expr.firma <- firma.expr(sub.firma.ps, probeset=2429277, option="transcript")

## get probeset expression levels for all probeset or for probeset=2429278
expr.firma <- firma.expr(sub.firma.ps, probeset=NULL, option="probeset")
expr.firma <- firma.expr(sub.firma.ps, probeset=2429278, option="probeset")

## get probeset expression levels for all probesets corresponding to transcript=2429277
expr.firma <- firma.expr(sub.firma.ps, probeset=2429277, option="probeset")
}
}
\keyword{manip}
