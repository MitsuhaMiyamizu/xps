\name{mas4}
\alias{mas4}
\alias{xpsMAS4-methods}
\alias{xpsMAS4}
\title{MAS 4.0 Expression Measure}
\description{
 This function converts a \code{\link{DataTreeSet}} into an \code{\link{ExprTreeSet}}
 using the XPS implementation of Affymetrix's MAS 4.0 expression measure.
}
\usage{
mas4(xps.data,
     filename   = character(0),
     filedir    = getwd(),
     tmpdir     = "",
     normalize  = FALSE,
     sc         = 500,
     option     = "transcript",
     exonlevel  = "",
     update     = FALSE,
     xps.scheme = NULL,
     add.data   = TRUE,
     verbose    = TRUE)

xpsMAS4(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{DataTreeSet}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{normalize}{logical. If \code{TRUE} scale normalization is used after an \code{ExprTreeSet} is obtained.}
  \item{sc}{value at which all arrays will be scaled to.}
  \item{option}{option determining the grouping of probes for summarization, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon arrays only.}
  \item{exonlevel}{exon annotation level determining which probes should be used for summarization; 
   exon/genome arrays only.}
  \item{update}{logical. If \code{TRUE} the existing ROOT data file \code{filename} will be updated.}
  \item{xps.scheme}{optional alternative \code{SchemeTreeSet}.}
  \item{add.data}{logical. If \code{TRUE} expression data will be included as slot \code{data}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{DataTreeSet}.}
  \item{\dots}{arguments \code{filename},\code{filedir},\code{tmpdir},\code{option},\code{exonlevel},\code{xps.scheme}.}
}
\details{
This function computes the Affymetrix MAS 4.0 expression measure, i.e. the 
\sQuote{Average Difference} expression level, as implemented in XPS.

If \code{normalize=TRUE} then the expression levels will be scaled to \code{sc}. 
For \code{sc=0} the expression levels will be scaled to the mean expression level.

\code{xpsMAS4} is the \code{DataTreeSet} method called by function \code{mas4}, however, 
expression levels will not be scaled to a common mean expression level.

For further details see \code{\link{mas5}}.
}
\value{
 An \code{\link{ExprTreeSet}}
}
\author{Christian Stratowa}
\note{
In contrast to function \code{mas4}, expression levels computed with \code{xpsMAS4} 
will not be scaled to a common mean expression level.
}
\references{
Affymetrix (1999) GeneChip Expression Analysis Algorithm Tutorial, Affymetrix Inc., 
Santa Clara, CA.
}
\seealso{\code{\link{xpsMAS4}}, \code{\link{express}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

data.mas4 <- mas4(data.test3,"tmp_Test3MAS4",tmpdir="",normalize=TRUE,sc=500, update=TRUE,verbose=FALSE)

## get data.frame
expr.mas4 <- validData(data.mas4)
head(expr.mas4)

## plot results (negative expression values!)
if (interactive()) {
boxplot(expr.mas4)
}

rm(scheme.test3, data.test3)
gc()
}
\keyword{manip}
