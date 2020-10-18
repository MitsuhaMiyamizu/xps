\name{attachMask-methods}
\docType{methods}
\alias{attachMask-methods}
\alias{attachMask}
\alias{removeMask-methods}
\alias{removeMask}
\title{Attach/Remove Scheme Mask}
\description{
 Attach/remove scheme mask to/from \code{\link{SchemeTreeSet}} or to slot \code{scheme}
 of \code{\link{DataTreeSet}}.
}
\usage{
attachMask(object)

removeMask(object)
}
\arguments{
  \item{object}{Object of class \code{"SchemeTreeSet"} or \code{"DataTreeSet"}.}
}
\details{
 \code{attachMask} exports mask from scheme tree from \code{\link{ROOT}} scheme file
 and and saves mask as data.frame \code{mask} of slot \code{scheme}.

 \code{removeMask} removes mask from \code{\link{SchemeTreeSet}} or from slot \code{scheme}
 of \code{\link{DataTreeSet}} and replaces data.frame \code{mask} with an empty data.frame
 of dim(0,0).
}
\value{
  A \code{\link{DataTreeSet}} object or \code{\link{SchemeTreeSet}}.
}
\author{Christian Stratowa}
\note{
 Do not use \code{attachMask} unless you know that your computer has sufficient RAM,
 especially for exon array schemes.
}
\seealso{\code{\link{import.expr.scheme}}, \code{\link{import.exon.scheme}}}
\examples{
## load existing ROOT scheme file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
dim(chipMask(scheme.test3))

scheme.test3 <- attachMask(scheme.test3)
dim(chipMask(scheme.test3))
head(chipMask(scheme.test3))

scheme.test3 <- removeMask(scheme.test3)
dim(chipMask(scheme.test3))
}
\keyword{methods}
