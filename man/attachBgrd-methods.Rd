\name{attachBgrd-methods}
\docType{methods}
\alias{attachBgrd-methods}
\alias{attachBgrd}
\alias{removeBgrd-methods}
\alias{removeBgrd}
\title{Attach/Remove Background Intensities}
\description{
 Attach/remove background intensities to/from \code{\link{DataTreeSet}}.

\emph{Usage}

\code{attachBgrd(object, treenames = "*")}

\code{removeBgrd(object)}
}
\arguments{
  \item{object}{Object of class \code{"DataTreeSet"}.}
  \item{treenames}{Object of class \code{"list"} representing the names of the 
   \code{\link{ROOT}} background trees.}
}
\details{
 Whenever one of the \code{\link{bgcorrect}} methods will be applied to raw CEL intensities,
 the background intensities will be stored in \code{\link{ROOT}} background trees.
 However, the background intensities will not be saved as data.frame \code{bgrd},
 thus avoiding memory problems. Function \code{attachBgrd} allows to fill slot \code{bgrd}
 on demand.

 \code{attachBgrd} exports intensities from background trees from \code{\link{ROOT}}
 data file and saves as data.frame \code{bgrd}. \code{treenames} is a vector of
 tree names to attach; for \code{treenames="*"} all trees from slot \code{treenames}
 will be exported and background intensities attached as data.frame \code{bgrd}.

 \code{removeBgrd} removes background intensities from \code{\link{DataTreeSet}} and
 replaces data.frame \code{bgrd} with an empty data.frame of dim(0,0).
}
\value{
  A \code{\link{DataTreeSet}} object.
}
\author{Christian Stratowa}
\note{
 Do not use \code{attachBgrd} unless you know that your computer has sufficient RAM,
 especially when using exon arrays. It may be advisible to use a subset of \code{treenames}
 only.
}
\seealso{\code{\link{attachInten}}, \code{\link{removeInten}}}
\keyword{methods}
