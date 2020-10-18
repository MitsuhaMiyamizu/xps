\name{attachInten-methods}
\docType{methods}
\alias{attachInten-methods}
\alias{attachInten}
\alias{removeInten-methods}
\alias{removeInten}
\title{Attach/Remove Intensities}
\description{
 Attach/remove raw CEL intensities to/from \code{\link{DataTreeSet}}.

\emph{Usage}

\code{attachInten(object, treenames = "*")}

\code{removeInten(object)}
}
\arguments{
  \item{object}{Object of class \code{"DataTreeSet"}.}
  \item{treenames}{Object of class \code{"list"} representing the names of the 
   \code{\link{ROOT}} data trees.}
}
\details{
 When CEL files will be imported using function \code{\link{import.data}},
 the raw intensities will be stored in \code{\link{ROOT}} data trees.
 However, the intensities will not be saved in class \code{\link{DataTreeSet}} as slot
 \code{data}, thus avoiding memory problems. Function \code{attachInten} allows to 
 fill slot \code{data} on demand.

 \code{attachInten} exports intensities from data trees from \code{\link{ROOT}} data
 file and and saves as data.frame \code{data}. \code{treenames} is a vector of
 tree names to attach; for \code{treenames="*"} all trees from slot \code{treenames}
 will be exported and intensities attached as data.frame \code{data}.

 \code{removeInten} removes intensities from \code{\link{DataTreeSet}} and replaces
 data.frame \code{data} with an empty data.frame of dim(0,0).
}
\value{
  A \code{\link{DataTreeSet}} object.
}
\author{Christian Stratowa}
\note{
 Do not use \code{attachInten} unless you know that your computer has sufficient RAM,
 especially when using exon arrays. It may be advisible to use a subset of \code{treenames}
 only.
}
\seealso{\code{\link{attachBgrd}}, \code{\link{removeBgrd}}}
\examples{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))
dim(intensity(data.test3))

data.test3 <- attachInten(data.test3)
dim(intensity(data.test3))
head(intensity(data.test3))

data.test3 <- removeInten(data.test3)
dim(intensity(data.test3))
}
\keyword{methods}
