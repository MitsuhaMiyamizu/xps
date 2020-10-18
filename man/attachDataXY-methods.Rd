\name{attachDataXY-methods}
\docType{methods}
\alias{attachDataXY-methods}
\alias{attachDataXY}
\alias{removeDataXY-methods}
\alias{removeDataXY}
\title{Attach/Remove (X,Y)-Coordinates}
\description{
 Attach/remove (x,y)-coordinates of raw CEL-files to/from \code{\link{DataTreeSet}}.

\emph{Usage}

\code{attachDataXY(object)}

\code{removeDataXY(object)}
}
\arguments{
  \item{object}{Object of class \code{"DataTreeSet"}.}
}
\details{
 \code{attachDataXY} exports (x,y)-coordinates only from data tree of \code{\link{ROOT}} data file
 and and saves it as \code{data.frame} in slot \code{data}.

 \code{removeDataXY} removes (x,y)-coordinates from slot \code{data} and replaces data.frame \code{data}
 with an empty data.frame of dim(0,0).
}
\value{
  A \code{\link{DataTreeSet}} object.
}
\author{Christian Stratowa}
\seealso{\code{\link{attachInten}}, \code{\link{removeInten}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## attach (x,y)-coordinates
data.test3 <- attachDataXY(data.test3)

## get data.frame
xy <- treeData(data.test3)
head(xy)

## remove (x,y)-coordinates
data.test3 <- removeDataXY(data.test3)

rm(scheme.test3, data.test3)
gc()
}
\keyword{methods}
