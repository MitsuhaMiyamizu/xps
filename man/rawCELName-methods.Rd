\name{rawCELName-methods}
\docType{methods}
\alias{rawCELName-methods}
\alias{rawCELName}
\title{Method for getting names of the raw CEL-files}
\description{
  Method for getting names (and full path) of the original CEL-files.

\emph{Usage}

\code{rawCELName(object, treename = "*", fullpath = TRUE)}
}
\arguments{
  \item{object}{object of class \code{DataTreeSet}.}
  \item{treename}{treename, for which the name of the original CEL-file should be returned.}
  \item{fullpath}{logical, if TRUE return full path.}
}
\details{
 Since CEL-files can be imported with \code{import.data} using alternative \code{celnames},
 method \code{rawCELName} allows to return the original name and optionally the full path for
 each CEL-file. 
}
\value{
  A character vector.
}
\author{Christian Stratowa}
\seealso{\code{\link{import.data}}}
\examples{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

rawCELName(data.test3)
rawCELName(data.test3, treename = "TestA2.cel", fullpath = FALSE)
}
\keyword{methods}
