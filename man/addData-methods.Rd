\name{addData-methods}
\docType{methods}
\alias{addData-methods}
\alias{addData}
\title{Import additional CEL files into a DataTreeSet}
\description{
 Import additional CEL files into a DataTreeSet and update \code{\link{ROOT}} data file.

\emph{Usage}

\code{
addData(object, celdir = NULL, celfiles = "", celnames = NULL, project = NULL, verbose = TRUE)
}
}
\arguments{
  \item{object}{object of class \code{DataTreeSet}.}
  \item{celdir}{system directory containing the CEL-files for corresponding scheme.}
  \item{celfiles}{optional vector of CEL-files to be imported.}
  \item{celnames}{optional vector of names which should replace the CEL-file names.}
  \item{project}{optional class \code{\link{ProjectInfo}}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
}
\details{
 Import additional CEL-files and update \code{\link{ROOT}} data file \code{rootfile}.

 To import CEL-files from different directories, vector \code{celfiles} must contain the
 full path for each CEL-file and \code{celdir} must be \code{celdir=NULL}. 
}
\value{
  A \code{DataTreeSet} object.
}
\author{Christian Stratowa}
\seealso{\code{\link{import.data}}, \code{\link{root.data}}}
\examples{
## get scheme and import subset of CEL-files from package
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- import.data(scheme.test3,"tmp_test3",celdir=paste(path.package("xps"),"raw",sep="/"),
              celfiles=c("TestA1.CEL","TestB2.CEL"),verbose=FALSE)

unlist(treeNames(data.test3))

## add further subset of CEL-files
data.test3 <- addData(data.test3,celdir=paste(path.package("xps"),"raw",sep="/"),
              celfiles=c("TestA2.CEL","TestB1.CEL"),verbose=FALSE)

unlist(treeNames(data.test3))
}
\keyword{methods}
