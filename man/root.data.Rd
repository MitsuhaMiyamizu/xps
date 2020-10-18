\name{root.data}
\alias{root.data}
\title{Create class DataTreeSet accessing ROOT data file}
\description{
  Create class DataTreeSet accessing ROOT data file.
}
\usage{
root.data(xps.scheme, rootfile = character(0), celnames = "*")
}
\arguments{
  \item{xps.scheme}{A \code{\link{SchemeTreeSet}} containing the correct scheme for the ROOT data file.}
  \item{rootfile}{name of ROOT data file, including full path.}
  \item{celnames}{optional \code{character} vector of tree names to get only subset of trees.}
}
\details{
An S4 class \code{\link{DataTreeSet}} will be created, serving as R wrapper to the 
existing \code{\link{ROOT}} data file \code{rootfile}.

If the \code{\link{DataTreeSet}} should only handle a subset of the trees stored in 
\code{rootfile}, the tree names must be supplied as vector \code{celnames}.

To get the names of all trees stored in \code{rootfile} you can call function 
\code{\link{getTreeNames}} first.
}
\value{
  A \code{DataTreeSet} object.
}
\author{Christian Stratowa}
\note{Use \code{\link{root.data}} to access the ROOT data file from new R sessions 
to avoid creating a new \code{\link{ROOT}} data file for every R session.
}
\seealso{\code{\link{import.data}}, \code{\link{DataTreeSet}}}
\examples{
## get scheme and import CEL-files from package
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- import.data(scheme.test3,"tmp_datatest3",celdir=paste(path.package("xps"),"raw",sep="/"),verbose=FALSE)

## use subset of CEL-files
subdata.test3 <- root.data(scheme.test3,"tmp_datatest3_cel.root", celnames=c("TestA1.cel","TestB2.cel"))
}
\keyword{manip}
