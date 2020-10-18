\name{root.merge.data}
\alias{root.merge.data}
\title{Create class DataTreeSet by merging ROOT data files}
\description{
  Create class DataTreeSet by merging different ROOT data files.
}
\usage{
root.merge.data(xps.scheme, rootfiles = list(), celnames = "*")
}
\arguments{
  \item{xps.scheme}{A \code{\link{SchemeTreeSet}} containing the correct scheme for the ROOT data file.}
  \item{rootfiles}{\code{list} of ROOT data file(s), including full path.}
  \item{celnames}{optional \code{character} vector of tree names to get only subset of trees.}
}
\details{
This function allows to merge data trees from different existing ROOT data files.

An S4 class \code{\link{DataTreeSet}} will be created, serving as R wrapper to the 
existing \code{\link{ROOT}} data file(s) \code{rootfiles}.

If the \code{\link{DataTreeSet}} should only handle a subset of the trees stored in 
\code{rootfiles}, the tree names must be supplied as vector \code{celnames}.

To get the names of all trees stored in separate \code{rootfiles} you can call function 
\code{\link{getTreeNames}} first.
}
\value{
  A \code{DataTreeSet} object.
}
\author{Christian Stratowa}
\seealso{\code{\link{root.data}}, \code{\link{DataTreeSet}}}
\examples{
## get scheme and import CEL-files from package
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- import.data(scheme.test3,"tmp_datatest3",celdir=paste(path.package("xps"),"raw",sep="/"),verbose=FALSE)

## get subset of CEL-files
subdataA <- root.data(scheme.test3,"tmp_datatest3_cel.root", celnames=c("TestA1.cel","TestA2.cel"))
subdataB <- root.data(scheme.test3,"tmp_datatest3_cel.root", celnames=c("TestB1.cel","TestB2.cel"))

## merge data
dataAB <- root.merge.data(scheme.test3,c(rootFile(subdataA),rootFile(subdataB)), celnames=c("TestB1.cel","TestA2.cel","TestB2.cel"))
}
\keyword{manip}
