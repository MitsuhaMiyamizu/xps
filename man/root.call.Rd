\name{root.call}
\alias{root.call}
\title{Create class CallTreeSet accessing ROOT detection call file}
\description{
  Create class CallTreeSet accessing ROOT detection call file.
}
\usage{
root.call(xps.scheme, rootfile = character(0), treetype = character(0), treenames = "*")
}
\arguments{
  \item{xps.scheme}{A \code{\link{SchemeTreeSet}} containing the correct scheme for the ROOT data file.}
  \item{rootfile}{name of ROOT data file, including full path.}
  \item{treetype}{tree type.}
  \item{treenames}{optional \code{character} vector of tree names to get only subset of trees.}
}
\details{
 An S4 class \code{\link{CallTreeSet}} will be created, serving as R wrapper to the 
 existing \code{\link{ROOT}} detection call file \code{rootfile}.

 Parameter \code{treetype} must be supplied to identify the ROOT trees for slots \code{data}
 and \code{detcall}. Valid tree types are listed in \code{\link{validTreetype}}.

 To get the names of all trees with their extensions \code{treetype}, which are stored in
 \code{rootfile}, you can call function \code{\link{getTreeNames}} first.

 If the \code{\link{CallTreeSet}} should only handle a subset of the trees stored in 
 \code{rootfile}, the tree names must be supplied as vector \code{treenames}.
}
\value{
  A \code{CallTreeSet} object.
}
\author{Christian Stratowa}
\seealso{\code{\link{root.data}}, \code{\link{root.expr}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## MAS5 detection call
detcall.mas5 <- mas5.call(data.test3,"tmp_Test3CallAll",tmpdir="",verbose=FALSE)

## use subset of trees
sub.call <- root.call(scheme.test3, "tmp_Test3CallAll.root", "dc5", c("TestA2", "TestB1"))
}
\keyword{manip}
