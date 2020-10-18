\name{root.expr}
\alias{root.expr}
\title{Create class ExprTreeSet accessing ROOT expression file}
\description{
  Create class ExprTreeSet accessing ROOT expression file.
}
\usage{
root.expr(xps.scheme, rootfile = character(0), treetype = character(0), treenames = "*")
}
\arguments{
  \item{xps.scheme}{A \code{\link{SchemeTreeSet}} containing the correct scheme for the ROOT data file.}
  \item{rootfile}{name of ROOT data file, including full path.}
  \item{treetype}{tree type.}
  \item{treenames}{optional \code{character} vector of tree names to get only subset of trees.}
}
\details{
 An S4 class \code{\link{ExprTreeSet}} will be created, serving as R wrapper to the 
 existing \code{\link{ROOT}} expression file \code{rootfile}.

 Parameter \code{treetype} must be supplied to identify the ROOT trees for slot \code{data}.
 Valid tree types are listed in \code{\link{validTreetype}}.

 To get the names of all trees with their extensions \code{treetype}, which are stored in
 \code{rootfile}, you can call function \code{\link{getTreeNames}} first.

 If the \code{\link{ExprTreeSet}} should only handle a subset of the trees stored in 
 \code{rootfile}, the tree names must be supplied as vector \code{treenames}.
}
\value{
  A \code{ExprTreeSet} object.
}
\author{Christian Stratowa}
\seealso{\code{\link{root.data}}, \code{\link{root.call}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

# rma
all.rma <- rma(data.test3,"tmp_Test3RMAAll",tmpdir="",background="pmonly",normalize=TRUE,verbose=FALSE)

## use subset of trees
sub.rma <- root.expr(scheme.test3, "tmp_Test3RMAAll.root", "mdp", c("TestA2.mdp", "TestB1"))
}
\keyword{manip}
