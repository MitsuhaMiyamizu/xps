\name{getTreeNames}
\alias{getTreeNames}
\title{Get Tree Names}
\description{
 Get tree names stored in a ROOT file.
}
\usage{
getTreeNames(rootfile, treetype = "*", setname = NULL, gettitle = FALSE)
}
\arguments{
  \item{rootfile}{name of ROOT file, including full path.}
  \item{treetype}{tree type.}
  \item{setname}{name of ROOT subdirectory containing trees.}
  \item{gettitle}{If \code{TRUE} the titles of the trees will be returned.}
}
\details{
 Extracts the tree names of \code{treetype} stored in \code{\link{ROOT}} file \code{rootfile}. \cr
 Valid tree types are listed in \code{\link{validTreetype}}. For \code{treetype="*"} names for all
 trees in \code{rootfile} are returned. \cr
 If \code{setname} is provided, only tree names in subdirectory \code{setname} are returned.
}
\value{
A \code{vector} of tree names. For \code{gettitle=TRUE} a \code{vector} of tree titles.
}
\author{Christian Stratowa}
\examples{
getTreeNames(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
getTreeNames(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"), "scm")
getTreeNames(paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))
}
\keyword{manip}
