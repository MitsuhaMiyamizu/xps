\name{getNumberTrees}
\alias{getNumberTrees}
\title{Get Number of Trees}
\description{
 Get number of trees stored in a ROOT file.
}
\usage{
getNumberTrees(rootfile, treetype = "*", setname = NULL)
}
\arguments{
  \item{rootfile}{name of ROOT file, including full path.}
  \item{treetype}{tree type.}
  \item{setname}{name of ROOT subdirectory containing trees.}
}
\details{
 Extracts the number of trees of \code{treetype} stored in \code{\link{ROOT}} file \code{rootfile}. \cr
 Valid tree types are listed in \code{\link{validTreetype}}. For \code{treetype="*"} the total number of trees in \code{rootfile}
 are returned. \cr
 If \code{setname} is provided, only trees in subdirectory \code{setname} are counted.
}
\value{
Number of trees.
}
\author{Christian Stratowa}
\examples{
getNumberTrees(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
getNumberTrees(paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))
}
\keyword{manip}
