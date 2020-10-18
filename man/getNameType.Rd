\name{getNameType}
\alias{getNameType}
\title{Get Chip Name and Type}
\description{
 Get chip name and type from ROOT scheme file.
}
\usage{
getNameType(rootfile)
}
\arguments{
  \item{rootfile}{name of ROOT scheme file, including full path.}
}
\details{
 Extracts the chip name and type directly from \code{\link{ROOT}} scheme file \code{rootfile}.
}
\value{
a \code{list} with parameters:
  \item{chipname}{chip name.}
  \item{chiptype}{chip type, either \sQuote{GeneChip} or \sQuote{ExonChip}.}
}
\author{Christian Stratowa}
\seealso{\code{\link{getChipName}}, \code{\link{getChipType}}}
\examples{
## correct usage
getNameType(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
## incorrect usage
getNameType(paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))
}
\keyword{manip}
