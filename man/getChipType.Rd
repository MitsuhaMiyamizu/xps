\name{getChipType}
\alias{getChipType}
\title{Get Chip Type}
\description{
 Get chip type from ROOT scheme file.
}
\usage{
getChipType(rootfile)
}
\arguments{
  \item{rootfile}{name of ROOT scheme file, including full path.}
}
\details{
 Extracts the chip type directly from \code{\link{ROOT}} scheme file \code{rootfile}.
}
\value{
a \code{character} with the chip type, either \sQuote{GeneChip} or \sQuote{ExonChip}.
}
\author{Christian Stratowa}
\seealso{\code{\link{getChipName}}, \code{\link{getNameType}}}
\examples{
## correct usage
getChipType(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
## incorrect usage
getChipType(paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))
}
\keyword{manip}
