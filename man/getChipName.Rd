\name{getChipName}
\alias{getChipName}
\title{Get Chip Name}
\description{
 Get chip name from ROOT scheme file.
}
\usage{
getChipName(rootfile)
}
\arguments{
  \item{rootfile}{name of ROOT scheme file, including full path.}
}
\details{
 Extracts the chip name directly from \code{\link{ROOT}} scheme file \code{rootfile}.
}
\value{
a \code{character} with the chip name.
}
\author{Christian Stratowa}
\seealso{\code{\link{getChipType}}, \code{\link{getNameType}}}
\examples{
## correct usage
getChipName(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
## incorrect usage
getChipName(paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))
}
\keyword{manip}
