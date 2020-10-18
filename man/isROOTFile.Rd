\name{isROOTFile}
\alias{isROOTFile}
\title{Test for ROOT File}
\description{
  Test if a file is a valid ROOT file.
}
\usage{
isROOTFile(filename)
}
\arguments{
  \item{filename}{name of ROOT file, including full path.}
}
\value{
 Return \code{TRUE} if file \code{filename} is a valid \code{\link{ROOT}} file.
}
\seealso{\code{\link{existsROOTFile}}}
\author{Christian Stratowa}
\examples{
isROOTFile(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
}
\keyword{manip}
