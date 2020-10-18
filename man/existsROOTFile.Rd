\name{existsROOTFile}
\alias{existsROOTFile}
\title{Test for Existing ROOT File}
\description{
  Test if a ROOT file does already exist.
}
\usage{
existsROOTFile(filename, tmp.rm = TRUE)
}
\arguments{
  \item{filename}{name of ROOT file, including full path.}
  \item{tmp.rm}{logical, if \code{TRUE} then exlude filenames beginning with dQuote(tmp\_).}
}
\value{
 Return \code{TRUE} if file \code{filename} is an already existing \code{\link{ROOT}} file.
}
\note{
 It is possible to create temporary \code{\link{ROOT}} files called \dQuote{tmp} or with
 \code{filename} starting with \dQuote{tmp\_} which can be overwritten. Thus by default
 temporary files will not be recognized by \code{existsROOTFile}. If you want to recognize
 temporary files, set \code{tmp.rm = TRUE}.
}
\seealso{\code{\link{isROOTFile}}}
\author{Christian Stratowa}
\examples{
existsROOTFile(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
}
\keyword{manip}
