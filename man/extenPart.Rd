\name{extenPart}
\alias{extenPart}
\title{Get Extension of Tree Names}
\description{
 Get the extension(s) of (tree) names.
}
\usage{
extenPart(names, as.unique=TRUE)
}
\arguments{
  \item{names}{vector of names.}
  \item{as.unique}{if TRUE return only unique extensions.}
}
\details{
 Extracts the extension part of \code{names}, e.g.of tree names of \code{treename.treetype} stored 
 in a \code{\link{ROOT}} file.
}
\value{
A \code{vector} of (unique) extensions.
}
\author{Christian Stratowa}
\seealso{\code{\link{namePart}}}
\examples{
names <- c("TestA1.int", "TestA2.int")
extenPart(names)
extenPart(names, as.unique=FALSE)
}
\keyword{manip}
