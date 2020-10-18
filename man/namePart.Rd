\name{namePart}
\alias{namePart}
\title{Get Tree Names w/o Extension}
\description{
 Get (tree) names w/o their extension.
}
\usage{
namePart(names)
}
\arguments{
  \item{names}{vector of names.}
}
\details{
 Extracts the name part of \code{names}, e.g.of tree names of \code{treename.treetype} stored in a 
\code{\link{ROOT}} file.
}
\value{
A \code{vector} of tree names w/o its extension.
}
\author{Christian Stratowa}
\seealso{\code{\link{extenPart}}}
\examples{
names <- c("TestA1.int", "TestA2.int")
namePart(names)
}
\keyword{manip}
