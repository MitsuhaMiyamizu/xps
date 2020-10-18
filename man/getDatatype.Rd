\name{getDatatype}
\alias{getDatatype}
\title{Get Data Type}
\description{
 Get data type corresponding to tree type.
}
\usage{
getDatatype(treetype)
}
\arguments{
  \item{treetype}{tree type.}
}
\details{
 Get data type corresponding to tree type.
 Valid tree types are described in \code{\link{validTreetype}}.
}
\value{
a \code{character} with the correct data type, 
i.e. \sQuote{rawdata}, \sQuote{preprocess} or \sQuote{normation}.
}
\author{Christian Stratowa}
\seealso{\code{\link{type2Exten}}, \code{\link{validTreetype}}}
\examples{
getDatatype("cel")
getDatatype("tbw")
}
\keyword{manip}
