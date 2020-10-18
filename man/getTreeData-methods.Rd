\name{getTreeData-methods}
\docType{methods}
\alias{getTreeData-methods}
\alias{getTreeData}
\title{Export Tree Data}
\description{
 Exports tree data from \code{\link{ROOT}} data file and and saves as \code{data.frame}.

\emph{Usage}

\code{
getTreeData(object, treetype = "cel", varlist = "fInten")
}
}
\arguments{
  \item{object}{Object of class \code{"ProcesSet"}.}
  \item{treetype}{type of tree to export, see \code{\link{validTreetype}}}
  \item{varlist}{names of tree leaves to export.}
}
\details{
Exports tree leaves from \code{\link{ROOT}} data file and and saves as \code{data.frame}.
}
\value{
  A \code{\link{data.frame}}.
}
\author{Christian Stratowa}
\seealso{\code{\link{export}}}
\keyword{methods}
