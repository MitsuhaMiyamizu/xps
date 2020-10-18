\name{validSE-methods}
\docType{methods}
\alias{validSE}
\title{Get Valid Standard Errors}
\description{
 Extracts valid standard errors with unit names as row names.

\emph{Usage}

\code{
validSE(object, which = "UnitName")
}
}
\arguments{
  \item{object}{object of class \code{ExprTreeSet}.}
  \item{which}{name of column containing unit name.}
}
\details{
Method \code{validSE} returns the standard errors (or standard deviations) from the expression trees and
uses column \code{which} as row names, usually the probeset IDs stored in column \dQuote{UnitName}.
}
\value{
  A \code{\link{data.frame}}.
}
\author{Christian Stratowa}
\seealso{\code{\link{validExpr}}}
\keyword{methods}
