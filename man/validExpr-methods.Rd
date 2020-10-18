\name{validExpr-methods}
\docType{methods}
\alias{validExpr}
\title{Get Valid Expression Levels}
\description{
 Extracts valid expression levels with unit names as row names from data.frame \code{data}.

\emph{Usage}

\code{
validExpr(object, which = "UnitName")
}
}
\arguments{
  \item{object}{object of class \code{ExprTreeSet}.}
  \item{which}{name of column containing unit name.}
}
\details{
Method \code{validExpr} returns the expression levels from slot \code{data} and
uses column \code{which} as row names, usually the probeset IDs stored in column \dQuote{UnitName}.
}
\value{
  A \code{\link{data.frame}}.
}
\author{Christian Stratowa}
\seealso{\code{\link{validData}}, \code{\link{validCall}}}
\keyword{methods}
