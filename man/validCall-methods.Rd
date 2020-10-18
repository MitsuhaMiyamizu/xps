\name{validCall-methods}
\docType{methods}
\alias{validCall}
\alias{validPVal}
\title{Get Valid Detection Call Values}
\description{
 Extracts valid present call values with unit names as row names.

\emph{Usage}

\code{validCall(object, which = "UnitName")} \cr
\code{validPVal(object, which = "UnitName")}
}
\arguments{
  \item{object}{object of class \code{CallTreeSet}.}
  \item{which}{name of column containing unit name.}
}
\details{
 Method \code{validCall} returns the present calls from slot \code{detcall} as \code{data.frame} 
 and uses column \code{which} as row names, usually the probeset IDs stored in column 
 \dQuote{UnitName}.

 Method \code{validPVal} returns the detection call p-values from slot \code{data} as 
 \code{data.frame} and uses column \code{which} as row names, usually the probeset IDs 
 stored in column \dQuote{UnitName}.
}
\value{
  A \code{\link{data.frame}}.
}
\author{Christian Stratowa}
\seealso{\code{\link{validData}}, \code{\link{validExpr}}}
\keyword{methods}
