\name{type2Exten}
\alias{type2Exten}
\title{Convert Method Type to Tree Extension}
\description{
  Convert Method Type to Tree Extension.
}
\usage{
type2Exten(type, datatype)
}
\arguments{
  \item{type}{method type.}
  \item{datatype}{data type.}
}
\details{
For every \code{datatype} different methods, i.e. algorithms exist which can be applied.
Valid \code{datatype}s are \sQuote{preprocess} and \sQuote{normation}.

For \code{datatype} \sQuote{preprocess} the following methods can be applied:\cr
\tabular{ll}{
  \code{mean}: \tab trimmed mean \cr
  \code{median}:  \tab median \cr
  \code{quantile}:  \tab quantile \cr
  \code{tukeybiweight}:  \tab tukey biweight \cr
  \code{medianpolish}:  \tab median polish 
}

For \code{datatype} \sQuote{normation} the following methods can be applied:\cr
\tabular{ll}{
  \code{mean}: \tab trimmed mean \cr
  \code{median}:  \tab median \cr
  \code{quantile}:  \tab quantile \cr
  \code{lowess}:  \tab lowess \cr
  \code{supsmu}:  \tab supsmu 
}

The tree extensions are described in \code{\link{validTreetype}}.
}
\value{
A \code{character} with the correct tree extension.
}
\author{Christian Stratowa}
\seealso{\code{\link{getDatatype}}, \code{\link{validTreetype}}}
\examples{
type2Exten("quantile","preprocess")
type2Exten("medianpolish","preprocess")
type2Exten("supsmu","normation")
}
\keyword{manip}
