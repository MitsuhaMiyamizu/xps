\name{attachData-methods}
\docType{methods}
\alias{attachData-methods}
\alias{attachData}
\alias{removeData-methods}
\alias{removeData}
\title{Attach/Remove Data}
\description{
 Attach/remove data from trees to/from \code{\link{ProcesSet}}.

\emph{Usage}

\code{attachData(object, treenames = character(0), varlist = character(0), outfile = "data.txt")}

\code{removeData(object)}
}
\arguments{
  \item{object}{Object of class \code{"ProcesSet"}.}
  \item{treename}{vector of tree names to export.}
  \item{varlist}{names of tree leaves to export}
  \item{outfile}{name of output file.}
}
\details{
 \code{attachData} exports \code{varlist} from tree(s) with \code{treenames} and and saves the result
 as \code{data.frame} in slot \code{data}. Possible values of parameter \code{varlist} are described
 in \code{\link{export}}.

 \code{removeData} removes data from slot \code{data} and replaces data.frame \code{data}
 with an empty data.frame of dim(0,0).
}
\value{
  A \code{\link{ProcesSet}} object.
}
\author{Christian Stratowa}
\seealso{\code{\link{attachDataXY}}, \code{\link{attachInten}}}
\keyword{methods}
