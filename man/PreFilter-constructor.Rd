\name{PreFilter-constructor}
\alias{PreFilter}
\title{Constructor for Class PreFilter}
\description{Constructor for class  PreFilter allows to apply different filters to class
 \code{\link{ExprTreeSet}}, i.e. to the expression level data.frame \code{data}.}
\usage{
PreFilter(mad         = character(),
          cv          = character(),
          variance    = character(),
          difference  = character(),
          ratio       = character(),
          gap         = character(),
          lothreshold = character(),
          hithreshold = character(),
          quantile    = character(),
          prescall    = character())
}
\arguments{
  \item{mad}{\code{"character"} vector describing parameters for \code{\link{madFilter}}.}
  \item{cv}{\code{"character"} vector describing parameters for \code{\link{cvFilter}}.}
  \item{variance}{\code{"character"} vector describing parameters for \code{\link{varFilter}}.}
  \item{difference}{\code{"character"} vector describing parameters for \code{\link{diffFilter}}.}
  \item{ratio}{\code{"character"} vector describing parameters for \code{\link{ratioFilter}}.}
  \item{gap}{\code{"character"} vector describing parameters for \code{\link{gapFilter}}.}
  \item{lothreshold}{\code{"character"} vector describing parameters for \code{\link{lowFilter}}.}
  \item{hithreshold}{\code{"character"} vector describing parameters for \code{\link{highFilter}}.}
  \item{quantile}{\code{"character"} vector describing parameters for \code{\link{quantileFilter}}.}
  \item{prescall}{\code{"character"} vector describing parameters for \code{\link{callFilter}}.}
}
\details{
 The PreFilter constructor allows to apply the following filters to class \code{\link{ExprTreeSet}}:
 \tabular{lll}{
  \tab \code{mad}:\tab character vector c(cutoff,epsilon). \cr
  \tab \code{cv}:\tab character vector c(cutoff,trim,epsilon). \cr
  \tab \code{variance}:\tab character vector c(cutoff,trim,epsilon). \cr
  \tab \code{difference}:\tab character vector c(cutoff,trim,epsilon). \cr
  \tab \code{ratio}:\tab character vector c(cutoff). \cr
  \tab \code{gap}:\tab character vector c(cutoff,window,trim,epsilon). \cr
  \tab \code{lothreshold}:\tab character vector c(cutoff,parameter,condition). \cr
  \tab \code{hithreshold}:\tab character vector c(cutoff,parameter,condition). \cr
  \tab \code{quantile}:\tab character vector c(cutoff,loquantile,hiquantile). \cr
  \tab \code{prescall}:\tab character vector c(cutoff,samples,condition). 
 }
}
\value{
 An object of type \code{"\linkS4class{PreFilter}"}
}
\author{Christian Stratowa}
\note{Function \code{PreFilter} is used as constructor for class \code{PreFilter} so that the user
 need not know details for creating S4 classes.
}
\seealso{\code{\linkS4class{Filter}}, \code{\link{UniFilter}}}
\examples{
## fill character vectors within constructor
prefltr <- PreFilter(mad=c(0.5,0.01), prescall=c(0.002, 6,"samples"),
                     lothreshold=c(6.0,0.02,"mean"), hithreshold=c(10.5,80.0,"percent"))
str(prefltr)

## alternatively add character vectors as methods after creation of constructor
prefltr <- PreFilter()
madFilter(prefltr) <- c(0.5,0.01)
gapFilter(prefltr) <- c(0.3,0.05,0.0,0.01)
lowFilter(prefltr) <- c(4.0,3,"samples")
highFilter(prefltr) <- c(14.5,75.0,"percent")
str(prefltr)
}
\keyword{manip}
