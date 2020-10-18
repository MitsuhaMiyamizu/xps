\name{UniFilter-constructor}
\alias{UniFilter}
\title{Constructor for Class UniFilter}
\description{Constructor for class UniFilter allows to apply different unitest filters to class
 \code{\link{ExprTreeSet}}, i.e. to the expression level data.frame \code{data}.}
\usage{
UniFilter(unitest    = "t.test",
          foldchange = character(),
          prescall   = character(),
          unifilter  = character())
}
\arguments{
  \item{unitest}{\code{"character"} vector describing parameters for \code{\link{uniTest}}.}
  \item{foldchange}{\code{"character"} vector describing parameters for \code{\link{fcFilter}}.}
  \item{prescall}{\code{"character"} vector describing parameters for \code{\link{callFilter}}.}
  \item{unifilter}{\code{"character"} vector describing parameters for \code{\link{unitestFilter}}.}
}
\details{
 The UniFilter constructor allows to apply the following unitest filters to class \code{\link{ExprTreeSet}}:
 \tabular{lll}{
  \tab \code{unitest}:\tab character vector c(type,alternative,correction.numperm,mu,paired,conflevel,varequ). \cr
  \tab \code{foldchange}:\tab character vector c(cutoff,direction). \cr
  \tab \code{prescall}:\tab character vector c(cutoff,samples,condition). \cr
  \tab \code{unifilter}:\tab character vector c(cutoff,variable). 
 }
}
\value{
 An object of type \code{"\linkS4class{UniFilter}"}
}
\author{Christian Stratowa}
\note{Function \code{UniFilter} is used as constructor for class \code{UniFilter} so that the user
 need not know details for creating S4 classes.
}
\seealso{\code{\linkS4class{UniFilter}}, \code{\link{PreFilter}}}
\examples{
## fill character vectors within constructor
unifltr <- UniFilter(unitest=c("t.test","two.sided","none",0,0.0,FALSE,0.95,TRUE),
                     foldchange=c(1.3,"both"),unifilter=c(0.1,"pval"))
str(unifltr)

## alternatively add character vectors as methods after creation of constructor
unifltr <- UniFilter()
fcFilter(unifltr) <- c(1.5,"both")
unitestFilter(unifltr) <- c(0.01,"pval")
str(unifltr)
}
\keyword{manip}
