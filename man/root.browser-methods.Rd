\name{root.browser-methods}
\docType{methods}
\alias{root.browser}
\title{Open the ROOT object browser}
\description{
Open the \code{\link{ROOT}} object browser to see all objects stored in a \code{\link{ROOT}} 
file including \code{\link{ROOT}} trees.
}
\usage{
root.browser(object)
}
\arguments{
  \item{object}{an object of type \code{\link{SchemeTreeSet}}, \code{\link{DataTreeSet}},
  \code{\link{ExprTreeSet}}, or \code{\link{CallTreeSet}}}
}
\author{Christian Stratowa}
\note{
Always select menu item \dQuote{Quit ROOT} from menu \dQuote{File} to close the ROOT browser,
otherwise you are in the CINT C/C++ interpreter from \code{\link{ROOT}}. To exit CINT,
you need to type \dQuote{.q}.
}
\keyword{methods}
