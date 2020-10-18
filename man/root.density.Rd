\name{root.density}
\alias{root.density}
\title{ROOT Density Plot}
\description{
  Creates a ROOT density plot for one or all ROOT tree(s).
}
\usage{
root.density(x, treename = "*", logbase = "log2", canvasname = "DensityPlot", save.as = "", w = 540, h = 540)
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}} or \code{\link{ExprTreeSet}}.}
  \item{treename}{name of tree, must be present in \code{rootfile} of object \code{x}.}
  \item{logbase}{usually \dQuote{log2}, or \dQuote{0}, determines if leaf data should be converted to log.}
  \item{canvasname}{name of ROOT canvas}
  \item{save.as}{graphics type for saving canvas, one of \dQuote{ps}, \dQuote{eps}, \dQuote{pdf},
   \dQuote{jpg}, \dQuote{gif}, \dQuote{png}, \dQuote{tiff}}
  \item{w}{the width of the canvas in pixels.}
  \item{h}{the height of the canvas in pixels.}
}
\details{
 Creates a ROOT density plot for one or all tree(s) present in \code{rootfile}.

 By selecting menu \dQuote{File->Save->canvasname.xxx} you can save the figure as e.g. *gif, *.jpg,
 *.pdf, *.ps or even as C++ macro.

 Alternatively, you can save the plot by setting \code{save.as}. However, this will close the
 canvas immediately after opening it.
}
\author{Christian Stratowa}
\note{
Always select menu item \dQuote{Quit ROOT} from menu \dQuote{File} to close the ROOT canvas,
otherwise you are in the CINT C/C++ interpreter from \code{\link{ROOT}}. To exit CINT,
you need to type \dQuote{.q}.
}
\seealso{\code{\link{root.hist1D}}}
\examples{
\dontrun{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

root.density(data.test3, "*")
root.density(data.test3, "TestA1.cel")
root.density(data.test3, "TestA1.cel", save.as="png")
}
}
\keyword{device}
