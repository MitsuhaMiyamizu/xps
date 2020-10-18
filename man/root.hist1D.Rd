\name{root.hist1D}
\alias{root.hist1D}
\title{ROOT 1D-Histogram}
\description{
  Creates a ROOT 1D-histogram for a ROOT tree.
}
\usage{
root.hist1D(x, treename = character(0), logbase = "log2", type = "hist", option = "HIST", canvasname = "Histogram1D", save.as = "", w = 540, h = 540)
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}} or \code{\link{ExprTreeSet}}.}
  \item{treename}{name of tree, must be present in \code{rootfile} of object \code{x}.}
  \item{logbase}{usually \dQuote{log2}, or \dQuote{0}, determines if leaf data should be converted to log.}
  \item{type}{ROOT 1D-hist or density, i.e. \dQuote{hist} or \dQuote{density}.}
  \item{option}{ROOT 1D-hist option only, usually one of \dQuote{HIST}, \dQuote{B}, \dQuote{C}, \dQuote{E}.}
  \item{canvasname}{name of ROOT canvas}
  \item{save.as}{graphics type for saving canvas, one of \dQuote{ps}, \dQuote{eps}, \dQuote{pdf},
   \dQuote{jpg}, \dQuote{gif}, \dQuote{png}, \dQuote{tiff}}
  \item{w}{the width of the canvas in pixels.}
  \item{h}{the height of the canvas in pixels.}
}
\details{
 Creates a ROOT 1D-histogram for tree \code{treename} present in \code{rootfile}.

 By selecting menu \dQuote{File->Save->canvasname.xxx} you can save the figure as e.g. *gif, *.jpg, *.pdf, *.ps or even as C++ macro.

 Alternatively, you can save the plot by setting \code{save.as}. However, this will close the
 canvas immediately after opening it.
}
\author{Christian Stratowa}
\note{
Always select menu item \dQuote{Quit ROOT} from menu \dQuote{File} to close the ROOT canvas,
otherwise you are in the CINT C/C++ interpreter from \code{\link{ROOT}}. To exit CINT,
you need to type \dQuote{.q}.
}
\seealso{\code{\link{root.hist2D}}, \code{\link{root.hist3D}}}
\examples{
\dontrun{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

root.hist1D(data.test3, "TestA1.cel")
root.hist1D(data.test3, "TestA1.cel", type="density")
}
}
\keyword{device}
