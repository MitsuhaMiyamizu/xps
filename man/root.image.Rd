\name{root.image}
\alias{root.image}
\title{ROOT Image}
\description{
  Creates a ROOT image for a ROOT tree.
}
\usage{
root.image(x, treename = character(0), leafname = "fInten", logbase = "log2", option = "COLZ", zlim = NULL, canvasname = "Image", save.as = "", w = 540, h = 540)
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}}.}
  \item{treename}{name of tree, must be present in \code{rootfile} of object \code{x}.}
  \item{leafname}{leaf name of tree, usual \dQuote{fInten} or \dQuote{fBg}.}
  \item{logbase}{usually \dQuote{log2}, or \dQuote{0}, determines if leaf data should be converted to log.}
  \item{option}{ROOT graph option, usually. one of \dQuote{COL}, \dQuote{COLZ}.}
  \item{zlim}{size limits c(min,max) of \code{leafname}.}
  \item{canvasname}{name of ROOT canvas}
  \item{save.as}{graphics type for saving canvas, one of \dQuote{ps}, \dQuote{eps}, \dQuote{pdf},
   \dQuote{jpg}, \dQuote{gif}, \dQuote{png}, \dQuote{tiff}}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
}
\details{
 Creates a ROOT image for tree \code{treename} present in \code{rootfile}.

 To zoom-in move the mouse cursor to the x-axis (y-axis) until it changes to a hand and click-drag
 to select an axis-range. 
 To unzoom move the mouse cursor to the x-axis (y-axis) until it changes to a hand and right-click
 to select \dQuote{Unzoom}.

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
\seealso{\code{\link{image-methods}}, \code{\link{image}}}
\examples{
\dontrun{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

root.image(data.test3, "TestA1.cel")
root.image(data.test3, "TestA1.cel", save.as="png")
}
}
\keyword{device}
