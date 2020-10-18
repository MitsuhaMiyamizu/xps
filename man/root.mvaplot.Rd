\name{root.mvaplot}
\alias{root.mvaplot}
\title{ROOT M vs A Plot}
\description{
  Creates a ROOT M vs A plot for a ROOT tree.
}
\usage{
root.mvaplot(x, treename1 = character(0), treename2 = character(0),logbase = "log2", option = "P", canvasname = "MvAPlot", save.as = "", w = 540, h = 540)
}
\arguments{
  \item{x}{object of class \code{\link{ExprTreeSet}} or \code{\link{DataTreeSet}}.}
  \item{treename1}{name of first tree, must be present in \code{rootfile} of object \code{x}.}
  \item{treename2}{name of second tree, must be present in \code{rootfile} of object \code{x}.}
  \item{logbase}{usually \dQuote{log2}, or \dQuote{0}, determines if leaf data should be converted to log.}
  \item{option}{ROOT TGraph::PaintGraph option, usually one of \dQuote{P}, \dQuote{*}.}
  \item{canvasname}{name of ROOT canvas}
  \item{save.as}{graphics type for saving canvas, one of \dQuote{ps}, \dQuote{eps}, \dQuote{pdf},
   \dQuote{jpg}, \dQuote{gif}, \dQuote{png}, \dQuote{tiff}}
  \item{w}{the width of the canvas in pixels.}
  \item{h}{the height of the canvas in pixels.}
}
\details{
 Creates a ROOT M vs A plot for trees \code{treename1} and \code{treename2} present in \code{rootfile}.

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
\seealso{\code{\link{root.graph1D}}}
\examples{
\dontrun{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

# compute RMA
data.rma <- rma(data.test3,"Test3RMA",tmpdir="",background="pmonly",normalize=TRUE)

root.mvaplot(data.rma, "TestA1.mdp", "TestB1.mdp")
}
}
\keyword{device}
