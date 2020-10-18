\name{root.profile}
\alias{root.profile}
\title{ROOT Profile Plot}
\description{
  Creates a ROOT profile plot, i.e. a plot of parallel coordinates
}
\usage{
root.profile(x, treename = "*", varlist = NULL, as.log = TRUE, globalscale = TRUE, boxes = TRUE, ylim = NULL, canvasname = "ProfilePlot", save.as = "", w = 800, h = 600)
}
\arguments{
  \item{x}{S4 object, usually of class \code{\link{DataTreeSet}} or \code{\link{ExprTreeSet}}.}
  \item{treename}{name of tree, usually all trees present in \code{rootfile} of object \code{x}.}
  \item{varlist}{leaf name of tree, usual \dQuote{fInten} or \dQuote{fLevel}.}
  \item{as.log}{logical indicating if \code{varlist} should be drawn as logarithmic data.}
  \item{globalscale}{logical indicating if all axes should be drawn at the same scale.}
  \item{boxes}{logical indicating if box-and-whisker plots should be drawn.}
  \item{ylim}{size limits c(min,max) of \code{varlist}.}
  \item{canvasname}{name of ROOT canvas}
  \item{save.as}{graphics type for saving canvas, one of \dQuote{ps}, \dQuote{eps}, \dQuote{pdf},
   \dQuote{jpg}, \dQuote{gif}, \dQuote{png}, \dQuote{tiff}}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
}
\details{
 Creates a ROOT profile plot for all trees \code{treename="*"} present in \code{rootfile}, or for 
 a subset of trees. In this case \code{varlist} must be the name of one tree leaf only;
 for \code{varlist=NULL} leaf "fInten" will be used for class \code{\link{DataTreeSet}} and
 leaf "fLevel" will be used for class \code{\link{ExprTreeSet}}.

 If \code{treename} is the name of one tree only then \code{varlist} can contain up to all leaves
 of the tree, separated by colons, e.g. \code{varlist="fLevel:fStdev"}.

 For \code{boxes=TRUE} the profile plot draws box-and-whisker plots and can thus be considered the
 equivalent of the  usual \code{boxplot}.

 A ROOT profile plot, i.e. a plot of parallel coordinates, is drawn in a \dQuote{TreeViewer}, a
 graphic user interface designed to handle \code{\link{ROOT}} trees. You can activate context menus
 by right-clicking on items or inside the right panel.

 The \dQuote{TreeViewer} is explained in \url{http://root.cern.ch/root/html/TTreeViewer.html}.

 By selecting menu \dQuote{File->Save->canvasname.xxx} you can save the figure as e.g. *gif, *.jpg,
 *.pdf, *.ps or even as C++ macro.

 Alternatively, you can save the plot by setting \code{save.as}. However, this will close the
 canvas immediately after opening it.
}
\author{Christian Stratowa}
\note{
Always select menu item \dQuote{Quit ROOT} from menu \dQuote{File} to close the ROOT tree viewer,
otherwise you are in the CINT C/C++ interpreter from \code{\link{ROOT}}. To exit CINT,
you need to type \dQuote{.q}.
}
\examples{
\dontrun{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

root.profile(data.test3)
}
}
\keyword{device}
