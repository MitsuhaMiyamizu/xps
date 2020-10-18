\name{nuseplot-methods}
\docType{methods}
\alias{nuseplot-methods}
\alias{nuseplot}
\title{Box Plots of Normalized Unscaled Standard Errors (NUSE)}
\description{
  Produce boxplot of Normalized Unscaled Standard Errors (NUSE) for the set of arrays.

\emph{Usage}

\code{
nuseplot(x,
         which    = "UnitName",
         size     = 0,
         range    = 0,
         names    = "namepart",
         main     = "NUSE Plot",
         ylim     = c(0.8, 1.2),
         las      = 2,
         add.line = TRUE,
         outline  = FALSE,
         ...) 
}
}
\arguments{
  \item{x}{object of class  or \code{\link{QualTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}} (only \code{ExprTreeSet}).}
  \item{size}{length of sequence to be generated as subset (only \code{ExprTreeSet}).}
  \item{range}{determines how far the plot whiskers extend out from the box.}
  \item{names}{optional vector of sample names.}
  \item{main}{the main title for the plot.}
  \item{ylim}{range for the plotted y values.}
  \item{las}{the style of axis labels.}
  \item{add.line}{\code{logical}, if TRUE a horizontal line is drawn.}
  \item{outline}{if \code{outline} is not true, the outliers are not drawn (only \code{ExprTreeSet}).}
  \item{\dots}{optional arguments to be passed to \code{boxplot}.}
}
\details{
 Create boxplots of Normalized Unscaled Standard Errors (NUSE) for the set of arrays.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as boxplot.

 If an object of class \code{\link{QualTreeSet}} was created by fitting a probe level model 
 with \code{qualopt="all"} then \code{nuseplot} will plot NUSE for "all" quality options. If 
 you want to plot NUSE for a certain quality option only, e.g. "normalized" data only, then 
 you can  use parameter \code{names} with \code{names="namepart:<qualopt>"}, 
 e.g. \code{names="namepart:normalized"}. 
}
\author{Christian Stratowa}
\seealso{\code{\link{NUSE}}, \code{\link{plotNUSE}}, \code{\link{rleplot}}}
\examples{
# load existing ROOT scheme file and ROOT expression file for rma
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.rma <- root.expr(scheme.test3, paste(path.package("xps"),"rootdata/tmp_Test3RMA.root",sep="/"), "mdp")

if (interactive()) {
nuseplot(data.rma)
}
}
\keyword{methods}
