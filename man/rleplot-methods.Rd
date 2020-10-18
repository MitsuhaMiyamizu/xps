\name{rleplot-methods}
\docType{methods}
\alias{rleplot-methods}
\alias{rleplot}
\title{Box Plots of Relative Log Expression (RLE)}
\description{
  Produce boxplots of Relative Log Expression (RLE) values for the set of arrays.

\emph{Usage}

\code{
rleplot(x,
        which    = "UnitName",
        size     = 0,
        range    = 0,
        names    = "namepart",
        main     = "RLE Plot",
        ylim     = c(-1.0, 1.0),
        las      = 2,
        add.line = TRUE,
        outline  = FALSE,
        ...) 
}
}
\arguments{
  \item{x}{object of class \code{\link{ExprTreeSet}} or \code{\link{QualTreeSet}}.}
  \item{which}{type of probes to be used, for details see \code{\link{validData}}.}
  \item{size}{length of sequence to be generated as subset.}
  \item{range}{determines how far the plot whiskers extend out from the box.}
  \item{names}{optional vector of sample names.}
  \item{main}{the main title for the plot.}
  \item{ylim}{range for the plotted y values.}
  \item{las}{the style of axis labels.}
  \item{add.line}{\code{logical}, if TRUE a horizontal line is drawn.}
  \item{outline}{if \code{outline} is not true, the outliers are not drawn.}
  \item{\dots}{optional arguments to be passed to \code{boxplot}.}
}
\details{
 Create boxplots of Relative Log Expression (RLE) values for the set of arrays, i.e. of M plots,
 where M is determined relative to a pseudo-median reference chip.

 For \code{names=NULL} full column names of slot \code{data} will be displayed while for 
 \code{names="namepart"} column names will be displayed without name extension. If \code{names}
 is a vector of column names, only these columns will displayed as boxplot.

 If an object of class \code{\link{QualTreeSet}} was created by fitting a probe level model 
 with \code{qualopt="all"} then \code{rleplot} will plot RLE for "all" quality options. If 
 you want to plot RLE for a certain quality option only, e.g. "normalized" data only, then 
 you can  use parameter \code{names} with \code{names="namepart:<qualopt>"}, 
 e.g. \code{names="namepart:normalized"}. 
}
\author{Christian Stratowa}
\seealso{\code{\link{RLE}}, \code{\link{plotRLE}}, \code{\link{mboxplot}}, \code{\link{nuseplot}}}
\examples{
# load existing ROOT scheme file and ROOT expression file for rma
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.rma <- root.expr(scheme.test3, paste(path.package("xps"),"rootdata/tmp_Test3RMA.root",sep="/"), "mdp")

if (interactive()) {
rleplot(data.rma)
}
}
\keyword{methods}
