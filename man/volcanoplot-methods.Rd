\name{volcanoplot-methods}
\docType{methods}
\alias{volcanoplot-methods}
\alias{volcanoplot}
\title{Volcano Plot}
\description{
  Produce a scatter plot of fold-change values vs p-values, called volcano plot.

\emph{Usage}

\code{
volcanoplot(x,
            labels      = "",
            p.value     = "pval",
            mask        = FALSE,
            show.cutoff = TRUE,
            cex.text    = 0.7,
            col.text    = "blue",
            col.cutoff  = "grey",
            xlim        = NULL,
            xlab        = "Log2(Fold-Change)",
            ylab        = "-Log10(P-Value)",
            pch         = '.',
            ...) 
}
}
\arguments{
  \item{x}{object of class \code{\link{AnalysisTreeSet}}.}
  \item{labels}{optional transcript labels to be drawn at plotting points.}
  \item{p.value}{type of p-value, 'pval' for p-value, 'padj' for adjusted p-value, or 'pcha' for p-chance.}
  \item{mask}{logical, if \code{TRUE} draw only points for transcripts satisfying the univariate test.}
  \item{show.cutoff}{logical, if \code{TRUE} draw lines indicating cutoff.}
  \item{cex.text}{magnification to be used for optional \code{labels}.}
  \item{col.text}{color to be used for optional \code{labels}.}
  \item{col.cutoff}{color to be used for lines indicating cutoff, if \code{show.cutoff=TRUE}.}
  \item{xlim}{optional range for the plotted fold-change values.}
  \item{xlab}{label of x-axis.}
  \item{ylab}{label of y-axis.}
  \item{pch}{either an integer specifying a symbol or a single character to be used as the default in plotting points.}
  \item{\dots}{optional arguments to be passed to \code{plot}.}
}
\details{
 Produces a volcano plot for slot \code{data} for an object of class \code{\link{AnalysisTreeSet}}.

 It is possible to label the points of the volcano plot, whereby the following
 \code{labels} parameters are valid:
 \tabular{lll}{
  \tab \code{fUnitName}:\tab unit name (probeset ID). \cr
  \tab \code{fName}:\tab gene name. \cr
  \tab \code{fSymbol}:\tab gene symbol. \cr
  \tab \code{fChromosome}:\tab chromosome. \cr
  \tab \code{fCytoBand}:\tab cytoband.
 }
}
\author{Christian Stratowa}
\keyword{methods}
