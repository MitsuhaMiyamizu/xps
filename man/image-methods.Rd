\name{image-methods}
\docType{methods}
\alias{image-methods}
\alias{image}
\title{Display an Image}
\description{
 Creates an image of intensities or residuals, respectively, for each sample.

\emph{Usage}

\code{image(x,
            bg        = FALSE,
            transfo    = log2,
            col        = NULL,
            names      = "namepart",
            xlab       = "",
            ylab       = "",
            add.legend = FALSE,
            ...)}

\code{image(x,
            type       = c("resids", "pos.resids", "neg.resids", "sign.resids", "weights"),
            qualopt    = c("raw", "adjusted", "normalized"),
            transfo    = log2,
            col        = NULL,
            names      = "namepart",
            xlab       = "",
            ylab       = "",
            add.legend = FALSE,
            ...)}
}
\arguments{
  \item{x}{object of class \code{\link{ProcesSet}}.}
  \item{bg}{logical. If \code{FALSE}, intensities from slot \code{data} will be used;
   if \code{TRUE}, background intensities from slot \code{bgrd} will be used.}
  \item{type}{character string specifying the type of residual image.}
  \item{qualopt}{character string specifying whether to draw residual image for \dQuote{raw},
  \dQuote{adjusted}, or \dQuote{normalized} intensities.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{col}{color range for intensities.}
  \item{names}{optional vector of sample names.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{	a label for the y axis.}
  \item{add.legend}{logical, if TRUE then a color bar will be drawn.}
  \item{\dots}{optional arguments to be passed to \code{image}.}
}
\details{
 Creates an image of intensities or residuals, respectively, for each array, i.e.
 \sQuote{pseudo chip images}.

 If \code{x} belongs to class \code{DataTreeSet} then images of raw intensities will be drawn.

 If \code{x} belongs to class \code{ExprTreeSet} and \code{bg=FALSE} then images of background
 corrected intensities will be drawn.

 If \code{x} belongs to class \code{ExprTreeSet} and \code{bg=TRUE} the distribution of the
 background intensities will be shown; this can be useful to see potential density gradients
 caused by hybridization conditions. For the computation of background intensities see function
 \code{\link{bgcorrect}}; it is suggested to use \code{\link{bgcorrect.mas4}} to identify
 density gradients.

 If \code{x} belongs to class \code{QualTreeSet} then images of the residuals or the probe weights,
 respectively, will be drawn. For \code{col=NULL} the same colors will be used as described in 
 vignette \dQuote{QualityAssess.pdf} of package \code{affyPLM}, using internally function
 \code{pseudoPalette} described in \code{affyPLM}.

 For \code{names=NULL} full tree names will be displayed while for  \code{names="namepart"} column
 names will be displayed without name extension. If \code{names} is a vector of tree names then
 data from these trees only will displayed as image(s).
}
\author{Christian Stratowa}
\seealso{\code{\link{plotImage}}}
\examples{
\dontrun{
## images of raw intensities as imported using import.data()
unlist(treeNames(data.test3)) # show available tree names
image(data.test3, names="TestA2.cel")
image(data.test3)

## images of background adjusted or background intensities, created by e.g. rma()
getTreeNames(rootFile(data.rma))
image(data.rma, names="TestA2.int")
image(data.rma, names="TestA2.rbg", bg=TRUE)

## residual images, created by e.g. rmaPLM()
getTreeNames(rootFile(rlm.all), treetype="res")
image(rlm.all, type="resids")
image(rlm.all, type="resids", names="TestA2_raw.res", add.legend=TRUE)
image(rlm.all, type="pos.resids", names="TestA2_raw.res", add.legend=TRUE)
image(rlm.all, type="neg.resids", names="TestA2_raw.res", add.legend=TRUE)
image(rlm.all, type="sign.resids", names="TestA2_raw.res", add.legend=TRUE)
image(rlm.all, type="weights", names="TestA2_raw.res", add.legend=TRUE)
image(rlm.all, type="resids", qualopt="adjusted", names="TestA2_adjusted.res", add.legend=TRUE)
}
}
\keyword{methods}
