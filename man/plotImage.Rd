\name{plotImage}
\alias{plotImage}
\title{Plot Image(s) for Device}
\description{
 Creates an image for each sample for the selected device.
}
\usage{
plotImage(x,
          type     = character(),
          qualopt  = c("raw", "adjusted", "normalized"),
          transfo  = log2,
          col      = NULL,
          names    = character(),
          dev      = "screen",
          outfile  = "Image",
          w        = 800,
          h        = 800,
          verbose  = TRUE,
          ...) 
}
\arguments{
  \item{x}{object of class \code{\link{DataTreeSet}} or \code{\link{QualTreeSet}}.}
  \item{type}{character string specifying the type of image.}
  \item{qualopt}{character string specifying whether to draw residual image for \dQuote{raw},
  \dQuote{adjusted}, or \dQuote{normalized} intensities.}
  \item{transfo}{a valid function to transform the data, usually \dQuote{log2}, or \dQuote{0}.}
  \item{col}{color range for intensities.}
  \item{names}{vector of sample names.}
  \item{dev}{graphics device to plot to, i.e. one of \dQuote{screen}, \dQuote{jpeg},\dQuote{png},
   \dQuote{pdf} or \dQuote{ps}.}
  \item{outfile}{the name of the output file.}
  \item{w}{the width of the device in pixels.}
  \item{h}{the height of the device in pixels.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{\dots}{optional arguments to be passed to \code{image}.}
}
\details{
 Creates intensity image(s) or residual image(s), respectively, for each array for the
 selected graphics device, see \code{\link{image}} for more details.

 For intensity image(s) \code{type} must be one of \dQuote{intensity}.

 For residual image(s) \code{type} must be one of \dQuote{resids},\dQuote{pos.resids},
 \dQuote{neg.resids}, \dQuote{sign.resids}, or \dQuote{weights}. Furthermore, \code{qualopt}
 determines if images should be drawn for \dQuote{raw}, \dQuote{adjusted}, or \dQuote{normalized}
 data.

 For \code{names="*"} names of all samples will be displayed as images. If \code{names}
 is a vector of column names, only these samples will displayed as image(s).
}
\author{Christian Stratowa}
\seealso{\code{\link{image-methods}}, \code{\link{image}}}
\examples{
\dontrun{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## qualification - rlm
rlm.all <- rmaPLM(data.test3, "tmp_Test3RLMall", filedir=getwd(), tmpdir="", qualopt="all", option="transcript", add.data=FALSE)

if (interactive()) {
## image(s) of raw data
plotImage(data.test3, type="intensity", names="*")
plotImage(data.test3, type="intensity", names="TestA2.cel")

## image(s) of residuals/weights
plotImage(rlm.all, type="weights", names="*")
plotImage(rlm.all, type="weights", qualopt="adjusted", names="*")
plotImage(rlm.all, type="resids", names="TestA2_raw.res")
}

## function image.dev() will be deprecated since it needs attachInten!!
## need to attach scheme mask and data
data.test3 <- attachMask(data.test3)
data.test3 <- attachInten(data.test3)
if (interactive()) {
image.dev(data.test3)
}
## to avoid memory comsumption of R remove data:
data.test3 <- removeInten(data.test3)
data.test3 <- removeMask(data.test3)
}
}
\keyword{device}
