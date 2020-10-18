\name{coiplot-methods}
\docType{methods}
\alias{coiplot-methods}
\alias{coiplot}
\title{Center-Of-Intensity QC Plots}
\description{
 Produce Center-Of-Intensity plot(s) of the positive and negative feature intensities.

\emph{Usage}

\code{
coiplot(x,
        type    = c("pos", "neg"),
        qualopt = "raw",
        radius  = 0.5,
        linecol = "gray70",
        visible = TRUE,
       ...)
}
}
\arguments{
  \item{x}{object of class \code{\link{QualTreeSet}}.}
  \item{type}{type of border elements to be used, one of \dQuote{pos}, \dQuote{neg}, or both.}
  \item{qualopt}{character string specifying whether to draw boxplots for \dQuote{raw},
  \dQuote{adjusted}, or \dQuote{normalized} border intensities.}
  \item{radius}{determines the radius within which the COI for each array should be located.}
  \item{linecol}{the color of the ablines and the circle to be drawn.}
  \item{visible}{logical, if TRUE then arrays outside the circle with \code{radius} will be
   flagged by labeling the data point with the array name.}
  \item{\dots}{optional arguments to be passed to \code{coiplot}.}
}
\details{
 Produces Center-Of-Intensity (COI) plot(s) of the positive and negative feature intensities for
 an object of class \code{\link{QualTreeSet}}. This plot is useful for detecting spatial biases
 in intensities on an array. 

 Mean intensities for the left, right, top and bottom border elements are calculated, separated
 into positive and negative controls, and the \dQuote{center of intensity} is calculated on a
 relative scale [-1,1]. Arrays with a COI outside a range with \code{radius} are considered to 
 be outliers. If \code{visible = TRUE} then outlier arrays will be flagged by labeling the data
 point(s) with the array name(s). 
}
\value{
 The names ot the outlier arrays, otherwise \code{NULL}.
}
\author{Christian Stratowa}
\seealso{\code{\link{plotCOI}}, \code{\link{borderplot}}}
\examples{
\dontrun{
## border intensities, created by e.g. rmaPLM()
coiplot(rlm.all)
coiplot(rlm.all, type="pos")
coiplot(rlm.all, type="neg", radius=0.1)
}
}
\keyword{methods}
