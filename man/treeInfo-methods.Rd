\name{treeInfo-methods}
\docType{methods}
\alias{treeInfo-methods}
\alias{treeInfo}
\title{Get UserInfo from ROOT Trees}
\description{
 Extract the UserInfo from \code{\link{ROOT}} trees, i.e. quality control information.

\emph{Usage}

\code{
 treeInfo(object,
          treename = "*",
          treetype = character(0),
          varlist  = "*",
          qualopt  = NULL,
          verbose  = FALSE,
          ...)
}
}
\arguments{
  \item{object}{Object of class \code{"TreeSet"}.}
  \item{treename}{Object of class \code{"list"} representing the names of the 
   \code{\link{ROOT}} trees.}
  \item{treetype}{type of tree to export, see \code{\link{validTreetype}}}
  \item{varlist}{names of tree leaves to export.}
  \item{qualopt}{option determining the data to which to apply qualification, one of 
   \sQuote{raw}, \sQuote{adjusted}, \sQuote{normalized}, \sQuote{all}.}
}
\details{
 \code{\link{ROOT}} trees have a pointer to a list fUserInfo where it is possible to store data
 which do not fit into the usual tree structure. Taking advantage of this feature \code{xps}
 stores certain pre-processed results of the tree(s) in this list. For example, data trees store
 the minimal/maximal intensities and the number of oligos with minimal/maximal intensities of
 the CEL-files in list fUserInfo, while call trees store the number and percentage of P/M/A
 calls.

 Function \code{treeInfo} allows to export this user information as a \code{data.frame},
 whereby the parameters of \code{varlist} depend on the \code{treetype}:

 Parameters for data trees with extensions \code{"cel"}, \code{"int"}, and background trees: \cr
  \code{fMinInten}: minimal intensity. \cr
  \code{fMaxInten}: maximal intensity. \cr
  \code{fNMinInten}: number of probes with minimal intensity. \cr
  \code{fNMaxInten}: number of probes with maximal intensity. \cr
  \code{fMaxNPixels}: maximal number of pixels. \cr
  \code{fNQuantiles}: number of precalculated quantiles. \cr
  \code{fQuantiles}: quantiles. \cr
  \code{fIntenQuant}: intensities at quantiles.

 Parameters for expression trees: \cr
  \code{fNUnits}: number of units, i.e. probesets. \cr
  \code{fMinLevel}: minimal expression level. \cr
  \code{fMaxLevel}: maximal expression level. \cr
  \code{fNQuantiles}: number of precalculated quantiles. \cr
  \code{fQuantiles}: quantiles. \cr
  \code{fLevelQuant}: expression levels at quantiles.

 Parameters for call trees: \cr
  \code{fNUnits}: number of units, i.e. probesets. \cr
  \code{fNAbsent}: number of units with absent call. \cr
  \code{fNMarginal}: number of units with marginal call. \cr
  \code{fNPresent}: number of units with present call. \cr
  \code{fPcAbsent}: percentage of units with absent call. \cr
  \code{fPcMarginal}: percentage of units with marginal call. \cr
  \code{fPcPresent}: percentage of units with present call. \cr
  \code{fMinPValue}: minimal p-value. \cr
  \code{fMaxPValue}: maximal p-value.

 Parameters for border trees with extension \code{"brd"}: \cr
  \code{fMeanLeft}: mean intensity of left border. \cr
  \code{fMeanRight}: mean intensity of right border. \cr
  \code{fMeanTop}: mean intensity of top border. \cr
  \code{fMeanBottom}: mean intensity of bottom border. \cr
  \code{fCOIXhi}: x-location of COI for the positive elements. \cr
  \code{fCOIYhi}: y-location of COI for the positive elements. \cr
  \code{fCOIXlo}: x-location of COI for the negative elements. \cr
  \code{fCOIYlo}: y-location of COI for the negative elements.

 Parameters for quality trees with extension \code{"rlm"}: \cr
  \code{fNUnits}: number of units, i.e. probesets. \cr
  \code{fMinLevel}: minimal expression level. \cr
  \code{fMaxLevel}: maximal expression level. \cr
  \code{fNQuantiles}: number of precalculated quantiles. \cr
  \code{fQuantiles}: quantiles. \cr
  \code{fLevelQuant}: expression levels at quantiles. \cr
  \code{fNUSEQuant}: NUSE at quantiles. \cr
  \code{fRLEQuant}: RLE at quantiles. \cr
  \code{fQualOption}: value of \code{qualopt}.

 Parameters for residual trees with extension \code{"res"}: \cr
  \code{fNQuantiles}: number of precalculated quantiles. \cr
  \code{fQuantiles}: quantiles. \cr
  \code{fResiduQuant}: residual at quantiles. \cr
  \code{fWeightQuant}: weight at quantiles. \cr
  \code{fQualOption}: value of \code{qualopt}. 
 }
\value{
  A \code{data.frame}.
}
\author{Christian Stratowa}
\note{
 Taking advantage of function \code{treeInfo} plotting methods \code{\link{boxplot}}, 
 \code{\link{callplot}}, \code{\link{coiplot}}, \code{\link{nuseplot}} and \code{\link{rleplot}}
 are able to display their results much faster, which is especially useful for large datasets.
}
\seealso{\code{\link{validTreetype}}}
\examples{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

userinfo <- treeInfo(data.test3, treetype="cel", varlist="*")
userinfo

userinfo <- treeInfo(data.test3, treename="TestB1", treetype="cel", varlist = "fNQuantiles:fIntenQuant")
userinfo

\dontrun{
userinfo <- treeInfo(rlm.all, treetype="rlm", varlist = "fNQuantiles:fNUSEQuant:fRLEQuant", qualopt = "raw")
userinfo

userinfo <- treeInfo(rlm.all, treetype="brd")
userinfo

userinfo <- treeInfo(rlm.all, treetype="res", qualopt = "raw")
userinfo

userinfo <- treeInfo(rlm.all, treetype="res", varlist = "fResiduQuant", qualopt = "raw")
userinfo
}
}
\keyword{methods}
