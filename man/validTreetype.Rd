\name{validTreetype}
\alias{validTreetype}
\title{Validate Tree Type}
\description{
  Validate tree type for corresponding data type.
}
\usage{
validTreetype(treetype, datatype)
}
\arguments{
  \item{treetype}{tree type.}
  \item{datatype}{data type.}
}
\details{
Every \code{\link{ROOT}} tree has an extension, which describes the type of data stored in this tree.
For example, \sQuote{TestA1.cel} is the tree name that stores the CEL-file data for \sQuote{TestA1.CEL}.

Trees with \code{datatype="scheme"} have the following extensions: \cr
  \code{scm}: scheme tree containing (x,y)-coordinates and mask for UNIT_ID. \cr
  \code{idx}: unit tree containing UnitName (i.e. probeset id), NumCells, NumAtoms, UnitType, for UNIT_ID. \cr
  \code{prb}: probe tree containing probe sequences. \cr
  \code{ann}: transcript annotation tree. \cr
  \code{anx}: exon annotation tree; exon arrays only. \cr
  \code{anp}: probeset annotation tree; exon arrays only. \cr
  \code{cxy}: coordinate tree containing CLF-file information; exon arrays only. \cr
  \code{exn}: exon tree; exon arrays only. \cr
  \code{pbs}: probeset tree; exon arrays only. 

Trees with \code{datatype="rawdata"} have the following extensions: \cr
  \code{cel}: data tree containing CEL-file data.

Trees with \code{datatype="preprocess"} have the following extensions: \cr
  \code{int}: intensity tree containing background-corrected intensities. \cr
  \code{sbg}: background tree containing MAS4 sector background levels. \cr
  \code{wbg}: background tree containing MAS5 weighted sector background levels. \cr
  \code{rbg}: background tree containing RMA background levels. \cr
  \code{gbg}: background tree containing GC-content background levels. \cr
  \code{cmn}: cell tree containing preprocessed intensities using algorithm \sQuote{mean}. \cr
  \code{cmd}: cell tree containing preprocessed intensities using algorithm \sQuote{median}. \cr
  \code{clw}: cell tree containing preprocessed intensities using algorithm \sQuote{lowess}. \cr
  \code{css}: cell tree containing preprocessed intensities using algorithm \sQuote{supsmu}. \cr
  \code{cqu}: cell tree containing preprocessed intensities using algorithm \sQuote{quantile}. \cr
  \code{dc5}: detection tree containing MAS5 detection call and p-value. \cr
  \code{dab}: detection tree containing DABG detection call and p-value. \cr
  \code{amn}: expression tree containing expression levels computed with \sQuote{arithmetic mean}. \cr
  \code{gmn}: expression tree containing expression levels computed with \sQuote{geometric mean}. \cr
  \code{wmn}: expression tree containing expression levels computed with \sQuote{weighted mean}. \cr
  \code{wdf}: expression tree containing expression levels computed with \sQuote{weighted difference}. \cr
  \code{adf}: expression tree containing expression levels computed with \sQuote{average difference}. \cr
  \code{tbw}: expression tree containing expression levels computed with \sQuote{tukey biweight}. \cr
  \code{mdp}: expression tree containing expression levels computed with \sQuote{median polish}. \cr
  \code{rlm}: quality tree containing expression levels, NUSE, RLE computed with \sQuote{median polish}. \cr
  \code{res}: residual tree containing the residual SE and the model fit weights. \cr
  \code{brd}: border tree containing border intensities, mean border intensities and COI.

Trees with \code{datatype="normation"} have the following extensions: \cr
  \code{tmn}: expression tree after normalization using algorithm \sQuote{trimmed mean}. \cr
  \code{med}: expression tree after normalization using algorithm \sQuote{median}. \cr
  \code{ksm}: expression tree after normalization using algorithm \sQuote{kernel smoother}. \cr
  \code{low}: expression tree after normalization using algorithm \sQuote{lowess}. \cr
  \code{sup}: expression tree after normalization using algorithm \sQuote{supsmu}. \cr
  \code{qua}: expression tree after normalization using algorithm \sQuote{quantile}. \cr
  \code{mdp}: expression tree after normalization using algorithm \sQuote{median polish}.
}
\value{
Returns the valid \code{treetype}, otherwise an error message is returned.
}
\author{Christian Stratowa}
\note{Not all tree types are used in the current package.
}
\seealso{\code{\link{getDatatype}}, \code{\link{type2Exten}}}
\examples{
validTreetype("prb", "scheme")
validTreetype("cel", "rawdata")
validTreetype("tbw", "preprocess")
}
\keyword{manip}

