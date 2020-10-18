\name{export}
\alias{export.scheme}
\alias{export.data}
\alias{export.expr}
\alias{export.call}
\alias{export-methods}
\alias{export}
\title{Export data as text files}
\description{
  Export data from classes \code{\link{SchemeTreeSet}}, \code{\link{DataTreeSet}},
 \code{\link{ExprTreeSet}}, or \code{\link{CallTreeSet}} to \code{outfile}.
}
\usage{
export.scheme(xps.scheme, treetype = character(0), varlist = "*", outfile = character(0), sep = "\t", as.dataframe = FALSE, verbose = TRUE)

export.data(xps.data, treename = "*", treetype = "cel", varlist = "*", outfile = character(0), sep = "\t", as.dataframe = FALSE, verbose = TRUE)

export.expr(xps.expr, treename = "*", treetype = character(0), varlist = "*", outfile = character(0), sep = "\t", as.dataframe = FALSE, verbose = TRUE)

export.call(xps.call, treename = "*", treetype = character(0), varlist = "*", outfile = character(0), sep = "\t", as.dataframe = FALSE, verbose = TRUE)

export(object, ...)
}
\arguments{
  \item{xps.scheme}{an object of type \code{\link{SchemeTreeSet}}.}
  \item{xps.data}{an object of type \code{\link{DataTreeSet}}.}
  \item{xps.expr}{an object of type \code{\link{ExprTreeSet}}.}
  \item{xps.call}{an object of type \code{\link{CallTreeSet}}.}
  \item{treename}{vector of tree names to export.}
  \item{treetype}{type of tree(s) to export, see \code{\link{validTreetype}}}
  \item{varlist}{names of tree leaves to export}
  \item{outfile}{name of output file.}
  \item{sep}{column separator}
  \item{as.dataframe}{if \code{TRUE} a data.frame will be returned.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{DataTreeSet}.}
  \item{\dots}{arguments \code{treenames},\code{treetype},\code{varlist},\code{outfile},\code{sep},\code{as.dataframe}.}
}
\details{
  Export data from classes \code{\link{SchemeTreeSet}}, \code{\link{DataTreeSet}},
 \code{\link{ExprTreeSet}}, or \code{\link{CallTreeSet}} to \code{outfile}.

 Parameter \code{varlist} lists the parameters to export: \cr
 - parameters are separated by ":", e.g. \code{varlist="fInten:fStdev"}. \cr
 - for \code{varlist="*"} all valid parameters will be exported.

 For class \code{DataTreeSet} the following \code{varlist} parameters are valid:
 \tabular{lll}{
  \tab \code{fInten}:\tab intensities from e.g. tree.cel. \cr
  \tab \code{fStdev}:\tab standard deviation from e.g. tree.cel. \cr
  \tab \code{fNPixels}:\tab number of pixels from e.g. tree.cel. \cr
  \tab \code{fBg}:\tab background values (background trees only).
 }
 For classes \code{ExprTreeSet} and \code{CallTreeSet} \code{varlist} can contain annotation
 parameters and parameters of the resulting data. \cr
  Following \code{varlist} annotation parameters are valid:
 \tabular{lll}{
  \tab \code{fUnitName}:\tab unit name (probeset ID). \cr
  \tab \code{fTranscriptID}:\tab transcript_id (probeset ID). \cr
  \tab \code{fName}:\tab gene name. \cr
  \tab \code{fSymbol}:\tab gene symbol. \cr
  \tab \code{fAccession}:\tab mRNA accession such as Refseq ID. \cr
  \tab \code{fEntrezID}:\tab entrez ID. \cr
  \tab \code{fChromosome}:\tab chromosome. \cr
  \tab \code{fStart}:\tab start position. \cr
  \tab \code{fStop}:\tab stop position. \cr
  \tab \code{fStrand}:\tab strand on chromosome. \cr
  \tab \code{fCytoBand}:\tab cytoband.
 }
  Following \code{varlist} parameters are valid for \code{ExprTreeSet}:
 \tabular{lll}{
  \tab \code{fLevel}:\tab expression level. \cr
  \tab \code{fStdev}:\tab standard deviation. \cr
  \tab \code{fNPairs}:\tab number of pairs.
 }
  Following \code{varlist} parameters are valid for \code{CallTreeSet}:
 \tabular{lll}{
  \tab \code{fCall}:\tab detection call. \cr
  \tab \code{fPValue}:\tab detection p-value.
 }
  Following \code{varlist} parameters are valid for \code{QualTreeSet}:
 \tabular{lll}{
  \tab \code{fLevel}:\tab expression level. \cr
  \tab \code{fStderr}:\tab standard error. \cr
  \tab \code{fNUSE}:\tab normalized unscaled standard error. \cr
  \tab \code{fRLE}:\tab relative log expression value.
 }
 An example: \code{varlist="fUnitName:fName:fSymbol:fLevel:fStdev:fEntrezID"}

 \code{export} is a generic method to export data from \code{\link{ROOT}} trees as text file.
}
\value{
If \code{as.dataframe} is \code{TRUE}, the data will be imported into the
current R session as \code{data.frame}. Otherwise, \code{NULL} will be returned.
}
\author{Christian Stratowa}
\seealso{\code{\link{export-methods}}}
\examples{
## load existing ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## export as table only
export(scheme.test3, treetype="idx", outfile="Test3_idx.txt",verbose=FALSE)

## export as table and import as data.frame
ann <- export.scheme(scheme.test3, treetype="ann", outfile="Test3_ann.txt",as.dataframe=TRUE,verbose=FALSE)
head(ann)
data <- export.data(data.test3, outfile="Test3_cel.txt",as.dataframe=TRUE,verbose=FALSE)
head(data)
}
\keyword{manip}
