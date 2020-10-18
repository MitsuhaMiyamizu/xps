\name{validData-methods}
\docType{methods}
\alias{validData}
\title{Extract Subset of Data}
\description{
 Extracts a subset of valid data from data.frame \code{data}.

\emph{Usage}

\code{
validData(object, which = "", unitID = NULL, unittype = "transcript")
}
}
\arguments{
  \item{object}{object of class \code{DataTreeSet}, \code{ExprTreeSet} or \code{CallTreeSet}.}
  \item{which}{type of probes to be returned for \code{DataTreeSet}, otherwise name of column
   containing unit name.}
  \item{unitID}{optional vector of UNIT_IDs.}
  \item{unittype}{\code{character} vector, \dQuote{transcript} or \dQuote{probeset}.}
}
\details{
For class \code{DataTreeSet} and expression arrays, \code{validData} returns all the 
perfect match or mismatch probes on the arrays the object represents as data.frame, 
i.e. \code{which} can have the following values:
\tabular{lll}{
  \tab \code{pm}:\tab perfect match probes. \cr
  \tab \code{mm}:\tab mismatch probes. \cr
  \tab \code{both}:\tab both perfect match and mismatch probes. 
}

For class \code{DataTreeSet} and exon arrays, \code{validData} returns the probes of 
the different exon levels as data.frame, i.e. \code{which} can have one of the following values:
\tabular{lll}{
  \tab \code{core}:\tab probesets supported by RefSeq and full-length GenBank transcripts. \cr
  \tab \code{metacore}:\tab core meta-probesets. \cr
  \tab \code{extended}:\tab probesets with other cDNA support. \cr
  \tab \code{metaextended}:\tab extended meta-probesets. \cr
  \tab \code{full}:\tab probesets supported by gene predictions only. \cr
  \tab \code{metafull}:\tab full meta-probesets. \cr
  \tab \code{affx}:\tab standard AFFX controls. \cr
  \tab \code{all}:\tab combination of above. \cr
  \tab \code{genomic}:\tab genomic background probes. \cr
  \tab \code{antigenomic}:\tab antigenomic background probes.
}

For class \code{ExprTreeSet} \code{validData} returns the valid expression levels from 
slot \code{data} with unit names as row names, usually the probeset IDs stored in 
column \code{which="UnitName"}.

For class \code{CallTreeSet} \code{validData} returns the valid detection call p-values 
from slot \code{data} with unit names as row names, usually the probeset IDs stored in 
column \code{which="UnitName"}.
}
\value{
  A \code{\link{data.frame}}.
}
\author{Christian Stratowa}
\seealso{\code{\link{pm}}, \code{\link{mm}}, \code{\link{validExpr}}, \code{\link{validCall}}}
\keyword{methods}
