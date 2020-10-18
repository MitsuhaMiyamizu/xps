\name{SchemeTreeSet-class}
\docType{class}
\alias{SchemeTreeSet-class}
\alias{SchemeTreeSet}
\alias{attachMask,SchemeTreeSet-method}
\alias{attachProbe,SchemeTreeSet-method}
\alias{attachProbeContentGC,SchemeTreeSet-method}
\alias{attachProbeSequence,SchemeTreeSet-method}
\alias{attachUnitNames,SchemeTreeSet-method}
\alias{chipMask,SchemeTreeSet-method}
\alias{chipMask}
\alias{chipMask<-,SchemeTreeSet,data.frame-method}
\alias{chipMask<-}
\alias{chipName,SchemeTreeSet-method}
\alias{chipName}
\alias{chipProbe,SchemeTreeSet-method}
\alias{chipProbe}
\alias{chipProbe<-,SchemeTreeSet,data.frame-method}
\alias{chipProbe<-}
\alias{chipType,SchemeTreeSet-method}
\alias{chipType}
\alias{chipType<-,SchemeTreeSet,character-method}
\alias{chipType<-}
\alias{export,SchemeTreeSet-method}
\alias{ncols,SchemeTreeSet-method}
\alias{ncols}
\alias{nrows,SchemeTreeSet-method}
\alias{nrows}
\alias{probeContentGC,SchemeTreeSet-method}
\alias{probeInfo,SchemeTreeSet-method}
\alias{probeInfo}
\alias{probeSequence,SchemeTreeSet-method}
\alias{probesetID2unitID,SchemeTreeSet-method}
\alias{removeMask,SchemeTreeSet-method}
\alias{removeProbe,SchemeTreeSet-method}
\alias{removeProbeContentGC,SchemeTreeSet-method}
\alias{removeProbeSequence,SchemeTreeSet-method}
\alias{removeUnitNames,SchemeTreeSet-method}
\alias{symbol2unitID,SchemeTreeSet-method}
\alias{transcriptID2unitID,SchemeTreeSet-method}
\alias{unitID2probesetID,SchemeTreeSet-method}
\alias{unitID2symbol,SchemeTreeSet-method}
\alias{unitID2transcriptID,SchemeTreeSet-method}
\alias{unitNames,SchemeTreeSet-method}
\alias{unitNames}
\alias{unitNames<-,SchemeTreeSet,data.frame-method}
\alias{unitNames<-}

\title{Class SchemeTreeSet}
\description{This class provides the link to the \code{\link{ROOT}} scheme file and 
the \code{\link{ROOT}} trees contained therein. It extends class \code{\link{TreeSet}}.}
\section{Objects from the Class}{
Objects can be created using the functions \code{\link{import.expr.scheme}},
\code{\link{import.exon.scheme}}, \code{\link{import.genome.scheme}} or \code{\link{root.scheme}}.
}
\section{Slots}{
	 \describe{
    \item{\code{chipname}:}{Object of class \code{"character"} representing the Affymetrix
     chip name.}
    \item{\code{chiptype}:}{Object of class \code{"character"} representing the chip tpye, 
     either \sQuote{GeneChip}, \sQuote{GenomeChip} or \sQuote{ExonChip}.}
    \item{\code{probeinfo}:}{Object of class \code{"list"} representing chip information, 
     including nrows, ncols, number of probes, etc.}
    \item{\code{unitname}:}{Object of class \code{"data.frame"}. The data.frame can contain the  
     mapping between the internal \code{UNIT_ID}s and the \code{UnitName}s, i.e. the probeset IDs.}
    \item{\code{mask}:}{Object of class \code{"data.frame"}. The data.frame can contain the mask 
     used to identify the probes as e.g. PM, MM or control probes.}
    \item{\code{probe}:}{Object of class \code{"data.frame"}. The data.frame can contain the probe 
     info for the oligos as e.g. probe sequence, G/C content.}
    \item{\code{setname}:}{Object of class \code{"character"} representing the name to the
     \code{\link{ROOT}} file subdirectoy where the \code{\link{ROOT}} scheme trees are stored;
     it is identical to \code{chipname}.}
    \item{\code{settype}:}{Object of class \code{"character"} describing the type of treeset
     stored in \code{setname}, i.e. \sQuote{scheme}.}
    \item{\code{rootfile}:}{Object of class \code{"character"} representing the name 
     of the \code{\link{ROOT}} scheme file, including full path.}
    \item{\code{filedir}:}{Object of class \code{"character"} describing the full path to
     the system directory where \code{rootfile} is stored.}
    \item{\code{numtrees}:}{Object of class \code{"numeric"} representing the number of 
     \code{\link{ROOT}} trees stored in subdirectoy \code{setname}.}
    \item{\code{treenames}:}{Object of class \code{"list"} representing the names of the 
     \code{\link{ROOT}} trees stored in subdirectoy \code{setname}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{TreeSet}"}, directly.
}
\section{Methods}{
  \describe{
    \item{attachMask}{\code{signature(object = "SchemeTreeSet")}: exports scheme tree from 
     \code{\link{ROOT}} scheme file and and saves as data.frame \code{mask}.}
    \item{attachProbe}{\code{signature(object = "SchemeTreeSet")}: exports probe tree from 
     \code{\link{ROOT}} scheme file and and saves \code{varlist} as data.frame \code{probe}.}
    \item{attachProbeContentGC}{\code{signature(object = "SchemeTreeSet")}: exports probe tree from 
     \code{\link{ROOT}} scheme file and and saves \code{fNumberGC} as data.frame \code{probe}.}
    \item{attachProbeSequence}{\code{signature(object = "SchemeTreeSet")}: exports probe tree from 
     \code{\link{ROOT}} scheme file and and saves \code{fSequence} as data.frame \code{probe}.}
    \item{attachUnitNames}{\code{signature(object = "SchemeTreeSet")}: exports unit tree from 
     \code{\link{ROOT}} scheme file and and saves as data.frame \code{unitname}.}
    \item{chipMask}{\code{signature(object = "SchemeTreeSet")}: extracts data.frame \code{mask}.}
    \item{chipMask<-}{\code{signature(object = "SchemeTreeSet", value = "data.frame")}: replaces 
     data.frame \code{mask}.}
    \item{chipName}{\code{signature(object = "SchemeTreeSet")}: extracts slot \code{chipname}.}
    \item{chipProbe}{\code{signature(object = "SchemeTreeSet")}: extracts data.frame \code{probe}.}
    \item{chipProbe<-}{\code{signature(object = "SchemeTreeSet", value = "data.frame")}: replaces 
     data.frame \code{probe}.}
    \item{chipType}{\code{signature(object = "SchemeTreeSet")}: extracts slot \code{chiptype}.}
    \item{chipType<-}{\code{signature(object = "SchemeTreeSet", value = "character")}: replaces 
     slot \code{chiptype}.}
    \item{export}{\code{signature(object = "SchemeTreeSet")}: exports \code{\link{ROOT}} trees 
     as text file, see \code{\link{export-methods}}.}
    \item{ncols}{\code{signature(object = "SchemeTreeSet")}: extracts the physical number of 
     array columns from slot \code{probeinfo}.}
    \item{nrows}{\code{signature(object = "SchemeTreeSet")}: extracts the physical number of 
     array rows from slot \code{probeinfo}.}
    \item{probeContentGC}{\code{signature(object = "SchemeTreeSet")}: extracts all or selected 
     GC contents from data.frame \code{probe}.}
    \item{probeInfo}{\code{signature(object = "SchemeTreeSet")}: extracts slot \code{probeinfo}.}
    \item{probeSequence}{\code{signature(object = "SchemeTreeSet")}: extracts all or selected 
     probe sequences from data.frame \code{probe}.}
    \item{probesetID2unitID}{\code{signature(object = "SchemeTreeSet")}: extracts all or selected 
     probesetIDs from data.frame \code{unitname} with UnitName, i.e. probeset ID, as (row)names.}
    \item{removeMask}{\code{signature(object = "SchemeTreeSet")}: replaces data.frame \code{mask} 
     with an empty data.frame of dim(0,0).}
    \item{removeProbe}{\code{signature(object = "SchemeTreeSet")}: replaces data.frame \code{probe} 
     with an empty data.frame of dim(0,0).}
    \item{removeProbeContentGC}{\code{signature(object = "SchemeTreeSet")}: replaces data.frame \code{probe} 
     with an empty data.frame of dim(0,0).}
    \item{removeProbeSequence}{\code{signature(object = "SchemeTreeSet")}: replaces data.frame \code{probe} 
     with an empty data.frame of dim(0,0).}
    \item{removeUnitNames}{\code{signature(object = "SchemeTreeSet")}: replaces data.frame 
     \code{unitname} with an empty data.frame of dim(0,0).}
    \item{symbol2unitID}{\code{signature(object = "SchemeTreeSet")}: extracts internal UNIT_ID(s) for one
     or more gene symbols.}
    \item{transcriptID2unitID}{\code{signature(object = "SchemeTreeSet")}: extracts all or selected 
     transcriptIDs from data.frame \code{unitname} with UnitName, i.e. transcript ID, as (row)names.}
    \item{unitID2probesetID}{\code{signature(object = "SchemeTreeSet")}: extracts all or selected unitIDs
     from data.frame \code{unitname} with UNIT_ID as (row)names.}
    \item{symbol2unitID}{\code{signature(object = "SchemeTreeSet")}: extracts gene symbols for one
     or more internal UNIT_ID(s).}
    \item{unitID2transcriptID}{\code{signature(object = "SchemeTreeSet")}: extracts all or selected unitIDs
     from data.frame \code{unitname} with UNIT_ID as (row)names.}
    \item{unitNames}{\code{signature(object = "SchemeTreeSet")}: extracts data.frame \code{unitname}.}
    \item{unitNames<-}{\code{signature(object = "SchemeTreeSet", value = "data.frame")}: replaces 
     data.frame \code{unitname}.}
	 }
}
\author{Christian Stratowa}
\examples{
showClass("SchemeTreeSet")
}
\keyword{classes}
