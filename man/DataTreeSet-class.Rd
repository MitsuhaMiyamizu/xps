\name{DataTreeSet-class}
\docType{class}
\alias{DataTreeSet-class}
\alias{DataTreeSet}
\alias{addData,DataTreeSet-method}
\alias{attachBgrd,DataTreeSet-method}
\alias{attachDataXY,DataTreeSet-method}
\alias{attachInten,DataTreeSet-method}
\alias{attachMask,DataTreeSet-method}
\alias{attachProbeContentGC,DataTreeSet-method}
\alias{attachUnitNames,DataTreeSet-method}
\alias{background,DataTreeSet-method}
\alias{background}
\alias{background<-,DataTreeSet,data.frame-method}
\alias{background<-}
\alias{bgtreeNames,DataTreeSet-method}
\alias{bgtreeNames}
\alias{hist,DataTreeSet-method}
\alias{indexUnits,DataTreeSet-method}
\alias{intensity,DataTreeSet-method}
\alias{intensity<-,DataTreeSet,data.frame-method}
\alias{intensity2GCplot,DataTreeSet-method}
\alias{mm,DataTreeSet-method}
\alias{mmindex,DataTreeSet-method}
\alias{ncols,DataTreeSet-method}
\alias{nrows,DataTreeSet-method}
\alias{pm,DataTreeSet-method}
\alias{pmindex,DataTreeSet-method}
\alias{pmplot,DataTreeSet-method}
\alias{probesetID2unitID,DataTreeSet-method}
\alias{probesetplot,DataTreeSet-method}
\alias{projectInfo,DataTreeSet-method}
\alias{projectInfo<-,DataTreeSet,ProjectInfo-method}
\alias{rawCELName,DataTreeSet-method}
\alias{removeBgrd,DataTreeSet-method}
\alias{removeDataXY,DataTreeSet-method}
\alias{removeInten,DataTreeSet-method}
\alias{removeMask,DataTreeSet-method}
\alias{removeProbeContentGC,DataTreeSet-method}
\alias{removeUnitNames,DataTreeSet-method}
\alias{symbol2unitID,DataTreeSet-method}
\alias{transcriptID2unitID,DataTreeSet-method}
\alias{unitID2probesetID,DataTreeSet-method}
\alias{unitID2symbol,DataTreeSet-method}
\alias{unitID2transcriptID,DataTreeSet-method}
\alias{validBgrd,DataTreeSet-method}
\alias{validBgrd}
\alias{validData,DataTreeSet-method}
\alias{xpsBgCorrect,DataTreeSet-method}
\alias{xpsDABGCall,DataTreeSet-method}
\alias{xpsFIRMA,DataTreeSet-method}
\alias{xpsINICall,DataTreeSet-method}
\alias{xpsMAS4,DataTreeSet-method}
\alias{xpsMAS5,DataTreeSet-method}
\alias{xpsMAS5Call,DataTreeSet-method}
\alias{xpsNormalize,DataTreeSet-method}
\alias{xpsPreprocess,DataTreeSet-method}
\alias{xpsQualify,DataTreeSet-method}
\alias{xpsQualityControl,DataTreeSet-method}
\alias{xpsRMA,DataTreeSet-method}
\alias{xpsSummarize,DataTreeSet-method}

\title{Class DataTreeSet}
\description{This class provides the link to the \code{\link{ROOT}} data file and 
the \code{\link{ROOT}} trees contained therein. It extends class \code{\link{ProcesSet}}.}
\section{Objects from the Class}{
Objects can be created using the functions \code{\link{import.data}} or \code{\link{root.data}}.
}
\section{Slots}{
	 \describe{
    \item{\code{bgtreenames}:}{Object of class \code{"list"} representing the names of  
     optional \code{\link{ROOT}} background trees.}
    \item{\code{bgrd}:}{Object of class \code{"data.frame"}. The data.frame can contain 
     background intensities stored in \code{\link{ROOT}} background trees.}
    \item{\code{projectinfo}:}{Object of class \code{"ProjectInfo"} containing information 
     about the project.}
    \item{\code{scheme}:}{Object of class \code{"SchemeTreeSet"} providing access to 
     \code{\link{ROOT}} scheme file.}
    \item{\code{data}:}{Object of class \code{"data.frame"}. The data.frame can contain 
     the data (e.g. intensities) stored in \code{\link{ROOT}} data trees.}
    \item{\code{params}:}{Object of class \code{"list"} representing relevant parameters.}
    \item{\code{setname}:}{Object of class \code{"character"} representing the name to the
     \code{\link{ROOT}} file subdirectoy where the \code{\link{ROOT}} data trees are stored,
     usually \sQuote{DataTreeSet}.}
    \item{\code{settype}:}{Object of class \code{"character"} describing the type of treeset
     stored in \code{setname}, usually \sQuote{rawdata}.}
    \item{\code{rootfile}:}{Object of class \code{"character"} representing the name 
     of the \code{\link{ROOT}} data file, including full path.}
    \item{\code{filedir}:}{Object of class \code{"character"} describing the full path to
     the system directory where \code{rootfile} is stored.}
    \item{\code{numtrees}:}{Object of class \code{"numeric"} representing the number of 
     \code{\link{ROOT}} trees stored in subdirectoy \code{setname}.}
    \item{\code{treenames}:}{Object of class \code{"list"} representing the names of the 
     \code{\link{ROOT}} trees stored in subdirectoy \code{setname}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{ProcesSet}"}, directly.
Class \code{"\linkS4class{TreeSet}"}, by class "ProcesSet", distance 2.
}
\section{Methods}{
  \describe{
    \item{addData}{\code{signature(object = "DataTreeSet")}: import additional CEL-files and 
     update \code{\link{ROOT}} data file \code{rootfile}.}
    \item{attachBgrd}{\code{signature(object = "DataTreeSet")}: exports background trees from 
     \code{\link{ROOT}} data file and and saves as data.frame \code{bgrd}.}
    \item{attachDataXY}{\code{signature(object = "DataTreeSet")}: exports (x,y)-coordinates from 
     \code{\link{ROOT}} data file and and saves as data.frame \code{data}.}
    \item{attachInten}{\code{signature(object = "DataTreeSet")}: exports intensity trees from 
     \code{\link{ROOT}} data file and and saves as data.frame \code{data}.}
    \item{attachMask}{\code{signature(object = "DataTreeSet")}: exports scheme tree from 
     \code{\link{ROOT}} scheme file and and saves as data.frame \code{mask} of slot \code{scheme}.}
    \item{attachProbeContentGC}{\code{signature(object = "DataTreeSet")}: exports probe tree from 
     \code{\link{ROOT}} scheme file and and saves \code{fNumberGC} as data.frame \code{probe}.}
    \item{attachUnitNames}{\code{signature(object = "DataTreeSet")}: exports unit tree from 
     \code{\link{ROOT}} scheme file and and saves as data.frame \code{unitname} of slot \code{scheme}.}
    \item{background}{\code{signature(object = "DataTreeSet")}: extracts slot \code{bgrd}.}
    \item{background<-}{\code{signature(object = "DataTreeSet", value = "data.frame")}: replaces slot \code{bgrd}.}
    \item{bgtreeNames}{\code{signature(object = "DataTreeSet")}: extracts slot \code{bgtreenames}.}
    \item{indexUnits}{\code{signature(object = "DataTreeSet")}: extracts (x,y)-coordinates and corresponding
     indices for all or selected unitIDs.}
    \item{intensity}{\code{signature(object = "DataTreeSet")}: extracts slot \code{data}.}
    \item{intensity<-}{\code{signature(object = "DataTreeSet", value = "data.frame")}: replaces slot \code{data}.}
    \item{intensity2GCplot}{\code{signature(x = "DataTreeSet")}: creates a boxplot of probe intensities  
     stratified by GC content.}
    \item{mm}{\code{signature(object = "DataTreeSet")}: extracts the mismatch intensities.}
    \item{mmindex}{\code{signature(object = "DataTreeSet")}: extracts (x,y)-coordinates and corresponding
     MM indices for all or selected unitIDs.}
    \item{ncols}{\code{signature(object = "DataTreeSet")}: extracts the physical number of 
     array columns from slot \code{scheme}.}
    \item{nrows}{\code{signature(object = "DataTreeSet")}: extracts the physical number of 
     array rows from slot \code{scheme}.}
    \item{pm}{\code{signature(object = "DataTreeSet")}: extracts the perfect match intensities.}
    \item{pmindex}{\code{signature(object = "DataTreeSet")}: extracts (x,y)-coordinates and corresponding
     PM indices for all or selected unitIDs.}
    \item{pmplot}{\code{signature(x = "DataTreeSet")}: creates a barplot of mean perfect match and
     mismatch intensities.}
    \item{probesetID2unitID}{\code{signature(object = "DataTreeSet")}: extracts all or selected 
     probesetIDs from data.frame \code{unitname} of slot \code{scheme} with UnitName, i.e. probeset ID, as (row)names.}
    \item{probesetplot}{\code{signature(x = "DataTreeSet")}: creates a line plot of probe intensities for 
     a probeset.}
    \item{projectInfo}{\code{signature(object = "DataTreeSet")}: extracts slot \code{projectinfo}.}
    \item{projectInfo<-}{\code{signature(object = "DataTreeSet", value = "ProjectInfo")}: replaces slot \code{projectinfo}.}
    \item{rawCELName}{\code{signature(object = "DataTreeSet")}: returns the name(s) of the 
     imported raw CEL-files.}
    \item{removeBgrd}{\code{signature(object = "DataTreeSet")}: replaces data.frame \code{bgrd} 
     with an empty data.frame of dim(0,0).}
    \item{removeDataXY}{\code{signature(object = "DataTreeSet")}: replaces data.frame \code{data} 
     with an empty data.frame of dim(0,0).}
    \item{removeInten}{\code{signature(object = "DataTreeSet")}: replaces data.frame \code{data} 
     with an empty data.frame of dim(0,0).}
    \item{removeMask}{\code{signature(object = "DataTreeSet")}: replaces data.frame \code{mask} 
     from slot \code{scheme} with an empty data.frame of dim(0,0).}
    \item{removeProbeContentGC}{\code{signature(object = "DataTreeSet")}: replaces data.frame \code{probe} 
     with an empty data.frame of dim(0,0).}
    \item{removeUnitNames}{\code{signature(object = "DataTreeSet")}: replaces data.frame 
     \code{unitname} from slot \code{scheme} with an empty data.frame of dim(0,0).}
    \item{symbol2unitID}{\code{signature(object = "DataTreeSet")}: extracts internal UNIT_ID(s) for one
     or more gene symbols.}
    \item{transcriptID2unitID}{\code{signature(object = "DataTreeSet")}: extracts all or selected 
     transcriptIDs from data.frame \code{unitname} of slot \code{scheme} with UnitName, i.e. transcript ID, as (row)names.}
    \item{unitID2probesetID}{\code{signature(object = "DataTreeSet")}: extracts all or selected unitIDs
     from data.frame \code{unitname} of slot \code{scheme} with UNIT_ID as (row)names.}
    \item{symbol2unitID}{\code{signature(object = "DataTreeSet")}: extracts gene symbols for one
     or more internal UNIT_ID(s).}
    \item{unitID2transcriptID}{\code{signature(object = "DataTreeSet")}: extracts all or selected unitIDs
     from data.frame \code{unitname} of slot \code{scheme} with UNIT_ID as (row)names.}
    \item{validBgrd}{\code{signature(object = "DataTreeSet")}: extracts the valid data 
     from data.frame \code{bgrd}.}
    \item{validData}{\code{signature(object = "DataTreeSet")}: extracts a subset of valid data 
     from data.frame \code{data}.}
    \item{xpsBgCorrect}{\code{signature(object = "DataTreeSet")}: applies background correction methods.
     See \code{\link{bgcorrect}}.}
    \item{xpsDABGCall}{\code{signature(object = "DataTreeSet")}: computes DABG call.}
    \item{xpsFIRMA}{\code{signature(object = "DataTreeSet")}: computes FIRMA expression level and splice score.}
    \item{xpsINICall}{\code{signature(object = "DataTreeSet")}: computes I/NI call.}
    \item{xpsMAS4}{\code{signature(object = "DataTreeSet")}: computes MAS4 expression levels.}
    \item{xpsMAS5}{\code{signature(object = "DataTreeSet")}: computes MAS5 expression levels.}
    \item{xpsMAS5Call}{\code{signature(object = "DataTreeSet")}: computes MAS5 detection call.}
    \item{xpsNormalize}{\code{signature(object = "DataTreeSet")}: applies normalization methods.}
    \item{xpsPreprocess}{\code{signature(object = "DataTreeSet")}: applies normalization methods.}
    \item{xpsQualify}{\code{signature(object = "DataTreeSet")}: applies quality control methods.}
    \item{xpsQualityControl}{\code{signature(object = "DataTreeSet")}: applies quality control methods.}
    \item{xpsRMA}{\code{signature(object = "DataTreeSet")}: computes RMA expression levels.}
    \item{xpsSummarize}{\code{signature(object = "DataTreeSet")}: applies summarization methods.}
	 }
}
\author{Christian Stratowa}
\seealso{
related classes \code{\linkS4class{ExprTreeSet}}, \code{\linkS4class{CallTreeSet}}.
}
\examples{
showClass("DataTreeSet")
}
\keyword{classes}
