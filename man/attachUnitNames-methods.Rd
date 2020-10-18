\name{attachUnitNames-methods}
\docType{methods}
\alias{attachUnitNames-methods}
\alias{attachUnitNames}
\alias{removeUnitNames-methods}
\alias{removeUnitNames}
\title{Attach/Remove Unit Names}
\description{
 Attach/remove unit names, i.e. the Affymetrix probeset IDs to/from \code{\link{SchemeTreeSet}}
 or to slot \code{scheme} of \code{\link{DataTreeSet}}.

\emph{Usage}

\code{attachUnitNames(object, treetype = "idx")}

\code{removeUnitNames(object)}
}
\arguments{
  \item{object}{Object of class \code{"SchemeTreeSet"} or \code{"DataTreeSet"}.}
  \item{treetype}{the unit tree type, i.e. \sQuote{idx} or \sQuote{pbs}.}
}
\details{
 \code{attachUnitNames} exports \dQuote{UnitName} from unit tree of \code{\link{ROOT}} scheme file
 and and saves it as \code{data.frame} in slot \code{unitname}.

 \code{removeUnitNames} removes unitname from slot \code{unitname} and replaces data.frame \code{unitname}
 with an empty data.frame of dim(0,0).

 For \code{treetype="idx"} the internal \dQuote{UNIT_ID} will be mapped to the Affymetrix probeset IDs of
 the expression arrays or to the transcript_cluster_ids of the exon arrays, respectively, as \dQuote{UnitName}.

 For \code{treetype="pbs"} the internal \dQuote{UNIT_ID} will be mapped to the Affymetrix probeset_ids of
 the exon arrays as \dQuote{UnitName}.
}
\value{
  A \code{\link{DataTreeSet}} object or \code{\link{SchemeTreeSet}}.
}
\author{Christian Stratowa}
\note{
 Do not use \code{attachUnitNames} unless you know that your computer has sufficient RAM,
 especially for exon array schemes.
}
\seealso{\code{\link{attachMask}}, \code{\link{removeMask}}}
\examples{
## first, load ROOT scheme file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))

## attach unitname
scheme.test3 <- attachUnitNames(scheme.test3)

## get data.frame
unitnames <- unitNames(scheme.test3)
head(unitnames)

## remove unitname
scheme.test3 <- removeUnitNames(scheme.test3)

rm(scheme.test3)
gc()
}
\keyword{methods}
