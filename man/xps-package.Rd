\name{xps-package}
\alias{xps-package}
\alias{xps}
\docType{package}
\title{xps Package Overview}
\description{
   xps Package Overview
}
\details{
  Important data classes: \code{\link{SchemeTreeSet}},
  \code{\link{DataTreeSet}},
  \code{\link{ExprTreeSet}},
  \code{\link{CallTreeSet}},
  \code{\link{FilterTreeSet}},
  \code{\link{AnalysisTreeSet}}.
  Full help on methods and associated functions is available from
  within class help pages.

  Additional data classes: \code{\link{ProjectInfo}},
  \code{\link{PreFilter}}, \code{\link{UniFilter}}.

  The package handles pre-processing, normalization, filtering and analysis
  of Affymetrix GeneChip expression arrays, including exon array systems 
  (Exon 1.0 ST: core, extended, full probesets), gene array systems (Gene 1.0 ST) 
  and plate array systems on computers with 1 GB RAM only. It imports Affymetrix 
  .CDF, .CLF, .PGF and .CEL as well as Affymetrix annotation files, and computes
  e.g. RMA, MAS5, FARMS, DFW, MAS5-calls, DABG-calls, I/NI-calls. It is an R wrapper to 
  XPS (eXpression Profiling System), which is based on ROOT, an object-oriented 
  framework developed at CERN. Thus, the prior installation of ROOT is a prerequisite
  for the usage of this package, see the README file. However, no knowledge of ROOT
  is required. ROOT is licensed under LGPL and can be downloaded from \url{http://root.cern.ch}.
}
\author{
   Christian Stratowa \email{cstrato@aon.at}
}
\keyword{package}
