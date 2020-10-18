\name{firma}
\alias{firma}
\alias{xpsFIRMA-methods}
\alias{xpsFIRMA}
\title{Finding Isoforms using Robust Multichip Analysis}
\description{
 This function converts a \code{\link{DataTreeSet}} for exon arrays into an \code{\link{ExprTreeSet}}
 using the Finding Isoforms using Robust Multichip Analysis (FIRMA).
}
\usage{
firma(xps.data,
      filename   = character(0),
      filedir    = getwd(),
      tmpdir     = "",
      background = "antigenomic",
      normalize  = TRUE,
      option     = "probeset",
      exonlevel  = "metacore",
      method     = "mdp",
      params     = list(16384, 0.0, 1.0, 10, 0.01, 1.0),
      xps.scheme = NULL,
      add.data   = TRUE,
      verbose    = TRUE)

xpsFIRMA(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{\link{DataTreeSet}}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{background}{probes used to compute background, one of \sQuote{genomic}, \sQuote{antigenomic}}
  \item{normalize}{logical. If \code{TRUE} normalize data using quantile normalization.}
  \item{option}{option determining the grouping of probes for summarization, one of 
   \sQuote{exon}, \sQuote{probeset}.}
  \item{exonlevel}{exon annotation level determining which probes should be used for summarization.}
  \item{method}{method to be used for summarization, currently \sQuote{mdp}.}
  \item{params}{list of (default) parameters for \code{rma}.}
  \item{xps.scheme}{optional alternative \code{SchemeTreeSet}.}
  \item{add.data}{logical. If \code{TRUE} expression data will be included as slot \code{data}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{DataTreeSet}.}
  \item{\dots}{the arguments described above.}
}
\details{
This function computes FIRMA (Finding Isoforms using Robust Multichip Analysis) for detecting 
differential alternative splicing for exon arrays, as described in Purdom et al. 

Following \code{option}s are valid for exon arrays:
\tabular{ll}{
  \code{probeset}: \tab expression levels are computed for individual probe sets,
   i.e. for each 'probeset_id'. \cr
  \code{exon}: \tab expression levels are computed for exon clusters,
   i.e. probe sets containing the same 'exon_id', where each exon cluster
   consists of one or more \code{probeset}s. \cr
}
Following \code{exonlevel} annotations are valid for exon arrays:
\tabular{lll}{
  \tab \code{core}:\tab probesets supported by RefSeq and full-length GenBank transcripts. \cr
  \tab \code{metacore}:\tab core meta-probesets. \cr
  \tab \code{extended}:\tab probesets with other cDNA support. \cr
  \tab \code{metaextended}:\tab extended meta-probesets. \cr
  \tab \code{full}:\tab probesets supported by gene predictions only. \cr
  \tab \code{metafull}:\tab full meta-probesets. \cr
  \tab \code{ambiguous}:\tab ambiguous probesets only. \cr
  \tab \code{affx}:\tab standard AFFX controls. \cr
  \tab \code{all}:\tab combination of above (including affx).
}
Exon levels can also be combined, with following combinations being most useful:
\tabular{ll}{
  \code{exonlevel="metacore+affx"}:  \tab core meta-probesets plus AFFX controls \cr
  \code{exonlevel="core+extended"}: \tab probesets with cDNA support \cr
  \code{exonlevel="core+extended+full"}:  \tab supported plus predicted probesets \cr
}

Exon level annotations are described in the Affymetrix whitepaper exon_probeset_trans_clust_whitepaper.pdf: \cr 
\dQuote{Exon Probeset Annotations and Transcript Cluster Groupings}.

Method \code{xpsFIRMA} is the \code{DataTreeSet} method called by function \code{firma}, 
containing the same parameters.
}
\value{
 An \code{\link{ExprTreeSet}}
}
\author{Christian Stratowa}
\note{
In contrary to other implementations of (FI)RMA the expression measure of FIRMA is given in linear scale,
analogously to the expression measures computed with \code{\link{mas5}} and \code{\link{mas4}}. 

Please note that the current implementation of FIRMA is based on median-polish only, see:
\url{http://www.aroma-project.org/node/81}

Please note that the default settings of \code{params} gives results which are identical to the 
results obtained with APT (Affymetrix Power Tools) and with package \code{affy_1.14.2} or earlier.
If you want to obtain results which are identical to the results obtained with \code{affy_1.16.0} 
or later then you need to set \code{params = list(16384, 0.0, 0.4, 10, 0.01, 1.0)}. 

By setting parameter \code{background="none"} it is possible to skip background correction .

For the analysis of many exon arrays it may be better to define a \code{tmpdir}, since this will
store only the results in the main file and not e.g. background and normalized intensities, and
thus will reduce the file size of the main file. For quantile normalization memory should not
be an issue, however medianpolish depends on RAM unless you are using a temporary file.

Parameter \code{exonlevel} determines not only which probes are used for medianpolish, but also
the probes used for background calculation and for quantile normalization. If you want to use
seperate probes for background calculation, quantile normalization and medianpolish summarization,
you can pass a numeric vector containing three integer values corresponding to the respective 
\code{exonlevel}, e.g. you can use \code{exonlevel=c(16316,8252,8252)}, see function 
\code{\link{exonLevel}} for more details.
}
\references{
 Purdom, E., Simpson K.M., Robinson M.D., Conboy J.G., Lapuk A.V. and Speed, T.P. (2008),
 FIRMA: a method for detection of alternative splicing from exon array data.
 Bioinformatics 24(15):1707-1714
}
\examples{
\dontrun{
## load ROOT scheme file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
scheme.exon <- root.scheme(paste(scmdir,"Scheme_HuEx10stv2r2_na27.root",sep="/"))

## load subset of ROOT data file
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
subnames <- c("HeartA","HeartB","HeartC", "MuscleA","MuscleB","MuscleC")
sub.exon <- root.data(scheme.exon, rootFile(data.exon), celnames=subnames)

## firma
outdir <- getwd()
sub.firma.ps <- firma(sub.exon,"HeartMuscleFIRMAcorePS",filedir=outdir,tmpdir="",background="antigenomic",
                      normalize=TRUE,option="probeset",exonlevel="core")

## get transcript expression levels for all transcripts or transcript=2429277
expr.firma <- firma.expr(sub.firma.ps, probeset=NULL, option="transcript")
expr.firma <- firma.expr(sub.firma.ps, probeset=2429277, option="transcript")

## get probeset expression levels for all probeset or probeset=2429278 or transcript=2429277
expr.firma <- firma.expr(sub.firma.ps, probeset=NULL, option="probeset")
expr.firma <- firma.expr(sub.firma.ps, probeset=2429278, option="probeset")
expr.firma <- firma.expr(sub.firma.ps, probeset=2429277, option="probeset")

## get probeset splice scores for all probeset or probeset=2429278 or transcript=2429277
score.firma <- firma.score(sub.firma.ps, probeset=NULL, option="probeset")
score.firma <- firma.score(sub.firma.ps, probeset=2429278, option="probeset")
score.firma <- firma.score(sub.firma.ps, probeset=2429277, option="probeset")

## different plots
boxplot(sub.firma.ps, which="UnitName:LEVEL_PS")
boxplot(sub.firma.ps, which="UnitName:LEVEL_TS")

hist(sub.firma.ps, which="UnitName:LEVEL_PS")
hist(sub.firma.ps, which="UnitName:LEVEL_TS")

rleplot(sub.firma.ps, which="UnitName:LEVEL_PS")
rleplot(sub.firma.ps, which="UnitName:LEVEL_TS")

nuseplot(sub.firma.ps, which="UnitName:STDEV_PS")
nuseplot(sub.firma.ps, which="UnitName:STDEV_TS")
}
}
\keyword{manip}
