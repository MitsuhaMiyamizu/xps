\name{dfw}
\alias{dfw}
\title{Distribution Free Weighted Expression Measure}
\description{
 This function converts a \code{\link{DataTreeSet}} into an \code{\link{ExprTreeSet}}
 using the Distribution Free Weighted Fold Change (DFW) method.
}
\usage{
dfw(xps.data,
    filename   = character(0),
    filedir    = getwd(),
    tmpdir     = "",
    normalize  = TRUE,
    m          = 3,
    n          = 1,
    c          = 0.01,
    option     = "transcript",
    exonlevel  = "",
    xps.scheme = NULL,
    add.data   = TRUE,
    verbose    = TRUE)
}
\arguments{
  \item{xps.data}{object of class \code{\link{DataTreeSet}}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{normalize}{logical. If \code{TRUE} normalize data using quantile normalization.}
  \item{m}{positive number as exponent of the weighted range WR.}
  \item{n}{positive number as exponent of the weighted standard deviation WSD.}
  \item{c}{scaling parameter.}
  \item{option}{option determining the grouping of probes for summarization, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon arrays only.}
  \item{exonlevel}{exon annotation level determining which probes should be used for summarization;
   exon/genome arrays only.}
  \item{xps.scheme}{optional alternative \code{SchemeTreeSet}.}
  \item{add.data}{logical. If \code{TRUE} expression data will be included as slot \code{data}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
}
\details{
This function computes the DFW (Distribution Free Weighted Fold Change) expression measure 
described in Chen et al. for both expression arrays and exon arrays. For exon arrays it is 
necessary to supply the requested \code{option} and \code{exonlevel}.

Following \code{option}s are valid for exon arrays:
\tabular{ll}{
  \code{transcript}: \tab expression levels are computed for transcript clusters,
   i.e. probe sets containing the same 'transcript_cluster_id'. \cr
  \code{exon}: \tab expression levels are computed for exon clusters,
   i.e. probe sets containing the same 'exon_id', where each exon cluster
   consists of one or more \code{probeset}s. \cr
  \code{probeset}: \tab expression levels are computed for individual probe sets,
   i.e. for each 'probeset_id'. \cr
}
Following \code{exonlevel} annotations are valid for exon arrays:
\tabular{lll}{
  \tab \code{core}:\tab probesets supported by RefSeq and full-length GenBank transcripts. \cr
  \tab \code{metacore}:\tab core meta-probesets. \cr
  \tab \code{extended}:\tab probesets with other cDNA support. \cr
  \tab \code{metaextended}:\tab extended meta-probesets. \cr
  \tab \code{full}:\tab probesets supported by gene predictions only. \cr
  \tab \code{metafull}:\tab full meta-probesets. \cr
  \tab \code{affx}:\tab standard AFFX controls. \cr
  \tab \code{all}:\tab combination of above (including affx).
}
Following \code{exonlevel} annotations are valid for whole genome arrays:
\tabular{lll}{
  \tab \code{core}:\tab probesets with category 'unique', 'similar' and 'mixed'. \cr
  \tab \code{metacore}:\tab probesets with category 'unique' only. \cr
  \tab \code{affx}:\tab standard AFFX controls. \cr
  \tab \code{all}:\tab combination of above (including affx).
}
Exon levels can also be combined, with following combinations being most useful:
\tabular{ll}{
  \code{exonlevel="metacore+affy"}:  \tab core meta-probesets plus AFFX controls \cr
  \code{exonlevel="core+extended"}: \tab probesets with cDNA support \cr
  \code{exonlevel="core+extended+full"}:  \tab supported plus predicted probesets \cr
}

Exon level annotations are described in the Affymetrix whitepaper exon_probeset_trans_clust_whitepaper.pdf: \cr 
\dQuote{Exon Probeset Annotations and Transcript Cluster Groupings}.

In order to use an alternative \code{\link{SchemeTreeSet}} set the corresponding SchemeSet \code{xps.scheme}.
}
\value{
 An \code{\link{ExprTreeSet}}
}
\author{Christian Stratowa}
\note{
The expression measure obtained with DFW is given in linear scale, analogously to the expression 
measures computed with \code{\link{mas5}} and \code{\link{rma}}. 

For the analysis of many exon arrays it may be better to define a \code{tmpdir}, since this will
store only the results in the main file and not e.g. background and normalized intensities, and
thus will reduce the file size of the main file. For quantile normalization memory should not
be an issue, however DFW depends on RAM unless you are using a temporary file.
}
\references{
 Chen, Z., McGee M., Liu Q., and Scheuermann, R.H. (2007), A distribution free summarization 
 method for Affymetrix GeneChip arrays. Bioinformatics 23(3):321-327
}
\seealso{\code{\link{express}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

data.dfw <- dfw(data.test3,"tmp_Test3DFW",verbose=FALSE)

## get data.frame
expr.dfw <- validData(data.dfw)
head(expr.dfw)
}
\keyword{manip}
