\name{dabg.call}
\alias{dabg.call}
\alias{xpsDABGCall-methods}
\alias{xpsDABGCall}
\title{Detection Above Background Call}
\description{
Computes the Detection Above Background Call first implemented for the Exon arrays.
}
\usage{
dabg.call(xps.data, filename = character(0), filedir = getwd(),
          alpha1 = 0.04, alpha2 = 0.06,
          option = "transcript", exonlevel = "", xps.scheme = NULL, add.data = TRUE, verbose = TRUE)

xpsDABGCall(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{DataTreeSet}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{alpha1}{a significance threshold in (0,alpha2).}
  \item{alpha2}{a significance threshold in (alpha1,0.5).}
  \item{option}{option determining the grouping of probes for summarization, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon arrays only.}
  \item{exonlevel}{exon annotation level determining which probes should be used for summarization; exon/genome arrays only.}
  \item{xps.scheme}{optional alternative \code{SchemeTreeSet}.}
  \item{add.data}{logical. If \code{TRUE} call data will be added to slots \code{data} and \code{detcall}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{DataTreeSet}.}
  \item{\dots}{the arguments described above.}
}
\details{
 This function generates a detection p-value based on comparing the perfect match probe
 intensity to the intensity distribution provided by background probes sharing the same
 GC-content as the PM probe under consideration. For exon/genome arrays special \sQuote{antigenomic}
 background probes of defined GC-content are used, while for expression arrays the Mismatch
 probes will be grouped by their GC-content. 

For exon/genome arrays it is necessary to supply \code{option} and \code{exonlevel}.

Following \code{option}s are valid for exon arrays only:
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
  \tab \code{ambiguous}:\tab ambiguous probesets only. \cr
  \tab \code{affx}:\tab standard AFFX controls. \cr
  \tab \code{all}:\tab combination of above.
}
Following \code{exonlevel} annotations are valid for whole genome arrays:
\tabular{lll}{
  \tab \code{core}:\tab probesets with category 'unique' and 'mixed'. \cr
  \tab \code{metacore}:\tab probesets with category 'unique' only. \cr
  \tab \code{affx}:\tab standard AFFX controls. \cr
  \tab \code{all}:\tab combination of above.
}
Exon levels can also be combined, with following combinations being most useful:
\tabular{ll}{
  \code{exonlevel="metacore+affx"}:  \tab core meta-probesets plus AFFX controls \cr
  \code{exonlevel="core+extended"}: \tab probesets with cDNA support \cr
  \code{exonlevel="core+extended+full"}:  \tab supported plus predicted probesets \cr
}

 Exon level annotations are described in the Affymetrix whitepaper
 'exon_probeset_trans_clust_whitepaper.pdf'.

 In order to use an alternative \code{\link{SchemeTreeSet}} set the corresponding
 SchemeTreeSet \code{xps.scheme}.

 \code{xpsDABGCall} is the \code{DataTreeSet} method called by function \code{dabg.call}, 
 containing the same parameters.
}
\value{
 A \code{\link{CallTreeSet}}
}
\author{Christian Stratowa}
\note{
 Yes, it is possible to compute DABG detection call for expression arrays, but it is
 very slow and thus not recommended.
}
\references{
Affymetrix (2005) Exon Probeset Annotations and Transcript Cluster Groupings, Affymetrix Inc., 
Santa Clara, CA, exon_probeset_trans_clust_whitepaper.pdf.
}
\seealso{\code{\link{mas5.call}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## DABG detection call
call.dabg <- dabg.call(data.test3,"tmp_Test3DABG",verbose=FALSE)

## get data.frames
pval.dabg <- pvalData(call.dabg)
pres.dabg <- presCall(call.dabg)
head(pval.dabg)
head(pres.dabg)

## plot results
if (interactive()) {
callplot(call.dabg)
}

rm(scheme.test3, data.test3)
gc()
}
\keyword{manip}
