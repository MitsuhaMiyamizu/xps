\name{ini.call}
\alias{ini.call}
\alias{xpsINICall-methods}
\alias{xpsINICall}
\title{Informative/Non-Informative Call}
\description{
Computes the Informative/Non-Informative Call for the exclusion of non-informative probe sets.
}
\usage{
ini.call(xps.data,
         filename   = character(0),
         filedir    = getwd(),
         tmpdir     = "",
         weight     = 0.5,
         mu         = 0.0,
         scale      = 1.0,
         tol        = 0.00001,
         cyc        = 100,
         alpha1     = 0.4,
         alpha2     = 0.6,
         version    = "1.3.1",
         option     = "transcript",
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)

xpsINICall(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{DataTreeSet}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{weight}{hyperparameter, usually set to 0.5 for \code{version="1.3.1"} and to 8.0
   for \code{version="1.3.0"}.}
  \item{mu}{hyperparameter allowing to correct for potential bias.}
  \item{scale}{scaling parameter, usually set to 1.0 for \code{version="1.3.1"} and to 2.0
   for \code{version="1.3.0"}.}
  \item{tol}{termination tolerance for EM algorithm.}
  \item{cyc}{maximum number of cycles of EM algorithm.}
  \item{alpha1}{a significance threshold in (0,alpha2).}
  \item{alpha2}{a significance threshold in (alpha1,1.0).}
  \item{version}{version of original farms package. Currently, \code{version="1.3.1"} and
   \code{version="1.3.0"} are implemented. Default is \code{version="1.3.1"}.}
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
 In contrast to \code{\link{mas5.call}} this function quantifies the signal-to-noise ratio for
 each probe set, as described in Talloen et al. Thus, the returned p-values and detection calls
 have a different meaning:

 The p-value that is returned estimates the signal-to-noise ratio (SNR): \cr
 P-values (SNR) of less than 0.5 indicate that there is more signal than noise and the
 corresponding genes are considered to be \sQuote{informative} for further analysis. In
 contrast, values greater than 0.5 indicate \sQuote{non-informative} genes.

 The informative call is computed by thresholding the p-value as in:

 call \dQuote{P} if p-value < alpha1 \cr
 call \dQuote{M} if alpha1 <= p-value < alpha2 \cr
 call \dQuote{A} if alpha2 <= p-value

 Here \dQuote{P} should be considered as informative \dQuote{I}, \dQuote{M} as marginally
 informative, and \dQuote{A} as non-informative \dQuote{NI}.

 The defaults for \code{alpha1=0.4} and \code{alpha2=0.6} are set to allow \dQuote{M} calls.
 In order to get the same results as package \sQuote{farms\_1.3.1}, you need to set
 \code{alpha1=0.5} and \code{alpha2=0.5}.

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
  \code{exonlevel="metacore+affy"}:  \tab core meta-probesets plus AFFX controls \cr
  \code{exonlevel="core+extended"}: \tab probesets with cDNA support \cr
  \code{exonlevel="core+extended+full"}:  \tab supported plus predicted probesets \cr
}

 Exon level annotations are described in the Affymetrix whitepaper
 'exon_probeset_trans_clust_whitepaper.pdf'.

 In order to use an alternative \code{\link{SchemeTreeSet}} set the corresponding
 SchemeTreeSet \code{xps.scheme}.

 \code{xpsINICall} is the \code{DataTreeSet} method called by function \code{ini.call}, 
 containing the same parameters.
}
\value{
 A \code{\link{CallTreeSet}}
}
\author{Christian Stratowa}
\note{
 Since I/NI-calls distinguish only between informative and non-informative genes, the
 calls are identical for all samples. 
}
\references{
 Talloen, W., Clevert D.-A., Hochreiter, S., Amaratunga, D., Bijnens, J., Kass, S., and
 Gohlmann, H.W.H. (2006), I/NI-calls for the exclusion of non-informative genes: a highly
 effective filtering tool for microarray data. Bioinformatics 23(21):2897-2902
}
\seealso{\code{\link{farms}}, \code{\link{mas5.call}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## I/NI call
call.ini <- ini.call(data.test3,"tmp_Test3INI",verbose=FALSE)

## get data.frames
snr.ini <- pvalData(call.ini)
inf.ini <- presCall(call.ini)
head(snr.ini)
head(inf.ini)

## plot results
if (interactive()) {
callplot(call.ini)
}

rm(scheme.test3, data.test3)
gc()
}
\keyword{manip}
