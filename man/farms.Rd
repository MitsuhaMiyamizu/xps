\name{farms}
\alias{farms}
\title{Factor Analysis for Robust Microarray Summarization Expression Measure}
\description{
 This function converts a \code{\link{DataTreeSet}} into an \code{\link{ExprTreeSet}}
 using the Factor Analysis for Robust Microarray Summarization (FARMS) method.
}
\usage{
farms(xps.data,
      filename   = character(0),
      filedir    = getwd(),
      tmpdir     = "",
      normalize  = TRUE,
      weight     = 0.5,
      mu         = 0.0,
      scale      = 1.0,
      tol        = 0.00001,
      cyc        = 100,
      weighted   = TRUE,
      version    = "1.3.1",
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
  \item{weight}{hyperparameter, usually set to 0.5 for \code{version="1.3.1"} and to 8.0
   for \code{version="1.3.0"}.}
  \item{mu}{hyperparameter allowing to correct for potential bias.}
  \item{scale}{scaling parameter, usually set to 1.0 for \code{version="1.3.1"} and to 2.0
   for \code{version="1.3.0"}.}
  \item{tol}{termination tolerance for EM algorithm.}
  \item{cyc}{maximum number of cycles of EM algorithm.}
  \item{weighted}{logical, used only with \code{version="1.3.1"}. Default is TRUE.}
  \item{version}{version of original farms package. Currently, \code{version="1.3.1"} and
   \code{version="1.3.0"} are implemented. Default is \code{version="1.3.1"}.}
  \item{option}{option determining the grouping of probes for summarization, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon arrays only.}
  \item{exonlevel}{exon annotation level determining which probes should be used for summarization;
   exon/genome arrays only.}
  \item{xps.scheme}{optional alternative \code{SchemeTreeSet}.}
  \item{add.data}{logical. If \code{TRUE} expression data will be included as slot \code{data}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
}
\details{
This function computes the FARMS (Factor Analysis for Robust Microarray Summarization) expression 
measure described in Hochreiter et al. for both expression arrays and exon arrays. 

Parameter \code{version} currently allows the user to choose between the original implementation 
of FARMS as implemented in package \sQuote{farms_1.3.0} or enhanced FARMS as implemented in 
package \sQuote{farms_1.3.1}. By default \code{version="1.3.1"} is used.

Parameter \code{weight} is a hyperparameter which determines the influence of the prior. For 
\code{version="1.3.1"} the value in the range of [0,1]. 

Parameter \code{mu} is a hyperparameter which allows to quantify different aspects of potential 
prior knowledge. Values near zero assume that most genes do not contain a signal and introduce 
a bias for loading matrix elements near zero. 

Parameter \code{weighted} is a logical and indicates whether a weighted mean or a least square 
fit is used to summarize the loading matrix. It is applicable only to \code{version="1.3.1"}. 

For exon arrays it is necessary to supply the requested \code{option} and \code{exonlevel}.

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
The expression measure obtained with FARMS is given in linear scale, analogously to the expression 
measures computed with \code{\link{mas5}} and \code{\link{rma}}. 

For the analysis of many exon arrays it may be better to define a \code{tmpdir}, since this will
store only the results in the main file and not e.g. background and normalized intensities, and
thus will reduce the file size of the main file. For quantile normalization memory should not
be an issue, however DFW depends on RAM unless you are using a temporary file.
}
\references{
 Hochreiter, S., Clevert D.-A., and Obermayer, K. (2006), A new summarization method for 
 Affymetrix probe level data. Bioinformatics 22(8):943-949
}
\seealso{\code{\link{express}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

data.farms <- farms(data.test3,"tmp_Test3FARMS",verbose=FALSE)

## get data.frame
expr.farms <- validData(data.farms)
head(expr.farms)
}
\keyword{manip}
