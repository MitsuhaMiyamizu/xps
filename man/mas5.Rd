\name{mas5}
\alias{mas5}
\alias{xpsMAS5-methods}
\alias{xpsMAS5}
\title{MAS 5.0 Expression Measure}
\description{
 This function converts a \code{\link{DataTreeSet}} into an \code{\link{ExprTreeSet}}
 using the XPS implementation of Affymetrix's MAS 5.0 expression measure.
}
\usage{
mas5(xps.data,
     filename   = character(0),
     filedir    = getwd(),
     tmpdir     = "",
     normalize  = FALSE,
     sc         = 500,
     option     = "transcript",
     exonlevel  = "",
     update     = FALSE,
     xps.scheme = NULL,
     add.data   = TRUE,
     verbose    = TRUE)


xpsMAS5(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{DataTreeSet}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{normalize}{logical. If \code{TRUE} scale normalization is used after an \code{ExprTreeSet} is obtained.}
  \item{sc}{value at which all arrays will be scaled to.}
  \item{option}{option determining the grouping of probes for summarization, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon arrays only.}
  \item{exonlevel}{exon annotation level determining which probes should be used for summarization;
   exon/genome arrays only.}
  \item{update}{logical. If \code{TRUE} the existing ROOT data file \code{filename} will be updated.}
  \item{xps.scheme}{optional alternative \code{SchemeTreeSet}.}
  \item{add.data}{logical. If \code{TRUE} expression data will be included as slot \code{data}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{DataTreeSet}.}
  \item{\dots}{arguments \code{filename},\code{filedir},\code{tmpdir},\code{option},\code{exonlevel},\code{xps.scheme}.}
}
\details{
This function computes the Affymetrix MAS 5.0 expression measure as implemented in XPS.
Although this implementation is based on the Affymetrix \sQuote{sadd_whitepaper.pdf}, it 
can be used to compute an expression level for both expression arrays and exon arrays. 
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
  \tab \code{ambiguous}:\tab ambiguous probesets only. \cr
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
  \code{exonlevel="metacore+affx"}:  \tab core meta-probesets plus AFFX controls \cr
  \code{exonlevel="core+extended"}: \tab probesets with cDNA support \cr
  \code{exonlevel="core+extended+full"}:  \tab supported plus predicted probesets \cr
}

Exon level annotations are described in the Affymetrix whitepaper 
\sQuote{exon_probeset_trans_clust_whitepaper.pdf}.

If \code{normalize=TRUE} then the expression levels will be scaled to \code{sc}. 
For \code{sc=0} the expression levels will be scaled to the mean expression level.

If \code{update=TRUE} then the existing \code{\link{ROOT}} file \code{filename} will be
updated, however, this is usually only recommended as option for function \code{\link{express}}.

In order to use an alternative \code{\link{SchemeTreeSet}} set the corresponding SchemeTreeSet 
\code{xps.scheme}.

\code{xpsMAS5} is the \code{DataTreeSet} method called by function \code{mas5}, however, 
expression levels will not be scaled to a common mean expression level.
}
\value{
 An \code{\link{ExprTreeSet}}
}
\author{Christian Stratowa}
\note{
In contrast to function \code{mas5}, expression levels computed with \code{xpsMAS5} 
will not be scaled to a common mean expression level.
}
\references{
Affymetrix (2002) Statistical Algorithms Description Document, Affymetrix Inc., 
Santa Clara, CA, whitepaper.
\url{http://www.affymetrix.com/support/technical/whitepapers/sadd_whitepaper.pdf}

Affymetrix (2005) Exon Probeset Annotations and Transcript Cluster Groupings, Affymetrix Inc., 
Santa Clara, CA, exon_probeset_trans_clust_whitepaper.pdf.
}
\seealso{\code{\link{express}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

data.mas5 <- mas5(data.test3,"tmp_Test3MAS5",tmpdir="",normalize=TRUE,sc=500,update=TRUE,verbose=FALSE)

## get data.frame
expr.mas5 <- validData(data.mas5)
head(expr.mas5)

## plot results
if (interactive()) {
boxplot(data.mas5)
boxplot(log2(expr.mas5))
}

rm(scheme.test3, data.test3)
gc()
}
\keyword{manip}
