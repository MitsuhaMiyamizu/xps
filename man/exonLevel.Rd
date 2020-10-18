\name{exonLevel}
\alias{exonLevel}
\title{Conversion of Parameter exonlevel to Integer}
\description{
 Conversion of parameter \code{exonlevel} to an integer vector.
}
\usage{
exonLevel(exonlevel = "", chiptype = "GeneChip", as.sum = TRUE)
}
\arguments{
  \item{exonlevel}{exon annotation level determining which probes should be used for summarization;
   exon/genome arrays only.}
  \item{chiptype}{chip tpye, one of \sQuote{GeneChip}, \sQuote{GenomeChip}, \sQuote{ExonChip}.}
  \item{as.sum}{logical, if \code{TRUE} an integer vector of size three will be returned, if 
   \code{FALSE} then the levels will be split into the basic integer representations.}
}
\details{
 Conversion of parameter \code{exonlevel} to an integer; this function is a utility function,
 which is usually only used internally.

Following \code{exonlevel} annotations are valid for exon arrays:
\tabular{lll}{
  \tab \code{core}:\tab \code{(=8192+1024)} probesets supported by RefSeq and full-length GenBank transcripts. \cr
  \tab \code{metacore}:\tab \code{(=8192)} core meta-probesets. \cr
  \tab \code{extended}:\tab \code{(=4096+512)} probesets with other cDNA support. \cr
  \tab \code{metaextended}:\tab \code{(=4096)} extended meta-probesets. \cr
  \tab \code{full}:\tab \code{(=2048+256)} probesets supported by gene predictions only. \cr
  \tab \code{metafull}:\tab \code{(=2048)} full meta-probesets. \cr
  \tab \code{ambiguous}:\tab \code{(=128)} probesets that fall within multiple genes. \cr
  \tab \code{affx}:\tab \code{(=60)} standard AFFX controls. \cr
  \tab \code{all}:\tab \code{(=16316)} combination of above (including affx).
}
Following \code{exonlevel} annotations are valid for whole genome arrays:
\tabular{lll}{
  \tab \code{core}:\tab \code{(=8192+1024)} probesets with category 'unique', 'similar' and 'mixed'. \cr
  \tab \code{metacore}:\tab \code{(=8192)} probesets with category \sQuote{unique} only. \cr
  \tab \code{affx}:\tab \code{(=60)} standard AFFX controls. \cr
  \tab \code{all}:\tab \code{(=9276)} combination of above (including affx).
}
Exon levels can also be combined, with following combinations being most useful:
\tabular{ll}{
  \code{exonlevel="metacore+affx"}:  \tab core meta-probesets plus AFFX controls \cr
  \code{exonlevel="core+extended"}: \tab probesets with cDNA support \cr
  \code{exonlevel="core+extended+full"}:  \tab supported plus predicted probesets \cr
}

Exon level annotations are described in the Affymetrix whitepaper exon_probeset_trans_clust_whitepaper.pdf: \cr 
\dQuote{Exon Probeset Annotations and Transcript Cluster Groupings}.

Parameter \code{exonlevel} determines not only which probes are used for medianpolish, but also
the probes used for background calculation and for quantile normalization. If you want to use
seperate probes for background calculation, quantile normalization and medianpolish summarization,
you can pass a numeric vector containing three integer values corresponding to the respective 
\code{exonlevel}. These integers must be the sum of the integers shown above, e.g. you can use 
\code{exonlevel=c(16316,8252,8252)}, where \code{8252=8192+60} for \code{"metacore+affx"}.
}
\value{
an \code{integer} vector.
}
\author{Christian Stratowa}
\note{
The following \code{exonlevel}s are unsupported:
\tabular{lll}{
  \tab \code{control->bgp->genomic}:\tab \code{(=32768)} genomic background probes. \cr
  \tab \code{control->bgp->antigenomic}:\tab \code{(=65536)} antigenomic background probes. \cr
  \tab \code{normgene->intron}:\tab \code{(=131072)} intronic controls. \cr
  \tab \code{normgene->exon}:\tab \code{(=262144)} exronic controls. \cr
  \tab \code{rescue->FLmRNA->unmapped}:\tab \code{(=524288)} unmapped mRNAs. \cr
}
For whole genome arrays it is possible (but not recommended) to use all probesets 
by using \code{exonlevel=c(992316,992316,992316)}. \cr
For exon arrays it is possible to use e.g. \code{exonlevel=c(1032124,1032124,631868)}. \cr
However, please note that these settings are not recommended and not supported.
}
\seealso{\code{\link{rma}}, \code{\link{mas5}}}
\examples{
exonLevel("core", "GenomeChip")
exonLevel("all", "GenomeChip")
exonLevel("core+extended+full", "ExonChip")
exonLevel("core+extended+full", "ExonChip", as.sum=FALSE)
exonLevel(c(16316,8252,8252), "ExonChip")
}
\keyword{manip}
