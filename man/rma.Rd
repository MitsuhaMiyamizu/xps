\name{rma}
\alias{rma}
\alias{xpsRMA-methods}
\alias{xpsRMA}
\title{Robust Multi-Array Average Expression Measure}
\description{
 This function converts a \code{\link{DataTreeSet}} into an \code{\link{ExprTreeSet}}
 using the robust multi-array average (RMA) expression measure.
}
\usage{
rma(xps.data,
    filename   = character(0),
    filedir    = getwd(),
    tmpdir     = "",
    background = "pmonly",
    normalize  = TRUE,
    option     = "transcript",
    exonlevel  = "",
    params     = list(16384, 0.0, 1.0, 10, 0.01, 1),
    xps.scheme = NULL,
    add.data   = TRUE,
    verbose    = TRUE)

xpsRMA(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{\link{DataTreeSet}}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{background}{probes used to compute background, one of \sQuote{pmonly}, \sQuote{mmonly},
   \sQuote{both}; for genome/exon arrays one of \sQuote{genomic}, \sQuote{antigenomic}}
  \item{normalize}{logical. If \code{TRUE} normalize data using quantile normalization.}
  \item{option}{option determining the grouping of probes for summarization, one of 
   \sQuote{transcript}, \sQuote{exon}, \sQuote{probeset}; exon arrays only.}
  \item{exonlevel}{exon annotation level determining which probes should be used for summarization;
   exon/genome arrays only.}
  \item{params}{list of (default) parameters for \code{rma}.}
  \item{xps.scheme}{optional alternative \code{SchemeTreeSet}.}
  \item{add.data}{logical. If \code{TRUE} expression data will be included as slot \code{data}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{DataTreeSet}.}
  \item{\dots}{the arguments described above.}
}
\details{
This function computes the RMA (Robust Multichip Average) expression measure described in 
Irizarry et al. for both expression arrays and exon arrays. For exon arrays it is necessary 
to supply the requested \code{option} and \code{exonlevel}.

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

Exon level annotations are described in the Affymetrix whitepaper exon_probeset_trans_clust_whitepaper.pdf: \cr 
\dQuote{Exon Probeset Annotations and Transcript Cluster Groupings}.

In order to use an alternative \code{\link{SchemeTreeSet}} set the corresponding SchemeSet \code{xps.scheme}.

\code{xpsRMA} is the \code{DataTreeSet} method called by function \code{rma}, 
containing the same parameters.
}
\value{
 An \code{\link{ExprTreeSet}}
}
\author{Christian Stratowa}
\note{
In contrary to other implementations of RMA the expression measure is given to you in linear scale,
analogously to the expression measures computed with \code{\link{mas5}} and \code{\link{mas4}}. 

Please note that the default settings of \code{params} gives results which are identical to the 
results obtained with APT (Affymetrix Power Tools) and with package \code{affy_1.14.2} or earlier.
If you want to obtain results which are identical to the results obtained with \code{affy_1.16.0} 
or later then you need to set \code{params = list(16384, 0.0, 0.4, 10, 0.01, 1)}. 

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
 Rafael. A. Irizarry, Benjamin M. Bolstad, Francois Collin, Leslie M. Cope, Bridget Hobbs
 and Terence P. Speed (2003), Summaries of Affymetrix GeneChip probe level data Nucleic
 Acids Research 31(4):e15

 Bolstad, B.M., Irizarry R. A., Astrand M., and Speed, T.P. (2003), A Comparison of
 Normalization Methods for High Density Oligonucleotide Array Data Based on Bias and
 Variance. Bioinformatics 19(2):185-193

 Irizarry, RA, Hobbs, B, Collin, F, Beazer-Barclay, YD, Antonellis, KJ, Scherf, U,
 Speed, TP (2003) Exploration, Normalization, and Summaries of High Density Oligonucleotide
 Array Probe Level Data. Biostatistics .Vol. 4, Number 2: 249-264 
}
\seealso{\code{\link{express}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

data.rma <- rma(data.test3,"tmp_Test3RMA",tmpdir="",background="pmonly",normalize=TRUE,verbose=FALSE)

## get data.frame
expr.rma <- validData(data.rma)
head(expr.rma)

## plot results
if (interactive()) {
boxplot(data.rma)
boxplot(log2(expr.rma))
}

rm(scheme.test3, data.test3)
gc()

\dontrun{
## examples using Affymetrix human tissue dataset (see also xps/examples/script4exon.R)
## first, load ROOT scheme file and ROOT data file from e.g.:
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"

## 1. example - expression array, e.g. HG-U133_Plus_2:
scheme.u133p2 <- root.scheme(paste(scmdir,"Scheme_HGU133p2_na25.root",sep="/"))
data.u133p2   <- root.data(scheme.u133p2, paste(datdir,"HuTissuesU133P2_cel.root",sep="/"))

workdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/u133p2"
data.rma <- rma(data.u133p2,"MixU133P2RMA",filedir=workdir,tmpdir="",
                background="pmonly",normalize=TRUE)

## 2. example - whole genome array, e.g. HuGene-1_0-st-v1:
scheme.genome <- root.scheme(paste(scmdir,"Scheme_HuGene10stv1r3_na25.root",sep="/"))
data.genome   <- root.data(scheme.genome, paste(datdir,"HuTissuesGenome_cel.root",sep="/"))

workdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/hugene"
data.g.rma <- rma(data.genome,"HuGeneMixRMAMetacore",filedir=workdir,tmpdir="",
                  background="antigenomic",normalize=T,exonlevel="metacore+affx")

## 3. example - exon array, e.g. HuEx-1_0-st-v2:
scheme.exon <- root.scheme(paste(scmdir,"Scheme_HuEx10stv2r2_na25.root",sep="/"))
data.exon   <- root.data(scheme.exon, paste(datdir,"HuTissuesExon_cel.root",sep="/"))

workdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/exon"
data.x.rma <- rma(data.exon,"MixRMAMetacore",filedir=workdir,tmpdir="",background="antigenomic",
                  normalize=T,option="transcript",exonlevel="metacore")
}
}
\keyword{manip}
