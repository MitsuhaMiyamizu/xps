\name{bgcorrect}
\alias{bgcorrect}
\alias{bgcorrect.gc}
\alias{bgcorrect.mas4}
\alias{bgcorrect.mas5}
\alias{bgcorrect.rma}
\alias{xpsBgCorrect-methods}
\alias{xpsBgCorrect}
\title{Background Correction}
\description{
 Background corrects probe intensities in an object of class \code{\link{DataTreeSet}}.
}
\usage{
bgcorrect(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, select = "none", method = character(0), option = character(0), exonlevel = "", params = list(), verbose = TRUE)

bgcorrect.gc(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, select = "antigenomic", exonlevel = "", verbose = TRUE)

bgcorrect.mas4(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, select = "all", exonlevel = "", verbose = TRUE)

bgcorrect.mas5(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, select = "both", exonlevel = "", verbose = TRUE)

bgcorrect.rma(xps.data, filename = character(0), filedir = getwd(), tmpdir = "", update = FALSE, select = "none", exonlevel = "", verbose = TRUE)

xpsBgCorrect(object, ...)
}
\arguments{
  \item{xps.data}{object of class \code{DataTreeSet}.}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{tmpdir}{optional temporary directory where temporary ROOT files should be stored.}
  \item{update}{logical. If \code{TRUE} the existing ROOT data file \code{filename} will be updated.}
  \item{select}{type of probes to select for background correction.}
  \item{method}{background method to use.}
  \item{option}{type of background correction to use.}
  \item{exonlevel}{exon annotation level determining which probes should be used
   for summarization; exon/genome arrays only.}
  \item{params}{vector of parameters for background method.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{object}{object of class \code{DataSet}.}
  \item{\dots}{the arguments described above.}
}
\details{
 Background corrects probe intensities in an object of class \code{\link{DataTreeSet}}.

 \code{xpsBgCorrect} is the \code{DataSet} method called by function \code{bgcorrect},
 containing the same parameters.
}
\value{
 An \code{\link{DataTreeSet}}
}
\author{Christian Stratowa}
\seealso{\code{\link{express}}}
\examples{
## first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## MAS4 sector background
data.bg.mas4 <- bgcorrect.mas4(data.test3,"tmp_Test3MAS4Bgrd",filedir=getwd(),tmpdir="",verbose=FALSE)

## need to attach background intensities
data.bg.mas4 <- attachBgrd(data.bg.mas4)

## get data.frame
bg.mas4 <- validBgrd(data.bg.mas4)
head(bg.mas4)

## plot images
if (interactive()) {
image.dev(data.bg.mas4,bg=TRUE,col=rainbow(32))
image(matrix(bg.mas4[,1], ncol=ncols(schemeSet(data.bg.mas4)), nrow=nrows(schemeSet(data.bg.mas4))))
}

\dontrun{
## examples using Affymetrix human tissue dataset (see also xps/examples/script4exon.R)

## example - exon array, e.g. HuEx-1_0-st-v2:
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
scheme.exon <- root.scheme(paste(scmdir,"Scheme_HuEx10stv2r2_na25.root",sep="/"))
data.exon   <- root.data(scheme.exon, paste(datdir,"HuTissuesExon_cel.root",sep="/"))

## compute rma background
workdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/exon"
data.bg.rma <- bgcorrect(data.exon, "HuExonRMABgrd", filedir=workdir, tmpdir="", 
               method="rma", select="antigenomic", option="pmonly:epanechnikov",
               params=c(16384), exonlevel="metacore+affx")

# or alternatively:
data.bg.rma <- bgcorrect.rma(data.exon, "HuExonRMABgrd", filedir=workdir, tmpdir="", 
               select="antigenomic", exonlevel="metacore+affx")
}
}
\keyword{manip}
