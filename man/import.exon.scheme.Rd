\name{import.exon.scheme}
\alias{import.exon.scheme}
\title{Import CLF, PGF and annotation files into a SchemeTreeSet}
\description{
Import the Affymetrix CLF, PGF, and probeset and transcript annotation files into a ROOT file and create S4 class SchemeTreeSet
}
\usage{
import.exon.scheme(filename   = character(0),
                   filedir    = getwd(),
                   layoutfile = character(0),
                   schemefile = character(0),
                   probeset   = character(0),
                   transcript = character(0),
                   control    = "",
                   add.mask   = FALSE,
                   verbose    = TRUE)
}
\arguments{
  \item{filename}{file name of ROOT scheme file.}
  \item{filedir}{system directory where ROOT scheme file should be stored.}
  \item{layoutfile}{name of CLF-file, including full path.}
  \item{schemefile}{name of PGF-file, including full path.}
  \item{probeset}{name of probeset annotation-file, including full path.}
  \item{transcript}{name of transcript annotation-file, including full path.}
  \item{control}{optional name of controls.ps-file, including full path.}
  \item{add.mask}{logical. If \code{TRUE} mask information will be included as slot \code{mask}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
}
\details{
 \code{import.exon.scheme} is used to import all information for an Affymetrix 
 exon array into a \code{\link{ROOT}} scheme file, including CLF and PGF-files, and the current 
 Afymetrix probeset and transcript annotation files. 

 An S4 class \code{\link{SchemeTreeSet}} will be created, serving as R wrapper to the 
 \code{\link{ROOT}} scheme file \code{filename}.

 Since a new \code{\link{ROOT}} scheme file needs only to be created when new annotation files
 are available from the Affymetrix website, it is recommended to store all \code{\link{ROOT}}
 scheme files in a commonly accessible system directory \code{filedir}.

 Use function \code{\link{root.scheme}} to access the \code{\link{ROOT}} scheme file from new
 R sessions to avoid creating a new \code{\link{ROOT}} scheme file for every session.
}
\value{
  A \code{SchemeTreeSet} object.
}
\author{Christian Stratowa}
\note{As mentioned above, use function \code{\link{root.scheme}} to access the 
 \code{\link{ROOT}} scheme file from new R sessions to avoid creating a new \code{\link{ROOT}} 
 scheme file for every R session.

 Do not separate \code{filename} of ROOT files with dots, use underscores, e.g. do not use 
 \code{filename="Scheme.HuEx10stv2r2.na32"} but use \code{filename="Scheme_HuEx10stv2r2_na32"}
 instead. Extension \dQuote{root} is added automatically, so that ROOT is able to recognize
 the file as ROOT file.

 Do not set \code{add.mask=TRUE} unless you know that your computer has sufficient RAM.

 Do not add item \code{control} unless you want to use one of the old annotation
 files where the probeset annotation file does not contain the AFFX controls.
}
\section{Warning }{
 The current version of \sQuote{xps} should be able to import all Affymetrix exon array
 annotation files up to September 2011. However, since Affymetrix is still changing the headers
 and/or columns of the annotation files, future annotation files may require adaptation
 of the source code, thus the current version of \sQuote{xps} may not be able to read
 these files.
}
\seealso{\code{\link{import.expr.scheme}}, \code{\link{root.scheme}}, \code{\link{SchemeTreeSet}}}
\examples{
\dontrun{
## define paths
scmdir <- "/common/path/schemes"
libdir <- "/my/path/Affy/libraryfiles"
anndir <- "/my/path/Affy/Annotation"

## create scheme for HuEx-1_0-st-v2.r2 Exon array
scheme.huex10stv2r2.na32 <- import.exon.scheme("Scheme_HuEx10stv2r2_na32", filedir=scmdir,
                            layoutfile=file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf"),
                            schemefile=file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf"),
                            probeset=file.path(anndir, "HuEx-1_0-st-v2.na32.hg19.probeset.csv"),
                            transcript=file.path(anndir, "HuEx-1_0-st-v2.na32.hg19.transcript.csv"))

## access ROOT scheme file from new R session
scheme.exon <- root.scheme(paste(scmdir,"Scheme_HuEx10stv2r2_na32.root",sep="/"))


## create scheme for HuGene-1_0-st-v1.r4 as exon array
scheme.hugene10stv1r4.na32 <- import.exon.scheme("Scheme_HuGene10stv1r4_na32",filedir=scmdir,
                              layoutfile=file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.clf"),
                              schemefile=file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.pgf"),
                              probeset=file.path(anndir, "HuGene-1_0-st-v1.na32.hg19.probeset.csv"),
                              transcript=file.path(anndir, "HuGene-1_0-st-v1.na32.hg19.transcript.csv"))

## access ROOT scheme file from new R session
scheme.gene <- root.scheme(file.path(scmdir, "Scheme_HuGene10stv1r4_na32.root"))


## create scheme for HuEx-1_0-st-v2.r2 Exon array with the old annotation file
scheme.huex10stv2r2.old <- import.exon.scheme("Scheme_HuEx10stv2r2_old",filedir=scmdir,
                           layoutfile=file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.clf"),
                           schemefile=file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.pgf"),
                           probeset=file.path(anndir, "HuEx-1_0-st-probeset-annot.csv"),
                           transcript=file.path(anndir, "HuEx-1_0-st-transcript-annot.csv"),
                           control=file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.controls.ps"))
}
}
\keyword{manip}
