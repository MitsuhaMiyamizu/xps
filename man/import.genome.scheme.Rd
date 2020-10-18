\name{import.genome.scheme}
\alias{import.genome.scheme}
\title{Import CLF, PGF and annotation files into a SchemeTreeSet}
\description{
Import the Affymetrix CLF, PGF and transcript annotation files into a ROOT file and create S4 class SchemeTreeSet
}
\usage{
import.genome.scheme(filename   = character(0),
                     filedir    = getwd(),
                     layoutfile = character(0),
                     schemefile = character(0),
                     transcript = character(0),
                     add.mask   = FALSE,
                     verbose    = TRUE)
}
\arguments{
  \item{filename}{file name of ROOT scheme file.}
  \item{filedir}{system directory where ROOT scheme file should be stored.}
  \item{layoutfile}{name of CLF-file, including full path.}
  \item{schemefile}{name of PGF-file, including full path.}
  \item{transcript}{name of transcript annotation-file, including full path.}
  \item{add.mask}{logical. If \code{TRUE} mask information will be included as slot \code{mask}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
}
\details{
 \code{import.genome.scheme} is used to import all information for an Affymetrix whole
 genome array into a \code{\link{ROOT}} scheme file, including CLF and PGF-files, and the  
 current Afymetrix  transcript annotation files. 

 An S4 class \code{\link{SchemeTreeSet}} will be created, serving as R wrapper to the 
 \code{\link{ROOT}} scheme file \code{filename}.

 Since a new \code{\link{ROOT}} scheme file needs only to be created when  new annotation files
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
 \code{filename="Scheme.HuGene10stv1.na27"} but use \code{filename="Scheme_HuGene10stv1_na27"}
 instead. Extension \dQuote{root} is added automatically, so that ROOT is able to recognize
 the file as ROOT file.

 Do not set \code{add.mask=TRUE} unless you know that your computer has sufficient RAM.

 Do not add item \code{control} unless you want to use one of the old annotation
 files where the probeset annotation file does not contain the AFFX controls.
}
\section{Warning }{
 The current version of \sQuote{xps} is able to import all Affymetrix genome array annotation
 files up to November 2008, i.e. all files of release 3 (r3) and earlier. However, in January 2009
 Affymetrix has updated all CLF, PGF and annotation files to release 4 (r4) and added a new
 probeset annotation file, thus in effect changing the whole genome arrays to exon arrays!

 Thus, for release 4 (r4) files, function \code{import.genome.scheme} can no longer be used,
 but you must use function \code{\link{import.exon.scheme}} instead (see examples). 
}
\seealso{\code{\link{import.exon.scheme}}, \code{\link{root.scheme}}, \code{\link{SchemeTreeSet}}}
\examples{
\dontrun{
## define paths
scmdir <- "/common/path/schemes"
libdir <- "/my/path/Affy/libraryfiles"
anndir <- "/my/path/Affy/Annotation"

## create scheme for HuGene-1_0-st-v1 whole genome array
scheme.hugene10stv1r3.na27 <- import.genome.scheme("Scheme_HuEx10stv1r3_na27", filedir=scmdir,
                              layoutfile=file.path(libdir, "HuGene-1_0-st-v1.r3.analysis_libraryfile", "HuGene-1_0-st-v1.r3.clf"),
                              schemefile=file.path(libdir, "HuGene-1_0-st-v1.r3.analysis_libraryfile", "HuGene-1_0-st-v1.r3.pgf"),
                              transcript=file.path(anndir, "HuGene-1_0-st-v1.na27.hg18.transcript.csv"))

## access ROOT scheme file from new R session
scheme.hugene10stv1r3 <- root.scheme(file.path(scmdir, "Scheme_HuEx10stv1r3_na27.root"))
}
}
\keyword{manip}
