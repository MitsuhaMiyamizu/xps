\name{import.data}
\alias{import.data}
\title{Import CEL files into a DataTreeSet}
\description{
Import the Affymetrix CEL files into a ROOT file and create S4 class DataTreeSet
}
\usage{
import.data(xps.scheme,
            filename = character(0),
            filedir  = getwd(),
            celdir   = NULL,
            celfiles = "*",
            celnames = NULL,
            project  = NULL,
            verbose  = TRUE)
}
\arguments{
  \item{xps.scheme}{a \code{\link{SchemeTreeSet}} containing the correct scheme for the CEL-files}
  \item{filename}{file name of ROOT data file.}
  \item{filedir}{system directory where ROOT data file should be stored.}
  \item{celdir}{system directory containing the CEL-files for corresponding scheme.}
  \item{celfiles}{optional vector of CEL-files to be imported.}
  \item{celnames}{optional vector of names which should replace the CEL-file names.}
  \item{project}{optional class \code{\link{ProjectInfo}}.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
}
\details{
 \code{import.data} is used to import CEL-files from directory \code{celdir} into 
 a \code{\link{ROOT}} data file. To import only a subset of CEL-files, list these CEL-files 
 as vector \code{celfiles}.

 To import CEL-files from different directories, vector \code{celfiles} must contain the
 full path for each CEL-file and \code{celdir} must be \code{celdir=NULL}. 

 The optional parameter \code{celnames} allows you to replace the original CEL-file names with
 names of your choice, otherwise the names of the CEL-files will be used as \code{celnames}.  

 Currently, the following types of Affymetrix CEL-files can be imported: 
 text (version 3), xml, binary (xda), generic (agcc,calvin)

 An S4 class \code{\link{DataTreeSet}} will be created, serving as R wrapper to the 
 \code{\link{ROOT}} data file \code{filename}.

 Use function \code{\link{root.data}} to access the \code{\link{ROOT}} data file from new
 R sessions to avoid creating a new \code{\link{ROOT}} data file for every session.
}
\value{
  A \code{DataTreeSet} object.
}
\author{Christian Stratowa}
\note{As mentioned above, use function \code{\link{root.data}} to access the 
 ROOT data file from new R sessions to avoid creating a new \code{\link{ROOT}} data file 
 for every R session.

 Do not separate \code{filename} of ROOT files with dots, use underscores, e.g. do not use 
 \code{filename="Data.Test3"} but use \code{filename="Data_Test3"} or \code{filename="DataTest3"}
 instead.

 To every ROOT data file the extension \dQuote{\_cel} is attached to \code{filename} to easily
 recognize ROOT data files containing the raw CEL data, e.g. for \code{filename="Data_Test3"}
 the final name is \dQuote{Data\_Test3\_cel.root}. Extension \dQuote{root} is added automatically,
 so that ROOT is able to recognize the file as ROOT file.

 Once a \code{\link{ROOT}} file is created it can not be overwritten, it must be deleted
 manually first. Only \code{\link{ROOT}} files called \dQuote{tmp} or with \code{filename} 
 starting with \dQuote{tmp\_} will be re-created automatically.

 If CEL-file names contain dots, colons, parenthesis, etc. as characters, these characters will be
 replaced by underscores. It is recommended to use parameter \code{celnames} to create shorter CEL
 names and to replace special characters. 
}
\seealso{\code{\link{root.data}}, \code{\link{DataTreeSet}}}
\examples{
## get scheme and import CEL-files from package
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- import.data(scheme.test3,"tmp_data_test3",celdir=paste(path.package("xps"),"raw",sep="/"))
unlist(treeNames(data.test3))

## import only subset of CEL-files
subdata.test3 <- import.data(scheme.test3,"tmpdt_data_test3",celdir=paste(path.package("xps"),"raw",sep="/"),
                 celfiles=c("TestA1.CEL","TestB2.CEL"),verbose=FALSE)
unlist(treeNames(subdata.test3))
}
\keyword{manip}
