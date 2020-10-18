\name{export.root}
\alias{export.root}
\title{Export data from ROOT file}
\description{
  Export data as text files directly from a \code{\link{ROOT}} file.
}
\usage{
export.root(datafile = character(0), schemefile = character(0), treeset = character(0), treename = "*", treetype = character(0), varlist = "*", outfile = character(0), sep = "\t", as.dataframe = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datafile}{name of ROOT data file including full path}
  \item{schemefile}{name of ROOT scheme file including full path}
  \item{treeset}{name of subdirectory in ROOT file where trees are stored}
  \item{treename}{name of ROOT tree to export.}
  \item{treetype}{type of tree(s) to export, see \code{\link{validTreetype}}.}
  \item{varlist}{names of tree leaves to export.}
  \item{outfile}{name of output file.}
  \item{sep}{column separator}
  \item{as.dataframe}{if \code{TRUE} a data.frame will be returned.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
}
\details{
 Export data as text files directly from a \code{\link{ROOT}} file.
}
\value{
If \code{as.dataframe} is \code{TRUE}, the data will be imported into the
current R session as \code{data.frame}. Otherwise, \code{NULL} will be returned.
}
\author{Christian Stratowa}
\seealso{\code{\link{export}}, \code{\link{export-methods}}}
\examples{
## export data directly from root file
schemefile <- paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/")
datafile   <- paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/")
data <- export.root(datafile, schemefile, "DataSet", "*", "cel", "*", "DataOutFile.txt", as.dataframe = TRUE,verbose=FALSE)
head(data)
}
\keyword{manip}
