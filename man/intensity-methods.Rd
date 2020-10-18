\name{intensity-methods}
\docType{methods}
\alias{intensity-methods}
\alias{intensity}
\alias{intensity<-}
\title{Get/Set Data Values}
\description{
  Get/set data values from/for class \code{DataTreeSet}.

\emph{Usage}

\code{intensity(object)} \cr
\code{intensity(object, filename = NULL, verbose = FALSE) <- value}
}
\arguments{
  \item{object}{object of class \code{\linkS4class{DataTreeSet}}.}
  \item{filename}{\code{character} vector containing optional ROOT file name.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{value}{\code{data.frame} containing expression values.}
}
\details{
 Get the intensity values from slot \code{data} or set slot \code{data} to \code{value}.

 Method \code{intensity} returns the data values from slot \code{data} as \code{data.frame},
 while replacement method \code{intensity<-} allows to replace slot \code{data} with
 a \code{data.frame}.

 Using replacement method \code{intensity<-} with default settings will not change the data
 stored in the ROOT data file, and thus will not have any effect on subsequent processing
 methods. If you really want to use the replacement data for further processing you must
 supply a new ROOT \code{filename}. This will export each intensity column of \code{value}
 as CEL-file (version 3), which will then be imported into the new ROOT data file \code{filename}.

 Warning: Do not use replacement method \code{intensity<-} until you really know what you
 are doing!

 Note: The first two columns of replacement \code{data.frame} \code{value} must be the
 (X,Y) coordinates, followed by the intensities whereby the number of intensity columns must
 be identical to the columns to be replaced.

 Note: If you do not want to replace your current object, create first a copy of type
 \code{DataTreeSet} by simply writing \code{newobj <- oldobj}, and use \code{newobj} for
 replacement. This is important since \code{intensity<-} does also update slots \code{rootfile},
 \code{filedir} and \code{treenames} when a new \code{filename} was chosen.

 Note: The CEL-files created are fully functional CEL-files (version 3), however some header
 rows such as GridCornerUL, AlgorithmParameters, and some of the data in DatHeader are 
 placeholders only.

 Warning: The CEL-files created WILL REPLACE THE ORIGINAL CEL-files, if they have identical names
 to the original CEL-files and the original CEL-files are located in the working directory.
 Thus the original CEL-files should preferable be located in directory \code{celdir} of function
 \code{\link{import.data}}.
}
\author{Christian Stratowa}
\seealso{\code{\link{validData}}}
\examples{
\dontrun{
## load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

## get intensity values
value <- intensity(data.test3)

## make a copy of your object if you do not want to replace it
newdata.test3 <- data.test3

## replace slot data with value
intensity(newdata.test3, "ReplacementData", FALSE) <- value
str(newdata.test3)

## now you can create an ExprTreeSet using the new intensity data
data.rma <- rma(newdata.test3,"ReplacementRMA",tmpdir="",background="none",normalize=TRUE,verbose=FALSE)
}
}
\keyword{methods}
