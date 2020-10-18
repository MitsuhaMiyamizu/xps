\name{xpsQAReport}
\alias{xpsQAReport}
\title{Create Quality Assessment Report.}
\description{
 Create a quality assessment report.
}
\usage{
xpsQAReport(xps.data,
            xps.expr    = NULL,
            xps.call    = NULL,
            xps.qual    = NULL,
            dataset     = character(0), 
            title       = "Quality Report",
            date        = "October, 2011",
            author      = "Christian Stratowa",
            outdir      = file.path(getwd(), "QAReport"),
            add.pseudo  = FALSE,
            overwrite   = FALSE,
            verbose     = TRUE,
            ...) 
}
\arguments{
  \item{xps.data}{object of class \code{\link{DataTreeSet}}.}
  \item{xps.expr}{object of class \code{\link{ExprTreeSet}}.}
  \item{xps.call}{object of class \code{\link{CallTreeSet}}.}
  \item{xps.qual}{object of class \code{\link{QualTreeSet}}.}
  \item{dataset}{name of the dataset.}
  \item{title}{title of quality report.}
  \item{date}{date of quality report.}
  \item{author}{author(s) of quality report.}
  \item{outdir}{name of directory where to create the quality report.}
  \item{add.pseudo}{logical, if \code{TRUE} add pseudo-images to the quality report.}
  \item{overwrite}{logical, if \code{TRUE} overwrite \code{outdir} and its contents.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
  \item{\dots}{optional arguments to be passed to \code{xpsQAReport}.}
}
\details{
 Function \code{xpsQAReport} creates a quality assessment report "QAReport.pdf" for all \code{TreeSet}s,
 which are passed as parameters to the function. It calls \code{library(tools)} and uses its function
 \code{buildVignettes} to create the report.  

 If parameter \code{xps.qual} is supplied, it is possible to create pseudo-images for every CEL-file 
 by setting parameter \code{add.pseudo=TRUE}. 
}
\value{
 None, the output is a pdf-file.
}
\author{Christian Stratowa, based on ideas of package affyQCReport.}
\note{
 Function \code{xpsQAReport} requires a working LaTeX implementation and so will only work on Windows
 platforms, and on OS X, if the user has installed the necessary LaTeX tools. 
}
\examples{
\dontrun{
 ## first, load ROOT scheme file and ROOT data file
 scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
 data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

 ## optional normalized expression levels
 data.rma <- rma(data.test3,"Test3RMA",tmpdir="",background="pmonly",normalize=TRUE,verbose=FALSE)

 ## optional MAS5 detection call
 call.mas5 <- mas5.call(data.test3,"Test3Call",tmpdir="",verbose=FALSE)

 ## optional quality measures
 rlm.all <- rmaPLM(data.test3, "tmp_Test3RLMall", filedir=getwd(), tmpdir="", qualopt="all", option="transcript", add.data=FALSE)

 ## quality assessment report
 xpsQAReport(data.test3, data.rma, call.mas5, rlm.all, dataset="My Dataset", add.pseudo=TRUE, overwrite=TRUE)
}
}
\keyword{manip}
