\name{root.scheme}
\alias{root.scheme}
\title{Create class SchemeTreeSet accessing ROOT scheme file}
\description{
  Create class SchemeTreeSet accessing ROOT scheme file.
}
\usage{
root.scheme(rootfile = character(0), add.mask = FALSE)
}
\arguments{
  \item{rootfile}{name of ROOT scheme file, including full path.}
  \item{add.mask}{if \code{TRUE} mask information will be included as slot \code{mask}.}
}
\details{
An S4 class \code{\link{SchemeTreeSet}} will be created, serving as R wrapper to the 
\code{\link{ROOT}} scheme file \code{rootfile}.
}
\value{
  A \code{SchemeTreeSet} object.
}
\author{Christian Stratowa}
\note{Use this function to access the \code{\link{ROOT}} scheme file from 
new R sessions to avoid creating a new \code{\link{ROOT}} scheme file for every R session.

Do not set \code{add.mask=TRUE} for exon arrays unless you know that your computer has 
sufficient RAM.
}
\seealso{\code{\link{import.expr.scheme}}, \code{\link{import.exon.scheme}}, \code{\link{SchemeTreeSet}}}
\examples{
## create class SchemeSet to access the ROOT scheme file for the Test3 GeneChip
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
str(scheme.test3)

\dontrun{
## scheme set for existing human root exon scheme file
scheme.huex10stv2r2.na22 <- root.scheme("/my/path/schemes/Scheme_HuEx10stv2r2_na22.root")
}
}
\keyword{manip}
