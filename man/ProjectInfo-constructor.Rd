\name{ProjectInfo-constructor}
\alias{ProjectInfo}
\title{Constructor for Class ProjectInfo}
\description{Constructor for class ProjectInfo class allows to save the relevant project 
information in the \code{\link{ROOT}} data file and in class \code{\link{DataTreeSet}}.}
\usage{
ProjectInfo(submitter      = character(),
            laboratory     = character(),
            contact        = character(),
            project        = character(),
            author         = character(),
            dataset        = character(),
            source         = character(),
            sample         = character(),
            celline        = character(),
            primarycell    = character(),
            tissue         = character(),
            biopsy         = character(),
            arraytype      = character(),
            hybridizations = character(),
            treatments     = character())
}
\arguments{
    \item{submitter}{\code{"character"} representing the name of the submitter.}
    \item{laboratory}{\code{"character"} representing the laboratory of the submitter.}
    \item{contact}{\code{"character"} representing the contact address of the submitter.}
    \item{project}{\code{"character"} vector representing the project information.}
    \item{author}{\code{"character"} vector representing the author information.}
    \item{dataset}{\code{"character"} vector representing the dataset information.}
    \item{source}{\code{"character"} vector representing the sample source information.}
    \item{sample}{\code{"character"} vector representing the sample  information.}
    \item{celline}{\code{"character"} vector representing the sample information for cell lines.}
    \item{primarycell}{\code{"character"} vector representing the sample information for primary cells.}
    \item{tissue}{\code{"character"} vector representing the sample information for tissues.}
    \item{biopsy}{\code{"character"} vector representing the sample information for biopsies.}
    \item{arraytype}{\code{"character"} vector representing the array information.}
    \item{hybridizations}{\code{"character"} vector representing the hybridization information for each hybridization.}
    \item{treatments}{\code{"character"} vector representing the treatment information for each hybridization.}
}
\details{
 The ProjectInfo constructor allows to save the following project information in the
 \code{\link{ROOT}} data file and in class \code{\link{DataTreeSet}}:
 \tabular{lll}{
  \tab \code{submitter}:\tab name of the submitter. \cr
  \tab \code{laboratory}:\tab laboratory of the submitter. \cr
  \tab \code{contact}:\tab contact address of the submitter. \cr
  \tab \code{project}:\tab character vector c(name,date,type,description,comments). \cr
  \tab \code{author}:\tab character vector c(lastname,firstname,type,company,department,email,
     phone,comments).. \cr
  \tab \code{dataset}:\tab character vector c(name,type,sample,submitter,date,description,comments). \cr
  \tab \code{source}:\tab character vector c(name,type,species,subspecies,description,comments). \cr
  \tab \code{sample}:\tab character vector c(name,type,sex,phenotype,genotype,extraction,
     isxenograft,xenostrain,xenosex,xenoage,xenoageunit,comments). \cr
  \tab \code{celline}:\tab character vector c(name,type,parent,atcc,modification,sex,phenotype,
     genotype,extraction,isxenograft,xenostrain,xenosex,xenoage,xenoageunit,comments). \cr
  \tab \code{primarycell}:\tab character vector c(name,type,date,description,sex,phenotype,
     genotype,extraction,isxenograft,xenostrain,xenosex,xenoage,xenoageunit,comments). \cr
  \tab \code{tissue}:\tab character vector c(name,type,development,morphology,disease,stage,
     donorage,ageunit,status,sex,phenotype,genotype,extraction,isxenograft,xenostrain,xenosex,
     xenoage,xenoageunit,comments). \cr
  \tab \code{biopsy}:\tab character vector c(name,type,morphology,disease,stage,donorage,ageunit,
     status,sex,phenotype,genotype,extraction,isxenograft,xenostrain,xenosex,xenoage,xenoageunit,comments). \cr
  \tab \code{arraytype}:\tab character vector c(chipname,chiptype,description,comments). \cr
  \tab \code{hybridizations}:\tab vector of character vectors with each containing
     c(name,type,inputname,date,preparation,protocol,repname,replica,comments). \cr
  \tab \code{treatments}:\tab vector of character vectors with each containing
     c(name,type,concentration,concentrationunit,time,timeunit,administration,comments).
 }
}
\value{
 An object of type \code{"\linkS4class{ProjectInfo}"}
}
\author{Christian Stratowa}
\note{Function \code{ProjectInfo} is used as constructor for class \code{ProjectInfo} so that the
 user need not know details for creating S4 classes.
}
\seealso{\code{\linkS4class{ProjectInfo}}}
\examples{
## fill character vectors within constructor
project <- ProjectInfo(submitter="Christian", laboratory="home",contact="email",
                       project=c("TestProject","20060106","Project Type","use Test3 data for testing","my comment"),
                       hybridizations=c(c("TestA1","hyb type","TestA1.CEL",20071117,"my prep1","standard protocol","A1",1,"my comment"),
                                        c("TestA2","hyb type","TestA2.CEL",20071117,"my prep2","standard protocol","A2",1,"my comment"),
                                        c("TestB1","hyb type","TestB1.CEL",20071117,"my prep1","standard protocol","B1",2,"my comment"),
                                        c("TestB2","hyb type","TestB2.CEL",20071117,"my prep2","standard protocol","B2",2,"my comment")))
str(project)

## alternatively add character vectors as methods after creation of constructor
authorInfo(project)    <- c("Stratowa","Christian","Project Leader","Company","Dept","cstrato.at.aon.at","++43-1-1234","my comment")
datasetInfo(project)   <- c("Test3Set","MC","Tissue","Stratowa","20060106","description","my comment")
treatmentInfo(project) <- c(c("TestA1","DMSO",4.3,"mM",1.0,"hours","intravenous","my comment"),
                            c("TestA2","DMSO",4.3,"mM",8.0,"hours","intravenous","my comment"),
                            c("TestB1","DrugA2",4.3,"mM",1.0,"hours","intravenous","my comment"),
                            c("TestB2","DrugA2",4.3,"mM",8.0,"hours","intravenous","my comment"))
str(project)
}
\keyword{manip}
