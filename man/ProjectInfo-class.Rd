\name{ProjectInfo-class}
\docType{class}
\alias{ProjectInfo-class}
\alias{arrayInfo,ProjectInfo-method}
\alias{arrayInfo}
\alias{arrayInfo<-,ProjectInfo,character-method}
\alias{arrayInfo<-}
\alias{authorInfo,ProjectInfo-method}
\alias{authorInfo}
\alias{authorInfo<-,ProjectInfo,character-method}
\alias{authorInfo<-}
\alias{biopsyInfo,ProjectInfo-method}
\alias{biopsyInfo}
\alias{biopsyInfo<-,ProjectInfo,character-method}
\alias{biopsyInfo<-}
\alias{cellineInfo,ProjectInfo-method}
\alias{cellineInfo}
\alias{cellineInfo<-,ProjectInfo,character-method}
\alias{cellineInfo<-}
\alias{datasetInfo,ProjectInfo-method}
\alias{datasetInfo}
\alias{datasetInfo<-,ProjectInfo,character-method}
\alias{datasetInfo<-}
\alias{hybridizInfo,ProjectInfo-method}
\alias{hybridizInfo}
\alias{hybridizInfo<-,ProjectInfo,character-method}
\alias{hybridizInfo<-}
\alias{primcellInfo,ProjectInfo-method}
\alias{primcellInfo}
\alias{primcellInfo<-,ProjectInfo,character-method}
\alias{primcellInfo<-}
\alias{projectInfo,ProjectInfo-method}
\alias{projectInfo}
\alias{projectInfo<-,ProjectInfo,character-method}
\alias{projectInfo<-}
\alias{sampleInfo,ProjectInfo-method}
\alias{sampleInfo}
\alias{sampleInfo<-,ProjectInfo,character-method}
\alias{sampleInfo<-}
\alias{show,ProjectInfo-method}
\alias{sourceInfo,ProjectInfo-method}
\alias{sourceInfo}
\alias{sourceInfo<-,ProjectInfo,character-method}
\alias{sourceInfo<-}
\alias{tissueInfo,ProjectInfo-method}
\alias{tissueInfo}
\alias{tissueInfo<-,ProjectInfo,character-method}
\alias{tissueInfo<-}
\alias{treatmentInfo,ProjectInfo-method}
\alias{treatmentInfo}
\alias{treatmentInfo<-,ProjectInfo,character-method}
\alias{treatmentInfo<-}

\title{Class ProjectInfo}
\description{This class allows to save the relevant project information in the 
\code{\link{ROOT}} data file and in class \code{\link{DataTreeSet}}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \cr
\code{new("ProjectInfo", submitter=[character], laboratory=[character], contact=[character], ...)}. \cr
Alternatively, the constructor \code{\link{ProjectInfo}} can be used.
}
\section{Slots}{
	 \describe{
    \item{\code{submitter}:}{Object of class \code{"character"} representing the name of the submitter.}
    \item{\code{laboratory}:}{Object of class \code{"character"} representing the laboratory of the submitter.}
    \item{\code{contact}:}{Object of class \code{"character"} representing the contact address of the submitter.}
    \item{\code{project}:}{Object of class \code{"list"} representing the project information.}
    \item{\code{author}:}{Object of class \code{"list"} representing the author information.}
    \item{\code{dataset}:}{Object of class \code{"list"} representing the dataset information.}
    \item{\code{source}:}{Object of class \code{"list"} representing the sample source information.}
    \item{\code{sample}:}{Object of class \code{"list"} representing the sample  information.}
    \item{\code{celline}:}{Object of class \code{"list"} representing the sample information for cell lines.}
    \item{\code{primarycell}:}{Object of class \code{"list"} representing the sample information for primary cells.}
    \item{\code{tissue}:}{Object of class \code{"list"} representing the sample information for tissues.}
    \item{\code{biopsy}:}{Object of class \code{"list"} representing the sample information for biopsies.}
    \item{\code{arraytype}:}{Object of class \code{"list"} representing the array information.}
    \item{\code{hybridizations}:}{Object of class \code{"data.frame"} representing the hybridization information for each hybridization.}
    \item{\code{treatments}:}{Object of class \code{"data.frame"} representing the treatment information for each hybridization.}
  }
}
\section{Methods}{
  \describe{
    \item{projectInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{project}.}
    \item{projectInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{project} with character vector c(name,date,type,description,comments).}
    \item{authorInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{author}.}
    \item{authorInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{author} with character vector c(lastname,firstname,type,company,department,email,
     phone,comments).}
    \item{datasetInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{dataset}.}
    \item{datasetInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{dataset} with character vector c(name,type,sample,submitter,date,description,comments).}
    \item{sourceInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{source}.}
    \item{sourceInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{source} with character vector c(name,type,species,subspecies,description,comments).}
    \item{sampleInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{sample}.}
    \item{sampleInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{sample} with character vector c(name,type,sex,phenotype,genotype,extraction,
     isxenograft,xenostrain,xenosex,xenoage,xenoageunit,comments).}
    \item{cellineInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{celline}.}
    \item{cellineInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{celline} with character vector c(name,type,parent,atcc,modification,sex,phenotype,
     genotype,extraction,isxenograft,xenostrain,xenosex,xenoage,xenoageunit,comments).}
    \item{primcellInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{primarycell}.}
    \item{primcellInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{primarycell} with character vector c(name,type,date,description,sex,phenotype,
     genotype,extraction,isxenograft,xenostrain,xenosex,xenoage,xenoageunit,comments).}
    \item{tissueInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{tissue}.}
    \item{tissueInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{tissue} with character vector c(name,type,development,morphology,disease,stage,
     donorage,ageunit,status,sex,phenotype,genotype,extraction,isxenograft,xenostrain,xenosex,
     xenoage,xenoageunit,comments).}
    \item{biopsyInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{biopsy}.}
    \item{biopsyInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{biopsy} with character vector c(name,type,morphology,disease,stage,donorage,ageunit,
     status,sex,phenotype,genotype,extraction,isxenograft,xenostrain,xenosex,xenoage,xenoageunit,comments).}
    \item{arrayInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{arraytype}.}
    \item{arrayInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{arraytype} with character vector c(chipname,chiptype,description,comments).}
    \item{hybridizInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{hybridizations}.}
    \item{hybridizInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{hybridizations} with vector of character vectors with each containing
     c(name,type,inputname,date,preparation,protocol,repname,replica,comments).}
    \item{treatmentInfo}{\code{signature(object = "ProjectInfo")}: extracts slot \code{treatments}.}
    \item{treatmentInfo<-}{\code{signature(object = "ProjectInfo", value = "character")}: replaces 
     slot \code{treatments} with vector of character vectors with each containing
     c(name,type,concentration,concentrationunit,time,timeunit,administration,comments).}
    \item{show}{\code{signature(object = "ProjectInfo")}: shows the content of \code{ProjectInfo}.}
	 }
}
\author{Christian Stratowa}
\examples{
project <- new("ProjectInfo",submitter="Christian", laboratory="home",contact="email")
projectInfo(project)    <- c("TestProject","20060106","Project Type","use Test3 data for testing","my comment")
authorInfo(project)     <- c("Stratowa","Christian","Project Leader","Company","Dept","cstrato.at.aon.at","++43-1-1234","my comment")
datasetInfo(project)    <- c("Test3Set","MC","Tissue","Stratowa","20060106","description","my comment")
sourceInfo(project)     <- c("Unknown","source type","Homo sapiens","caucasian","description","my comment")
primcellInfo(project)   <- c("Mel31","primary cell",20071123,"extracted from patient","male","my pheno","my genotype","RNA extraction",TRUE,"NMRI","female",7.0,"months", "my comment")
arrayInfo(project)      <- c("Test3","GeneChip","description","my comment")
hybridizInfo(project)   <- c(c("TestA1","hyb type","TestA1.CEL",20071117,"my prep1","standard protocol","A1",1,"my comment"),
                             c("TestA2","hyb type","TestA2.CEL",20071117,"my prep2","standard protocol","A2",1,"my comment"),
                             c("TestB1","hyb type","TestB1.CEL",20071117,"my prep1","standard protocol","B1",2,"my comment"),
                             c("TestB2","hyb type","TestB2.CEL",20071117,"my prep2","standard protocol","B2",2,"my comment"))
treatmentInfo(project)  <- c(c("TestA1","DMSO",4.3,"mM",1.0,"hours","intravenous","my comment"),
                             c("TestA2","DMSO",4.3,"mM",8.0,"hours","intravenous","my comment"),
                             c("TestB1","DrugA2",4.3,"mM",1.0,"hours","intravenous","my comment"),
                             c("TestB2","DrugA2",4.3,"mM",8.0,"hours","intravenous","my comment"))
show(project)
}
\keyword{classes}
