#------------------------------------------------------------------------------#
# Script: step-by-step functions to demonstrate how to use package "xps"
#
# Note: please feel free to copy-paste the examples of interest and adapt the
#       examples to your own needs
#
# Copyright (c) 2007-2012 Christian Stratowa, Vienna, Austria.
# All rights reserved.
#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# 1. step: import Affymetrix chip definition and annotation files into
#          ROOT scheme files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: ROOT scheme files need only be created when new annotation files
#          are available. In order to allow all users to access the ROOT 
#          scheme files it is recommended to create all ROOT scheme files
#          in a common system directory and use function "root.scheme()"
#          to load the desired scheme into the current R session to analyse
#          a project.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: See script "script4schemes.R" how to create root scheme files
#          for the latest Affymetrix annotation files.
#------------------------------------------------------------------------------#

### new R session: load library xps
library(xps)

### define directories:
# directory containing Affymetrix library files
libdir <- "/Volumes/GigaDrive/Affy/libraryfiles"
# directory containing Affymetrix annotation files
anndir <- "/Volumes/GigaDrive/Affy/Annotation"
# directory to store ROOT scheme files
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for expression arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# note: see script "script4schemes.R" how to create the newest scheme files
#       for many Affymetrix expression arrays
# note: do not separate name of ROOT files with dots, use underscores,
#       e.g. do not use "Scheme.Test3.na32" but "Scheme_Test3_na32"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt expression arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Test3:
scheme.test3.na32 <- import.expr.scheme("test3", filedir = file.path(scmdir, "na32"),
                     schemefile = file.path(libdir, "Test3.CDF"), 
                     probefile  = file.path(libdir, "Test3_probe.tab"), 
                     annotfile  = file.path(anndir, "Version11Jul", "Test3.na32.annot.csv"))

# HG-U133_Plus_2:
scheme.hgu133plus2.na32 <- import.expr.scheme("hgu133plus2", filedir = file.path(scmdir, "na32"),
                           schemefile = file.path(libdir, "HG-U133_Plus_2.CDF"), 
                           probefile  = file.path(libdir, "HG-U133-PLUS_probe.tab"), 
                           annotfile  = file.path(anndir, "Version11Jul", "HG-U133_Plus_2.na32.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for ivt plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HT_HG-U133_Plus_PM
scheme.hthgu133pluspm.na32 <- import.expr.scheme("hthgu133pluspm", filedir = file.path(scmdir, "na32"),
                              schemefile = file.path(libdir, "HT_HG-U133_Plus_PM.CDF"), 
                              probefile  = file.path(libdir, "HT_HG-U133_Plus_PM.probe.tab"), 
                              annotfile  = file.path(anndir, "Version11Jul", "HT_HG-U133_Plus_PM.na32.annot.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_0-st-v1.r4: used as exon array
scheme.hugene10stv1.na32 <- import.exon.scheme("hugene10stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "HuGene-1_0-st-v1.na32.hg19.probeset.csv"),
                            file.path(anndir, "Version11Jul", "HuGene-1_0-st-v1.na32.hg19.transcript.csv"))

# HuEx-1_0-st-v2.r2:
scheme.huex10stv2.na32 <- import.exon.scheme("huex10stv2", filedir = file.path(scmdir, "na32"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.clf"),
                          file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.pgf"),
                          file.path(anndir, "Version11Jul", "HuEx-1_0-st-v2.na32.hg19.probeset.csv"),
                          file.path(anndir, "Version11Jul", "HuEx-1_0-st-v2.na32.hg19.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for whole genome and exon plate arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# HuGene-1_1-st-v1:
scheme.hugene11stv1.na32 <- import.exon.scheme("hugene11stv1", filedir = file.path(scmdir, "na32"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.clf"),
                            file.path(libdir, "HuGene-1_1-st-v1.r4.analysis-lib-files", "HuGene-1_1-st-v1.r4.pgf"),
                            file.path(anndir, "Version11Jul", "HuGene-1_1-st-v1.na32.hg19.probeset.csv"),
                            file.path(anndir, "Version11Jul", "HuGene-1_1-st-v1.na32.hg19.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create ROOT scheme files for alternative CDF-files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# alternative CDF-files, e.g.:
#AffyProbeMiner
apmdir <- "/Volumes/GigaDrive/Affy/CDF_alternative/AffyProbeMiner"
scheme.hgu133p2.apm <- import.expr.scheme("Scheme_HGU133p2_apm_refseq",filedir=scmdir,paste(apmdir,"HG-U133_Plus_2_transcript_refseq/HG-U133_Plus_2_transcript_refseq.cdf",sep="/"),paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),"","HG-U133_Plus_2_apm_refseq")
#BrainArray (UniMichigan)
umidir <- "/Volumes/GigaDrive/Affy/CDF_alternative/UniMichigan"
scheme.hgu133p2.umi <- import.expr.scheme("Scheme_HGU133p2_umi_refseq_v10",filedir=scmdir,paste(umidir,"Version10/Hs133P_Hs_REFSEQ/Hs133P_Hs_REFSEQ.cdf",sep="/"),paste(umidir,"Version10/Hs133P_Hs_REFSEQ/Hs133P_Hs_REFSEQ_probe.tab",sep="/"),"","HG-U133_Plus_2_umi_refseq")



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# utility functions to demonstrate how to access scheme files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### export different trees from ROOT scheme file
# Test3: export as table only
export(scheme.test3.na32, treetype="idx", outfile="Test3_idx.txt")
export(scheme.test3.na32, treetype="scm", outfile="Test3_scm.txt")
export(scheme.test3.na32, treetype="prb", outfile="Test3_prb.txt")
export(scheme.test3.na32, treetype="ann", outfile="Test3_ann.txt")

# export as table and import as data.frame
idx <- export(scheme.test3.na32, treetype="idx", outfile="Test3_idx.txt",as.dataframe=TRUE)
ann <- export(scheme.test3.na32, treetype="ann", outfile="Test3_ann.txt",as.dataframe=TRUE)

### attach mask later: if import parameter was: as.dataframe=FALSE
scheme.test3.na32 <- attachMask(scheme.test3.na32)
str(scheme.test3.na32)
### export scheme mask
msk <- chipMask(scheme.test3.na32)
scheme.test3.na32 <- removeMask(scheme.test3.na32)
str(scheme.test3.na32)

### scheme accessors
rootFile(scheme.test3.na32)
chipName(scheme.test3.na32)
chipType(scheme.test3.na32)
probeInfo(scheme.test3.na32)

### browse ROOT scheme files
root.browser(scheme.test3.na32)



#------------------------------------------------------------------------------#
# 2. step: import CEL-files into ROOT data files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: ROOT raw data files need only be created once for each new project.
#          They can be updated later to add additional CEL-files to a project.
#          In order to allow all users to access the ROOT data files it is
#          recommended to create all ROOT data files in a common system
#          directory and use function "root.data()" to load the desired data
#          into the current R session to analyse a project.
#          To distinguish ROOT raw data files from other ROOT files, "_cel" is
#          added to each raw data file, e.g. "DataTest3_cel.root".
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 1: Test3 samples from package xps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### define directories:
# directory containing Test3 CEL files
celdir <- paste(path.package("xps"),"raw",sep="/")

### import raw data
# first, import ROOT scheme file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
# import CEL files
data.test3 <- import.data(scheme.test3, "DataTest3", celdir=celdir)
str(data.test3)

# alternatively import CEL files and information about the project
project <- new("ProjectInfo",submitter="Christian", laboratory="home",contact="email")
projectInfo(project)    <- c("TestProject","20060106","Project Type","use Test3 data for testing","my comment")
authorInfo(project)     <- c("Stratowa","Christian","Project Leader","Company","Dept","cstrato.at.aon.at","++43-1-1234","my comment")
datasetInfo(project)    <- c("Test3Set","MC","Tissue","Stratowa","20060106","description","my comment")
sourceInfo(project)     <- c("Unknown","source type","Homo sapiens","caucasian","description","my comment")
cellineInfo(project)    <- c("HeLa-S3","cell type","HeLa","ATCC-12.3","pCSV transfected","female","my pheno","my genotype","RNA extraction",FALSE,"","",0.0,"", "my comment")
arrayInfo(project)      <- c("Test3","GeneChip","description","my comment")
hybridizInfo(project)   <- c(c("TestA1","hyb type","TestA1.CEL",20071117,"my prep1","standard protocol","A1",1,"my comment"),
                             c("TestA2","hyb type","TestA2.CEL",20071117,"my prep2","standard protocol","A2",1,"my comment"),
                             c("TestB1","hyb type","TestB1.CEL",20071117,"my prep1","standard protocol","B1",2,"my comment"),
                             c("TestB2","hyb type","TestB2.CEL",20071117,"my prep2","standard protocol","B2",2,"my comment"))
treatmentInfo(project)  <- c(c("TestA1","DMSO",4.3,"mM",1.0,"hours","intravenous","my comment"),
                             c("TestA2","DMSO",4.3,"mM",8.0,"hours","intravenous","my comment"),
                             c("TestB1","DrugA2",4.3,"mM",1.0,"hours","intravenous","my comment"),
                             c("TestB2","DrugA2",4.3,"mM",8.0,"hours","intravenous","my comment"))
# need to delete old ROOT file "DataTest3_cel.root" first!
data.test3 <- import.data(scheme.test3,"DataTest3",filedir="./data",celdir=celdir,project=project)
str(data.test3)

### update raw data files
# import first set of CEL files
data.test3 <- import.data(scheme.test3, "NewDataTest3", celdir=celdir, celfiles=c("TestA1.CEL","TestA2.CEL"))
str(data.test3)
# import additional CEL-files and update ROOT data file
data.test3 <- addData(data.test3, celdir=celdir, celfiles=c("TestB1.CEL","TestB2.CEL"))
str(data.test3)

### to get only subset of trees from an exisiting ROOT data file
subdata.test3 <- root.data(scheme.test3, rootFile(data.test3), c("TestA1","TestB2","TestA2"))
str(subdata.test3)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# demonstration how to load raw data, access the data, and plot the data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### 1.load existing ROOT scheme file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))

### 2.load existing ROOT data file
rootfile <- paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/")
data.test3 <- root.data(scheme.test3, rootfile=rootfile)

# alternatively load only subset from a ROOT data file into a different R session
subdata.test3 <- root.data(scheme.test3, rootfile=rootfile, c("TestA1","TestB2","TestA2"))

### two possibilities to get tree data as data.frame
# 1.possibility: need to attach data first to "data.test3"
# attach intensity data
data.test3 <- attachInten(data.test3)
# get data
tmp <- intensity(data.test3)
head(tmp)
# also possible to attach only subset
data.test3 <- attachInten(data.test3, c("TestB1.cel","TestA2"))
# get data
tmp <- intensity(data.test3)
head(tmp)
# to avoid memory comsumption of R do:
data.test3 <- removeInten(data.test3)

# 2.possibility:
tmp <- getTreeData(data.test3)
head(tmp)
tmp <- getTreeData(data.test3, varlist="fStdev")
head(tmp)
tmp <- getTreeData(data.test3, varlist="fNPixels")
head(tmp)
tmp <- getTreeData(data.test3, varlist="fInten:fStdev")
head(tmp)

### export different trees from ROOT data file
# export tree data
export(data.test3, treetype="cel", outfile="TestAB_cel.txt")
export(data.test3, treetype="cel", varlist = "fInten", outfile="TestAB_int_cel.txt")
export(data.test3, treename="TestB1", treetype="cel", varlist = "fInten", outfile="TestB1_int_cel.txt")
export(data.test3, treename="TestB1", treetype="cel", varlist = "fInten:fStdev", outfile="TestB1_int_sd__cel.txt")
export(data.test3, treename="TestB1", treetype="cel", varlist = "*", outfile="TestB1_cel.txt")
export(data.test3, treename=c("TestB1","TestA2"), treetype="cel", varlist = "fInten", outfile="TestB1A2_int_cel.txt")

# export data directly from root file
schemefile <- paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/")
datafile   <- paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/")
export.root(datafile, schemefile, "DataSet", "*", "cel", "*", "DataOutFile")

# inspect ROOT file with ROOT browser (to quit ROOT, type ".q")
root.browser(data.test3)

### Note: no longer needed to attachInten()
# plots
hist(data.test3)
image(data.test3)
image(data.test3, names="TestA2.cel", add.legend=TRUE)
boxplot(data.test3, which="userinfo:fIntenQuant")
mboxplot(data.test3, ylim=c(-6,6))

# plots to export
plotImage(data.test3, type="intensity", names="*")
plotImage(data.test3, type="intensity", names="TestA2.cel")
plotImage(data.test3, type="intensity", names="*", dev="png", outfile="Image_data_inten1")
plotImage(data.test3, type="intensity", names="TestA2.cel", dev="png", outfile="Image_TestA2")

plotBoxplot(data.test3, which="userinfo:fIntenQuant")
plotBoxplot(data.test3, which="userinfo:fIntenQuant", dev="png", outfile="Boxplot_DataTest3")
plotBoxplot(data.test3, which="userinfo:fIntenQuant", dev="jpeg", outfile="Boxplot_DataTest3")
plotBoxplot(data.test3, which="userinfo:fIntenQuant", dev="pdf", outfile="Boxplot_DataTest3")


## PM-plot raw data: need to attachInten()
# need to attach scheme mask, since it was not attached to scheme
data.test3 <- attachMask(data.test3)
# need to attach data first
data.test3 <- attachInten(data.test3)
str(data.test3)

pmplot(data.test3)

# to avoid memory comsumption of R remove data:
data.test3 <- removeInten(data.test3)
data.test3 <- removeMask(data.test3)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 2: Tissues from Affymetrix Exon Array Dataset for HG-U133_Plus_2 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### define directories:
# directory of ROOT scheme files
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
# directory containing Tissues CEL files
celdir <- "/Volumes/GigaDrive/ChipData/Exon/HuMixture"
# directory to store ROOT raw data files
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"

### HG-U133_Plus_2 data: import raw data
# first, import ROOT scheme file
scheme.u133p2 <- root.scheme(paste(scmdir,"hgu133plus2.root",sep="/"))

# subset of CEL files to import
celfiles <- c("u1332plus_ivt_breast_A.CEL","u1332plus_ivt_breast_B.CEL","u1332plus_ivt_breast_C.CEL",
              "u1332plus_ivt_prostate_A.CEL","u1332plus_ivt_prostate_B.CEL","u1332plus_ivt_prostate_C.CEL")
# rename CEL files
celnames <- c("BreastA","BreastB","BreastC","ProstateA","ProstateB","ProstateC")
# import CEL files
data.mix.u133p2 <- import.data(scheme.u133p2, "HuTissuesU133P2", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 3: Tissues from Affymetrix Exon Array Dataset for HuEx-1_0-st-v2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# same R session for example

### HuEx-1_0-st-v2 data: import raw data
# first, import ROOT scheme file
scheme.exon <- root.scheme(paste(scmdir,"huex10stv2.root",sep="/"))

# subset of CEL files to import
celfiles <- c("huex_wta_breast_A.CEL","huex_wta_breast_B.CEL","huex_wta_breast_C.CEL",
              "huex_wta_prostate_A.CEL","huex_wta_prostate_B.CEL","huex_wta_prostate_C.CEL")
# rename CEL files
celnames <- c("BreastA","BreastB","BreastC","ProstateA","ProstateB","ProstateC")
# import CEL files
data.mix.exon <- import.data(scheme.exon, "HuTissuesExon", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 4: Tissues from Affymetrix Exon Array Dataset for HuGene-1_0-st-v1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# same R session for example

# directory containing Tissues CEL files
celdir <- "/Volumes/GigaDrive/ChipData/Exon/HuGene"

### HuGene-1_0-st-v1 data: import raw data
# first, import ROOT scheme file
scheme.genome <- root.scheme(paste(scmdir,"hugene10stv1.root",sep="/"))

# subset of CEL files to import
celfiles <- c("TisMap_Breast_01_v1_WTGene1.CEL","TisMap_Breast_02_v1_WTGene1.CEL","TisMap_Breast_03_v1_WTGene1.CEL",
              "TisMap_Prostate_01_v1_WTGene1.CEL","TisMap_Prostate_02_v1_WTGene1.CEL","TisMap_Prostate_03_v1_WTGene1.CEL")
# rename CEL files
celnames <- c("Breast01","Breast02","Breast03","Prostate01","Prostate02","Prostate03")
# import CEL files
data.mix.genome <- import.data(scheme.genome, "HuTissuesGenome", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 4a: Tissues from Affymetrix Human Gene 2.1 ST Plate Data Set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# same R session for example

# directory containing Tissues CEL files
celdir <- "/Volumes/GigaDrive/ChipData/Exon/HuGene2/human2.0/HuGene2.1_Plate"

### HuGene-2_1-st data: import raw data
# first, import ROOT scheme file
scheme.genome <- root.scheme(file.path(scmdir, "hugene21stv1.root"))

# subset of CEL files to import
celfiles <- c("Liver_HuGene-2_1_GT_Rep1_A03_MC.CEL","Liver_HuGene-2_1_GT_Rep2_D06_MC.CEL","Liver_HuGene-2_1_GT_Rep3_F02_MC.CEL",
              "Spleen_HuGene-2_1_GT_Rep1_A11_MC.CEL","Spleen_HuGene-2_1_GT_Rep2_C07_MC.CEL","Spleen_HuGene-2_1_GT_Rep3_F04_MC.CEL")
# rename CEL files
celnames <- c("LiverRep1","LiverRep2","LiverRep3","SpleenRep1","SpleenRep2","SpleenRep3")
# import CEL files
data.genome <- import.data(scheme.genome, "HuTissuesGenome21", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 5: Test data from Affymetrix "HT_PM_human_tissue_panel" Dataset for HT_HG-U133_Plus_PM 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### define directories:
# directory of ROOT scheme files
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
# directory containing Tissues CEL files
celdir <- "/Volumes/GigaDrive/ChipData/Plate/HT_PM_human_tissue_panel"
# directory to store ROOT raw data files
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"

### HT_HG-U133_Plus_PM data: import raw data
# first, import ROOT scheme file
scheme.u133ppm <- root.scheme(paste(scmdir,"hthgu133pluspm.root",sep="/"))

# subset of CEL files to import
celfiles <- c("Human_PM_TestData.A01.CEL","Human_PM_TestData.A02.CEL","Human_PM_TestData.A03.CEL",
              "Human_PM_TestData.B01.CEL","Human_PM_TestData.B02.CEL","Human_PM_TestData.B03.CEL")
# rename CEL files
celnames <- c("TestDataA01","TestDataA02","TestDataA03","TestDataB01","TestDataB02","TestDataB03")
# import CEL files
data.mix.u133ppm <- import.data(scheme.u133ppm, "TestDataHTU133PPM", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 6: Test data from Affymetrix "genetitan_plate_1_sample_data" for HT_HG-U133_Plus_PM 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### define directories:
# directory of ROOT scheme files
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
# directory containing Tissues CEL files
celdir <- "/Volumes/GigaDrive/ChipData/Plate/genetitan_plate_1_sample_data/Plate_1_9628"
# directory to store ROOT raw data files
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"

### HT_HG-U133_Plus_PM data: import raw data
# first, import ROOT scheme file
scheme.u133ppm <- root.scheme(paste(scmdir,"hthgu133pluspm.root",sep="/"))

# subset of CEL files to import
celfiles <- c("plate_1_id9628_A01.CEL","plate_1_id9628_B02.CEL","plate_1_id9628_C03.CEL",
              "plate_1_id9628_A07.CEL","plate_1_id9628_B08.CEL","plate_1_id9628_C09.CEL")
# rename CEL files
celnames <- c("Hela_A01","Hela_B02","Hela_C03","maqc_A_A07","maqc_A_B08","maqc_A_C09")
# import CEL files
data.genetitan <- import.data(scheme.u133ppm, "TestDataGeneTitan", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 7: Sample data from Affymetrix "MAQC A and B" Dataset for HG-U219 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### define directories:
# directory of ROOT scheme files
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
# directory containing sample CEL files
celdir <- "/Volumes/GigaDrive/ChipData/Plate/hg-u219-ap-sampledata"
# directory to store ROOT raw data files
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"

### HG-U219 data: import raw data
# first, import ROOT scheme file
scheme.u219 <- root.scheme(paste(scmdir,"hgu219.root",sep="/"))

# subset of CEL files to import
celfiles <- c("HG-U219_MaqcA_1.CEL","HG-U219_MaqcA_2.CEL","HG-U219_MaqcA_3.CEL",
              "HG-U219_MaqcB_1.CEL","HG-U219_MaqcB_2.CEL","HG-U219_MaqcB_3.CEL")
# rename CEL files
celnames <- c("MaqcA1","MaqcA2","MaqcA3","MaqcB1","MaqcB2","MaqcB3")
# import CEL files
data.u219 <- import.data(scheme.u219, "MaqcDataHGU219", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# demonstration how to access the data, and plot the data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

# import ROOT scheme files
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
scheme.u133p2 <- root.scheme(paste(scmdir,"hgu133plus2.root",sep="/"))
scheme.exon   <- root.scheme(paste(scmdir,"huex10stv2.root",sep="/"))
scheme.genome <- root.scheme(paste(scmdir,"hugene10stv1.root",sep="/"))

# import ROOT data files
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.u133p2 <- root.data(scheme.u133p2, paste(datdir,"HuTissuesU133P2_cel.root",sep="/"))
data.exon   <- root.data(scheme.exon, paste(datdir,"HuTissuesExon_cel.root",sep="/"))
data.genome <- root.data(scheme.genome, paste(datdir,"HuTissuesGenome_cel.root",sep="/"))

# inspect ROOT file with ROOT browser
root.browser(data.exon)


### plot raw data for HG-U133_Plus_2
# plots
hist(data.u133p2)
image(data.u133p2)
boxplot(data.u133p2, which="userinfo:fIntenQuant")

# plots to export
plotImage(data.u133p2, type="intensity", dev="png", col=rainbow(32), outfile="Image_DataMixU133P2_BrA",names="BreastA.cel")
plotBoxplot(data.u133p2, which="userinfo:fIntenQuant", dev="png", w=600, h=480, outfile="Boxplot_DataMixU133P2")


### plot raw data for HuEx-1_0-st-v2
# plots
names <- unlist(treeNames(data.exon))
image(data.exon, names=names[1], add.legend=TRUE)
hist(data.exon)
hist(data.exon, which="core")
boxplot(data.exon, which="userinfo:fIntenQuant")

# plots to export
outfile <- paste("Image_", names[1], sep="")
plotImage(data.exon, type="intensity", names=names[1], dev="png", outfile=outfile)
plotDensity(data.exon, add.legend=TRUE, dev="png", outfile="Density")
plotBoxplot(data.exon, which="userinfo:fIntenQuant", mar=NULL, dev="png", outfile="Boxplot")

# Note: You can also use the corresponding methods "root.drawxxx()", such as:
root.density(data.exon, "*")
root.image(data.exon, "BreastA.cel")
root.hist2D(data.test3, "BreastA.cel", "BreastB.cel", option="COLZ")


### plot raw data for HuGene-1_0-st-v1
# plots
hist(data.genome, which="core")
image(data.genome, names="Breast01.cel_MEAN", add.legend=TRUE)
boxplot(data.genome, which="userinfo:fIntenQuant")

# plots to export
plotImage(data.genome, type="intensity", names=names"Breast01.cel", dev="png", outfile="Image_DataMixHuGene_Br01")
plotBoxplot(data.genome, which="userinfo:fIntenQuant", dev="png", w=600, h=480, outfile="Boxplot_DataMixHuGene")


#------------------------------------------------------------------------------#
# 3. step: convert raw data to expression levels
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: ROOT scheme files and ROOT raw data files are usually already stored
#          in special system directories. When a new R session is created for the
#          first time, they must fist be loaded using "root.scheme()" and "root.data()".
#          However, this is not necessary when re-opening a saved R session later.
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 1: Test3 samples from package xps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))


### preprocess raw data ###

# 1. RMA
data.rma <- rma(data.test3,"Test3RMA",tmpdir="",background="pmonly",normalize=TRUE)

# 2. MAS5
data.mas5 <- mas5(data.test3,"Test3MAS5",,tmpdir="",normalize=TRUE,sc=500)
# to store all trees (including e.g. background trees) in same ROOT file, use "update=TRUE"
data.mas5 <- mas5(data.test3,"Test3MAS5All",,tmpdir="",normalize=TRUE,sc=500, update=TRUE)

# 3. MAS5 detection call
call.mas5 <- mas5.call(data.test3,"Test3Call",tmpdir="")

# 4. DABG detection call (yes, this is possible for expression arrays)
call.dabg <- dabg.call(data.test3,"Test3DABG")

# get data.frames
expr.rma <- validData(data.rma)
expr.mas5 <- validData(data.mas5)
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)


### plot results ###

# compare mas5 to rma
plot(expr.rma[,1],expr.mas5[,1])
plot(expr.rma[,1],expr.mas5[,1],log="xy",xlim=c(1,20000),ylim=c(1,20000))

# density plots
hist(data.rma)
hist(data.mas5)

# boxplots
boxplot(data.rma)
boxplot.dev(data.rma)
boxplot.dev(data.rma, dev="png", mar=c(4,4,1,1), w=480, h=480, outfile="BoxPlot_Test3_rma")

boxplot(data.mas5)
boxplot.dev(data.mas5)
boxplot.dev(data.mas5, dev="png", mar=c(4,4,1,1), w=480, h=480, outfile="BoxPlot_Test3_mas5")

# relative boxplots
mboxplot(data.rma)
mboxplot(data.rma, ylim=c(-2,3))
mboxplot(data.mas5, ylim=c(-4,5))

# M vs A plots
mvaplot(data.rma, pch=20, ylim=c(-2,2))
mvaplot(data.mas5, pch=20)

# present call plots
callplot(call.mas5)
callplot(call.mas5, beside=F, ylim=c(0,125))
callplot(call.dabg)
callplot(call.dabg, beside=F, ylim=c(0,125))

# image for rma background intensities
# 1. find background treenames for data.rma
getTreeNames(rootFile(data.rma))
# 2. get background intensity for e.g. tree "TestA2.rbg"
bgrd <- export.root(rootFile(data.rma),schemeFile(data.rma),"PreprocesSet","TestA2","rbg","fBg","BgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
rbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.rma)), nrow=nrows(schemeSet(data.rma)))
# 4. create image
image(rbg)
image(log2(rbg))

# image for mas5 background intensities
# 1. find background treenames for data.mas5
getTreeNames(rootFile(data.mas5))
# 2. get background intensity for e.g. tree "TestA2.wbg"
bgrd <- export.root(rootFile(data.mas5),schemeFile(data.mas5),"PreprocesSet","TestA2","wbg","fBg","BgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
wbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.mas5)), nrow=nrows(schemeSet(data.mas5)))
# 4. create image
image(wbg)
image(log2(wbg))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 2: Tissues from Affymetrix Exon Array Dataset for HG-U133_Plus_2 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
scheme.u133p2 <- root.scheme(paste(scmdir,"Scheme_HGU133p2_na31.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.u133p2 <- root.data(scheme.u133p2, paste(datdir,"HuTissuesU133P2_cel.root",sep="/"))


### preprocess raw data ###

# 1. RMA
data.rma <- rma(data.u133p2,"MixU133P2RMA",tmpdir="",background="pmonly",normalize=TRUE)

# 2. MAS5
data.mas5 <- mas5(data.u133p2,"MixU133P2MAS5",,tmpdir="",normalize=TRUE,sc=500)
# to store all trees (including e.g. background trees) in same ROOT file, use "update=TRUE"
data.mas5 <- mas5(data.u133p2,"MixU133P2MAS5All",,tmpdir="",normalize=TRUE,sc=500, update=TRUE)

# 3. MAS5 detection call
call.mas5 <- mas5.call(data.u133p2,"MixU133P2Call",tmpdir="")

# 4. DABG detection call
# Note: yes, it is possible to compute DABG call for expression arrays, but:
#       on my MacBook Pro with 2.3 GHz Intel Core 2 Duo it takes 44 min per array!!!
#call.dabg <- dabg.call(data.u133p2,"MixU133P2DABG")

# get data.frames
expr.rma <- validData(data.rma)
expr.mas5 <- validData(data.mas5)
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)


### plot results ###

# compare mas5 to rma
plot(expr.rma[,1],expr.mas5[,1])
plot(expr.rma[,1],expr.mas5[,1],log="xy",xlim=c(1,20000),ylim=c(1,20000))

# density plots
hist(data.rma)
hist(data.mas5)

# boxplots
boxplot(data.rma)
boxplot(data.mas5)

# relative boxplots
mboxplot(data.rma)
mboxplot(data.rma, ylim=c(-2,3))
mboxplot(data.mas5, ylim=c(-4,5))

# M vs A plots
mvaplot(data.rma, pch=20, ylim=c(-4,4))
mvaplot(data.mas5, pch=20)

# present call plots
callplot(call.mas5)
callplot(call.mas5, beside=FALSE, ylim=c(0,125))

# image for rma background intensities
# 1. find background treenames for data.rma
getTreeNames(rootFile(data.rma))
# 2. get background intensity for e.g. tree "BreastA.rbg"
bgrd <- export.root(rootFile(data.rma),schemeFile(data.rma),"PreprocesSet","BreastA","rbg","fBg","BreastAU133BgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
rbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.rma)), nrow=nrows(schemeSet(data.rma)))
# 4. create image
image(rbg)
image(log2(rbg))

# image for mas5 background intensities
# 1. find background treenames for data.mas5
getTreeNames(rootFile(data.mas5))
# 2. get background intensity for e.g. tree "BreastA.wbg"
bgrd <- export.root(rootFile(data.mas5),schemeFile(data.mas5),"PreprocesSet","BreastA","wbg","fBg","BreastAU133BgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
wbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.mas5)), nrow=nrows(schemeSet(data.mas5)))
# 4. create image
image(wbg)
image(log2(wbg))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 3: Tissues from Affymetrix Exon Array Dataset for HuEx-1_0-st-v2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
scheme.exon <- root.scheme(paste(scmdir,"Scheme_HuEx10stv2r2_na31.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.exon <- root.data(scheme.exon, paste(datdir,"HuTissuesExon_cel.root",sep="/"))


### preprocess raw data ###
datdir <- getwd()

# 1. RMA
# transcript: metacore
# ok for 6 exon arrays in RAM
data.rma <- rma(data.exon,"HuExonMixRMAMetacore",filedir=datdir,tmpdir="",background="antigenomic",
                normalize=TRUE,option="transcript",exonlevel="metacore+affx")
# for many exon arrays you may decide to use tmpdir (see helpfile for more information)
tmpdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/temp"
data.rma.tmp <- rma(data.exon,"HuExonMixRMAtmpMetacore",filedir=datdir,tmpdir=tmpdir,background="antigenomic",
                    normalize=TRUE,option="transcript",exonlevel="metacore+affx")
# probeset: metacore
data.rma.ps <- rma(data.exon,"HuExonMixRMAMetacorePS",filedir=datdir,tmpdir="",background="antigenomic",
                   normalize=TRUE,option="probeset",exonlevel="metacore+affx")

# 2. MAS5
data.mas5 <- mas5(data.exon,"HuExonMixMAS5Metacore",filedir=datdir,tmpdir="",
                  normalize=TRUE,sc=500,option="transcript",exonlevel="metacore+affx")
# to store all trees (including e.g. background trees) in same ROOT file, use "update=TRUE"
data.mas5 <- mas5(data.exon,"HuExonMixExonMAS5MetacoreAll",filedir=datdir,tmpdir="",
                  normalize=TRUE,sc=500,option="transcript",exonlevel="metacore+affx", update=TRUE)

# 3. MAS5 detection call (yes, this is possible for exon arrays)
# note: alpha1 and alpha2 need to be adjusted to get usable P/M/A calls
call.mas5 <- mas5.call(data.exon,"HuExonMixCallMetacore",filedir=datdir,tmpdir="",
                       option="transcript",exonlevel="metacore+affx")

# 4. DABG detection call
# transcript: metacore
call.dabg <- dabg.call(data.exon,"HuExonMixDABGMetacore",filedir=datdir,
                       option="transcript",exonlevel="metacore+affx")
# note: alpha1 and alpha2 need to be adjusted to get usable P/M/A calls for transcripts
#call.dabg <- dabg.call(data.exon,"HuExonMixDABGMetacore", alpha1=???,alpha2=???)
# probeset: metacore
call.dabg.ps <- dabg.call(data.exon,"HuExonMixDABGMetacorePS",filedir=datdir,
                          option="probeset",exonlevel="metacore", alpha1=0.01,alpha2=0.015)

# get data.frames
expr.rma <- validData(data.rma)
expr.mas5 <- validData(data.mas5)
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)
pval.dabg <- pvalData(call.dabg)
pres.dabg <- presCall(call.dabg)


### plot results ###

# compare mas5 to rma
plot(expr.rma[,1],expr.mas5[,1])
plot(expr.rma[,1],expr.mas5[,1],log="xy",xlim=c(1,20000),ylim=c(1,20000))

# density plots
hist(data.rma)
hist(data.mas5)

# boxplots
boxplot(data.rma)
boxplot(data.mas5)

# relative boxplots
mboxplot(data.rma)
mboxplot(data.rma, ylim=c(-3,3))
mboxplot(data.mas5, ylim=c(-4,5))

# M vs A plots
mvaplot(data.rma, pch=20, ylim=c(-4,4))
mvaplot(data.mas5, pch=20)

# present call plots
callplot(call.mas5)
callplot(call.mas5, beside=F, ylim=c(0,125))
callplot(call.dabg)
callplot(call.dabg, beside=F, ylim=c(0,125))
callplot(call.dabg.ps)
callplot(call.dabg.ps, beside=F, ylim=c(0,125))

# image for rma background intensities
# 1. find background treenames for data.rma
getTreeNames(rootFile(data.rma))
# 2. get background intensity for e.g. tree "BreastA.rbg"
bgrd <- export.root(rootFile(data.rma),schemeFile(data.rma),"PreprocesSet","BreastA","rbg","fBg","BreastAExonBgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
rbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.rma)), nrow=nrows(schemeSet(data.rma)))
# 4. create image
image(rbg)
image(log2(rbg))

# image for mas5 background intensities
# 1. find background treenames for data.mas5
getTreeNames(rootFile(data.mas5))
# 2. get background intensity for e.g. tree "BreastA.wbg"
bgrd <- export.root(rootFile(data.mas5),schemeFile(data.mas5),"PreprocesSet","BreastA","wbg","fBg","BreastAExonBgrdROOTOut.txt",as.dataframe=T)
rootfile <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/HuMixture/MixExon/MixMAS5Metacore_tbw.root"
bgrd <- export.root(rootfile,schemeFile(data.mas5),"PreprocesSet","BreastA","wbg","fBg","BreastAExonBgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
wbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.mas5)), nrow=nrows(schemeSet(data.mas5)))
# 4. create image
image(wbg)
image(log2(wbg))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 4: Tissues from Affymetrix Exon Array Dataset for HuGene-1_0-st-v1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
scheme.genome <- root.scheme(paste(scmdir,"Scheme_HuGene10stv1r3_na31.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.genome <- root.data(scheme.genome, paste(datdir,"HuTissuesGenome_cel.root",sep="/"))


### preprocess raw data ###
datdir <- getwd()

# 1. RMA
data.rma <- rma(data.genome,"HuGeneMixRMAMetacore",filedir=datdir,tmpdir="",
                background="antigenomic",normalize=TRUE,exonlevel="metacore+affx")

# 2. MAS5
data.mas5 <- mas5(data.genome,"HuGeneMixMAS5Metacore",filedir=datdir,tmpdir="",
                  normalize=TRUE,sc=500,exonlevel="metacore+affx")
# to store all trees (including e.g. background trees) in same ROOT file, use "update=TRUE"
data.mas5 <- mas5(data.genome,"HuGeneMixExonMAS5MetacoreAll",filedir=datdir,tmpdir="",
                  normalize=TRUE,sc=500,exonlevel="metacore+affx", update=TRUE)

# 3. MAS5 detection call (yes, this is possible for exon arrays)
# note: alpha1 and alpha2 need to be adjusted to get usable P/M/A calls
call.mas5 <- mas5.call(data.genome,"HuGeneMixCallMetacore",filedir=datdir,tmpdir="",
                       exonlevel="metacore+affx")

# 4. DABG detection call
# note: alpha1 and alpha2 need to be adjusted to get usable P/M/A calls for transcripts
call.dabg <- dabg.call(data.genome,"HuGeneMixDABGMetacore",filedir=datdir,
                       exonlevel="metacore+affx")
#call.dabg <- dabg.call(data.genome,"HuGeneMixDABGMetacore", alpha1=???,alpha2=???)

# get data.frames
expr.rma <- validData(data.rma)
expr.mas5 <- validData(data.mas5)
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)
pval.dabg <- pvalData(call.dabg)
pres.dabg <- presCall(call.dabg)


### plot results ###

# compare mas5 to rma
plot(expr.rma[,1],expr.mas5[,1])
plot(expr.rma[,1],expr.mas5[,1],log="xy",xlim=c(1,20000),ylim=c(1,20000))

# density plots
hist(data.rma)
hist(data.mas5)

# boxplots
boxplot(data.rma)
boxplot(data.mas5)

# relative boxplots
mboxplot(data.rma)
mboxplot(data.rma, ylim=c(-3,3))
mboxplot(data.mas5, ylim=c(-4,5))

# M vs A plots
mvaplot(data.rma, pch=20, ylim=c(-4,4))
mvaplot(data.mas5, pch=20)

# present call plots
callplot(call.mas5)
callplot(call.mas5, beside=F, ylim=c(0,125))
callplot(call.dabg)
callplot(call.dabg, beside=F, ylim=c(0,125))

# image for rma background intensities
# 1. find background treenames for data.rma
getTreeNames(rootFile(data.rma))
# 2. get background intensity for e.g. tree "Breast1.rbg"
bgrd <- export.root(rootFile(data.rma),schemeFile(data.rma),"PreprocesSet","Breast01","rbg","fBg","Breast01HuGeneBgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
rbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.rma)), nrow=nrows(schemeSet(data.rma)))
# 4. create image
image(rbg)
image(log2(rbg))

# image for mas5 background intensities
# 1. find background treenames for data.mas5
getTreeNames(rootFile(data.mas5))
# 2. get background intensity for e.g. tree "Breast01.wbg"
bgrd <- export.root(rootFile(data.mas5),schemeFile(data.mas5),"PreprocesSet","Breast01","wbg","fBg","Breast01HuGeneBgrdROOTOut.txt",as.dataframe=T)
# 3. convert bgrd column to matrix
wbg <- matrix(bgrd[,"BGRD"], ncol=ncols(schemeSet(data.mas5)), nrow=nrows(schemeSet(data.mas5)))
# 4. create image
image(wbg)
image(log2(wbg))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 4a: Tissues from Affymetrix Human Gene 2.1 ST Plate Data Set 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
scheme.genome <- root.scheme(file.path(scmdir, "hugene21stv1.root"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.genome <- root.data(scheme.genome, paste(datdir,"HuTissuesGenome21_cel.root",sep="/"))


### preprocess raw data ###
datdir <- getwd()

# 1. RMA
data.rma <- rma(data.genome,"HuGene21RMAcore",filedir=datdir,tmpdir="",
                background="antigenomic",normalize=TRUE,exonlevel="core+affx")

# 2. MAS5
data.mas5 <- mas5(data.genome,"HuGene21MAS5core",filedir=datdir,tmpdir="",
                  normalize=TRUE,sc=500,exonlevel="core+affx")

# 3. MAS5 detection call (yes, this is possible for exon arrays)
call.mas5 <- mas5.call(data.genome,"HuGene21Callcore",filedir=datdir,tmpdir="",
                       exonlevel="core+affx")

# 4. DABG detection call
call.dabg <- dabg.call(data.genome,"HuGene21DABGcore",filedir=datdir,
                       exonlevel="core+affx")

# get data.frames
expr.rma <- validData(data.rma)
expr.mas5 <- validData(data.mas5)
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)
pval.dabg <- pvalData(call.dabg)
pres.dabg <- presCall(call.dabg)

# export expression data
export.expr(data.rma, treename = "*", treetype = "mdp", varlist = "fUnitName:fSymbol:fLevel", outfile = "HuGene21RMAcoreSymbols.txt", sep = "\t", as.dataframe = FALSE, verbose = TRUE)
export.expr(data.rma, treename = "*", treetype = "mdp", varlist = "fUnitName:fName:fSymbol:fLevel", outfile = "HuGene21RMAcoreNamesSymbols.txt", sep = "\t", as.dataframe = FALSE, verbose = TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 5: Test data from Affymetrix HT_PM_human_tissue_panel Dataset for HT_HG-U133_Plus_PM 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
scheme.u133ppm <- root.scheme(paste(scmdir,"Scheme_HTHGU133pPM_na31.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.u133ppm <- root.data(scheme.u133ppm, paste(datdir,"TestDataHTU133PPM_cel.root",sep="/"))


### preprocess raw data ###

# 1. RMA
data.rma <- rma(data.u133ppm,"TestDataHTU133PPM_RMA",tmpdir="",background="pmonly",normalize=TRUE)

# 2. MAS5: 
data.mas5 <- mas5(data.u133ppm,"TestDataHTU133PPM_MAS5",tmpdir="",normalize=TRUE,sc=500)

# get data.frames
expr.rma <- validData(data.rma)

# density plots
hist(data.rma)

# boxplots
boxplot(data.rma)

root.density(data.u133ppm)
root.density(data.rma)
root.profile(data.u133ppm)
root.profile(data.rma)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 6: Test data from Affymetrix "genetitan_plate_1_sample_data" for HT_HG-U133_Plus_PM 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
scheme.u133ppm <- root.scheme(paste(scmdir,"Scheme_HTHGU133pPM_na31.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.u133ppm <- root.data(scheme.u133ppm, paste(datdir,"TestDataGeneTitan_cel.root",sep="/"))

### preprocess raw data ###
# 1. RMA
data.rma <- rma(data.u133ppm,"TestDataGeneTitan_RMA",tmpdir="",background="pmonly",normalize=TRUE)

# 2. MAS5: 
data.mas5 <- mas5(data.u133ppm,"TestDataHTU133PPM_MAS5",tmpdir="",normalize=TRUE,sc=500)

# 3. MAS5 detection call: need to use background-corrected data!!
call.mas5 <- mas5.call(data.u133ppm,"TestDataHTU133PPM_Call",tmpdir="",bgcorrect.option="correctbg")

# get data.frames
expr.rma <- validData(data.rma)
expr.mas5 <- validData(data.mas5)
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)

# compare mas5 to rma
plot(expr.rma[,1],expr.mas5[,1])
plot(expr.rma[,1],expr.mas5[,1],log="xy",xlim=c(1,20000),ylim=c(1,20000))

plot(expr.rma[,1],expr.rma[,2],log="xy")
plot(expr.mas5[,1],expr.mas5[,2],log="xy")

# density plots
hist(data.rma)
hist(data.mas5)

# boxplots
boxplot(data.rma)
boxplot(data.mas5)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 7: Sample data from Affymetrix "MAQC A and B" Dataset for HG-U219 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
scheme.u219 <- root.scheme(paste(scmdir,"Scheme_HGU219_na31.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.u219 <- root.data(scheme.u219, paste(datdir,"MaqcDataHGU219_cel.root",sep="/"))

### preprocess raw data ###
# 1. RMA
data.rma <- rma(data.u219,"MaqcDataHGU219_RMA",tmpdir="",background="pmonly",normalize=TRUE)

# 2. MAS5: not possible since no MM
data.mas5 <- mas5(data.u219,"MaqcDataHGU219_MAS5",tmpdir="",normalize=TRUE,sc=500)

# 3. MAS5 detection call: need to use background-corrected data!!
call.mas5 <- mas5.call(data.u219,"MaqcDataHGU219_Call",tmpdir="",bgcorrect.option="correctbg")

# get data.frames
expr.rma <- validData(data.rma)
expr.mas5 <- validData(data.mas5)
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)

# density plots
hist(data.rma)
hist(data.mas5)

# boxplots
boxplot(data.rma)
boxplot(data.mas5)


#------------------------------------------------------------------------------#
# 4. step: apply filters to expression levels
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: ROOT scheme files and ROOT raw data files are usually already stored
#          in special system directories. When a new R session is created for the
#          first time, they must fist be loaded using "root.scheme()" and "root.data()".
#          However, this is not necessary when re-opening a saved R session later.
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# example 1: Test3 samples from package xps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

### second, preprocess raw data if not already done
# e.g. RMA and MAS5 detection call
data.rma <- rma(data.test3,"Test3RMA",tmpdir="",background="pmonly",normalize=TRUE)
call.mas5 <- mas5.call(data.test3,"Test3Call",tmpdir="")


### apply non-specific filters
# create PreFilter
prefltr <- PreFilter(mad=c(0.5,0.01), prescall=c(0.002, 6,"samples"),
                     lothreshold=c(6.0,0.02,"mean"), hithreshold=c(10.5,80.0,"percent"))
# apply prefilter to data.rma
rma.pfr <- prefilter(data.rma,"Test3Prefilter",getwd(),prefltr,2,"log2","PreFilter",call.mas5)

### apply univariate filters
# create UniFilter
unifltr <- UniFilter(unitest=c("t.test","two.sided","none",0,0.0,FALSE,0.95,TRUE),
                     foldchange=c(1.3,"both"), unifilter=c(0.1,"pval"))
# apply unifilter to pre-filtered data
rma.ufr <- unifilter(data.rma,"Test3Unifilter",getwd(),unifltr,group=c("GrpA","GrpA","GrpB","GrpB"),
                     xps.fltr=rma.pfr)

### get data.frame of result
# get results only for genes satisfying unifltr (default):
ds.ufr <- validData(rma.ufr)
dim(ds.ufr)
head(ds.ufr)

# get results for all genes
ds.all <- validData(rma.ufr,which="UnitName")
dim(ds.all)
head(ds.all)

# alternatively use export.filter to export selected variables only
ds.ufr <- export.filter(rma.ufr,treetype="stt",varlist="fUnitName:fName:fSymbol:mn1:mn2:fc:pval:mask",as.dataframe=TRUE)
dim(ds.ufr)
head(ds.ufr)

ds.all <- export.filter(rma.ufr,treetype="stt",varlist="fUnitName:fName:fSymbol:mn1:mn2:fc:pval:flag",as.dataframe=TRUE)
dim(ds.all)
head(ds.all)



#------------------------------------------------------------------------------#
# Demonstrations of advanced methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: ROOT scheme files and ROOT raw data files are usually already stored
#          in special system directories. When a new R session is created for the
#          first time, they must fist be loaded using "root.scheme()" and "root.data()".
#          However, this is not necessary when re-opening a saved R session later.
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# demonstration 1: compute RMA stepwise for Test3 samples
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

# directory
datdir <- getwd()

### first, load ROOT scheme file and ROOT data file
scheme.test3 <- root.scheme(paste(path.package("xps"),"schemes/SchemeTest3.root",sep="/"))
data.test3 <- root.data(scheme.test3, paste(path.package("xps"),"rootdata/DataTest3_cel.root",sep="/"))

### 1.step: background - rma
data.bg.rma <- bgcorrect.rma(data.test3,"Test3RMABgrd",filedir=datdir)

# attach data
data.bg.rma <- attachMask(data.bg.rma)
data.bg.rma <- attachInten(data.bg.rma)
data.bg.rma <- attachBgrd(data.bg.rma)

# plot intensities
hist(data.bg.rma)
mboxplot(data.bg.rma, ylim=c(-6,6))
pmplot(data.bg.rma)
image(data.bg.rma,col=rainbow(32))

# remove data
data.bg.rma <- removeInten(data.bg.rma)
data.bg.rma <- removeBgrd(data.bg.rma)

### 2step: normalization - quantile
data.qu.rma <- normalize.quantiles(data.bg.rma,"Test3RMANorm",filedir=datdir)

# plot intensities
data.qu.rma <- attachInten(data.qu.rma)

hist(data.qu.rma)
mboxplot(data.qu.rma, ylim=c(-6,6))

data.qu.rma <- removeInten(data.qu.rma)

### 3.step: summarization - medpol
data.mp.rma <- summarize.rma(data.qu.rma,"Test3RMAExpr",filedir=datdir,tmpdir="")

# plot expression levels
hist(data.mp.rma)
boxplot(data.mp.rma)
mboxplot(data.mp.rma)
mvaplot(data.mp.rma, pch=20, ylim=c(-4,4))

### alternatively save all data in same ROOT file using "update=TRUE"
data.bg.rmall <- bgcorrect.rma(data.test3,"Test3RMAall",filedir=datdir,tmpdir="")
data.qu.rmall <- normalize.quantiles(data.bg.rmall,"Test3RMAall",filedir=datdir,tmpdir="",update=TRUE)
data.mp.rmall <- summarize.rma(data.qu.rmall,"Test3RMAall",filedir=datdir,tmpdir="",update=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# demonstration 2: compute RMA for Test3 samples using function "express()"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### compute rma stepwise:
#   important: for stepwise computation tmpdir must be tmpdir="" otherwise the root file will be empty!
expr.bg.rma <- express(data.test3,"Test3ExprsBgrd",filedir=datdir,tmpdir="",update=FALSE,
               bgcorrect.method="rma",bgcorrect.select="none",bgcorrect.option="pmonly:epanechnikov",bgcorrect.params=c(16384))

#   important: for stepwise computation tmpdir must be tmpdir="" otherwise the root file will be empty!
expr.qu.rma <- express(expr.bg.rma,"Test3ExprsNorm",filedir=datdir,tmpdir="",update=FALSE,
               normalize.method="quantile",normalize.select="pmonly",normalize.option="transcript:together:none",normalize.logbase="0",normalize.params=c(0.0))

#   important: only for summarization step can tmpdir be defined!
expr.mp.rma <- express(expr.qu.rma,"Test3ExprsSum",filedir=datdir,tmpdir="",update=FALSE,
               summarize.method="medianpolish",summarize.select="pmonly",summarize.option="transcript",summarize.logbase="log2",summarize.params=c(10, 0.01, 1.0))

# rma
data.rma <- rma(data.test3,"tmp_Test3RMA",tmpdir="",background="pmonly",normalize=TRUE)

# compare results
expr <- exprs(data.rma)
expr.mp <- exprs(expr.mp.rma)
# plot differences
plot((expr.mp[,1] - expr[,1])/expr[,1], ylim=c(-0.0001,0.0001))

### compute rma with a single call to express()
#   important: for single call to express() tmpdir can be defined!
expr.rma <- express(data.test3,"Test3Exprs",filedir=datdir,tmpdir="",update=FALSE,
            bgcorrect.method="rma",bgcorrect.select="none",bgcorrect.option="pmonly:epanechnikov",bgcorrect.params=c(16384),
            normalize.method="quantile",normalize.select="pmonly",normalize.option="transcript:together:none",normalize.logbase="0",normalize.params=c(0.0),
            summarize.method="medianpolish",summarize.select="pmonly",summarize.option="transcript",summarize.logbase="log2",summarize.params=c(10, 0.01, 1.0))

# compare results
expr <- exprs(data.rma)
expr.mp <- exprs(expr.rma)
# plot differences
plot((expr.mp[,1] - expr[,1])/expr[,1], ylim=c(-0.0001,0.0001))











