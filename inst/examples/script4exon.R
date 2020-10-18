#------------------------------------------------------------------------------#
# Script: step-by-step functions to demonstrate how to analyze exon arrays
#         using the Affymetrix human tissue-mixture exon array dataset
#
# Note: please feel free to copy-paste the examples of interest and adapt the
#       examples to your own needs
#
# Copyright (c) 2007-2011 Christian Stratowa, Vienna, Austria.
# All rights reserved.
#
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

#------------------------------------------------------------------------------#
# 1. step: import Affymetrix chip definition and annotation files into
#          ROOT scheme files
#------------------------------------------------------------------------------#

# HG-U133_Plus_2:
scheme.hgu133plus2 <- import.expr.scheme("hgu133plus2", filedir = file.path(scmdir, "na32"),
                      schemefile = file.path(libdir, "HG-U133_Plus_2.CDF"), 
                      probefile  = file.path(libdir, "HG-U133-PLUS_probe.tab"), 
                      annotfile  = file.path(anndir, "Version11Jul", "HG-U133_Plus_2.na32.annot.csv"))

# HuGene-1_0-st-v1.r4: used as exon array
scheme.hugene10stv1 <- import.exon.scheme("hugene10stv1", filedir = file.path(scmdir, "na32"),
                       file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.clf"),
                       file.path(libdir, "HuGene-1_0-st-v1.r4.analysis-lib-files", "HuGene-1_0-st-v1.r4.pgf"),
                       file.path(anndir, "Version11Jul", "HuGene-1_0-st-v1.na32.hg19.probeset.csv"),
                       file.path(anndir, "Version11Jul", "HuGene-1_0-st-v1.na32.hg19.transcript.csv"))

# HuEx-1_0-st-v2.r2:
scheme.huex10stv2 <- import.exon.scheme("huex10stv2", filedir = file.path(scmdir, "na32"),
                     file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.clf"),
                     file.path(libdir, "HuEx-1_0-st-v2_libraryfile", "HuEx-1_0-st-r2", "HuEx-1_0-st-v2.r2.pgf"),
                     file.path(anndir, "Version11Jul", "HuEx-1_0-st-v2.na32.hg19.probeset.csv"),
                     file.path(anndir, "Version11Jul", "HuEx-1_0-st-v2.na32.hg19.transcript.csv"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# utility functions to demonstrate how to access scheme files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### scheme accessors
rootFile(scheme.huex10stv2r2)
chipName(scheme.huex10stv2r2)
chipType(scheme.huex10stv2r2)
probeInfo(scheme.huex10stv2r2)

### get tree names
unlist(treeNames(scheme.hgu133p2))
getTreeNames(rootFile(scheme.huex10stv2r2))

### browse ROOT scheme files
root.browser(scheme.huex10stv2r2)

### HG-U133_Plus_2: export trees from ROOT scheme file
# export as table only
export(scheme.hgu133p2, treetype="scm", outfile="HGU133Plus2_scm.txt")
export(scheme.hgu133p2, treetype="prb", outfile="HGU133Plus2_prb.txt")
# export as table and import as data.frame
idx <- export(scheme.hgu133p2, treetype="idx", outfile="HGU133Plus2_idx.txt",as.dataframe=TRUE)
head(idx)
ann <- export(scheme.hgu133p2, treetype="ann", outfile="HGU133Plus2_ann.txt",as.dataframe=TRUE)
head(ann)

### HuEx-1_0-st-v2.r2: export trees from ROOT scheme file
# export as table only
export(scheme.huex10stv2r2, treetype="scm", outfile="HuEx10stv2r2_scm.txt")
export(scheme.huex10stv2r2, treetype="idx", outfile="HuEx10stv2r2_idx.txt")
export(scheme.huex10stv2r2, treetype="ann", outfile="HuEx10stv2r2_ann.txt")
export(scheme.huex10stv2r2, treetype="anx", outfile="HuEx10stv2r2_anx.txt")
export(scheme.huex10stv2r2, treetype="anp", outfile="HuEx10stv2r2_anp.txt")


#------------------------------------------------------------------------------#
# 2. step: import CEL-files into ROOT data files
#------------------------------------------------------------------------------#

### new R session: load library xps
library(xps)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tissues from Affymetrix Exon Array Dataset for HG-U133_Plus_2 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
celnames <- c("BreastA","BreastB","BreastC",
              "ProstateA","ProstateB","ProstateC")
# import CEL files
data.mix.u133p2 <- import.data(scheme.u133p2, "HuTissuesU133P2", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tissues from Affymetrix Exon Array Dataset for HuGene-1_0-st-v1
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
# Tissues from Affymetrix Exon Array Dataset for HuEx-1_0-st-v2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# same R session

# directory containing Tissues CEL files
celdir <- "/Volumes/GigaDrive/ChipData/Exon/HuMixture"

### HuEx-1_0-st-v2 data: import raw data
# first, import ROOT scheme file
scheme.exon <- root.scheme(paste(scmdir,"huex10stv2.root",sep="/"))

# subset of CEL files to import
celfiles <- c("huex_wta_breast_A.CEL","huex_wta_breast_B.CEL","huex_wta_breast_C.CEL",
              "huex_wta_prostate_A.CEL","huex_wta_prostate_B.CEL","huex_wta_prostate_C.CEL")
# rename CEL files
celnames <- c("BreastA","BreastB","BreastC",
              "ProstateA","ProstateB","ProstateC")
# import CEL files
data.mix.exon <- import.data(scheme.exon, "HuTissuesExon", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# demonstration how to access the data, and plot the data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

# import ROOT scheme files
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
scheme.u133p2 <- root.scheme(paste(scmdir,"hgu133plus2.root",sep="/"))
scheme.genome <- root.scheme(paste(scmdir,"hugene10stv1.root",sep="/"))
scheme.exon   <- root.scheme(paste(scmdir,"huex10stv2.root",sep="/"))

# import ROOT data files
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.u133p2 <- root.data(scheme.u133p2, paste(datdir,"HuTissuesU133P2_cel.root",sep="/"))
data.genome <- root.data(scheme.genome, paste(datdir,"HuTissuesGenome_cel.root",sep="/"))
data.exon   <- root.data(scheme.exon, paste(datdir,"HuTissuesExon_cel.root",sep="/"))


### plot raw data for HG-U133_Plus_2
# plots
hist(data.u133p2)
image(data.u133p2)
boxplot(data.u133p2, which="userinfo:fIntenQuant")

# plots - alternative possibility:
# select File->SaveAs...->nnn.png to save to disk
root.density(data.u133p2)
root.image(data.u133p2, treename="BreastA.cel")


### plot raw data for HuGene-1_0-st-v1
# plots
hist(data.genome)
boxplot(data.genome, which="userinfo:fIntenQuant")
image(data.genome, col=rainbow(32), names="Breast01.cel")

# plots - alternative possibility:
# select File->SaveAs...->nnn.png to save to disk
root.density(data.genome)
root.image(data.genome, treename="Breast01.cel")


### plot raw data for HuEx-1_0-st-v2
# plots
names <- unlist(treeNames(data.exon))
image(data.exon, names=names[1], add.legend=TRUE)
hist(data.exon)
hist(data.exon, which="core")
boxplot(data.exon, which="userinfo:fIntenQuant")

# plots - alternative possibility:
# select File->SaveAs...->nnn.png to save to disk
root.density(data.exon)
root.image(data.exon, treename="BreastA.cel")


#------------------------------------------------------------------------------#
# 3. step: convert raw data to expression levels
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tissues from Affymetrix Exon Array Dataset for HG-U133_Plus_2 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
scheme.u133p2 <- root.scheme(paste(scmdir,"hgu133plus2.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.u133p2 <- root.data(scheme.u133p2, paste(datdir,"HuTissuesU133P2_cel.root",sep="/"))

### preprocess raw data ###
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/u133p2"

# 1. RMA
data.rma <- rma(data.u133p2,"MixU133P2RMA",filedir=datdir,tmpdir="",
                background="pmonly",normalize=TRUE)
tmpdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/temp"
data.rma <- rma(data.u133p2,"MixU133P2RMAtmp",filedir=datdir,tmpdir=tmpdir,
                background="pmonly",normalize=TRUE)

# 2. MAS5
# to store all trees (including e.g. background trees) in same ROOT file, use "update=TRUE"
data.mas5 <- mas5(data.u133p2,"MixU133P2MAS5All",filedir=datdir,tmpdir="",
                  normalize=TRUE,sc=500, update=TRUE)

# 3. MAS5 detection call
call.mas5 <- mas5.call(data.u133p2,"MixU133P2Call",filedir=datdir,tmpdir="")

# get data.frames
expr.rma <- validData(data.rma)
expr.mas5 <- validData(data.mas5)
pval.mas5 <- pvalData(call.mas5)
pres.mas5 <- presCall(call.mas5)

### plot results ###

# compare mas5 to rma
plot(expr.rma[,1],expr.mas5[,1],log="xy",xlim=c(1,20000),ylim=c(1,20000))

# density plots
hist(data.rma)
hist(data.mas5)

# boxplots
boxplot(data.rma)
boxplot(data.mas5)

# relative boxplots
mboxplot(data.rma, ylim=c(-4,5))
mboxplot(data.mas5, ylim=c(-4,5))

# M vs A plots
mvaplot.dev(data.rma, pch=20, ylim=c(-6,6), names="BreastA.mdp_LEVEL")
mvaplot.dev(data.mas5, pch=20, names="BreastA.tmn_LEVEL")

# present call plots
callplot(call.mas5)
callplot(call.mas5, beside=FALSE, ylim=c(0,125))

# save image
rm(data.u133p2,datdir,scheme.u133p2,scmdir)
save.image(file="HuTissues.U133P2.Rdata");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tissues from Affymetrix Exon Array Dataset for HuGene-1_0-st-v1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
scheme.genome <- root.scheme(paste(scmdir,"hugene10stv1.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.genome <- root.data(scheme.genome, paste(datdir,"HuTissuesGenome_cel.root",sep="/"))


### preprocess raw data ###
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/hugene"

# 1. RMA
data.g.rma <- rma(data.genome,"HuGeneMixRMAMetacore",filedir=datdir,tmpdir="",
                  background="antigenomic",normalize=TRUE,exonlevel="metacore+affx")

# 2. MAS5
data.g.mas5 <- mas5(data.genome,"HuGeneMixMAS5Metacore",filedir=datdir,tmpdir="",
                    normalize=TRUE,sc=500,exonlevel="metacore+affx")
# to store all trees (including e.g. background trees) in same ROOT file, use "update=T"
data.g.mas5 <- mas5(data.genome,"HuGeneMixMAS5MetacoreAll",filedir=datdir,tmpdir="",
                    normalize=TRUE,sc=500,exonlevel="metacore+affx", update=TRUE)

# 3. MAS5 detection call (yes, this is possible for genome/exon arrays)
# note: alpha1 and alpha2 need to be adjusted to get usable P/M/A calls
call.g.mas5 <- mas5.call(data.genome,"HuGeneMixCallMetacore",filedir=datdir,tmpdir="",
                         exonlevel="metacore+affx")

# 4. DABG detection call
# note: alpha1 and alpha2 need to be adjusted to get usable P/M/A calls for transcripts
call.g.dabg <- dabg.call(data.genome,"HuGeneMixDABGMetacore",filedir=datdir,
                         exonlevel="metacore+affx")
#call.g.dabg <- dabg.call(data.genome,"HuGeneMixDABGMetacore", alpha1=???,alpha2=???)

# get data.frames
expr.g.rma <- validData(data.g.rma)
expr.g.mas5 <- validData(data.g.mas5)
pval.g.mas5 <- pvalData(call.g.mas5)
pres.g.mas5 <- presCall(call.g.mas5)
pval.g.dabg <- pvalData(call.g.dabg)
pres.g.dabg <- presCall(call.g.dabg)


### plot results ###

# compare mas5 to rma
plot(expr.g.rma[,1],expr.g.mas5[,1])
plot(expr.g.rma[,1],expr.g.mas5[,1],log="xy",xlim=c(1,20000),ylim=c(1,20000))

# density plots
hist(data.g.rma)
hist(data.g.mas5)

# boxplots
boxplot(data.g.rma)
boxplot(data.g.mas5)

# relative boxplots
mboxplot(data.g.rma, ylim=c(-3,3))
mboxplot(data.g.mas5, ylim=c(-4,5))

# M vs A plots
mvaplot(data.g.rma, pch=20, ylim=c(-4,4))
mvaplot.dev(data.g.rma, pch=20, ylim=c(-6,6))
mvaplot.dev(data.g.rma, pch=20, ylim=c(-6,6), names="Breast01.mdp_LEVEL")

mvaplot(data.g.mas5, pch=20)
mvaplot.dev(data.g.mas5, pch=20)
mvaplot.dev(data.g.mas5, pch=20, names="Breast01.tmn_LEVEL")

# present call plots
callplot(call.g.mas5)
callplot(call.g.dabg)

# save image
rm(data.genome,datdir,scheme.genome,scmdir)
save.image(file="HuTissues.Genome.Rdata");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tissues from Affymetrix Exon Array Dataset for HuEx-1_0-st-v2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
scheme.exon <- root.scheme(paste(scmdir,"huex10stv2.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.exon <- root.data(scheme.exon, paste(datdir,"HuTissuesExon_cel.root",sep="/"))

### preprocess raw data ###
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/exon"

# 1. RMA
# transcript: metacore
# ok for 6 exon arrays in RAM
data.x.rma <- rma(data.exon,"MixRMAMetacore",filedir=datdir,tmpdir="",background="antigenomic",
                  normalize=TRUE,option="transcript",exonlevel="metacore")
# for many exon arrays you may decide to use tmpdir (see helpfile for more information)
tmpdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/temp"
data.x.rma.tmp <- rma(data.exon,"MixRMAtmpMetacore",filedir=datdir,tmpdir=tmpdir,background="antigenomic",
                      normalize=TRUE,option="transcript",exonlevel="metacore")
# probeset: metacore
data.x.rma.ps <- rma(data.exon,"MixRMAMetacorePS",filedir=datdir,tmpdir="",background="antigenomic",
                     normalize=TRUE,option="probeset",exonlevel="metacore")

# 2. MAS5
# to store all trees (including e.g. background trees) in same ROOT file, use "update=T"
data.x.mas5 <- mas5(data.exon,"MixExonMAS5MetacoreAll",filedir=datdir,tmpdir="",
                    normalize=TRUE,sc=500,option="transcript",exonlevel="metacore", update=TRUE)

# 3. MAS5 detection call (yes, this is possible for exon arrays)
# note: alpha1 and alpha2 need to be adjusted to get usable P/M/A calls
call.x.mas5 <- mas5.call(data.exon,"MixCallMetacore",filedir=datdir,tmpdir="",
                         option="transcript",exonlevel="metacore")

# 4. DABG detection call
# transcript: metacore
call.x.dabg <- dabg.call(data.exon,"MixDABGMetacore",filedir=datdir,option="transcript",exonlevel="metacore")
# note: alpha1 and alpha2 need to be adjusted to get usable P/M/A calls for transcripts
#call.x.dabg <- dabg.call(data.exon,"MixDABGMetacore", alpha1=???,alpha2=???)
# probeset: metacore
call.x.dabg.ps <- dabg.call(data.exon,"MixDABGMetacorePS",filedir=datdir,option="probeset",exonlevel="metacore")

# get data.frames
expr.x.rma <- validData(data.x.rma)
expr.x.mas5 <- validData(data.x.mas5)
pval.x.mas5 <- pvalData(call.x.mas5)
pres.x.mas5 <- presCall(call.x.mas5)
pval.x.dabg <- pvalData(call.x.dabg)
pres.x.dabg <- presCall(call.x.dabg)

### plot results ###

# compare mas5 to rma
plot(expr.x.rma[,1],expr.x.mas5[,1],log="xy",xlim=c(1,20000),ylim=c(1,20000))

# density plots
hist(data.x.rma)
hist(data.x.mas5)

# boxplots
boxplot(data.x.rma)
boxplot(data.x.mas5)

# relative boxplots
mboxplot(data.x.rma, ylim=c(-4,5))
mboxplot(data.x.mas5, ylim=c(-4,5))

# M vs A plots
mvaplot.dev(data.x.rma, pch=20, ylim=c(-6,6), names="BreastA.mdp_LEVEL")
mvaplot.dev(data.x.mas5, pch=20, names="BreastA.tmn_LEVEL")

# present call plots
callplot(call.x.mas5)
callplot(call.x.mas5, beside=F, ylim=c(0,125))
callplot(call.x.dabg)
callplot(call.x.dabg, beside=F, ylim=c(0,125))

# save image
rm(data.exon,datdir,scheme.exon,scmdir)
save.image(file="HuTissues.Exon.Rdata");


#------------------------------------------------------------------------------#
# UPDATE: import more CEL-files into ROOT data files
#------------------------------------------------------------------------------#

### new R session: load library xps
library(xps)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tissues from Affymetrix Exon Array Dataset for HG-U133_Plus_2 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### define directories:
# directory of ROOT scheme files
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
# directory containing Tissues CEL files
celdir <- "/Volumes/GigaDrive/ChipData/Exon/HuMixture"
# directory to store ROOT raw data files
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"

### HG-U133_Plus_2 data: import raw data
# import ROOT scheme file
scheme.u133p2 <- root.scheme(paste(scmdir,"hgu133plus2.root",sep="/"))

# import current ROOT data file
data.u133p2 <- root.data(scheme.u133p2, paste(datdir,"HuTissuesU133P2_cel.root",sep="/"))

# subset of CEL files to import
celfiles <- c("u1332plus_ivt_cerebellum_A.CEL","u1332plus_ivt_cerebellum_B.CEL","u1332plus_ivt_cerebellum_C.CEL",
              "u1332plus_ivt_heart_A.CEL","u1332plus_ivt_heart_B.CEL","u1332plus_ivt_heart_C.CEL",
              "u1332plus_ivt_kidney_A.CEL","u1332plus_ivt_kidney_B.CEL","u1332plus_ivt_kidney_C.CEL",
              "u1332plus_ivt_liver_A.CEL","u1332plus_ivt_liver_B.CEL","u1332plus_ivt_liver_C.CEL",
              "u1332plus_ivt_muscle_A.CEL","u1332plus_ivt_muscle_B.CEL","u1332plus_ivt_muscle_C.CEL",
              "u1332plus_ivt_pancreas_A.CEL","u1332plus_ivt_pancreas_B.CEL","u1332plus_ivt_pancreas_C.CEL",
              "u1332plus_ivt_spleen_A.CEL","u1332plus_ivt_spleen_B.CEL","u1332plus_ivt_spleen_C.CEL",
              "u1332plus_ivt_testes_A.CEL","u1332plus_ivt_testes_B.CEL","u1332plus_ivt_testes_C.CEL",
              "u1332plus_ivt_thyroid_A.CEL","u1332plus_ivt_thyroid_B.CEL","u1332plus_ivt_thyroid_C.CEL")
# rename CEL files
celnames <- c("CerebellumA","CerebellumB","CerebellumC",
              "HeartA","HeartB","HeartC",
              "KidneyA","KidneyB","KidneyC",
              "LiverA","LiverB","LiverC",
              "MuscleA","MuscleB","MuscleC",
              "PancreasA","PancreasB","PancreasC",
              "SpleenA","SpleenB","SpleenC",
              "TestesA","TestesB","TestesC",
              "ThyroidA","ThyroidB","ThyroidC")
# import additional CEL-files and update ROOT data file
data.u133p2 <- addData(data.u133p2, celdir=celdir, celfiles=celfiles, celnames=celnames)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tissues from Affymetrix Exon Array Dataset for HuEx-1_0-st-v2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# same R session

### HuEx-1_0-st-v2 data: import raw data
# import ROOT scheme file
scheme.exon <- root.scheme(paste(scmdir,"huex10stv2.root",sep="/"))

# import current ROOT data file
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.exon <- root.data(scheme.exon, paste(datdir,"HuTissuesExon_cel.root",sep="/"))

# subset of CEL files to import
celfiles <- c("huex_wta_cerebellum_A.CEL","huex_wta_cerebellum_B.CEL","huex_wta_cerebellum_C.CEL",
              "huex_wta_heart_A.CEL","huex_wta_heart_B.CEL","huex_wta_heart_C.CEL",
              "huex_wta_kidney_A.CEL","huex_wta_kidney_B.CEL","huex_wta_kidney_C.CEL",
              "huex_wta_liver_A.CEL","huex_wta_liver_B.CEL","huex_wta_liver_C.CEL",
              "huex_wta_muscle_A.CEL","huex_wta_muscle_B.CEL","huex_wta_muscle_C.CEL",
              "huex_wta_pancreas_A.CEL","huex_wta_pancreas_B.CEL","huex_wta_pancreas_C.CEL",
              "huex_wta_spleen_A.CEL","huex_wta_spleen_B.CEL","huex_wta_spleen_C.CEL",
              "huex_wta_testes_A.CEL","huex_wta_testes_B.CEL","huex_wta_testes_C.CEL",
              "huex_wta_thyroid_A.CEL","huex_wta_thyroid_B.CEL","huex_wta_thyroid_C.CEL")
# rename CEL files
celnames <- c("CerebellumA","CerebellumB","CerebellumC",
              "HeartA","HeartB","HeartC",
              "KidneyA","KidneyB","KidneyC",
              "LiverA","LiverB","LiverC",
              "MuscleA","MuscleB","MuscleC",
              "PancreasA","PancreasB","PancreasC",
              "SpleenA","SpleenB","SpleenC",
              "TestesA","TestesB","TestesC",
              "ThyroidA","ThyroidB","ThyroidC")
# import additional CEL-files and update ROOT data file
data.exon <- addData(data.exon, celdir=celdir, celfiles=celfiles, celnames=celnames)


# plots - alternative possibility to avoid memory problems:
# select File->SaveAs...->nnn.png to save to disk
root.density(data.exon, "*", w=400, h=400)
root.image(data.exon, treename="BreastA.cel", w=400, h=400)
root.hist1D(data.exon, "BreastA.cel", type="density", option="COLZ", w=400, h=400)
root.hist2D(data.exon, "BreastA.cel", "BreastB.cel", option="COLZ", w=400, h=400)
root.hist3D(data.exon, "BreastA.cel", "BreastB.cel", "BreastC.cel", option="SCAT", w=400, h=400)
root.graph2D(data.exon, "BreastA.cel", "BreastB.cel", w=400, h=400)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# RMA for Tissues from Affymetrix Exon Array Dataset for HG-U133_Plus_2 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
scheme.u133p2 <- root.scheme(paste(scmdir,"hgu133plus2.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.u133p2 <- root.data(scheme.u133p2, paste(datdir,"HuTissuesU133P2_cel.root",sep="/"))

### RMA
data.rma <- rma(data.u133p2,"MixU133P2RMA",tmpdir="",background="pmonly",normalize=TRUE)

# get data.frames
expr.rma <- validData(data.rma)
exprs.rma <- exprs(data.rma)

### plot results ###
# density plots
hist(data.rma)

# boxplots
boxplot(data.rma,las=2)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# RMA for Tissues from Affymetrix Exon Array Dataset for HuEx-1_0-st-v2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
scheme.exon <- root.scheme(paste(scmdir,"huex10stv2.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.exon <- root.data(scheme.exon, paste(datdir,"HuTissuesExon_cel.root",sep="/"))

### preprocess raw data ###
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/exon"

# 1. RMA
# transcript: metacore
data.x.rma <- rma(data.exon,"MixRMAMetacore",filedir=datdir,tmpdir="",background="antigenomic",
                  normalize=TRUE,option="transcript",exonlevel="metacore")

# root plots
root.density(data.x.rma, "*", w=400, h=400)
root.hist1D(data.x.rma, "BreastA.mdp", type="density", option="COLZ", w=400, h=400)
root.hist2D(data.x.rma, "BreastA.mdp", "BreastB.mdp", option="COLZ", w=400, h=400)
root.hist2D(data.x.rma, "BreastA.mdp", "BreastB.mdp", option="SURF2", w=400, h=400)
root.hist2D(data.x.rma, "BreastA.mdp", "BreastB.mdp", option="SURF3", w=400, h=400)
root.hist3D(data.x.rma, "BreastA.mdp", "BreastB.mdp", "BreastC.mdp", option="SCAT", w=400, h=400)
root.graph2D(data.x.rma, "BreastA.mdp", "BreastB.mdp", w=400, h=400)
root.mvaplot(data.x.rma, "BreastA.mdp", "BreastB.mdp", w=400, h=400)


#------------------------------------------------------------------------------#
# 4. step: apply filters to expression levels
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tissues from Affymetrix Exon Array Dataset for HuEx-1_0-st-v2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
scheme.exon <- root.scheme(file.path(scmdir, "huex10stv2.root"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.exon   <- root.data(scheme.exon, file.path(datdir, "HuTissuesExon_cel.root"))

### RMA for transcription and probeset
rma.tc.core <- rma(data.exon, "MixExonRMAcore", filedir=outdir, tmpdir="", background="antigenomic",
                  normalize=TRUE, option="transcript", exonlevel="affx+core")
rma.ps.core <- rma(data.exon, "MixExonRMAcorePS", filedir=outdir, tmpdir="", background="antigenomic",
                  normalize=TRUE, option="probeset", exonlevel="affx+core")

### apply non-specific filters, e.g. mad
# create PreFilter
prefltr <- PreFilter(mad=c(0.5,0.01))
# apply prefilter to rma
pfr.tc.core <- prefilter(rma.tc.core, filename="PrefilterExonCore", filedir=getwd(), filter=prefltr)
pfr.ps.core <- prefilter(rma.ps.core, filename="PrefilterExonCorePS", filedir=getwd(), filter=prefltr)

# test
pre <- validData(pfr.tc.core)
pre <- pre[pre[, "FLAG"] == 1,]

### apply univariate filters
# create UniFilter
unifltr <- UniFilter(unitest=c("t.test","two.sided","none",0,0.0,FALSE,0.95,TRUE),
                     foldchange=c(1.3,"both"), unifilter=c(0.1,"pval"))

# a, apply unifilter only
ufr.tc.core <- unifilter(rma.tc.core, filename="UnifilterExonCore", filedir=getwd(), filter=unifltr,
                         group=c("Breast","Breast","Breast","Prostate","Prostate","Prostate"))
ufr.ps.core <- unifilter(rma.ps.core, filename="UnifilterExonCorePS", filedir=getwd(), filter=unifltr,
                         group=c("Breast","Breast","Breast","Prostate","Prostate","Prostate"))

# b, apply unifilter to pre-filtered data
puf.tc.core <- unifilter(rma.tc.core, filename="tmp_PreUnifilterExonCoreTC", filedir=getwd(), filter=unifltr,
                         group=c("Breast","Breast","Breast","Prostate","Prostate","Prostate"),xps.fltr=pfr.tc.core)
puf.s.core <- unifilter(rma.ps.core, filename="tmp_PreUnifilterExonCorePS", filedir=getwd(), filter=unifltr,
                         group=c("Breast","Breast","Breast","Prostate","Prostate","Prostate"),xps.fltr=pfr.ps.core)


### advanced example: include present calls
## DABG detection call
dabg.tc <- dabg.call(data.exon, "MixDABGcoreTC", filedir=outdir, option="transcript", exonlevel="affx+core")
dabg.ps <- dabg.call(data.exon, "MixDABGcorePS", filedir=outdir, option="probeset", exonlevel="affx+core")

### prefilter: mad and dabg.call
prefltr <- PreFilter(mad=c(0.5,0.01), prescall=c(0.002, 6,"samples"))
pfr.tc.dabg <- prefilter(rma.tc.core, filename="tmp_PrefilterExonCore", filedir=getwd(), filter=prefltr, xps.call=dabg.tc)

pre <- validData(pfr.tc.dabg)
pre <- pre[pre[, "FLAG"] == 1,]
dim(pre)

flag <- validData(pfr.tc.dabg)
flag <- flag[, "FLAG"] == 1

pval <- validData(callTreeset(pfr.tc.dabg))
pval <- pval[flag,]
range(pval)

pres <- presCall(callTreeset(pfr.tc.dabg))
pres <- pres[flag,]
rownames(pres) <- pres[,"UnitName"]
pres <- pres[,3:ncol(pres)]

### unifilter for pre-filtered data
unifltr <- UniFilter(unitest=c("t.test","two.sided","none",0,0.0,FALSE,0.95,TRUE),
                     foldchange=c(1.5,"both"), unifilter=c(0.001,"pval"))
puf.tc.dabg <- unifilter(rma.tc.core, filename="tmp_PreUnifilterExonDABG", filedir=getwd(), filter=unifltr,
                         group=c("Breast","Breast","Breast","Prostate","Prostate","Prostate"),xps.fltr=pfr.tc.dabg)

puf <- validData(puf.tc.dabg, which = "UnitName")
dim(puf)

ufr <- validFilter(puf.tc.dabg)
tmp <- puf[ufr[,"FLAG"] == 1,]
dim(tmp)

### export all data
export.filter(puf.tc.dabg, treetype="stt", outfile="PreUnifilterExonDABG_stt.txt")

### export only data with flag=1
export.filter(puf.tc.dabg, treetype="stt", varlist="fUnitName:fTranscriptID:fSymbol:mn1:mn2:fc:pval:mask", outfile="PreUnifilterExonDABG_stt_mask.txt")



#------------------------------------------------------------------------------#
# Demonstrations of advanced methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    Note: ROOT scheme files and ROOT raw data files are usually already stored
#          in special system directories. When a new R session is created for the
#          first time, they must fist be loaded using "root.scheme()" and "root.data()".
#          However, this is not necessary when re-opening a saved R session later.
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tissues from Affymetrix Exon Array Dataset for HuEx-1_0-st-v2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### new R session: load library xps
library(xps)

### first, load ROOT scheme file and ROOT data file
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes/na32"
scheme.exon <- root.scheme(paste(scmdir,"huex10stv2.root",sep="/"))
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
data.exon <- root.data(scheme.exon, paste(datdir,"HuTissuesExon_cel.root",sep="/"))

### preprocess raw data ###
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/exon"


### RMA transcript: metacore
data.x.rma <- rma(data.exon,"HuExonRMAMetacore",filedir=datdir,tmpdir="",background="antigenomic",
                  normalize=T,option="transcript",exonlevel="metacore+affx")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# demonstration 1: compute RMA stepwise 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### 1.step: background - rma
data.bg.rma <- bgcorrect.rma(data.exon, "HuExonRMABgrd", filedir=datdir, 
               select="antigenomic", exonlevel="metacore+affx")
# or:
data.bg.rma <- bgcorrect(data.exon, "HuExonRMABgrd", filedir=datdir,  
               method="rma", select="antigenomic", option="pmonly:epanechnikov",
               params=c(16384), exonlevel="metacore+affx")

### 2step: normalization - quantile
data.qu.rma <- normalize.quantiles(data.bg.rma, "HuExonRMANorm", filedir=datdir , 
               exonlevel="metacore+affx")
# or:
data.qu.rma <- normalize(data.bg.rma, "HuExonRMANorm", filedir=datdir,
               method="quantile", select="pmonly", option="transcript:together:none", 
               logbase="0", params=c(0.0), exonlevel="metacore+affx")

### 3.step: summarization - medpol
data.mp.rma <- summarize.rma(data.qu.rma, "HuExonRMASum", filedir=datdir, tmpdir="", 
               exonlevel="metacore+affx")
# or:
data.mp.rma <- summarize(data.qu.rma, "HuExonRMASum", filedir=datdir, tmpdir="", 
               method="medianpolish", select="pmonly", option="transcript", 
               logbase="log2", params=c(10, 0.01, 1.0), exonlevel="metacore+affx")

# compare results
expr.x  <- exprs(data.x.rma)
expr.mp <- exprs(data.mp.rma)
# plot differences
plot((expr.mp[,1] - expr.x[,1])/expr.x[,1], ylim=c(-0.0001,0.0001))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# demonstration 2: compute RMA using function "express()"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### compute rma stepwise
#   important: for stepwise computation tmpdir must be tmpdir="" otherwise the root file will be empty!
expr.bg.rma <- express(data.exon, "HuExonExprsBgrd", filedir=datdir, tmpdir="", update=FALSE,
               bgcorrect.method="rma", bgcorrect.select="antigenomic",
               bgcorrect.option="pmonly:epanechnikov", bgcorrect.params=c(16384),
               exonlevel="metacore+affx")

#   important: for stepwise computation tmpdir must be tmpdir="" otherwise the root file will be empty!
expr.qu.rma <- express(expr.bg.rma, "HuExonExprsNorm", filedir=datdir, tmpdir="", update=FALSE,
               normalize.method="quantile", normalize.select="pmonly",
               normalize.option="transcript:together:none", normalize.logbase="0",
               normalize.params=c(0.0), exonlevel="metacore+affx")

#   important: only for summarization step can tmpdir be defined!
expr.mp.rma <- express(expr.qu.rma, "HuExonExprsSum", filedir=datdir, tmpdir="", update=FALSE,
               summarize.method="medianpolish", summarize.select="pmonly", 
               summarize.option="transcript", summarize.logbase="log2", 
               summarize.params=c(10, 0.01, 1.0), exonlevel="metacore+affx")

# compare results
expr.x  <- exprs(data.x.rma)
expr.mp <- exprs(expr.mp.rma)
# plot differences
plot((expr.mp[,1] - expr.x[,1])/expr.x[,1], ylim=c(-0.0001,0.0001))

### compute rma with a single call to express()
expr.rma <- express(data.exon,"HuExonExprs",filedir=datdir,tmpdir="",update=FALSE,
            bgcorrect.method="rma",bgcorrect.select="antigenomic",bgcorrect.option="pmonly:epanechnikov",bgcorrect.params=c(16384),
            normalize.method="quantile",normalize.select="pmonly",normalize.option="transcript:together:none",normalize.logbase="0",normalize.params=c(0.0),
            summarize.method="medianpolish",summarize.select="pmonly",summarize.option="transcript",summarize.logbase="log2",summarize.params=c(10, 0.01, 1.0),
            exonlevel="metacore+affx")

# compare results
expr.x <- exprs(data.x.rma)
expr   <- exprs(expr.rma)
# plot differences
plot((expr[,1] - expr.x[,1])/expr.x[,1], ylim=c(-0.0001,0.0001))










# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# demonstration 3: diverse tests
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# compute background for antigenomic probes
# a, using bgcorrect()
data.bg.gc <- bgcorrect.gc(data.exon,"tmp_HuExonGCBgrd",filedir=datdir,tmpdir="", exonlevel="metacore+affx")

# b, using express()
expr.bg.gc <- express(data.exon, "HuExonExprGCBgrd", filedir=datdir, tmpdir="", update=F,
              bgcorrect.method="gccontent", bgcorrect.select="antigenomic",
              bgcorrect.option="attenuatebg", bgcorrect.params=c(0.4, 0.005, -1.0),
              exonlevel="metacore+affx")




##############################

# 1. RMA
# transcript
data.rma.ts8192 <- rma(data.exon,"MixRMAMetacore",tmpdir="",background="antigenomic",
                       normalize=T,option="transcript",exonlevel="metacore")
data.rma.ts9216 <- rma(data.exon,"MixRMACore",tmpdir="",background="antigenomic",
                       normalize=T,option="transcript",exonlevel="core")
data.rma.ts8252 <- rma(data.exon,"MixRMAMetacoreAffx",tmpdir="",background="antigenomic",
                       normalize=T,option="transcript",exonlevel="metacore+affx")

# probeset
data.rma.ps8192 <- rma(data.exon,"MixRMAMetacorePS",tmpdir="",background="antigenomic",
                       normalize=T,option="probeset",exonlevel="metacore")
data.rma.ps9216 <- rma(data.exon,"MixRMACorePS",tmpdir="",background="antigenomic",
                       normalize=T,option="probeset",exonlevel="core")
data.rma.ps8252 <- rma(data.exon,"MixRMAMetacoreAffxPS",tmpdir="",background="antigenomic",
                       normalize=T,option="probeset",exonlevel="metacore+affx")

# 2. MAS5
data.mas5.ts9216T <- mas5(data.exon,"MixExonMAS5Core500",tmpdir="",
                         normalize=T,sc=500,option="transcript",exonlevel="core", update=T)
data.mas5.ts9216F <- mas5(data.exon,"MixExonMAS5Core",tmpdir="",
                         normalize=F,sc=500,option="transcript",exonlevel="core", update=T)

# 3. MAS5 detection call (yes, this is possible for exon arrays)
call.mas5.ts9216 <- mas5.call(data.exon,"MixCallCore",tmpdir="",option="transcript",exonlevel="core")

# 4. DABG detection call
# transcript: metacore
call.dabg.ts8192 <- dabg.call(data.exon,"MixDABGMetacore",option="transcript",exonlevel="metacore")
call.dabg.ts9216 <- dabg.call(data.exon,"MixDABGCore",option="transcript",exonlevel="core")
# probeset: metacore
call.dabg.ps9216 <- dabg.call(data.exon,"MixDABGCorePS",option="probeset",exonlevel="core")



#------------------------------------------------------------------------------#


