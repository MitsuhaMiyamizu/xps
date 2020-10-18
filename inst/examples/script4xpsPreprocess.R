#------------------------------------------------------------------------------#
# Script: step-by-step functions to demonstrate how to use function express()
#
# Note: please feel free to copy-paste the examples of interest and adapt the
#       examples to your own needs
#
# Copyright (c) 2009-2010 Christian Stratowa, Vienna, Austria.
# All rights reserved.
#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Tissues from Affymetrix Exon Array Dataset for HG-U133_Plus_2 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Software versions used:
# R:    R-2.9.1
# xps:  xps_1.5.14
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# create scheme and data files 
#------------------------------------------------------------------------------#

### define directories:
libdir <- "/Volumes/GigaDrive/Affy/libraryfiles"
anndir <- "/Volumes/GigaDrive/Affy/Annotation/Version09Mar"
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
celdir <- "/Volumes/GigaDrive/ChipData/Exon/hutissues"
outdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/u133p2"

### new R session: import scheme and data
library(xps)

## scheme for HG-U133_Plus_2:
scheme.u133p2 <- import.expr.scheme("Scheme_HGU133p2_na28", filedir=scmdir,
                 schemefile=paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),
                 probefile=paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),
                 annotfile=paste(anndir,"HG-U133_Plus_2.na28.annot.csv",sep="/"))

## subset of CEL files to import
celfiles <- c("u1332plus_ivt_breast_A.CEL", "u1332plus_ivt_breast_B.CEL", "u1332plus_ivt_breast_C.CEL",
              "u1332plus_ivt_prostate_A.CEL", "u1332plus_ivt_prostate_B.CEL", "u1332plus_ivt_prostate_C.CEL")
celnames <- c("BreastA", "BreastB", "BreastC", "ProstateA", "ProstateB", "ProstateC")
data.u133p2 <- import.data(scheme.hgu133p2, "HuTissuesU133P2", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)

## some plots: need to attach data first

## attach mask and data
data.u133p2 <- attachMask(data.u133p2)
data.u133p2 <- attachInten(data.u133p2)

## plots
hist(data.u133p2)
hist(data.u133p2, which="pm")
boxplot(data.u133p2)
boxplot(data.u133p2, which="pm")
boxplot.dev(data.u133p2, which="pm", dev="png", outfile="BoxPlot_DataU133P2", mar=c(6,3,1,1), w=360, h=320)

png(file="DensityPlot_DataU133P2.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(data.u133p2, which="pm")
dev.off()

## remove all data
data.u133p2 <- removeMask(data.u133p2)
data.u133p2 <- removeInten(data.u133p2)
gc()


#------------------------------------------------------------------------------#
# compute RMA using rma() or express() 
#------------------------------------------------------------------------------#

### new R session: rma
library(xps)

## import ROOT scheme and data files
scheme.u133p2 <- root.scheme(paste(scmdir,"Scheme_HGU133p2_na28.root",sep="/"))
data.u133p2   <- root.data(scheme.u133p2, paste(datdir,"HuTissuesU133P2_cel.root",sep="/"))

## RMA: 
data.rma <- rma(data.u133p2, "tmp_MixU133P2RMA", filedir=outdir, tmpdir="",
                background="pmonly", normalize=TRUE)

## get data.frames for rma
xps.rma  <- validData(data.rma)
head(xps.rma)


## express: rma 
expr.rma <- express(data.u133p2, "tmp_U133P2Exprs", filedir=outdir, tmpdir="", update=FALSE,
            bgcorrect.method="rma", bgcorrect.select="none", bgcorrect.option="pmonly:epanechnikov", bgcorrect.params=c(16384),
            normalize.method="quantile", normalize.select="pmonly", normalize.option="transcript:together:none", normalize.logbase="0", normalize.params=c(0.0),
            summarize.method="medianpolish", summarize.select="pmonly", summarize.option="transcript", summarize.logbase="log2", summarize.params=c(10, 0.01, 1.0))

## xpsPreprocess: since express is only wrapper function
expr.rma <- xpsPreprocess(data.u133p2,"tmp_U133P2Exprs",filedir=outdir,tmpdir="",update=FALSE,
            bgcorrect.method="rma", bgcorrect.select="none", bgcorrect.option="pmonly:epanechnikov", bgcorrect.params=c(16384),
            normalize.method="quantile", normalize.select="pmonly", normalize.option="transcript:together:none", normalize.logbase="0", normalize.params=c(0.0),
            summarize.method="medianpolish", summarize.select="pmonly", summarize.option="transcript", summarize.logbase="log2", summarize.params=c(10, 0.01, 1.0))

## get data.frames for express
xprs.rma  <- validData(expr.rma)
head(xprs.rma)


#------------------------------------------------------------------------------#
# compute RMA using express() stepwise
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. step: background correction: rma
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## rma background: numpoints=16384
bgrd.rma <- express(data.u133p2, "tmp_BgrdRMA", filedir=outdir, tmpdir="", update=FALSE,
            bgcorrect.method="rma", bgcorrect.select="none", bgcorrect.option="pmonly:epanechnikov",
            bgcorrect.params=c(16384))

## attach mask, bgrd and background-corrected intensities
bgrd.rma <- attachMask(bgrd.rma)
bgrd.rma <- attachBgrd(bgrd.rma)
bgrd.rma <- attachInten(bgrd.rma)

## get colnames of bgrd trees
bgrdnames <- colnames(validBgrd(bgrd.rma))
bgrdnames

## images for rma bgrd
image(bgrd.rma, bg=TRUE, transfo=NULL, col=heat.colors(12), names=bgrdnames[1])
image(bgrd.rma, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[1])
image.dev(bgrd.rma, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[1], mar=c(1, 1, 2, 1), dev="screen", w=360, h=360)
image.dev(bgrd.rma, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[1], mar=c(1, 1, 2, 1), dev="png", outfile="Image_BgrdRMA_bgrd1", w=360, h=360)
image.dev(bgrd.rma, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[2], mar=c(1, 1, 2, 1), dev="png", outfile="Image_BgrdRMA_bgrd2", w=300, h=300)
image.dev(bgrd.rma, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[4], mar=c(1, 1, 2, 1), dev="png", outfile="Image_BgrdRMA_bgrd4", w=300, h=300)

## get colnames of bgrd-corrected trees
datanames <- colnames(validData(bgrd.rma))
datanames

## images for background-corrected intensities
image(bgrd.rma, transfo=NULL, col=heat.colors(12), names=datanames[1])
image(bgrd.rma, transfo=log2, col=heat.colors(12), names=datanames[1])
image.dev(bgrd.rma, bg=FALSE, transfo=log2, col=heat.colors(12), names=datanames[1], mar=c(1, 1, 2, 1), dev="png", outfile="Image_BgrdRMA_data1", w=360, h=360)
image.dev(bgrd.rma, bg=FALSE, transfo=log2, col=heat.colors(12), names=datanames[2], mar=c(1, 1, 2, 1), dev="png", outfile="Image_BgrdRMA_data2", w=300, h=300)
image.dev(bgrd.rma, bg=FALSE, transfo=log2, col=heat.colors(12), names=datanames[4], mar=c(1, 1, 2, 1), dev="png", outfile="Image_BgrdRMA_data4", w=300, h=300)

## plots
hist(bgrd.rma)
boxplot(bgrd.rma)
# use pmonly since only these were corrected
hist(bgrd.rma, which="pm")
boxplot(bgrd.rma, which="pm")

boxplot.dev(bgrd.rma, which="pm", dev="png", outfile="BoxPlot_BgrdRMA", mar=c(6,3,1,1), w=360, h=320)

png(file="DensityPlot_BgrdRMA.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(bgrd.rma, which="pm")
dev.off()

## bgrd-corrected intensity: use pmonly since only these were corrected
data <- validData(bgrd.rma, which="pm")
colnames(data) <- namePart(colnames(data))
head(data)

png(file="BoxPlot_BgrdRMA_data.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(data), las=2)
dev.off()

## background: use pmonly since only these were corrected
bgrd <- validBgrd(bgrd.rma, which="pm")
colnames(bgrd) <- namePart(colnames(bgrd))
head(bgrd)

png(file="BoxPlot_BgrdRMA_bgrd.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(bgrd), las=2)
dev.off()

## to avoid memory comsumption of R remove data:
bgrd.rma <- removeMask(bgrd.rma)
bgrd.rma <- removeBgrd(bgrd.rma)
bgrd.rma <- removeInten(bgrd.rma)
gc()

# Note: To avoid the memory problems when plotting data, you can use the corresponding
#       methods "root.drawxxx()", such as:
root.density(bgrd.rma, "*")
root.image(bgrd.rma, "BreastA.int")
root.image(bgrd.rma, "BreastA.rbg", leafname="fBg")
root.profile(bgrd.rma)
root.hist2D(bgrd.rma, "BreastA.int", "BreastB.int", option="COLZ")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. step: probe-level normalization: quantile
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## quantile normalization: trim=0.0
norm.qu <- express(bgrd.rma, "tmp_NormQuan", filedir=outdir, tmpdir="", update=FALSE, 
           normalize.method="quantile", normalize.select="pmonly", normalize.option="transcript:together:none", normalize.logbase="0", normalize.params=c(0.0), 
           verbose=TRUE)

## attach mask and intensity
norm.qu <- attachMask(norm.qu)
norm.qu <- attachInten(norm.qu)

## get colnames of trees
names <- colnames(validData(norm.qu))
names

## plots
hist(norm.qu)
hist(norm.qu, which="pm")

png(file="DensityPlot_NormQuant.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(norm.qu, which="pm")
dev.off()

boxplot(norm.qu)  #too many zeros
boxplot(norm.qu, which="pm")
boxplot.dev(norm.qu, which="pm", dev="screen", mar=c(6,3,1,1))
boxplot.dev(norm.qu, which="pm", dev="png", outfile="BoxPlot_NormQuant", mar=c(6,3,1,1), w=360, h=320)
boxplot.dev(norm.qu, which="pm", dev="jpeg", outfile="BoxPlot_NormQuant", mar=c(6,3,1,1), w=360, h=320)

## use pmonly since only these were normalized
data <- validData(norm.qu, which="pm")
colnames(data) <- namePart(colnames(data))
head(data)

png(file="BoxPlot_NormQuant_data.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(data), las=2)
dev.off()

## scatterplot
png(file="ScatterPlot_NormQuant_PM_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_NormQuant_PM_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()

## remove data
norm.qu <- removeInten(norm.qu)
gc()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 3. step: summarization: medianpolish
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## medpol summarization: maxiter=10, eps=0.01, neglog=1.0
expr.mp <- express(norm.qu, "tmp_ExprMedpol", filedir=outdir, update=FALSE, 
           summarize.method="medianpolish", summarize.select="pmonly", summarize.option="transcript", summarize.logbase="log2", summarize.params=c(10, 0.01, 1.0))

## plots
hist(expr.mp)
boxplot(expr.mp) 

png(file="DensityPlot_ExprMedpol.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(expr.mp)
dev.off()

boxplot.dev(expr.mp, dev="png", outfile="BoxPlot_ExprMedpol", mar=c(6,3,1,1), w=360, h=320)

## get expression level
expr  <- validData(expr.mp)

## scatterplot
png(file="ScatterPlot_ExprMedpol_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(expr[,1]), log2(expr[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_ExprMedpol_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(expr[,1]), log2(expr[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()


#------------------------------------------------------------------------------#
# background correction using express()
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# RMA background correction
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## rma background: numpoints=16384
bgrd.rma <- express(data.u133p2, "tmp_BgrdRMA", filedir=outdir, tmpdir="", update=FALSE,
            bgcorrect.method="rma", bgcorrect.select="none", bgcorrect.option="pmonly:epanechnikov",
            bgcorrect.params=c(16384))

## see above


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAS4 background correction: sector 4x4
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## mas4 background: pcntcells=0.02, secrows=4, seccols=4, smoothiter=0
bgrd.mas4 <- express(data.u133p2, "tmp_BgrdMAS4", filedir=outdir, tmpdir="", update=FALSE, 
             bgcorrect.method="sector", bgcorrect.select="all", bgcorrect.option="subtractbg",
             bgcorrect.params=c(0.02, 4, 4, 0))

## attach mask, bgrd and background-corrected intensities
bgrd.mas4 <- attachMask(bgrd.mas4)
bgrd.mas4 <- attachBgrd(bgrd.mas4)
bgrd.mas4 <- attachInten(bgrd.mas4)

## get colnames of bgrd trees
bgrdnames <- colnames(validBgrd(bgrd.mas4))
bgrdnames

## images for mas4 bgrd
image(bgrd.mas4, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[4])
image.dev(bgrd.mas4, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[4], mar=c(1, 1, 2, 1), dev="png", outfile="Image_BgrdMAS4_bgrd4", w=360, h=360)

## plots
#hist(bgrd.mas4)  #negative values => NA
boxplot(bgrd.mas4)

boxplot.dev(bgrd.mas4, dev="png", outfile="BoxPlot_BgrdMAS4", mar=c(6,3,1,1), w=360, h=320)
boxplot.dev(bgrd.mas4, which="both", dev="png", outfile="BoxPlot_BgrdMAS4_both", mar=c(6,3,1,1), w=360, h=320)

## bgrd-corrected intensity: use both(PM+MM) 
data <- validData(bgrd.mas4, which="both")
colnames(data) <- namePart(colnames(data))
head(data)

png(file="BoxPlot_BgrdMAS4_data.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(data), las=2)
dev.off()

## background
bgrd <- validBgrd(bgrd.mas4)
colnames(bgrd) <- namePart(colnames(bgrd))
head(bgrd)

png(file="BoxPlot_BgrdMAS4_bgrd.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(bgrd), las=2)
dev.off()

## to avoid memory comsumption of R remove data:
bgrd.mas4 <- removeMask(bgrd.mas4)
bgrd.mas4 <- removeBgrd(bgrd.mas4)
bgrd.mas4 <- removeInten(bgrd.mas4)
gc()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAS5 background correction: weightedsector 4x4
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## mas5 background: pcntcells=0.02, secrows=4, seccols=4, smoothiter=0, smooth=100, noisefrac=0.5
bgrd.mas5 <- express(data.u133p2, "tmp_BgrdMAS5", filedir=outdir, tmpdir="", update=FALSE, 
             bgcorrect.method="weightedsector", bgcorrect.select="both", bgcorrect.option="correctbg",
             bgcorrect.params=c(0.02, 4, 4, 0, 100, 0.5))
#             bgcorrect.params=c(0.005, 4, 4, 0, 100, 0.5))

## attach mask, bgrd and background-corrected intensities
bgrd.mas5 <- attachMask(bgrd.mas5)
bgrd.mas5 <- attachBgrd(bgrd.mas5)
bgrd.mas5 <- attachInten(bgrd.mas5)

## get colnames of bgrd trees
bgrdnames <- colnames(validBgrd(bgrd.mas5))
bgrdnames

## images for mas5 bgrd
image(bgrd.mas5, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[4])
image.dev(bgrd.mas5, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[4], mar=c(1, 1, 2, 1), dev="png", outfile="Image_BgrdMAS5_bgrd4", w=360, h=360)

## plots
hist(bgrd.mas5)
boxplot(bgrd.mas5)
hist(bgrd.mas5, which="both")
boxplot(bgrd.mas5, which="both")

boxplot.dev(bgrd.mas5, dev="png", outfile="BoxPlot_BgrdMAS5", mar=c(6,3,1,1), w=360, h=320)
boxplot.dev(bgrd.mas5, which="both", dev="png", outfile="BoxPlot_BgrdMAS5_both", mar=c(6,3,1,1), w=360, h=320)

png(file="DensityPlot_BgrdMAS5.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(bgrd.mas5)
dev.off()

## use pcntcells=0.005
png(file="DensityPlot_BgrdMAS5_pc005.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(bgrd.mas5)
dev.off()

png(file="DensityPlot_BgrdMAS5_both.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(bgrd.mas5, which="both")
dev.off()

## bgrd-corrected intensity: use both(PM+MM) 
data <- validData(bgrd.mas5, which="both")
colnames(data) <- namePart(colnames(data))
head(data)

png(file="BoxPlot_BgrdMAS5_data.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(data), las=2)
dev.off()

## background
bgrd <- validBgrd(bgrd.mas5)
colnames(bgrd) <- namePart(colnames(bgrd))
head(bgrd)

png(file="BoxPlot_BgrdMAS5_bgrd.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(bgrd), las=2)
dev.off()

## to avoid memory comsumption of R remove data:
bgrd.mas5 <- removeMask(bgrd.mas5)
bgrd.mas5 <- removeBgrd(bgrd.mas5)
bgrd.mas5 <- removeInten(bgrd.mas5)
gc()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAS-like background correction: sector 8x8, smoothiter=3, correctbg
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## mas background: pcntcells=0.02, secrows=8, seccols=8, smoothiter=3
bgrd.mas <- express(data.u133p2, "tmp_BgrdMAS8x8", filedir=outdir, tmpdir="", update=FALSE, 
            bgcorrect.method="sector", bgcorrect.select="both", bgcorrect.option="correctbg",
            bgcorrect.params=c(0.02, 8, 8, 3))

## attach mask, bgrd and background-corrected intensities
bgrd.mas <- attachMask(bgrd.mas)
bgrd.mas <- attachBgrd(bgrd.mas)
bgrd.mas <- attachInten(bgrd.mas)

## get colnames of bgrd trees
bgrdnames <- colnames(validBgrd(bgrd.mas))
bgrdnames

## images for mas bgrd
image(bgrd.mas, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[4])
image.dev(bgrd.mas, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[4], mar=c(1, 1, 2, 1), dev="png", outfile="Image_BgrdMAS8x8x_bgrd4", w=360, h=360)

## plots
hist(bgrd.mas)
boxplot(bgrd.mas)
hist(bgrd.mas, which="both")
boxplot(bgrd.mas, which="both")

boxplot.dev(bgrd.mas, dev="png", outfile="BoxPlot_BgrdMAS8x8", mar=c(6,3,1,1), w=360, h=320)
boxplot.dev(bgrd.mas, which="both", dev="png", outfile="BoxPlot_BgrdMAS8x8_both", mar=c(6,3,1,1), w=360, h=320)

png(file="DensityPlot_BgrdMAS8x8.png",width=360,height=320)
par(mar=c(6,4,1,1));
hist(bgrd.mas)
dev.off()

png(file="DensityPlot_BgrdMAS8x8_both.png",width=360,height=320)
par(mar=c(6,4,1,1));
hist(bgrd.mas, which="both")
dev.off()

## bgrd-corrected intensity: use both(PM+MM) 
data <- validData(bgrd.mas, which="both")
colnames(data) <- namePart(colnames(data))
head(data)

png(file="BoxPlot_BgrdMAS8x8_data.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(data), las=2)
dev.off()

## background
bgrd <- validBgrd(bgrd.mas)
colnames(bgrd) <- namePart(colnames(bgrd))
head(bgrd)

png(file="BoxPlot_BgrdMAS8x8_bgrd.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(bgrd), las=2)
dev.off()

## to avoid memory comsumption of R remove data:
bgrd.mas <- removeMask(bgrd.mas)
bgrd.mas <- removeBgrd(bgrd.mas)
bgrd.mas <- removeInten(bgrd.mas)
gc()


## mas5 background: pcntcells=0.02, secrows=8, seccols=8, smoothiter=3, smooth=100, noisefrac=0.5
bgrd.mas <- express(data.u133p2, "tmp_BgrdMAS5_8x8", filedir=outdir, tmpdir="", update=FALSE, 
            bgcorrect.method="weightedsector", bgcorrect.select="both", bgcorrect.option="correctbg",
            bgcorrect.params=c(0.02, 8, 8, 3, 100, 0.5))

## attach mask, bgrd and background-corrected intensities
bgrd.mas <- attachMask(bgrd.mas)
bgrd.mas <- attachBgrd(bgrd.mas)

## get colnames of bgrd trees
bgrdnames <- colnames(validBgrd(bgrd.mas))
bgrdnames

## images for mas bgrd
image(bgrd.mas, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[4])
image.dev(bgrd.mas, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[4], mar=c(1, 1, 2, 1), dev="png", outfile="Image_BgrdMAS5_8x8_bgrd4", w=360, h=360)

## to avoid memory comsumption of R remove data:
bgrd.mas <- removeMask(bgrd.mas)
bgrd.mas <- removeBgrd(bgrd.mas)
gc()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GC-content background correction
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## gc background: trim=0.4, l=0.005, h=-1.0
bgrd.gc <- express(data.u133p2, "tmp_BgrdGC", filedir=outdir, tmpdir="", update=FALSE, 
           bgcorrect.method="gccontent", bgcorrect.select="none", bgcorrect.option="attenuatebg",
           bgcorrect.params=c(0.4, 0.005, -1.0))

## attach mask, bgrd and background-corrected intensities
bgrd.gc <- attachMask(bgrd.gc)
bgrd.gc <- attachBgrd(bgrd.gc)
bgrd.gc <- attachInten(bgrd.gc)

## get colnames of bgrd trees
bgrdnames <- colnames(validBgrd(bgrd.gc))
bgrdnames

## images for gc bgrd
image(bgrd.gc, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[4])
image.dev(bgrd.gc, bg=TRUE, transfo=log2, col=heat.colors(12), names=bgrdnames[4], mar=c(1, 1, 2, 1), dev="png", outfile="Image_BgrdGC_bgrd4", w=360, h=360)

## get colnames of bgrd-corrected trees
datanames <- colnames(validData(bgrd.gc))
datanames

image.dev(bgrd.gc, bg=FALSE, transfo=log2, col=heat.colors(12), names=datanames[4], mar=c(1, 1, 2, 1), dev="png", outfile="Image_BgrdGC_data4", w=360, h=360)

## plots
hist(bgrd.gc)
boxplot(bgrd.gc)
hist(bgrd.gc, which="pm")
boxplot(bgrd.gc, which="pm")

boxplot.dev(bgrd.gc, dev="png", outfile="BoxPlot_BgrdGC", mar=c(6,3,1,1), w=360, h=320)
boxplot.dev(bgrd.gc, which="pm", dev="png", outfile="BoxPlot_BgrdGC_pm", mar=c(6,3,1,1), w=360, h=320)

png(file="DensityPlot_BgrdGC.png",width=360,height=320)
par(mar=c(6,4,1,1));
hist(bgrd.gc)
dev.off()

png(file="DensityPlot_BgrdGC_pm.png",width=360,height=320)
par(mar=c(6,4,1,1));
hist(bgrd.gc, which="pm")
dev.off()

## bgrd-corrected intensity: use pmonly 
data <- validData(bgrd.gc, which="pm")
colnames(data) <- namePart(colnames(data))
head(data)

png(file="BoxPlot_BgrdGC_data.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(data), las=2)
dev.off()

## background
bgrd <- validBgrd(bgrd.gc, which="pm")
colnames(bgrd) <- namePart(colnames(bgrd))
head(bgrd)

png(file="BoxPlot_BgrdGC_bgrd.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(bgrd), las=2)
dev.off()

## to avoid memory comsumption of R remove data:
bgrd.gc <- removeMask(bgrd.gc)
bgrd.gc <- removeBgrd(bgrd.gc)
bgrd.gc <- removeInten(bgrd.gc)
gc()


#------------------------------------------------------------------------------#
# probe-level normalization using express()
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# RMA quantile: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## quantile normalization: bgrd.rma, trim=0.0
norm.qu <- express(bgrd.rma, "tmp_NormQuan", filedir=outdir, tmpdir="", update=FALSE, 
           normalize.method="quantile", normalize.select="pmonly", normalize.option="transcript:together:none",
           normalize.logbase="0", normalize.params=c(0.0, 1.0))

## see above

## quantile normalization: bgrd.mas, trim=0.0
norm.qu <- express(bgrd.mas, "tmp_NormQuan", filedir=outdir, tmpdir="", update=FALSE, 
           normalize.method="quantile", normalize.select="all", normalize.option="transcript:together:none",
           normalize.logbase="0", normalize.params=c(0.0, 1.0))

## attach mask and intensity
norm.qu <- attachMask(norm.qu)
norm.qu <- attachInten(norm.qu)

## plots
hist(norm.qu)
hist(norm.qu, which="pm")

png(file="DensityPlot_NormQuant_BgrdMAS.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(norm.qu)
dev.off()

boxplot(norm.qu) #ok since "all"
boxplot(norm.qu, which="pm")
boxplot.dev(norm.qu, dev="screen", mar=c(6,3,1,1))
boxplot.dev(norm.qu, dev="png", outfile="BoxPlot_NormQuant_BgrdMAS", mar=c(6,3,1,1), w=360, h=320)

## normalized intensity
data <- validData(norm.qu)
colnames(data) <- namePart(colnames(data))
head(data)

png(file="BoxPlot_NormQuant_BgrdMAS_data.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(data), las=2)
dev.off()

## scatterplot
png(file="ScatterPlot_NormQuant_BgrdMAS_all_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_NormQuant_BgrdMAS_all_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()

## normalized intensity: pm only
data <- validData(norm.qu, which="pm")

## scatterplot
png(file="ScatterPlot_NormQuant_BgrdMAS_PM_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_NormQuant_BgrdMAS_PM_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()

## remove data
norm.qu <- removeInten(norm.qu)
gc()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mean: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## mean: trim=0.0, targetinten=-1
norm.mn <- express(bgrd.mas, "tmp_NormMean", filedir=outdir, tmpdir="", update=FALSE, 
           normalize.method="mean", normalize.select="both", normalize.option="transcript:all",
           normalize.logbase="0", normalize.params=c(0.0, -1))

## attach mask and intensity
norm.mn <- attachMask(norm.mn)
norm.mn <- attachInten(norm.mn)

## plots: include reference tree - difference caused by scaling to 500
hist(norm.mn)
hist(norm.mn, which="both")
hist(norm.mn, which="pm")

png(file="DensityPlot_NormMean_BgrdMAS.png",width=360,height=320)
par(mar=c(6,4,1,1));
hist(norm.mn, which="both")
dev.off()

boxplot(norm.mn) 
boxplot(norm.mn, which="both")
boxplot.dev(norm.mn, which="both", dev="png", outfile="BoxPlot_NormMeanLg2_BgrdMAS", mar=c(6,3,1,1), w=360, h=320)

## normalized intensity: PM+MM
data <- validData(norm.mn, which="both")
colnames(data) <- namePart(colnames(data))
head(data)

png(file="BoxPlot_NormMean_BgrdMAS_data.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(data), las=2)
dev.off()

## scatterplot
png(file="ScatterPlot_NormMean_BgrdMAS_both_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_NormMean_BgrdMAS_both_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()

## remove data
norm.qu <- removeInten(norm.qu)
gc()


## mean: trim=0.0, targetinten=500
norm.mn <- express(bgrd.mas, "tmp_NormMean500", filedir=outdir, tmpdir="", update=FALSE, 
           normalize.method="mean", normalize.select="both", normalize.option="transcript:all",
           normalize.logbase="0", normalize.params=c(0.0, 500))

## attach mask and intensity
norm.mn <- attachMask(norm.mn)
norm.mn <- attachInten(norm.mn)

## plots: include reference tree - difference caused by scaling to 500
hist(norm.mn)
hist(norm.mn, which="both")
hist(norm.mn, which="pm")

png(file="DensityPlot_NormMean500_BgrdMAS.png",width=360,height=320)
par(mar=c(6,4,1,1));
hist(norm.mn, which="both")
dev.off()

boxplot(norm.mn) 
boxplot(norm.mn, which="both")
boxplot.dev(norm.mn, which="both", dev="png", outfile="BoxPlot_NormMean500_BgrdMAS", mar=c(6,3,1,1), w=360, h=320)

## normalized intensity: PM+MM
data <- validData(norm.mn, which="both")
colnames(data) <- namePart(colnames(data))
head(data)

png(file="BoxPlot_NormMean500_BgrdMAS_data.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(data), las=2)
dev.off()

## scatterplot
png(file="ScatterPlot_NormMean500_BgrdMAS_both_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_NormMean500_BgrdMAS_both_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()

## remove data
norm.qu <- removeInten(norm.qu)
gc()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# median: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## median (log2): targetinten=500
norm.md <- express(bgrd.mas, "tmp_NormMedian", filedir=outdir, tmpdir="", update=FALSE, 
           normalize.method="median", normalize.select="pmonly", normalize.option="transcript:all",
           normalize.logbase="log2", normalize.params=c(500), reference.index=1)

## attach mask and intensity
norm.md <- attachMask(norm.md)
norm.md <- attachInten(norm.md)

## plots: include reference tree
hist(norm.md)
hist(norm.md, which="pm")

png(file="DensityPlot_NormMedianLg2_BgrdMAS.png",width=360,height=320)
par(mar=c(6,4,1,1));
hist(norm.md, which="pm")
dev.off()

boxplot(norm.md) 
boxplot(norm.md, which="pm")
boxplot.dev(norm.md, which="pm", dev="png", outfile="BoxPlot_NormMedian_BgrdMAS", mar=c(6,3,1,1), w=360, h=320)

## normalized intensity: PM only
data <- validData(norm.md, which="pm")
colnames(data) <- namePart(colnames(data))
head(data)

png(file="BoxPlot_NormMedianLg2_BgrdMAS_data.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(data), las=2)
dev.off()

## scatterplot
png(file="ScatterPlot_NormMedianLg2_BgrdMAS_PM_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_NormMedianLg2_BgrdMAS_PM_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()

## remove data
norm.md <- removeInten(norm.md)
gc()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# supsmu: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## supsmu: bass=0.0, span=0.0
norm.sup <- express(bgrd.mas, "tmp_NormSupsmu", filedir=outdir, tmpdir="", update=FALSE, 
            normalize.method="supsmu", normalize.select="pmonly", normalize.option="transcript:all",
            normalize.logbase="log2", normalize.params=c(0.0, 0.0, 0.0, 0.0))

## attach mask and intensity
norm.sup <- attachMask(norm.sup)
norm.sup <- attachInten(norm.sup)

## plots: include reference tree
hist(norm.sup)
hist(norm.sup, which="pm")

png(file="DensityPlot_NormSup_BgrdMAS.png",width=360,height=320)
par(mar=c(6,4,1,1));
hist(norm.sup, which="pm")
dev.off()

boxplot(norm.sup) 
boxplot(norm.sup, which="pm")
boxplot.dev(norm.sup, which="pm", dev="png", outfile="BoxPlot_NormSup_BgrdMAS", mar=c(6,3,1,1), w=360, h=320)

## normalized intensity: PM only
data <- validData(norm.sup, which="pm")
colnames(data) <- namePart(colnames(data))
head(data)

png(file="BoxPlot_NormSup_BgrdMAS_data.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(data), las=2)
dev.off()

## scatterplot
png(file="ScatterPlot_NormSup_BgrdMAS_PM_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_NormSup_BgrdMAS_PM_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()

## remove data
norm.sup <- removeInten(norm.sup)
gc()


## supsmu: bass=0.0, span=0.0, rule=2, f=0.0
norm.sup <- express(bgrd.mas, "tmp_NormSupsmuR2", filedir=outdir, tmpdir="", update=FALSE, 
            normalize.method="supsmu", normalize.select="pmonly", normalize.option="transcript:all",
            normalize.logbase="log2", normalize.params=c(0.0, 0.0, 2.0,0.0))

## attach mask and intensity
norm.sup <- attachMask(norm.sup)
norm.sup <- attachInten(norm.sup)

## plots: include reference tree
hist(norm.sup)
hist(norm.sup, which="pm")

png(file="DensityPlot_NormSupR2_BgrdMAS.png",width=360,height=320)
par(mar=c(6,4,1,1));
hist(norm.sup, which="pm")
dev.off()

boxplot(norm.sup) 
boxplot(norm.sup, which="pm")
boxplot.dev(norm.sup, which="pm", dev="png", outfile="BoxPlot_NormSupR2_BgrdMAS", mar=c(6,3,1,1), w=360, h=320)

## normalized intensity: PM only
data <- validData(norm.sup, which="pm")
colnames(data) <- namePart(colnames(data))
head(data)

png(file="BoxPlot_NormSupR2_BgrdMAS_data.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(data), las=2)
dev.off()

## scatterplot
png(file="ScatterPlot_NormSupR2_BgrdMAS_PM_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_NormSupR2_BgrdMAS_PM_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()

## remove data
norm.sup <- removeInten(norm.sup)
gc()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# lowess: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## lowess: span=0.67, iter=3
norm.low <- express(bgrd.mas, "tmp_NormLowess", filedir=outdir, tmpdir="", update=FALSE, 
            normalize.method="lowess", normalize.select="pmonly", normalize.option="transcript:all",
            normalize.logbase="log2", normalize.params=c(0.67, 3.0, 0.0, 0.0))

## attach mask and intensity
norm.low <- attachMask(norm.low)
norm.low <- attachInten(norm.low)

## plots: include reference tree
hist(norm.low)
hist(norm.low, which="pm")

png(file="DensityPlot_NormLow_BgrdMAS.png",width=360,height=320)
par(mar=c(6,4,1,1));
hist(norm.low, which="pm")
dev.off()

boxplot(norm.low) 
boxplot(norm.low, which="pm")
boxplot.dev(norm.low, which="pm", dev="png", outfile="BoxPlot_NormLow_BgrdMAS", mar=c(6,3,1,1), w=360, h=320)

## normalized intensity: PM only
data <- validData(norm.low, which="pm")
colnames(data) <- namePart(colnames(data))
head(data)

png(file="BoxPlot_NormLow_BgrdMAS_data.png",width=360,height=320)
par(mar=c(6,3,1,1));
boxplot(log2(data), las=2)
dev.off()

## scatterplot
png(file="ScatterPlot_NormLow_BgrdMAS_PM_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_NormLow_BgrdMAS_PM_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(data[,1]), log2(data[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()

## remove data
norm.low <- removeInten(norm.low)
gc()


#------------------------------------------------------------------------------#
# summarization using express()
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# RMA median-polish: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## medpol summarization: maxiter=10, eps=0.01, neglog=1.0
expr.mp <- express(norm.qu, "tmp_ExprMedpol", filedir=outdir, update=FALSE, 
           summarize.method="medianpolish", summarize.select="pmonly", summarize.option="transcript",
           summarize.logbase="log2", summarize.params=c(10, 0.01, 1.0))

## see above


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAS4 Average difference: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## AvgDiff summarization: STP=3
expr.adf <- express(bgrd.mas, "tmp_ExprAvgDif", filedir=outdir, update=FALSE, 
            summarize.method="avgdiff", summarize.select="none", summarize.option="transcript",
            summarize.logbase="0", summarize.params=c(3.0))

## plots
# hist(expr.adf) #negative values => NA
boxplot(expr.adf) 

boxplot.dev(expr.adf, dev="png", outfile="BoxPlot_ExprAvgDif", mar=c(6,3,1,1), w=360, h=320)

## get expression level
expr  <- validData(expr.adf)

## scatterplot
png(file="ScatterPlot_ExprAvgDif_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(expr[,1]), log2(expr[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_ExprAvgDif_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(expr[,1]), log2(expr[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAS5 Tukey-biweight: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## TukeyBiweight summarization: tau=0.03, scaletau=10, delta=2.0e-20, c=5, eps=0.0001, neglog=1.0, nfrac=0.5
expr.tbw <- express(bgrd.mas, "tmp_ExprTukey", filedir=outdir, update=FALSE, 
            summarize.method="tukeybiweight", summarize.select="none", summarize.option="transcript",
            summarize.logbase="log2", summarize.params=c(0.03, 10.0, 2.0e-20, 5.0, 0.0001, 1.0, 0.5))

## plots
hist(expr.tbw)
boxplot(expr.tbw) 

png(file="DensityPlot_ExprTukey.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(expr.tbw)
dev.off()

boxplot.dev(expr.tbw, dev="png", outfile="BoxPlot_ExprTukey", mar=c(6,3,1,1), w=360, h=320)

## get expression level
expr  <- validData(expr.tbw)

## scatterplot
png(file="ScatterPlot_ExprTukey_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(expr[,1]), log2(expr[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_ExprTukey_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(expr[,1]), log2(expr[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# FARMS: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## default for FARMS and DFW :use quantile normalization w/o background-correction
norm.qu <- express(data.u133p2, "tmp_NormQuanU133", filedir=outdir, update=FALSE, 
           normalize.method="quantile", normalize.select="pmonly", normalize.option="transcript:together:none",
           normalize.logbase="0", normalize.params=c(0.0))


## FARMS summarization: version=131, weight=0.5, mu=0.0, scale=1.0, tol=0.00001, cyc=100, weighted=1
expr.frm <- express(norm.qu, "tmp_ExprFARMS", filedir=outdir, update=FALSE, 
            summarize.method="farms", summarize.select="pmonly", summarize.option="transcript",
            summarize.logbase="log2", summarize.params=c(131, 0.5, 0.0, 1.0, 0.00001, 100, 1))

## plots
hist(expr.frm)
boxplot(expr.frm) 

png(file="DensityPlot_ExprFARMS.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(expr.frm)
dev.off()

boxplot.dev(expr.frm, dev="png", outfile="BoxPlot_ExprFARMS", mar=c(6,3,1,1), w=360, h=320)

## get expression level
expr  <- validData(expr.frm)

## scatterplot
png(file="ScatterPlot_ExprFARMS_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(expr[,1]), log2(expr[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_ExprFARMS_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(expr[,1]), log2(expr[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DFW: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## DFW summarization: m=3, n=1, c=0.001
expr.dfw <- express(norm.qu, "tmp_ExprDFW", filedir=outdir, update=FALSE, 
            summarize.method="dfw", summarize.select="pmonly", summarize.option="transcript",
            summarize.logbase="log2", summarize.params=c(3.0, 1.0, 0.01))

## plots
hist(expr.dfw)
boxplot(expr.dfw) 

png(file="DensityPlot_ExprDFW.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(expr.dfw)
dev.off()

boxplot.dev(expr.dfw, dev="png", outfile="BoxPlot_ExprDFW", mar=c(6,3,1,1), w=360, h=320)

## get expression level
expr  <- validData(expr.dfw)

## scatterplot
png(file="ScatterPlot_ExprDFW_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(expr[,1]), log2(expr[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_ExprDFW_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(expr[,1]), log2(expr[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()


#------------------------------------------------------------------------------#
# probeset-level normalization using normalize()
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# normalize.constant: mean/median
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## mean normalization: trim=0.02, targetinten=500 (refindex=1 to prevent refmethod=mean)
nrxp.mn <- normalize(expr.tbw,  "tmp_NrmExpMean", filedir=outdir, update=FALSE, 
           select="separate", method="mean", option="transcript:all", logbase="0", 
           refindex=1, refmethod="mean", params=c(0.02, 500))
## or
nrxp.mn <- normalize.constant(expr.tbw, "tmp_NrmExpMean", filedir=outdir,
           refindex=1, refmethod="mean", params=c(0.02, 500))

## plots
hist(nrxp.mn)
boxplot(nrxp.mn) 

png(file="DensityPlot_NrmExpMean.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(nrxp.mn)
dev.off()

boxplot.dev(nrxp.mn, dev="png", outfile="BoxPlot_NrmExpMean", mar=c(6,3,1,1), w=360, h=320)

## get expression level
nrxp  <- validData(nrxp.mn)

## scatterplot
png(file="ScatterPlot_NrmExpMean_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(nrxp[,1]), log2(nrxp[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_NrmExpMean_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(nrxp[,1]), log2(nrxp[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# normalize.supsmu: supsmu
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## supsmu normalization: bass=0.0, span=0.0
nrxp.sup <- normalize(expr.tbw,  "tmp_NrmExpSup", filedir=outdir, update=FALSE, 
            select="separate", method="supsmu", option="transcript:all", logbase="log2", 
            refindex=0, refmethod="mean", params=c(0.0, 0.0, 0.0, 0.0))
## or
nrxp.sup <- normalize.supsmu(expr.tbw, "tmp_NrmExpSup", filedir=outdir)

## plots
hist(nrxp.sup)
boxplot(nrxp.sup) 

png(file="DensityPlot_NrmExpSup.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(nrxp.sup)
dev.off()

boxplot.dev(nrxp.sup, dev="png", outfile="BoxPlot_NrmExpSup", mar=c(6,3,1,1), w=360, h=320)

## get expression level
nrxp  <- validData(nrxp.sup)

## scatterplot
png(file="ScatterPlot_NrmExpSup_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(nrxp[,1]), log2(nrxp[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_NrmExpSup_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(nrxp[,1]), log2(nrxp[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()


## supsmu normalization: bass=0.0, span=0.0, rule=2, f=0.0
nrxp.sup <- normalize(expr.tbw,  "tmp_NrmExpSupR2", filedir=outdir, update=FALSE, 
            select="separate", method="supsmu", option="transcript:all", logbase="log2", 
            refindex=0, refmethod="mean", params=c(0.0, 0.0, 2, 0.0))

## plots
hist(nrxp.sup)
boxplot(nrxp.sup) 

png(file="DensityPlot_NrmExpSupR2.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(nrxp.sup)
dev.off()

boxplot.dev(nrxp.sup, dev="png", outfile="BoxPlot_NrmExpSupR2", mar=c(6,3,1,1), w=360, h=320)

## get expression level
nrxp  <- validData(nrxp.sup)

## scatterplot
png(file="ScatterPlot_NrmExpSupR2_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(nrxp[,1]), log2(nrxp[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_NrmExpSupR2_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(nrxp[,1]), log2(nrxp[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# normalize.lowess: lowess
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## lowess normalization: span=0.67, iter=3
nrxp.low <- normalize(expr.tbw,  "tmp_NrmExpLow", filedir=outdir, update=FALSE, 
            select="separate", method="lowess", option="transcript:all", logbase="log2", 
            refindex=0, refmethod="mean", params=c(0.67, 3, 0.0, 0.0))
## or
nrxp.low <- normalize.lowess(expr.tbw, "tmp_NrmExpLow", filedir=outdir)

## plots
hist(nrxp.low)
boxplot(nrxp.low) 

png(file="DensityPlot_NrmExpLow.png",width=360,height=320)
par(mar=c(6,3,1,1));
hist(nrxp.low)
dev.off()

boxplot.dev(nrxp.low, dev="png", outfile="BoxPlot_NrmExpLow", mar=c(6,3,1,1), w=360, h=320)

## get expression level
nrxp  <- validData(nrxp.low)

## scatterplot
png(file="ScatterPlot_NrmExpLow_BrABrB.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(nrxp[,1]), log2(nrxp[,2]), xlab="BreastA", ylab="BreastB")
dev.off()

png(file="ScatterPlot_NrmExpLow_BrAPrA.png",width=360,height=360)
par(mar=c(5,5,1,1));
plot(log2(nrxp[,1]), log2(nrxp[,4]), xlab="BreastA", ylab="ProstateA")
dev.off()







#------------------------------------------------------------------------------#

