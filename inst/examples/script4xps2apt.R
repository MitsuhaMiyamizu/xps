#------------------------------------------------------------------------------#
# Script: step-by-step functions to demonstrate how to compare results obtained
#         with xps to results obtained with apt (Affymetrix Power Tools)
#         using the Affymetrix human tissue-mixture exon array dataset
#
# Copyright (c) 2008-2009 Christian Stratowa, Vienna, Austria.
# All rights reserved.
#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Tissues from Affymetrix Exon Array Dataset for HG-U133_Plus_2 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Software versions used:
# R:    R-2.9.0
# xps:  xps_1.5.7
# affy: affy_1.22.0
# APT:  apt-1.10.2
#------------------------------------------------------------------------------#


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. XPS: compute rma, mas5, mas5calls
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### define directories:
libdir <- "/Volumes/GigaDrive/Affy/libraryfiles"
anndir <- "/Volumes/GigaDrive/Affy/Annotation/Version09Mar"
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
celdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/u133p2_apt"
outdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/u133p2_apt"

### new R session: import scheme and data
library(xps)

## scheme for HG-U133_Plus_2:
scheme.hgu133p2 <- import.expr.scheme("Scheme_HGU133p2_na28", filedir=scmdir,
                   schemefile=paste(libdir,"HG-U133_Plus_2.cdf",sep="/"),
                   probefile=paste(libdir,"HG-U133-PLUS_probe.tab",sep="/"),
                   annotfile=paste(anndir,"HG-U133_Plus_2.na28.annot.csv",sep="/"))

## subset of CEL files to import
celfiles <- c("u1332plus_ivt_breast_A.CEL", "u1332plus_ivt_breast_B.CEL", "u1332plus_ivt_breast_C.CEL",
              "u1332plus_ivt_prostate_A.CEL", "u1332plus_ivt_prostate_B.CEL", "u1332plus_ivt_prostate_C.CEL")
celnames <- c("BreastA", "BreastB", "BreastC", "ProstateA", "ProstateB", "ProstateC")
data.mix.u133p2 <- import.data(scheme.hgu133p2, "HuTissuesU133P2", filedir=datdir,celdir=celdir,celfiles=celfiles,celnames=celnames)


### new R session: rma, mas5, mas5call
library(xps)

## import ROOT scheme and data files
scheme.u133p2 <- root.scheme(paste(scmdir,"Scheme_HGU133p2_na28.root",sep="/"))
data.u133p2   <- root.data(scheme.u133p2, paste(datdir,"HuTissuesU133P2_cel.root",sep="/"))

## RMA
data.rma <- rma(data.u133p2,"MixU133P2RMA",filedir=outdir,tmpdir="",
                background="pmonly",normalize=TRUE)
## MAS5
data.mas5 <- mas5(data.u133p2,"MixU133P2MAS5All",filedir=outdir,tmpdir="",
                  normalize=TRUE,sc=500, update=TRUE)
## MAS5 detection call
call.mas5 <- mas5.call(data.u133p2,"MixU133P2Call",filedir=outdir,tmpdir="")

## get data.frames
xps.rma  <- validData(data.rma)
xps.rma  <- xps.rma[order(rownames(xps.rma)),]

xps.mas5 <- validData(data.mas5)
xps.mas5 <- xps.mas5[order(rownames(xps.mas5)),]

xps.pval <- validData(call.mas5)
xps.pval <- xps.pval[order(rownames(xps.pval)),]

xps.pres <- validCall(call.mas5)
xps.pres <- xps.pres[order(rownames(xps.pres)),]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. APT: compute rma, mas5, mas5calls
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### run on command line

## apt rma: sketch-quantile
apt-probeset-summarize -a rma-bg,quant-norm.sketch=-1.usepm=true.bioc=true,pm-only,med-polish.expon=true -d HG-U133_Plus_2.cdf -o ./rma_sk *.CEL

## apt rma: use quantile and not sketch-quantile
apt-probeset-summarize -a rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish.expon=true -d HG-U133_Plus_2.cdf -o ./rma *.CEL

## apt mas5 signal
apt-probeset-summarize -a mas5-bg,pm-mm,mas5-signal -d HG-U133_Plus_2.cdf -o ./mas5 *.CEL

### apt mas5 detection call
apt-probeset-summarize -a pm-mm,mas5-detect.calls=1.pairs=1 -d HG-U133_Plus_2.cdf -x 10 -o ./mas5call *.CEL


### import results into R

## apt rma: sketch-quantile
apt.rma.sk <- read.delim("./rma_sk/rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=52)
apt.rma.sk <- apt.rma.sk[order(rownames(apt.rma.sk)),]

## apt rma
apt.rma <- read.delim("./rma/rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=52)
apt.rma <- apt.rma[order(rownames(apt.rma)),]

## apt mas5
apt.mas5 <- read.delim("./mas5/mas5-bg.pm-mm.mas5-signal.summary.txt", row.names=1, comment.char="", skip=52)
# note: as mentioned in FAQ, apt does not allow for signal level normalization across the chip.
#       since xps mas5 scales to sc=500, we need to normalize the apt results first
apt.mas5 <- apply(apt.mas5, 2, function(x){x*(500/mean(x, trim=0.02))})
apt.mas5 <- apt.mas5[order(rownames(apt.mas5)),]

# apt mas5call
apt.pval <- read.delim("./mas5call/pm-mm.mas5-detect.summary.txt", row.names=1, comment.char="", skip=52)
apt.pval <- apt.pval[order(rownames(apt.pval)),]

## import ExpressionConsole mas5:
ec.mas5 <- read.delim("EC_MAS500_signal.txt", row.names=1, comment.char="")
ec.mas5 <- ec.mas5[order(rownames(ec.mas5)),]

## import ExpressionConsole mas5 call:
ec.pval <- read.delim("EC_MAS500_pval.txt", row.names=1, comment.char="")
ec.pval <- ec.pval[order(rownames(ec.pval)),]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 3. affy: compute rma, mas5, mas5calls
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

detach(package:xps)
library(affy)

## rma:
affy.rma <- justRMA()
affy.rma <- 2^exprs(affy.rma)

## mas5:
affy <- ReadAffy()
affy.mas5 <- mas5(affy, normalize=TRUE, sc=500)
affy.mas5 <- exprs(affy.mas5)

## mas5calls
affy.dc5  <- mas5calls(affy)
affy.pval <- assayData(affy.dc5)[["se.exprs"]]

rm(affy); gc()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 4. compare apt vs xps vs affy
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### compare RMA

## compare expression levels obtained with quantile-sketch vs quantile
png(file="U133P2_APT_RMA_RMAsketch_MvA.png",width=400,height=400)
plot(log2(apt.rma[,1] * apt.rma.sk[,1])/2, log2(apt.rma.sk[,1]/apt.rma[,1]), main="APT: RMA vs RMA-sketch", xlab="A = Log2(SketchRMA*RMA)", ylab="M = Log2(SketchRMA/RMA)",log="",ylim=c(-0.1,0.1))
dev.off()

## compare expression levels obtained with xps vs apt
png(file="U133P2_XPS_APT_RMA_MvA.png",width=400,height=400)
plot(log2(apt.rma[,1] * xps.rma[,1])/2, log2(xps.rma[,1]/apt.rma[,1]), main="RMA: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.1,0.1))
dev.off()

## plot xps vs affy vs apt
tmp <- cbind(xps.rma[,1],affy.rma[,1],apt.rma[,1])
colnames(tmp) <- c("xps.rma","affy.rma","apt.rma")
png(file="U133P2_XPS_AFFY_APT_RMA.png",width=400,height=400)
pairs(log2(tmp), labels=colnames(tmp))
dev.off()

# plot percent difference
plot(1:54675, 100*(apt.rma[,1] - xps.rma[,1])/apt.rma[,1])
plot(1:54675, 100*(affy.rma[,1] - xps.rma[,1])/affy.rma[,1])
plot(1:54675, 100*(affy.rma[,1] - apt.rma[,1])/affy.rma[,1])


### compare MAS5

## compare expression levels obtained with xps, apt, ExpressionConsole, affy
png(file="U133P2_XPS_APT_EC_AFFY_MAS5_MvA_3x2.png",width=400,height=600)
oldpar <- par(mfrow=c(3,2), oma=c(0,0,0,0), mar=c(4,4,2,1))
plot(log2(xps.mas5[,1] * apt.mas5[,1])/2, log2(xps.mas5[,1]/apt.mas5[,1]), main="MAS5: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-1,1))
plot(log2(affy.mas5[,1] * apt.mas5[,1])/2, log2(affy.mas5[,1]/apt.mas5[,1]), main="MAS5: AFFY vs APT", xlab="A = Log2(AFFY*APT)", ylab="M = Log2(AFFY/APT)",log="",ylim=c(-1,1))
plot(log2(xps.mas5[,1] * ec.mas5[,1])/2, log2(xps.mas5[,1]/ec.mas5[,1]), main="MAS5: XPS vs EC", xlab="A = Log2(XPS*EC)", ylab="M = Log2(XPS/EC)",log="",ylim=c(-1,1))
plot(log2(affy.mas5[,1] * ec.mas5[,1])/2, log2(affy.mas5[,1]/ec.mas5[,1]), main="MAS5: AFFY vs EC", xlab="A = Log2(AFFY*EC)", ylab="M = Log2(AFFY/EC)",log="",ylim=c(-1,1))
plot(log2(apt.mas5[,1] * ec.mas5[,1])/2, log2(apt.mas5[,1]/ec.mas5[,1]), main="MAS5: APT vs EC", xlab="A = Log2(APT*EC)", ylab="M = Log2(APT/EC)",log="",ylim=c(-1,1))
plot(log2(xps.mas5[,1] * affy.mas5[,1])/2, log2(xps.mas5[,1]/affy.mas5[,1]), main="MAS5: XPS vs AFFY", xlab="A = Log2(XPS*AFFY)", ylab="M = Log2(XPS/AFFY)",log="",ylim=c(-1,1))
par(oldpar);
dev.off()

# compare expression levels obtained with xps, apt, ExpressionConsole, affy
png(file="U133P2_Diff_XPS_APT_EC_AFFY_MAS5_2x3.png",width=600,height=440)
oldpar <- par(mfrow=c(2,3), oma=c(0,0,0,0), mar=c(4,4,2,1))
plot(xps.mas5[,1], 100*(apt.mas5[,1] - xps.mas5[,1])/apt.mas5[,1], xlab="Intensity XPS", ylab="(APT - XPS) / APT [%]", log="x", ylim=c(-100,100))
plot(affy.mas5[,1], 100*(apt.mas5[,1] - affy.mas5[,1])/apt.mas5[,1], xlab="Intensity AFFY", ylab="(APT - AFFY) / APT [%]", log="x", ylim=c(-100,100))
plot(ec.mas5[,1], 100*(apt.mas5[,1] - ec.mas5[,1])/apt.mas5[,1], xlab="Intensity EC", ylab="(APT - EC) / APT [%]", log="x", ylim=c(-100,100))
plot(xps.mas5[,1], 100*(ec.mas5[,1] - xps.mas5[,1])/ec.mas5[,1], xlab="Intensity XPS", ylab="(EC - XPS) / EC [%]", log="x", ylim=c(-100,100))
plot(affy.mas5[,1], 100*(ec.mas5[,1] - affy.mas5[,1])/ec.mas5[,1], xlab="Intensity AFFY", ylab="(AFFY - EC) / EC [%]", log="x", ylim=c(-100,100))
plot(xps.mas5[,1], 100*(affy.mas5[,1] - xps.mas5[,1])/affy.mas5[,1], xlab="Intensity XPS", ylab="(AFFY - XPS) / AFFY [%]", log="x", ylim=c(-100,100))
par(oldpar);
dev.off()

## plot xps vs affy vs EC vs apt
tmp <- cbind(xps.mas5[,1],affy.mas5[,1],ec.mas5[,1],apt.mas5[,1])
colnames(tmp) <- c("xps.mas5","affy.mas5","ec.mas5","apt.mas5")
png(file="U133P2_XPS_AFFY_EC_APT_MAS5.png",width=440,height=440)
pairs(log2(tmp), labels=colnames(tmp))
dev.off()


### compare MAS5 detection calls

# compare detection calls obtained with xps vs affy
png(file="U133P2_Diff_XPS_AFFY_MAS5_pval.png",width=400,height=400)
plot(xps.pval[,1], 100*(affy.pval[,1] - xps.pval[,1])/affy.pval[,1], main="MAS5 P-Value: XPS vs AFFY", xlab="XPS P-Value", ylab="Difference [%]", ylim=c(-1,1))
dev.off()

# plot xps vs affy vs apt
tmp <- cbind(xps.pval[,1],affy.pval[,1],ec.pval[,1],apt.pval[,1])
colnames(tmp) <- c("xps.pval","affy.pval","ec.pval","apt.pval")
png(file="U133P2_XPS_AFFY_EC_APT_MAS5_pval.png",width=440,height=440)
pairs(log2(tmp), labels=colnames(tmp),xlim=c(-15,0),ylim=c(-15,0))
dev.off()



#------------------------------------------------------------------------------#
# Tissues from Affymetrix Exon Array Dataset for HuGene-1_0-st-v1 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Software versions used:
# R:   R-2.9.0
# xps: xps_1.5.7
# APT: apt-1.10.2
#------------------------------------------------------------------------------#


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. XPS: compute rma, dabg calls
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### define directories:
libdir <- "/Volumes/GigaDrive/Affy/libraryfiles"
anndir <- "/Volumes/GigaDrive/Affy/Annotation/Version09Mar"
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
celdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/hugene_apt"
outdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/hugene_apt"

### new R session: import scheme and data
library(xps)

## scheme for HuGene-1_0-st-v1:
scheme.genome <- import.exon.scheme("Scheme_HuGene10stv1r4_na28",filedir=scmdir,
                 paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.clf",sep="/"),
                 paste(libdir,"HuGene-1_0-st-v1.r4.analysis-lib-files/HuGene-1_0-st-v1.r4.pgf",sep="/"),
                 paste(anndir,"HuGene-1_0-st-v1.na28.hg18.probeset.csv",sep="/"),
                 paste(anndir,"HuGene-1_0-st-v1.na28.hg18.transcript.csv",sep="/"))

# subset of CEL files to import
celfiles <- c("TisMap_Breast_01_v1_WTGene1.CEL","TisMap_Breast_02_v1_WTGene1.CEL","TisMap_Breast_03_v1_WTGene1.CEL",
              "TisMap_Prostate_01_v1_WTGene1.CEL","TisMap_Prostate_02_v1_WTGene1.CEL","TisMap_Prostate_03_v1_WTGene1.CEL")
celnames <- c("Breast01","Breast02","Breast03","Prostate01","Prostate02","Prostate03")
data.genome <- import.data(scheme.genome, "HuTissuesGenome", filedir=datdir, celdir=celdir, celfiles=celfiles, celnames=celnames)


### new R session: rma, dabg call
library(xps)

## import ROOT scheme and data files
scheme.genome <- root.scheme(paste(scmdir, "Scheme_HuGene10stv1r4_na28.root", sep="/"))
data.genome   <- root.data(scheme.genome, paste(datdir,"HuTissuesGenome_cel.root", sep="/"))

## RMA for transcripts
data.rma.tc <- rma(data.genome,"MixHuGeneRMAcore_tc", filedir=outdir, tmpdir="", background="antigenomic",
                   normalize=T, option="transcript", exonlevel="core")
xps.rma.tc  <- validData(data.rma.tc)

## RMA for probesets
data.rma.ps <- rma(data.genome, "MixHuGeneRMAcore_ps", filedir=outdir, tmpdir="", background="antigenomic",
                   normalize=T, option="probeset", exonlevel="core")
xps.rma.ps  <- validData(data.rma.ps)

## xps medpol only for transcripts
expr.mp.tc <- express(data.genome,"MixHuGeneMedPolcore_tc",summarize.method="medianpolish",
              summarize.select="pmonly",summarize.option="transcript",summarize.logbase="log2",
              summarize.params=c(10, 0.01, 1.0),exonlevel="core")
xps.mp.tc <- validData(expr.mp.tc)

## xps medpol only for probesets
expr.mp.ps <- express(data.genome,"MixHuGeneMedPolcore_ps",summarize.method="medianpolish",
              summarize.select="pmonly",summarize.option="probeset",summarize.logbase="log2",
              summarize.params=c(10, 0.01, 1.0),exonlevel="core")
xps.mp.ps <- validData(expr.mp.ps)

## xps rma: bg="all", quantile="all", medpol="core"
data.rma.tc.bq16 <- rma(data.genome,"HuGeneMixRMAbgqu16mp9",filedir=outdir,tmpdir="",
                        background="antigenomic",normalize=T,exonlevel=c(16316,16316,9216))
xps.rma.tc.bq16 <- validData(data.rma.tc.bq16)

## xps rma: bg=quantile="core+affx+unmapped+exon+intron+antigenomic", medpol="core"
data.rma.tc.bq99 <- rma(data.genome,"HuGeneMixRMAbgqu99mp9",filedir=outdir,tmpdir="",
                        background="antigenomic",normalize=T,exonlevel=c(992316,992316,9216))
xps.rma.tc.bq99 <- validData(data.rma.tc.bq99)

## xps rma: bg=quantile="core+affx+unmapped+exon+intron+antigenomic", medpol="core"
data.rma.ps.bq99 <- rma(data.genome,"HuGeneMixRMAbgqu99mp9_ps", filedir=outdir, tmpdir="",
                        background="antigenomic", normalize=T, option="probeset", exonlevel=c(992316,992316,9216))
xps.rma.ps.bq99 <- validData(data.rma.ps.bq99)


## need to create probeset list for APT so that APT can use same probesets as XPS
# transcriptList.txt
tmp <- as.data.frame(rownames(xps.rma.tc))
colnames(tmp) <- "probeset_id"
write.table(tmp, "transcriptList.txt", quote=FALSE, row.names=FALSE)
# probesetList.txt
tmp <- as.data.frame(rownames(xps.rma.ps))
colnames(tmp) <- "probeset_id"
write.table(tmp, "probesetList.txt", quote=FALSE, row.names=FALSE)
# coreList.mps
writeLines(rownames(xps.rma.tc), "core.txt")
metaProbesets(scheme.genome, "core.txt", "coreList.mps", "core")


## DABG for transcripts
call.dabg.tc <- dabg.call(data.genome,"HuGeneMixDABGcore_tc",filedir=datdir,
                          option="transcript", exonlevel="core")
xps.pval.tc <- validData(call.dabg.tc)

## DABG for probesets
call.dabg.ps <- dabg.call(data.genome,"HuGeneMixDABGcore_ps",filedir=datdir,
                          option="probeset", exonlevel="core")
xps.pval.ps <- validData(call.dabg.ps)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. APT: compute rma, dabg calls
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### run on command line

## apt rma: use quantile; use transcriptList.txt
apt-probeset-summarize -a rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish.expon=true -p HuGene-1_0-st-v1.r3.pgf -c HuGene-1_0-st-v1.r3.clf -s transcriptList.txt -o ./rma_r3 *.CEL

## apt rma: use quantile; use probesetList.txt
apt-probeset-summarize -a rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish.expon=true -p HuGene-1_0-st-v1.r4.pgf -c HuGene-1_0-st-v1.r4.clf -s probesetList.txt -o ./rma_r4_ps *.CEL

## apt rma: use quantile; use coreList.mps
apt-probeset-summarize -a rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish.expon=true -p HuGene-1_0-st-v1.r4.pgf -c HuGene-1_0-st-v1.r4.clf -m coreList.mps -o ./rma_r4_tc *.CEL


## apt medpol only; use transcriptList.txt
apt-probeset-summarize -a pm-only,med-polish.expon=true -p HuGene-1_0-st-v1.r3.pgf -c HuGene-1_0-st-v1.r3.clf -s transcriptList.txt -o ./rma_r3 *.CEL

## apt medpol only; use probesetList.txt
apt-probeset-summarize -a pm-only,med-polish.expon=true -p HuGene-1_0-st-v1.r4.pgf -c HuGene-1_0-st-v1.r4.clf -s probesetList.txt -o ./rma_r4_ps *.CEL

## apt medpol only; use coreList.mps
apt-probeset-summarize -a pm-only,med-polish.expon=true -p HuGene-1_0-st-v1.r4.pgf -c HuGene-1_0-st-v1.r4.clf -m coreList.mps -o ./rma_r4_tc *.CEL


### apt dabg call; transcriptList.txt
apt-probeset-summarize -a dabg -p HuGene-1_0-st-v1.r3.pgf -c HuGene-1_0-st-v1.r3.clf -b HuGene-1_0-st-v1.r3.bgp -s transcriptList.txt  -x 8 -o ./dabg_r3 *.CEL

## apt dabg call; probesetList.txt
apt-probeset-summarize -a dabg -p HuGene-1_0-st-v1.r4.pgf -c HuGene-1_0-st-v1.r4.clf -b HuGene-1_0-st-v1.r4.bgp -s probesetList.txt -x 8 -o ./dabg_r4_ps *.CEL

## apt dabg call; coreList.mps
apt-probeset-summarize -a dabg -p HuGene-1_0-st-v1.r4.pgf -c HuGene-1_0-st-v1.r4.clf -b HuGene-1_0-st-v1.r4.bgp -m coreList.mps -x 8 -o ./dabg_r4_tc *.CEL


### import results into R

## apt rma
apt.rma.r3.tc <- read.delim("./rma_r3/rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=52)

## apt rma
apt.rma.r4.ps <- read.delim("./rma_r4_ps/rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=52)

## apt rma
apt.rma.r4.tc <- read.delim("./rma_r4_tc/rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=52)


## apt medpol
apt.mp.r3.tc <- read.delim("./rma_r3/pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=52)

## apt medpol
apt.mp.r4.ps <- read.delim("./rma_r4_ps/pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=52)

## apt medpol
apt.mp.r4.tc <- read.delim("./rma_r4_tc/pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=52)


## apt dabg
apt.dabg.r3.tc <- read.delim("./dabg_r3/dabg.summary.txt", row.names=1, comment.char="", skip=52)

## apt dabg
apt.dabg.r4.ps <- read.delim("./dabg_r4_ps/dabg.summary.txt", row.names=1, comment.char="", skip=52)

## apt dabg
apt.dabg.r4.tc <- read.delim("./dabg_r4_tc/dabg.summary.txt", row.names=1, comment.char="", skip=52)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 3. compare apt vs xps vs affy
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### compare RMA

# compare RMA obtained with xps vs apt
png(file="HuGene_XPS_APT_RMA_MvA_tc.png",width=400,height=400)
plot(log2(apt.rma.r3.tc[,1] * xps.rma.tc[,1])/2, log2(xps.rma.tc[,1]/apt.rma.r3.tc[,1]), main="RMA TC: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,0.2))
dev.off()

png(file="HuGene_XPS_APT_RMA_MvA_tc_r4.png",width=400,height=400)
plot(log2(apt.rma.r4.tc[,1] * xps.rma.tc[,1])/2, log2(xps.rma.tc[,1]/apt.rma.r4.tc[,1]), main="RMA TC: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,0.2))
dev.off()

png(file="HuGene_XPS_APT_RMA_MvA_ps_r4.png",width=400,height=400)
plot(log2(apt.rma.r4.ps[,1] * xps.rma.ps[,1])/2, log2(xps.rma.ps[,1]/apt.rma.r4.ps[,1]), main="RMA PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,0.2))
dev.off()


## compare medpol obtained with xps vs apt
png(file="HuGene_XPS_APT_MedPol_MvA_tc.png",width=400,height=400)
plot(log2(apt.mp.r3.tc[,1] * xps.mp.tc[,1])/2, log2(xps.mp.tc[,1]/apt.mp.r3.tc[,1]), main="MedPol TC: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.001,0.001))
dev.off()

png(file="HuGene_XPS_APT_MedPol_MvA_tc_r4.png",width=400,height=400)
plot(log2(apt.mp.r4.tc[,1] * xps.mp.tc[,1])/2, log2(xps.mp.tc[,1]/apt.mp.r4.tc[,1]), main="MedPol TC: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.001,0.001))
dev.off()

png(file="HuGene_XPS_APT_MedPol_MvA_ps_r4.png",width=400,height=400)
plot(log2(apt.mp.r4.ps[,1] * xps.mp.ps[,1])/2, log2(xps.mp.ps[,1]/apt.mp.r4.ps[,1]), main="MedPol PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.001,0.001))
dev.off()


## compare rma: bg="all", quantile="all", medpol="core" obtained with xps vs apt
png(file="HuGene_XPS_APT_RMA_bg16qu16_MvA_tc.png",width=400,height=400)
plot(log2(apt.rma.r3.tc[,1] * xps.rma.tc.bq16[,1])/2, log2(xps.rma.tc.bq16[,1]/apt.rma.r3.tc[,1]), main="RMA TC: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,0.2))
dev.off()

## compare rma: bg=quantile="core+affx+unmapped+exon+intron+antigenomic", medpol="core" obtained with xps vs apt
png(file="HuGene_XPS_APT_RMA_bg99qu99_MvA_tc.png",width=400,height=400)
plot(log2(apt.rma.r3.tc[,1] * xps.rma.tc.bq99[,1])/2, log2(xps.rma.tc.bq99[,1]/apt.rma.r3.tc[,1]), main="RMA TC: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.001,0.001))
dev.off()

png(file="HuGene_XPS_APT_RMA_bg99qu99_MvA_tc_r4.png",width=400,height=400)
plot(log2(apt.rma.r4.tc[,1] * xps.rma.tc.bq99[,1])/2, log2(xps.rma.tc.bq99[,1]/apt.rma.r4.tc[,1]), main="RMA TC: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.001,0.001))
dev.off()

png(file="HuGene_XPS_APT_RMA_bg99qu99_MvA_ps_r4.png",width=400,height=400)
plot(log2(apt.rma.r4.ps[,1] * xps.rma.ps.bq99[,1])/2, log2(xps.rma.ps.bq99[,1]/apt.rma.r4.ps[,1]), main="RMA PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.001,0.001))
dev.off()


### compare DABG

# compare dabg obtained with xps vs apt
png(file="HuGeneDiff_XPS_APT_DABG_pval_tc.png",width=400,height=400)
plot(xps.pval.tc[,1], (xps.pval.tc[,1] - apt.dabg.r3.tc[,1]), main="DABG TC P-Value: XPS vs APT", xlab="XPS P-Value", ylab="Difference", log="x", ylim=c(-0.00001,0.00001))
dev.off()

png(file="HuGeneDiff_XPS_APT_DABG_pval_ps_r4.png",width=400,height=400)
plot(xps.pval.ps[,1], (xps.pval.ps[,1] - apt.dabg.r4.ps[,1]), main="DABG PS P-Value: XPS vs APT", xlab="XPS P-Value", ylab="Difference", log="x", ylim=c(-0.00001,0.00001))
dev.off()

png(file="HuGeneDiff_XPS_APT_DABG_pval_tc_r4.png",width=400,height=400)
plot(xps.pval.tc[,1], (xps.pval.tc[,1] - apt.dabg.r4.tc[,1]), main="DABG TC P-Value: XPS vs APT", xlab="XPS P-Value", ylab="Difference", log="x", ylim=c(-0.00001,0.00001))
dev.off()



#------------------------------------------------------------------------------#
# Tissues from Affymetrix Exon Array Dataset for HuEx-1_0-st-v2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Software versions used:
# R:   R-2.9.0
# xps: xps_1.5.7
# APT: apt-1.10.2
#------------------------------------------------------------------------------#


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. XPS: compute rma, dabg calls
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### define directories:
libdir <- "/Volumes/GigaDrive/Affy/libraryfiles"
anndir <- "/Volumes/GigaDrive/Affy/Annotation/Version09Mar"
scmdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Schemes"
datdir <- "/Volumes/GigaDrive/CRAN/Workspaces/ROOTData"
celdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/exon_apt"
outdir <- "/Volumes/GigaDrive/CRAN/Workspaces/Exon/hutissues/exon_apt"

### new R session: import scheme and data
library(xps)

## scheme for HuEx-1_0-st-v2:
scheme.exon <- import.exon.scheme("Scheme_HuEx10stv2r2_na28",filedir=scmdir,
               paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.clf",sep="/"),
               paste(libdir,"HuEx-1_0-st-v2_libraryfile/HuEx-1_0-st-r2/HuEx-1_0-st-v2.r2.pgf",sep="/"),
               paste(anndir,"HuEx-1_0-st-v2.na28.hg18.probeset.csv",sep="/"),
               paste(anndir,"HuEx-1_0-st-v2.na28.hg18.transcript.csv",sep="/"))

# subset of CEL files to import
celfiles <- c("huex_wta_breast_A.CEL", "huex_wta_breast_B.CEL", "huex_wta_breast_C.CEL",
              "huex_wta_prostate_A.CEL", "huex_wta_prostate_B.CEL", "huex_wta_prostate_C.CEL")
celnames <- c("BreastA", "BreastB", "BreastC", "ProstateA", "ProstateB", "ProstateC")
data.exon <- import.data(scheme.exon, "HuTissuesExon", filedir=datdir, celdir=celdir, celfiles=celfiles, celnames=celnames)


### new R session: rma, dabg call
library(xps)

## import ROOT scheme and data files
scheme.exon <- root.scheme(paste(scmdir, "Scheme_HuEx10stv2r2_na28.root", sep="/"))
data.exon   <- root.data(scheme.exon, paste(datdir,"HuTissuesExon_cel.root", sep="/"))

## RMA for transcripts
data.rma.tc <- rma(data.exon,"MixHuExonRMAcore_tc", filedir=outdir, tmpdir="", background="antigenomic",
                   normalize=T, option="transcript", exonlevel="core")
xps.rma.tc  <- validData(data.rma.tc)

## medpol only for transcripts
expr.mp.tc <- express(data.exon, "HuExonMedPolcore", summarize.method="medianpolish",
              summarize.select="pmonly", summarize.option="transcript", summarize.logbase="log2",
              summarize.params=c(10, 0.01, 1.0), exonlevel="core")
xps.mp.tc <- validData(expr.mp.tc)

## xps rma: bg=quantile="all+xxx", medpol="core"
data.rma.tc.bq103 <- rma(data.exon,"HuExonMixRMAbgqu103mp9_tc",filedir=outdir,tmpdir="",
                     background="antigenomic",normalize=T,option="transcript",exonlevel=c(1032124,1032124,9216))
xps.rma.tc.bq103 <- validData(data.rma.tc.bq103)

## xps rma: bg=quantile="all+xxx", medpol="core"
data.rma.tc.bq24 <- rma(data.exon,"HuExonMixRMAbgqu24mp9_tc",filedir=outdir,tmpdir="",
                     background="antigenomic",normalize=T,option="transcript",exonlevel=c(245692,245692,9216))
xps.rma.tc.bq24 <- validData(data.rma.tc.bq24)

## xps rma: bg=quantile="all+xxx+unmapped", medpol="core"
data.rma.tc.bq313 <- rma(data.exon,"HuExonMixRMAbgqu313mp9_tc",filedir=outdir,tmpdir="",
                     background="antigenomic",normalize=T,option="transcript",exonlevel=c(3129276,3129276,9216))
xps.rma.tc.bq313 <- validData(data.rma.tc.bq313)

## RMA for probesets
data.rma.ps <- rma(data.exon, "MixHuExonRMAcore_ps", filedir=outdir, tmpdir="", background="antigenomic",
                   normalize=T, option="probeset", exonlevel="core")
xps.rma.ps  <- validData(data.rma.ps)

## medpol only for probesets
expr.mp.ps <- express(data.exon,"HuExonMedPolcorePS",summarize.method="medianpolish",
              summarize.select="pmonly",summarize.option="probeset",summarize.logbase="log2",
              summarize.params=c(10, 0.01, 1.0),exonlevel="core")
xps.mp.ps <- validData(expr.mp.ps)

## xps rma: bg="all", quantile="all", medpol="core"
data.rma.ps.bq16 <- rma(data.exon,"HuExonMixRMAbgqu16mp8",filedir=outdir,tmpdir="",
                    background="antigenomic",normalize=T,option="probeset",exonlevel=c(16316,16316,9216))
xps.rma.ps.bq16 <- validData(data.rma.ps.bq16)

## xps rma: bg=quantile="all+genomic+antigenomic+intron+exon+unmapped", medpol="core"
data.rma.ps.bq103 <- rma(data.exon,"HuExonMixRMAbgqu103mp9_ps",filedir=outdir,tmpdir="",
                     background="antigenomic",normalize=T,option="probeset",exonlevel=c(1032124,1032124,9216))
xps.rma.ps.bq103 <- validData(data.rma.ps.bq103)

## xps rma: bg=quantile="core+genomic+antigenomic", medpol="core"
data.rma.ps.bq10 <- rma(data.exon,"HuExonMixRMAbgqu10mp9_ps",filedir=outdir,tmpdir="",
                     background="antigenomic",normalize=T,option="probeset",exonlevel=c(107520,107520,9216))
xps.rma.ps.bq10 <- validData(data.rma.ps.bq10)

## xps rma: bg=quantile="core+genomic+antigenomic+intron", medpol="core"
data.rma.ps.bq23 <- rma(data.exon,"HuExonMixRMAbgqu23mp9_ps",filedir=outdir,tmpdir="",
                     background="antigenomic",normalize=T,option="probeset",exonlevel=c(238592,238592,9216))
xps.rma.ps.bq23 <- validData(data.rma.ps.bq23)

## xps rma: bg=quantile="all+genomic+antigenomic+intron", medpol="core"
data.rma.ps.bq24 <- rma(data.exon,"HuExonMixRMAbgqu24mp9_ps",filedir=outdir,tmpdir="",
                     background="antigenomic",normalize=T,option="probeset",exonlevel=c(245692,245692,9216))
xps.rma.ps.bq24 <- validData(data.rma.ps.bq24)

## xps rma: bg=quantile="all+genomic+antigenomic+intron+exon+unmapped+unknowntype", medpol="core"
data.rma.ps.bq313 <- rma(data.exon,"HuExonMixRMAbgqu313mp9_ps",filedir=outdir,tmpdir="",
                     background="antigenomic",normalize=T,option="probeset",exonlevel=c(3129276,3129276,9216))
xps.rma.ps.bq313 <- validData(data.rma.ps.bq313)

## xps rma: bg=quantile="all+genomic+antigenomic+intron+exon+unmapped+unknowntype+free", medpol="core"
data.rma.ps.bq31 <- rma(data.exon,"HuExonMixRMAbgqu31mp9_ps",filedir=outdir,tmpdir="",
                     background="antigenomic",normalize=T,option="probeset",exonlevel=c(3129340,3129340,9216))
xps.rma.ps.bq31 <- validData(data.rma.ps.bq31)

## xps rma: bg=quantile="all+antigenomic+intron+unmapped", medpol="core"
data.rma.ps.bq73 <- rma(data.exon,"HuExonMixRMAbgqu73mp9_ps",filedir=outdir,tmpdir="",
                     background="antigenomic",normalize=T,option="probeset",exonlevel=c(737212,737212,9216))
xps.rma.ps.bq73 <- validData(data.rma.ps.bq73)

## xps rma: bg=quantile="all+antigenomic+intron+unmapped+unknowntype", medpol="core"
data.rma.ps.bq283 <- rma(data.exon,"HuExonMixRMAbgqu283mp9_ps",filedir=outdir,tmpdir="",
                     background="antigenomic",normalize=T,option="probeset",exonlevel=c(2834364,2834364,9216))
xps.rma.ps.bq283 <- validData(data.rma.ps.bq283)

## xps rma: bg=quantile="all+antigenomic+intron+unmapped+unknowntype+free", medpol="core"
data.rma.ps.bq28 <- rma(data.exon,"HuExonMixRMAbgqu28mp9_ps",filedir=outdir,tmpdir="",
                     background="antigenomic",normalize=T,option="probeset",exonlevel=c(2834428,2834428,9216))
xps.rma.ps.bq28 <- validData(data.rma.ps.bq28)

## xps rma: bg=quantile="all+genomic+antigenomic+intron+unmapped+unknowntype+free", medpol="core"
data.rma.ps.bq286 <- rma(data.exon,"HuExonMixRMAbgqu286mp9_ps",filedir=outdir,tmpdir="",
                     background="antigenomic",normalize=T,option="probeset",exonlevel=c(2867196,2867196,9216))
xps.rma.ps.bq286 <- validData(data.rma.ps.bq286)

## xps rma: bg="all", quantile="all", medpol="all"
data.rma.ps.all <- rma(data.exon,"HuExonMixRMAAllPS",filedir=outdir,tmpdir="",background="antigenomic",
                       normalize=T,option="probeset",exonlevel="all")
xps.rma.ps.all <- validData(data.rma.ps.all)

## need to create probeset list for APT so that APT can use same probesets as XPS
# probesetList.txt
tmp <- as.data.frame(rownames(xps.rma.ps))
colnames(tmp) <- "probeset_id"
write.table(tmp, "corePSList.txt", quote=FALSE, row.names=FALSE)
# coreList.mps
writeLines(rownames(xps.rma.tc), "core.txt")
metaProbesets(scheme.exon, "core.txt", "coreList.mps", "core")
# allPSList.txt
tmp <- as.data.frame(rownames(xps.rma.ps.all))
colnames(tmp) <- "probeset_id"
write.table(tmp, "allPSList.txt", quote=FALSE, row.names=FALSE)


## DABG for transcripts
call.dabg.tc <- dabg.call(data.exon,"HuExonDABGcore",filedir=datdir,option="transcript",exonlevel="core")
# note: alpha1 and alpha2 need to be adjusted to get usable P/M/A calls for transcripts
#call.dabg.tc <- dabg.call(data.exon,"MixDABGMetacore", alpha1=???,alpha2=???)
xps.pval.tc <- validData(call.dabg.tc)

## DABG for probesets
call.dabg.ps <- dabg.call(data.exon,"HuExonDABGcorePS",filedir=outdir,option="probeset",exonlevel="core")
xps.pval.ps <- validData(call.dabg.ps)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. APT: compute rma, dabg calls
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### run on command line

## apt rma: use quantile; use corePSList.txt
apt-probeset-summarize -a rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish.expon=true -p HuEx-1_0-st-v2.r2.pgf -c HuEx-1_0-st-v2.r2.clf -b HuEx-1_0-st-v2.r2.antigenomic.bgp -s corePSList.txt -o ./rma *.CEL

## apt rma: use quantile; use coreList.mps
apt-probeset-summarize -a rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish.expon=true -p HuEx-1_0-st-v2.r2.pgf -c HuEx-1_0-st-v2.r2.clf -m coreList.mps -o ./rma_ps *.CEL

## apt rma: use allPSList.txt
apt-probeset-summarize -a rma-bg,quant-norm.sketch=0.usepm=true.bioc=true,pm-only,med-polish.expon=true -p HuEx-1_0-st-v2.r2.pgf -c HuEx-1_0-st-v2.r2.clf -b HuEx-1_0-st-v2.r2.antigenomic.bgp -s allPSList.txt -o ./rma_all *.CEL

## apt medpol only; use corePSList.txt
apt-probeset-summarize -a pm-only,med-polish.expon=true -p HuEx-1_0-st-v2.r2.pgf -c HuEx-1_0-st-v2.r2.clf -s corePSList.txt -o ./rma *.CEL

## apt medpol only; use coreList.mps
apt-probeset-summarize -a pm-only,med-polish.expon=true -p HuEx-1_0-st-v2.r2.pgf -c HuEx-1_0-st-v2.r2.clf -m coreList.mps -o ./rma_ps *.CEL

### apt detection above background call; corePSList.txt
apt-probeset-summarize -a dabg -p HuEx-1_0-st-v2.r2.pgf -c HuEx-1_0-st-v2.r2.clf -b HuEx-1_0-st-v2.r2.antigenomic.bgp -s corePSList.txt  -x 12 -o ./dabg *.CEL


### import results into R

## apt rma
apt.rma.tc <- read.delim("./rma/rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=52)
apt.rma.ps <- read.delim("./rma_ps/rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=52)
apt.rma.ps.all <- read.delim("./rma_all/rma-bg.quant-norm.pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=52)

## apt medpol
apt.mp.tc <- read.delim("./rma/pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=52)
apt.mp.ps <- read.delim("./rma_ps/pm-only.med-polish.summary.txt", row.names=1, comment.char="", skip=52)

## apt dabg
apt.dabg.ps <- read.delim("./dabg/dabg.summary.txt", row.names=1, comment.char="", skip=52)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 3. compare apt vs xps vs affy
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### compare RMA

# compare RMA obtained with xps vs apt
png(file="HuExon_XPS_APT_RMA_MvA_tc.png",width=400,height=400)
plot(log2(apt.rma.tc[,1] * xps.rma.tc[,1])/2, log2(xps.rma.tc[,1]/apt.rma.tc[,1]), main="RMA TC: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

png(file="HuExon_XPS_APT_RMA_MvA_ps.png",width=400,height=400)
plot(log2(apt.rma.ps[,1] * xps.rma.ps[,1])/2, log2(xps.rma.ps[,1]/apt.rma.ps[,1]), main="RMA PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

png(file="HuExon_XPS_APT_RMA_MvA_ps_all.png",width=400,height=400)
plot(log2(apt.rma.ps.all[,1] * xps.rma.ps.all[,1])/2, log2(xps.rma.ps.all[,1]/apt.rma.ps.all[,1]), main="RMA PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

## compare medpol obtained with xps vs apt
png(file="HuExon_XPS_APT_MedPol_MvA_tc.png",width=400,height=400)
plot(log2(apt.mp[,1] * xps.mp[,1])/2, log2(xps.mp[,1]/apt.mp[,1]), main="MedPol TC: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.001,0.001))
dev.off()

png(file="HuExon_XPS_APT_MedPol_MvA_ps.png",width=400,height=400)
plot(log2(apt.mp.ps[,1] * xps.mp.ps[,1])/2, log2(xps.mp.ps[,1]/apt.mp.ps[,1]), main="MedPol PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.001,0.001))
dev.off()

## compare rma: bg="all", quantile="all", medpol="core" obtained with xps vs apt
png(file="HuExon_XPS_APT_RMA_bg16qu16_MvA_ps.png",width=400,height=400)
plot(log2(apt.rma.ps[,1] * xps.rma.ps.bq16[,1])/2, log2(xps.rma.ps.bq16[,1]/apt.rma.ps[,1]), main="RMA: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

## compare rma: bg=quantile="all+xxx", medpol="core" obtained with xps vs apt
png(file="HuExon_XPS_APT_RMA_bg103qu103_MvA_tc.png",width=400,height=400)
plot(log2(apt.rma.tc[,1] * xps.rma.tc.bq103[,1])/2, log2(xps.rma.tc.bq103[,1]/apt.rma.tc[,1]), main="RMA TC: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

png(file="HuExon_XPS_APT_RMA_bg24qu24_MvA_tc.png",width=400,height=400)
plot(log2(apt.rma.tc[,1] * xps.rma.tc.bq24[,1])/2, log2(xps.rma.tc.bq24[,1]/apt.rma.tc[,1]), main="RMA TC: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

png(file="HuExon_XPS_APT_RMA_bg103qu103_MvA_ps.png",width=400,height=400)
plot(log2(apt.rma.ps[,1] * xps.rma.ps.bq103[,1])/2, log2(xps.rma.ps.bq103[,1]/apt.rma.ps[,1]), main="RMA PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

png(file="HuExon_XPS_APT_RMA_bg10qu10_MvA_ps.png",width=400,height=400)
plot(log2(apt.rma.ps[,1] * xps.rma.ps.bq10[,1])/2, log2(xps.rma.ps.bq10[,1]/apt.rma.ps[,1]), main="RMA PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

png(file="HuExon_XPS_APT_RMA_bg23qu23_MvA_ps.png",width=400,height=400)
plot(log2(apt.rma.ps[,1] * xps.rma.ps.bq23[,1])/2, log2(xps.rma.ps.bq23[,1]/apt.rma.ps[,1]), main="RMA PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

## compare rma: bg=quantile="all+xxx+unmapped", medpol="core" obtained with xps vs apt
png(file="HuExon_XPS_APT_RMA_bg313qu313_MvA_tc.png",width=400,height=400)
plot(log2(apt.rma.tc[,1] * xps.rma.tc.bq313[,1])/2, log2(xps.rma.tc.bq313[,1]/apt.rma.tc[,1]), main="RMA TC: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

png(file="HuExon_XPS_APT_RMA_bg313qu313_MvA_ps.png",width=400,height=400)
plot(log2(apt.rma.ps[,1] * xps.rma.ps.bq313[,1])/2, log2(xps.rma.ps.bq313[,1]/apt.rma.ps[,1]), main="RMA PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

png(file="HuExon_XPS_APT_RMA_bg31qu31_MvA_ps.png",width=400,height=400)
plot(log2(apt.rma.ps[,1] * xps.rma.ps.bq31[,1])/2, log2(xps.rma.ps.bq31[,1]/apt.rma.ps[,1]), main="RMA PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

png(file="HuExon_XPS_APT_RMA_bg73qu73_MvA_ps.png",width=400,height=400)
plot(log2(apt.rma.ps[,1] * xps.rma.ps.bq73[,1])/2, log2(xps.rma.ps.bq73[,1]/apt.rma.ps[,1]), main="RMA PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

png(file="HuExon_XPS_APT_RMA_bg283qu283_MvA_ps.png",width=400,height=400)
plot(log2(apt.rma.ps[,1] * xps.rma.ps.bq283[,1])/2, log2(xps.rma.ps.bq283[,1]/apt.rma.ps[,1]), main="RMA PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

png(file="HuExon_XPS_APT_RMA_bg28qu28_MvA_ps.png",width=400,height=400)
plot(log2(apt.rma.ps[,1] * xps.rma.ps.bq28[,1])/2, log2(xps.rma.ps.bq28[,1]/apt.rma.ps[,1]), main="RMA PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()

png(file="HuExon_XPS_APT_RMA_bg286qu286_MvA_ps.png",width=400,height=400)
plot(log2(apt.rma.ps[,1] * xps.rma.ps.bq286[,1])/2, log2(xps.rma.ps.bq286[,1]/apt.rma.ps[,1]), main="RMA PS: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",log="",ylim=c(-0.2,1))
dev.off()


### compare DABG

# compare detection calls obtained with xps vs apt
png(file="HuExon_Diff_XPS_APT_DABG_pval_ps.png",width=400,height=400)
plot(xps.pval.ps[,1], (xps.pval.ps[,1] - apt.dabg.ps[,1]), main="DABG P-Value: XPS vs APT", xlab="XPS P-Value", ylab="Difference", log="x", ylim=c(-0.00001,0.00001))
dev.off()



#------------------------------------------------------------------------------#
# Appendix: APT PLIER
#------------------------------------------------------------------------------#

## apt plier: use quantile and not sketch-quantile
apt-probeset-summarize -a plier-mm -d HG-U133_Plus_2.cdf *.CEL

# import data.frame
apt.plier <- read.delim("plier-mm.summary.txt", row.names=1, comment.char="", skip=50)
apt.plier <- apt.plier[order(rownames(apt.plier)),]

## Scatter-plot for Breast_A vs Breast_B 
png(file="U133P2_APT_PLIER_Scatter.png",width=400,height=400)
plot(apt.plier[,1], apt.plier[,2], main="APT: PLIER", xlab="Breast_A", ylab="Breast_B",log="xy")
dev.off()

## MvA-plot for Breast_A vs Breast_B 
png(file="U133P2_APT_PLIER_MvA.png",width=400,height=400)
plot(log2(apt.plier[,1] * apt.plier[,2])/2, log2(apt.plier[,1]/apt.plier[,2]), main="APT: PLIER", xlab="A = Log2(BrA*BrB)", ylab="M = Log2(BrA/BrB)",log="")
dev.off()

## MvA-plot for Breast_A vs Breast_B using RMA
png(file="U133P2_APT_RMA_MvA.png",width=400,height=400)
plot(log2(apt.rma[,1] * apt.rma[,2])/2, log2(apt.rma[,1]/apt.rma[,2]), main="APT: RMA", xlab="A = Log2(BrA*BrB)", ylab="M = Log2(BrA/BrB)",log="",ylim=c(-10,10))
dev.off()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# create metacoreList.mps using R (too slow: about 6-8hrs)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

### load library xps
library(xps)

### metacoreList.mps
#   need to create probeset list for APT so that APT can use same probesets as XPS
xps.rma <- validData(data.x.rma)
tc.id <- rownames(xps.rma)
writeLines(tc.id, "metacore.txt")

# get probeset_id, transcript_cluster_id, probe_count from probeset annotation tree
anp <- export.scheme(scheme.exon, treetype="anp", varlist="fProbesetID:fTranscriptID:fNProbes:fCrossHybType:fLevelID:fBounded", outfile="HuExon_anp.txt",as.dataframe=TRUE,verbose=FALSE)


# 1. trial: create metacoreList.mps using for-loop
mps <- vector(mode="character",length=(length(tc.id) + 1))
mps[1] <- paste("probeset_id", "transcript_cluster_id", "probeset_list", "probe_count", sep="\t")
Sys.time()
for (i in 1:100) {
   tmp <- anp[anp[,"TranscriptClusterID"] == tc.id[i],]
   tmp <- tmp[tmp[,"Level"] == "core",]
   tmp <- tmp[tmp[,"CrossHybType"] == 1,]
   tmp <- tmp[tmp[,"Bounded_NoBoundedEvidence"] == 0,]
   txt <- paste(tmp[,"ProbesetID"], collapse=" ")
   txt <- paste(tmp[1,"TranscriptClusterID"], tmp[1,"TranscriptClusterID"], txt, sum(tmp[,"ProbeCount"]), sep="\t")
   mps[i+1] <- txt
}
Sys.time()
writeLines(mps, "metacoreList.mps")


# 2. trial: create metacoreList.mps using sapply()
"mclist" <- function(x) {
   tmp <- subset(anp, anp[,"TranscriptClusterID"] == x)
   tmp <- subset(tmp, (tmp[,"Level"] == "core" & tmp[,"CrossHybType"] == 1 & tmp[,"Bounded_NoBoundedEvidence"] == 0))
   txt <- paste(tmp[,"ProbesetID"], collapse="\t")
   txt <- paste(tmp[1,"TranscriptClusterID"], tmp[1,"TranscriptClusterID"], txt, sum(tmp[,"ProbeCount"]), sep="\t")
   return(txt)
}

Sys.time()
mps <- sapply(tc.id, mclist)
Sys.time()
mps <- append(paste("probeset_id", "transcript_cluster_id", "probeset_list", "probe_count", sep="\t"), mps)
writeLines(mps, "metacoreList.mps")


#------------------------------------------------------------------------------#
