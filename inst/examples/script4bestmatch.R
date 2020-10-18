#------------------------------------------------------------------------------#
# Scripts for Exon Arrays:  best match to HuGene and HGU133P2
#
# Copyright (c) 2008-2009 Dr. Christian Stratowa, Vienna, Austria.
# All rights reserved.
#
# Note: this script is only a trial, it may have mistakes and there
#       may be better solutions, thus use it at your own risk.
#       Bug reports and refined solutions are welcome
#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Functions 
#------------------------------------------------------------------------------#

## remove rows with duplicate probesets: keep only probeset with highest percentage
uniqueframe <- function(ma) {
   maxunique <- function(id, m) {
      m <- m[which(m[,1] == id),];
      m <- m[which(m[,2] == max(m[,2])),];
      return(m[1,]);
   }

   dup <- duplicated(ma[,1])
   uni <- unique(ma[dup,1])

   ds <- NULL
   for (i in uni) {ds <- rbind(ds, maxunique(i,ma))}

   tmp <- ma[dup==F,]
   tmp <- tmp[!tmp[,1] %in% uni,]

   ds <- rbind(ds, tmp)
   ds <- ds[order(rownames(ds)),]
   return(ds)
}#uniqueframe

## return subframe ma1 containing only rownames common to ma1 and ma2
intersectframes <- function(ma1, ma2) {
   id  <- intersect(rownames(ma1), rownames(ma2))
   ok  <- match(id, rownames(ma1))
   ok  <- ok[!is.na(ok)]
   ma1 <- ma1[ok,]
   ma1 <- ma1[order(rownames(ma1)),]
   return(ma1)
}#intersectframes

## return subframe ma1 containing only rownames common to ma1 and ma2
intersectrows <- function(ma1, ma2, col1=NULL, col2=NULL, ord=NULL) {
   if (is.null(col1)) {
      colA <- rownames(ma1)
   } else {
      colA <- ma1[,col1]
   }#if

   if (is.null(col2)) {
      col2 <- rownames(ma2)
   } else {
      col2 <- ma2[,col2]
   }#if

   id  <- intersect(colA, col2)
   ok  <- match(id, colA)
   ok  <- ok[!is.na(ok)]
   ma1 <- ma1[ok,]

   if (is.null(ord)) {
      ma1 <- ma1[order(rownames(ma1)),]
   } else if (ord > 0) {
      ma1 <- ma1[order(ma1[,ord]),]
   }#if

   return(ma1)
}#intersectrowss


#------------------------------------------------------------------------------#
# Affymetrix BestMatch tables 
#------------------------------------------------------------------------------#

## read BestMatch tables
hx2hg <- read.delim("HuExVsHuGene_BestMatch.txt", row.names=3, comment.char="")
hx2hg <- hx2hg[,5:6]
hx2hg <- uniqueframe(hx2hg)
colnames(hx2hg) <- c("HuGene", "PercentX2G")

up2hx <- read.delim("U133PlusVsHuEx_BestMatch.txt", row.names=3, comment.char="")
up2hx <- up2hx[,5:6]
up2hx <- uniqueframe(up2hx)
colnames(up2hx) <- c("HuExon", "PercentU2X")

up2hg <- read.delim("U133PlusVsHuGene_BestMatch.txt", row.names=3, comment.char="")
up2hg <- up2hg[,5:6]
up2hg <- uniqueframe(up2hg)
colnames(up2hg) <- c("HuGene", "PercentU2G")


## get unique exon and genome identical to u133p2
u2x  <- intersectframes(up2hx, up2hg)
u2g  <- intersectframes(up2hg, up2hx)
u2gx <- cbind(u2g[,1,drop=F],u2x[,1,drop=F])
write.table(u2x,"U133P2_Exon.txt",sep="\t",col.names=NA)
write.table(u2g,"U133P2_Genome.txt",sep="\t",col.names=NA)
write.table(u2gx,"U133P2_Genome_Exon.txt",sep="\t",col.names=NA)


#------------------------------------------------------------------------------#
# RMA normalized expression levels 
#------------------------------------------------------------------------------#

## read RMA data from xps
xps.urma <- read.delim("MixU133P2RMA.txt", row.names=2, comment.char="")
xps.grma <- read.delim("MixHuGeneRMAcore_tc.txt", row.names=2, comment.char="")
xps.xrma <- read.delim("MixHuExonRMAcore_tc.txt", row.names=2, comment.char="")
xps.grma.bq99  <- read.delim("HuGeneMixRMAbgqu99mp9.txt", row.names=2, comment.char="")
xps.xrma.bq103 <- read.delim("HuExonMixRMAbgqu103mp9_tc.txt", row.names=2, comment.char="")
xps.urma <- xps.urma[,-1]
xps.grma <- xps.grma[,-1]
xps.xrma <- xps.xrma[,-1]
xps.grma.bq99  <- xps.grma.bq99[,-1]
xps.xrma.bq103 <- xps.xrma.bq103[,-1]

## read RMA data from apt
apt.urma <- read.delim("apt-u133p2-rma-summary.txt", row.names=1, comment.char="", skip=52)
apt.grma <- read.delim("apt-hugene-core-rma-summary.txt", row.names=1, comment.char="", skip=52)
apt.xrma <- read.delim("apt-huexon-core-rma-summary.txt", row.names=1, comment.char="", skip=52)


#------------------------------------------------------------------------------#
# Common unique core probesets 
#------------------------------------------------------------------------------#

## get core unique exon and genome identical to u133p2
tmp  <- intersectrows(u2gx, xps.grma, 1, NULL)
core <- intersectrows(tmp,  xps.xrma, 2, NULL)


#------------------------------------------------------------------------------#
# RMA subsets with common unique core probesets 
#------------------------------------------------------------------------------#

## U133P2 core
xpsu <- intersectrows(xps.urma, core, NULL, NULL)


## HuGene vs U133P2 xps core
xpsg <- intersectrows(xps.grma, core, NULL, 1, -1)
xpsg <- xpsg[match(core[,1],rownames(xpsg)),]

## HuGene vs U133P2 xps-bq99 core
bq99g <- intersectrows(xps.grma.bq99, core, NULL, 1, -1)
bq99g <- bq99g[match(core[,1],rownames(bq99g)),]

## HuGene vs U133P2 apt core
aptg <- intersectrows(apt.grma, core, NULL, 1, -1)
aptg <- aptg[match(core[,1],rownames(aptg)),]


## HuExon vs U133P2 xps core
xpsx <- intersectrows(xps.xrma, core, NULL, 2, -1)
xpsx <- xpsx[match(core[,2],rownames(xpsx)),]

## HuExon vs U133P2 xps-bq103 core
bq103x <- intersectrows(xps.xrma.bq103, core, NULL, 2, -1)
bq103x <- bq103x[match(core[,2],rownames(bq103x)),]

## HuExon vs U133P2 apt core
aptx <- intersectrows(apt.xrma, core, NULL, 2, -1)
aptx <- aptx[match(core[,2],rownames(aptx)),]


## export common core probesets
tmp <- cbind(rownames(xpsu),rownames(xpsx),rownames(xpsg))
colnames(tmp) <- c("U133P2","HuExon","HuGene")
write.table(tmp,"U133P2_Genome_Exon_core.txt",sep="\t",col.names=NA)


#------------------------------------------------------------------------------#
# Plots 
#------------------------------------------------------------------------------#

## HuGene vs U133P2 xps core
plot(xpsu[,1], xpsg[,1], log="xy")
plot(log2(xpsu[,1] * xpsg[,1])/2, log2(xpsg[,1]/xpsu[,1]), main="XPS: U133P2 vs HuGene", xlab="A = Log2(Gene*U133)", ylab="M = Log2(Gene/U133)")

## HuGene vs U133P2 xps-bq99 core
plot(xpsu[,1], bq99g[,1], log="xy")
plot(log2(xpsu[,1] * bq99g[,1])/2, log2(bq99g[,1]/xpsu[,1]), main="XPS-BQ99: U133P2 vs HuGene", xlab="A = Log2(Gene*U133)", ylab="M = Log2(Gene/U133)")

## HuGene vs U133P2 apt core
plot(xpsu[,1], aptg[,1], log="xy")
plot(log2(xpsu[,1] * aptg[,1])/2, log2(aptg[,1]/xpsu[,1]), main="APT: U133P2 vs HuGene", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",)


## HuExon vs U133P2 xps core
plot(xpsu[,1], xpsx[,1], log="xy")
plot(log2(xpsu[,1] * xpsx[,1])/2, log2(xpsx[,1]/xpsu[,1]), main="XPS: U133P2 vs HuExon", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)")

## HuExon vs U133P2 xps-bq103 core
plot(xpsu[,1], bq103x[,1], log="xy")
plot(log2(xpsu[,1] * bq103x[,1])/2, log2(bq103x[,1]/xpsu[,1]), main="XPS-BQ103: U133P2 vs HuExon", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)")

## HuExon vs U133P2 apt core
plot(xpsu[,1], aptx[,1], log="xy")
plot(log2(xpsu[,1] * aptx[,1])/2, log2(aptx[,1]/xpsu[,1]), main="APT: U133P2 vs HuExon", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)",)

plot(aptx[,1], xpsx[,1], log="xy")
plot(log2(aptx[,1] * xpsx[,1])/2, log2(xpsx[,1]/aptx[,1]), main="RMA: XPS vs APT", xlab="A = Log2(XPS*APT)", ylab="M = Log2(XPS/APT)")


## HuExon vs HuGene xps core
plot(xpsx[,1], xpsg[,1], log="xy")
plot(log2(xpsx[,1] * xpsg[,1])/2, log2(xpsg[,1]/xpsx[,1]), main="XPS: HuExon vs HuGene", xlab="A = Log2(Exon*Gene)", ylab="M = Log2(Gene/Exon)")

plot(log2(bq103x[,1] * bg99g[,1])/2, log2(bg99g[,1]/bq103x[,1]), main="XPS-BQ: HuExon vs HuGene", xlab="A = Log2(Exon*Gene)", ylab="M = Log2(Gene/Exon)")

plot(log2(xpsg[,1] * bq99g[,1])/2, log2(bq99g[,1]/xpsg[,1]), main="HuGene: XPS-BQ99 vs XPS", xlab="A = Log2(XPS-BQ16*XPS)", ylab="M = Log2(XPS-BQ16/XPS)")


## use mean of triplicates!!! see paper T.Speed
plot(xpsu[,1], bg99g[,1], log="xy")
plot(rowMeans(bg99g), log="xy")

## pairs xps
tmp <- cbind(rowMeans(xpsu), rowMeans(xpsg), rowMeans(xpsx))
colnames(tmp) <- c("XPS.U133P2","XPS.HuGene","XPS.HuExon")
png(file="XPS_U133P2_Gene_Exon_mean.png",width=440,height=440)
pairs(log2(tmp), labels=colnames(tmp))
dev.off()

tmp <- cbind(rowMeans(xpsu), rowMeans(xpsg), rowMeans(xpsx))
colnames(tmp) <- c("XPS.U133P2","XPS.HuGene","XPS.HuExon")
png(file="XPS_U133P2_Gene_Exon_mean_dot.png",width=440,height=440)
pairs(log2(tmp), labels=colnames(tmp), pch=".")
dev.off()

## pairs apt  (xpsu is identical to aptu)
tmp <- cbind(rowMeans(xpsu), rowMeans(aptg), rowMeans(aptx))
colnames(tmp) <- c("APT.U133P2","APT.HuGene","APT.HuExon")
png(file="APT_U133P2_Gene_Exon_mean.png",width=440,height=440)
pairs(log2(tmp), labels=colnames(tmp))
dev.off()

## pairs xps vs apt
tmp <- cbind(rowMeans(xpsg), rowMeans(aptg), rowMeans(xpsx), rowMeans(aptx))
colnames(tmp) <- c("XPS.HuGene","APT.HuGene","XPS.HuExon","APT.HuExon")
png(file="XPS_APT_Gene_Exon_mean.png",width=440,height=440)
pairs(log2(tmp), labels=colnames(tmp))
dev.off()


#------------------------------------------------------------------------------#
# Limma 
#------------------------------------------------------------------------------#

library(limma)

## create design matrix
tissue <- c("Breast","Breast","Breast","Prostate","Prostate","Prostate") 
design <- model.matrix(~factor(tissue)) 
colnames(design) <- c("Breast","BreastvsProstate") 
design


## limma U133P2: xpsu
tmp <- as.matrix(log2(xpsu))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
xpsu.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
xpsu.lm <- xpsu.lm[order(xpsu.lm[,"ID"]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(xpsu.lm) <- c("xpsu.ID","xpsu.logFC","xpsu.P.Value","xpsu.adj.P.Val")


## limma HuGene vs U133P2 xps core: xpsg
tmp <- as.matrix(log2(xpsg))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
xpsg.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
xpsg.lm <- xpsg.lm[match(core[,1],xpsg.lm[,1]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(xpsg.lm) <- c("xpsg.ID","xpsg.logFC","xpsg.P.Value","xpsg.adj.P.Val")

## limma HuGene vs U133P2 xps-bq99 core: bq16g
tmp <- as.matrix(log2(bq99g))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
bq99g.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
bq99g.lm <- bq99g.lm[match(core[,1],bq99g.lm[,1]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(bq99g.lm) <- c("bq99g.ID","bq99g.logFC","bq99g.P.Value","bq99g.adj.P.Val")

## limma HuGene vs U133P2 apt core: aptg
tmp <- as.matrix(log2(aptg))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
aptg.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
aptg.lm <- aptg.lm[match(core[,1],aptg.lm[,1]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(aptg.lm) <- c("aptg.ID","aptg.logFC","aptg.P.Value","aptg.adj.P.Val")


## limma HuExon vs U133P2 xps core: xpsx
tmp <- as.matrix(log2(xpsx))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
xpsx.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
xpsx.lm <- xpsx.lm[match(core[,2],xpsx.lm[,1]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(xpsx.lm) <- c("xpsx.ID","xpsx.logFC","xpsx.P.Value","xpsx.adj.P.Val")

## limma HuExon vs U133P2 xps-bq103 core: bq16x
tmp <- as.matrix(log2(bq103x))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
bq103x.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
bq103x.lm <- bq103x.lm[match(core[,2],bq103x.lm[,1]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(bq103x.lm) <- c("bq103x.ID","bq103x.logFC","bq103x.P.Value","bq103x.adj.P.Val")

## limma HuExon vs U133P2 apt core: aptx
tmp <- as.matrix(log2(aptx))
fit <- lmFit(tmp, design) 
fit <- eBayes(fit) 
aptx.lm <- topTable(fit, coef=2, n=length(rownames(tmp)), adjust="BH")
aptx.lm <- aptx.lm[match(core[,2],aptx.lm[,1]),c("ID","logFC","P.Value","adj.P.Val")]
colnames(aptx.lm) <- c("aptx.ID","aptx.logFC","aptx.P.Value","aptx.adj.P.Val")


## export all
core.lm <- cbind(xpsu.lm, xpsg.lm, bq99g.lm, aptg.lm, xpsx.lm, bq103x.lm, aptx.lm)
write.table(core.lm,"Limma_U133P2_Genome_Exon_core.txt",sep="\t",col.names=NA)


## plots HuGene
plot(core.lm[,"xpsu.logFC"], core.lm[,"xpsg.logFC"])
plot(core.lm[,"xpsu.logFC"], core.lm[,"bq99g.logFC"])
plot(core.lm[,"xpsu.logFC"], core.lm[,"aptg.logFC"])


## plots HuExon
plot(core.lm[,"xpsu.logFC"], core.lm[,"xpsx.logFC"])
plot(core.lm[,"xpsu.logFC"], core.lm[,"bq103x.logFC"])
plot(core.lm[,"xpsu.logFC"], core.lm[,"aptx.logFC"])


# plot HuGene P-Value for xps vs apt
png(file="XPS_APT_HuGene_PValue.png",width=440,height=240)
oldpar <- par(mfrow=c(1,2), oma=c(0,0,0,0), mar=c(4,4,2,1))
plot(core.lm[,"xpsg.P.Value"], core.lm[,"aptg.P.Value"])
plot(core.lm[,"xpsg.P.Value"], core.lm[,"aptg.P.Value"],xlim=c(0,0.001),ylim=c(0,0.001))
par(oldpar);
dev.off()

## plot HuGene P-Value for xps vs apt
png(file="XPS_APT_HuGene_PValue.png",width=440,height=440)
oldpar <- par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,4,2,1))
plot(core.lm[,"xpsg.P.Value"], core.lm[,"aptg.P.Value"],xlab="XPS: P-Value",ylab="APT: P-Value")
plot(core.lm[,"xpsg.P.Value"], core.lm[,"aptg.P.Value"],xlab="XPS: P-Value",ylab="APT: P-Value",xlim=c(0,0.001),ylim=c(0,0.001))
plot(core.lm[,"bq99g.P.Value"], core.lm[,"aptg.P.Value"],xlab="XPS-BQ99: P-Value",ylab="APT: P-Value")
plot(core.lm[,"bq99g.P.Value"], core.lm[,"aptg.P.Value"],xlab="XPS-BQ99: P-Value",ylab="APT: P-Value",xlim=c(0,0.001),ylim=c(0,0.001))
par(oldpar);
dev.off()

# plot HuExon P-Value for xps vs apt
png(file="XPS_APT_HuExon_PValue.png",width=440,height=240)
oldpar <- par(mfrow=c(1,2), oma=c(0,0,0,0), mar=c(4,4,2,1))
plot(core.lm[,"xpsx.P.Value"], core.lm[,"aptx.P.Value"])
plot(core.lm[,"xpsx.P.Value"], core.lm[,"aptx.P.Value"],xlim=c(0,0.000001),ylim=c(0,0.000001))
par(oldpar);
dev.off()

# plot HuExon P-Value for xps vs apt
png(file="XPS_APT_HuExon_PValue.png",width=440,height=440)
oldpar <- par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,4,2,1))
plot(core.lm[,"xpsx.P.Value"], core.lm[,"aptx.P.Value"])
plot(core.lm[,"xpsx.P.Value"], core.lm[,"aptx.P.Value"],xlim=c(0,0.000001),ylim=c(0,0.000001))
plot(core.lm[,"bq103x.P.Value"], core.lm[,"aptx.P.Value"])
plot(core.lm[,"bq103x.P.Value"], core.lm[,"aptx.P.Value"],xlim=c(0,0.000001),ylim=c(0,0.000001))
par(oldpar);
dev.off()

## plot HuExon adjusted P-Value for xps vs apt
png(file="XPS_APT_HuExon_AdjPValue.png",width=440,height=440)
oldpar <- par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,4,2,1))
plot(core.lm[,"xpsx.adj.P.Val"], core.lm[,"aptx.adj.P.Val"],xlab="XPS: Adjusted P-Value",ylab="APT: Adjusted P-Value")
plot(core.lm[,"xpsx.adj.P.Val"], core.lm[,"aptx.adj.P.Val"],xlab="XPS: Adjusted P-Value",ylab="APT: Adjusted P-Value",xlim=c(0,0.0001),ylim=c(0,0.0001))
plot(core.lm[,"bq103x.adj.P.Val"], core.lm[,"aptx.adj.P.Val"],xlab="XPS-BQ103: Adjusted P-Value",ylab="APT: Adjusted P-Value")
plot(core.lm[,"bq103x.adj.P.Val"], core.lm[,"aptx.adj.P.Val"],xlab="XPS-BQ103: Adjusted P-Value",ylab="APT: Adjusted P-Value",xlim=c(0,0.0001),ylim=c(0,0.0001))
par(oldpar);
dev.off()

# plot logFC for xps
tmp <- cbind(core.lm[,"xpsu.logFC"], core.lm[,"xpsg.logFC"], core.lm[,"xpsx.logFC"])
colnames(tmp) <- c("U133P2.logFC","HuGene.logFC","HuExon.logFC")
png(file="XPS_U133P2_HuGene_HuExon_logFC.png",width=440,height=440)
pairs(tmp, labels=colnames(tmp))
dev.off()

# plot logFC for apt
tmp <- cbind(core.lm[,"xpsu.logFC"], core.lm[,"aptg.logFC"], core.lm[,"aptx.logFC"])
colnames(tmp) <- c("U133P2.logFC","HuGene.logFC","HuExon.logFC")
png(file="APT_U133P2_HuGene_HuExon_logFC.png",width=440,height=440)
pairs(tmp, labels=colnames(tmp))
dev.off()






