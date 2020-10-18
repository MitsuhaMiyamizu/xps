#==============================================================================#
# methods.AnalysisTreeSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# filterTreeset:
# validFilter:
# getTreeData:
# validData:
# volcanoplot:
#==============================================================================#


#------------------------------------------------------------------------------#
# AnalysisTreeSet initialization:
#------------------------------------------------------------------------------#

setMethod("initialize", "AnalysisTreeSet", 
   function(.Object, ...) {
      if (debug.xps()) print("------initialize:AnalysisTreeSet------")

      .Object <- callNextMethod(.Object, ...);
      .Object;
   }
)#initialize

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("AnalysisTreeSet",
   function(object) {
      if (debug.xps()) print("------setValidity:AnalysisTreeSet------")
      msg <- NULL;

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# AnalysisTreeSet accessors:
#------------------------------------------------------------------------------#

setMethod("filterTreeset", signature(object="AnalysisTreeSet"),
   function(object) object@fltrset
)#filterTreeset

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("validFilter", signature(object="AnalysisTreeSet"),
   function(object) validData(object@fltrset)
)#validFilter


#------------------------------------------------------------------------------#
# AnalysisTreeSet methods:
#------------------------------------------------------------------------------#

setMethod("getTreeData", signature(object="AnalysisTreeSet"),
   function(object,
            treename = "UniTest",
            treetype = "stt",
            varlist  = "fUnitName:mn1:mn2:fc:pval")
   {
      if (debug.xps()) print("------getTreeData.AnalysisTreeSet------")

      ds <- export(object,
                   treenames    = treename,
                   treetype     = treetype,
                   varlist      = varlist,
                   as.dataframe = TRUE);
      return(ds);
   }
)#getTreeData

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"dataAnalysisTreeSet" <-
function(object,
         which = "UnitName:mask") 
{
   if (debug.xps()) print("------dataAnalysisTreeSet------")

   ## check for presence of data
   data  <- object@data;
   if (min(dim(data)) == 0) {
      stop(paste("slot", sQuote("data"), "has no data"));
   }#if

   strg <- unlist(strsplit(which,":"));

   ## for which="mask" export only data with mask=1
   if (strg[1] == "mask" || (length(strg) == 2 && strg[2] == "mask")) {
      mask <- validFilter(object);
      data <- cbind(data, mask[,"FLAG",drop=F]);
      data <- data[data[,"FLAG"]==1,];
      data <- data[,-which(colnames(data) == "FLAG")];
   }#if

   ## use UnitName from column "which" as rownames
   if (!is.na(match(strg[1], colnames(data)))) {
      rownames(data) <- data[,strg[1]];
      return(data[,-which(colnames(data) == strg[1])]);
   }#if

   return(data);
}#dataAnalysisTreeSet

setMethod("validData", "AnalysisTreeSet", dataAnalysisTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("volcanoplot", signature(x="AnalysisTreeSet"),
   function(x,
            labels      = "",
            p.value     = "pval",
            mask        = FALSE,
            show.cutoff = TRUE,
            cex.text    = 0.7,
            col.text    = "blue",
            col.cutoff  = "grey",
            xlim        = NULL,
            xlab        = "Log2(Fold-Change)",
            ylab        = "-Log10(P-Value)",
            pch         = '.',
            ...) 
   {
      if (debug.xps()) print("------volcanoplot.AnalysisTreeSet------")

      PTYPE <- c("P-Value", "P-Adjusted", "P-Chance");
      ptype <- c("pval", "padj", "pcha");
      if (is.na(match(p.value, ptype))) {
         stop(paste(sQuote("p.value"), "must be <pval,padj,pcha>"));
      }#if

      show.labels <- !(is.null(labels) || labels == "");
      if (show.labels) {
         varlist <- paste("fUnitName", labels, p.value, "fc:flag", sep=":");
         if (mask) varlist <- paste(varlist, "mask", sep=":");
         ds <- export.filter(x, treetype="stt", varlist=varlist, as.dataframe=TRUE, verbose=FALSE);
      } else {
         if (mask) {
            ds <- validData(x);
         } else {
            ds <- validData(x, "UnitName");
         }#if
      }#if

      pv <- -log10(ds[,PTYPE[match(p.value, ptype)]]);
      fc <- log2(ds[,"FoldChange"]);

      if (is.null(xlim)) {
         xmax <- max(abs(fc));
         xlim <- c(-xmax, xmax);
      }#if

      ## plot 
      plot(x    = fc,
           y    = pv,
           xlim = xlim,
           xlab = xlab,
           ylab = ylab,
           pch  = pch,
           ...);

      ## plot labels
      if (show.labels) {
         LABEL <- c("fUnitName", "fName", "fSymbol", "fChromosome", "fCytoband");
         ANNOT <- c("UnitName", "GeneName", "GeneSymbol", "Chromosome", "Cytoband");
         label <- ANNOT[match(labels, LABEL)];
         text(fc, pv, labels=ds[,label], cex=cex.text, col=col.text);
      }#if

      unifltr   <- filterTreeset(x);
      cutoff.fc <- fcFilter(unifltr@filter)$cutoff;
      cutoff.pv <- unitestFilter(unifltr@filter)$cutoff;
#      cutoff.fc <- fcFilter(theFilter(unifltr))$cutoff;
#      cutoff.pv <- unitestFilter(theFilter(unifltr))$cutoff;

      if (show.cutoff & !is.null(cutoff.fc)) {
         abline(v = log2(cutoff.fc), col=col.cutoff);
         abline(v = -log2(cutoff.fc), col=col.cutoff);
      }#if
      if (show.cutoff & !is.null(cutoff.pv)) {
         abline(h = -log(cutoff.pv), col=col.cutoff);
      }#if
   }
)#volcanoplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
