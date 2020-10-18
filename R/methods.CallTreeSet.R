#==============================================================================#
# methods.CallTreeSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# pvalData:
# pvalData<-:
# presCall:
# presCall<-:
# validPVal:
# validCall:
# attachPVal:
# removePVal:
# attachCall:
# removeCall:
# callplot:
#==============================================================================#


#------------------------------------------------------------------------------#
# CallTreeSet initialization:
#------------------------------------------------------------------------------#

setMethod("initialize", "CallTreeSet", 
   function(.Object,
            calltype = "none",
            ...) 
   {
      if (debug.xps()) print("------initialize:CallTreeSet------")

      ## set default calltype
      if (calltype == "") {
         calltype <- "mas5";
      }#if

      .Object <- callNextMethod(.Object, calltype=calltype, ...);
      .Object@calltype = calltype;
      .Object;
   }
)#initialize

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("CallTreeSet",
   function(object) {
      if (debug.xps()) print("------setValidity:CallTreeSet------")

      msg <- NULL;

      ## check for correct settype
      if (object@settype != "preprocess") {
         msg <- validMsg(msg, paste(sQuote("settype"),
                        "must be of type <preprocess>"));
      }#if

      ## check calltype
      TYPE <- c("mas5", "dabg", "ini");
      if (is.na(match(object@calltype, TYPE))) {
         msg <- validMsg(msg, paste(sQuote("normation"), "must be one of",
                         "<mas5,dabg,ini>"));
      }#if

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# CallTreeSet accessors:
#------------------------------------------------------------------------------#

setMethod("pvalData", signature(object="CallTreeSet"),
   function(object) object@data
)#pvalData

setReplaceMethod("pvalData", signature(object="CallTreeSet", value="data.frame"),
   function(object, treenames = NULL, value) {
      ## keep columns "UNIT_ID" and "UnitName"
      ds <- value[, c("UNIT_ID", "UnitName")];

      ## remove columns "UNIT_ID" and "UnitName"
      value <- value[, is.na(match(colnames(value), "UNIT_ID")),  drop=FALSE];
      value <- value[, is.na(match(colnames(value), "UnitName")), drop=FALSE];

      ## result dataframe ds
      if (is.null(treenames)) {
         if (is.null(ds)) {ds <- value} else {ds <-cbind(ds, value)}
      } else if (length(treenames) == 1) {
         value <- value[, "PVALUE", drop=FALSE];
         colnames(value) <- paste(treenames, colnames(value), sep="_");

         if (is.null(ds)) {ds <- value} else {ds <-cbind(ds, value)}
      } else {
         treenames <- namePart(treenames);
         datanames <- namePart(colnames(value));

         pos <- match(datanames, treenames);
         if (length(pos[!is.na(pos)]) != length(treenames)) {
            stop("at least one treename is not present in value");
         }#if

         value <- value[,!is.na(match(datanames, treenames))];

         if (is.null(ds)) {ds <- value} else {ds <-cbind(ds, value)}
      }#if

      ## correct extension
      treenames <- colnames(value);
      treenames <- sub("_PVALUE", "", treenames);

      object@data      <- ds;
      object@treenames <- as.list(treenames);
      object@numtrees  <- length(treenames);
      return(object);
   }
)#pvalData<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("presCall", signature(object="CallTreeSet"),
   function(object) object@detcall
)#presCall

setReplaceMethod("presCall", signature(object="CallTreeSet", value="data.frame"),
   function(object, treenames = NULL, value) {
      ## keep columns "UNIT_ID" and "UnitName"
      ds <- value[, c("UNIT_ID", "UnitName")];

      ## remove columns "UNIT_ID" and "UnitName"
      value <- value[, is.na(match(colnames(value), "UNIT_ID")),  drop=FALSE];
      value <- value[, is.na(match(colnames(value), "UnitName")), drop=FALSE];

      ## result dataframe ds
      if (is.null(treenames)) {
         if (is.null(ds)) {ds <- value} else {ds <-cbind(ds, value)}
      } else if (length(treenames) == 1) {
         value <- value[, "CALL", drop=FALSE];
         colnames(value) <- paste(treenames, colnames(value), sep="_");

         if (is.null(ds)) {ds <- value} else {ds <-cbind(ds, value)}
      } else {
         treenames <- namePart(treenames);
         datanames <- namePart(colnames(value));

         pos <- match(datanames, treenames);
         if (length(pos[!is.na(pos)]) != length(treenames)) {
            stop("at least one treename is not present in value");
         }#if

         value <- value[,!is.na(match(datanames, treenames))];

         if (is.null(ds)) {ds <- value} else {ds <-cbind(ds, value)}
      }#if

      object@detcall <- ds;
      return(object);
   }
)#presCall<-


#------------------------------------------------------------------------------#
# CallTreeSet methods:
#------------------------------------------------------------------------------#

"pvalCallTreeSet" <-
function(object,
         which = "UnitName") 
{
   if (debug.xps()) print("------pvalCallTreeSet------")

   return(validData(object, which));
}#pvalCallTreeSet

setMethod("validPVal", "CallTreeSet", pvalCallTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"callCallTreeSet" <-
function(object,
         which = "UnitName") 
{
   if (debug.xps()) print("------callCallTreeSet------")

   ## check for presence of detcall
   dcall <- object@detcall;
   ntree <- object@numtrees;
   if (min(dim(dcall)) == 0 || ncol(dcall) < ntree) {
      stop(paste("slot", sQuote("detcall"), "has no data"));
   }#if

   ## use names from column "which" as rownames
   if (!is.na(match(which, colnames(dcall)))) {
      len <- length(which(duplicated(dcall[,which])==TRUE));
      if (len == 0) {
         rownames(dcall) <- dcall[, which];
      } else {
         warning(paste("cannot use ", sQuote(which), "as row.names since it has <",
                       len, "> non-unique values.", sep=""));

         ## use names from column "UNIT_ID" as rownames
         if (!is.na(match("UNIT_ID", colnames(dcall)))) {
            rownames(dcall) <- dcall[, "UNIT_ID"];
         }#if
      }#if
   }#if

   ## get name part for matching columns
   treenames <- namePart(object@treenames);
   callnames <- namePart(colnames(dcall));

   return(dcall[,!is.na(match(callnames, treenames))]);
}#callCallTreeSet

setMethod("validCall", "CallTreeSet", callCallTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("attachPVal", signature(object="CallTreeSet"),
   function(object, treenames="*") {
      if (debug.xps()) print("------attachPVal.CallTreeSet------")

      oldtrees <- object@treenames;
      treetype <- extenPart(oldtrees);
      if (treenames[1] == "*") treenames <- oldtrees;
      if (length(treenames) > 0) {
         pvalData(object, treenames) <- export(object,
                                               treenames    = treenames,
                                               treetype     = treetype,
                                               varlist      = "fUnitName:fPValue",
                                               as.dataframe = TRUE,
                                               verbose      = FALSE);
         ## necessary since "pvalData<-" updates slots treenames, numtrees
         object@treenames <- as.list(oldtrees);
         object@numtrees  <- length(oldtrees);
      } else {
         warning("missing data tree names, data will not be added.");
      }#if
      return(object);
   }
)#attachPVal

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("removePVal", signature(object="CallTreeSet"),
   function(object) {
      if (debug.xps()) print("------removePVal.CallTreeSet------")

      return(removeData(object));
   }
)#removePVal

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("attachCall", signature(object="CallTreeSet"),
   function(object, treenames="*") {
      if (debug.xps()) print("------attachCall.CallTreeSet------")

      oldtrees <- object@treenames;
      treetype <- extenPart(oldtrees);
      if (treenames[1] == "*") treenames <- oldtrees;
      if (length(treenames) > 0) {
         presCall(object, treenames) <- export(object,
                                               treenames    = treenames,
                                               treetype     = treetype,
                                               varlist      = "fUnitName:fCall",
                                               as.dataframe = TRUE,
                                               verbose      = FALSE);
      } else {
         warning("missing data tree names, data will not be added.");
      }#if
      return(object);
   }
)#attachCall

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("removeCall", signature(object="CallTreeSet"),
   function(object) {
      if (debug.xps()) print("------removeCall.CallTreeSet------")

      object@detcall <- data.frame(matrix(nrow=0,ncol=0));
      gc(); #????
      return(object);
   }
)#removeCall

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("callplot", signature(x="CallTreeSet"),
   function(x,
            beside = TRUE,
            names  = "namepart",
            col    = c("red","green","blue"),
            legend = c("P","M","A"),
            ylim   = c(0,100),
            ylab   = "detection call [%]",
            las    = 2,
            ...) 
   {
      if (debug.xps()) print("------callplot.CallTreeSet------")

      dc <- treeInfo(x,
                     treetype = extenPart(treeNames(x)),
                     varlist  = "userinfo:fPcAbsent:fPcMarginal:fPcPresent",
                     verbose  = FALSE
                    );

      if (is.null(names))              names <- colnames(dc)
      else if (names[1] == "namepart") names <- namePart(colnames(dc))
      else                             dc    <- dc[, names, drop=F];

      barplot(as.matrix(dc[nrow(dc):1,]),
              beside    = beside,
              col       = col,
              names.arg = names,
              legend    = legend,
              ylim      = ylim,
              ylab      = ylab,
              las       = las,
              ...);
   }
)#callplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
