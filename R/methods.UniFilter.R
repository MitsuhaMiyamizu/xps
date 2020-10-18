#==============================================================================#
# methods.UniFilter.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# fcFilter:
# fcFilter<-:
# unitestFilter:
# unitestFilter<-:
# callFilter:
# callFilter<-:
#==============================================================================#


#------------------------------------------------------------------------------#
# UniFilter initialization:
#------------------------------------------------------------------------------#

setMethod("initialize", "UniFilter", 
   function(.Object, 
            foldchange = list(),
            prescall   = list(),
            unifilter  = list(),
            unitest    = list(),
            ...) 
   {
      if (debug.xps()) print("------initialize:UniFilter------")

      .Object@numfilters <- 0;

      if (length(foldchange)) fcFilter(.Object)      <- unlist(foldchange);
      if (length(prescall))   callFilter(.Object)    <- unlist(prescall);
      if (length(unifilter))  unitestFilter(.Object) <- unlist(unifilter);
      if (length(unitest))    uniTest(.Object)       <- unlist(unitest);

      .Object <- callNextMethod(.Object,
                                foldchange = foldchange,
                                prescall   = prescall,
                                unifilter  = unifilter,
                                unitest    = unitest,
                                ...);

      if (length(foldchange)) fcFilter(.Object)      <- unlist(foldchange);
      if (length(prescall))   callFilter(.Object)    <- unlist(prescall);
      if (length(unifilter))  unitestFilter(.Object) <- unlist(unifilter);
      if (length(unitest))    uniTest(.Object)       <- unlist(unitest);

      .Object;
   }
)#initialize

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("UniFilter",
   function(object) {
      if (debug.xps()) print("------setValidity:UniFilter------")
      msg <- NULL;

      ## check unitest
      if (!length(object@unitest)) {
         msg <- validMsg(msg, paste(sQuote("UniFilter"),
                                   "must at least contain <unitest>"));
      }#if

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# UniFilter accessors:
#------------------------------------------------------------------------------#

setMethod("fcFilter", signature(object="UniFilter"),
   function(object) object@foldchange
)#fcFilter

setReplaceMethod("fcFilter", signature(object="UniFilter", value="character"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:fcFilter------")

      if (length(value) == 1) {
         value[2] <- "both";  #default direction
      } else if (length(value) != 2) {
         stop(paste(sQuote("foldchange"), "must have <cutoff,direction>"));
      }#if

      DIR <- c("both", "up", "down");
      if (is.na(match(value[2], DIR))) {
         stop(paste(sQuote("foldchange@direction"), "must be <up,down,both>"));
      }#if

      if (length(object@foldchange) == 0) {
         object@numfilters <- object@numfilters + 1;
      }#if
      object@foldchange <- list(cutoff    = as.double(value[1]),
                                direction = value[2]);
      return(object);
   }
)#fcFilter<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("callFilter", signature(object="UniFilter"),
   function(object) object@prescall
)#callFilter

setReplaceMethod("callFilter", signature(object="UniFilter", value="character"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:callFilter------")

       if (length(value) == 2) {
         value[3] <- "samples";  #default condition
      } else if (length(value) != 3) {
         stop(paste(sQuote("prescall"), "must have <cutoff,samples,condition>"));
      }#if

      if (!isFilterCondition(value[length(value)])) {
         stop(paste(sQuote("prescall@condition"), "must be <percent,samples>"));
      }#if

      if (length(object@prescall) == 0) {
         object@numfilters <- object@numfilters + 1;
      }#if
      object@prescall <- list(cutoff  = as.double(value[1]),
                              samples = as.double(value[2]),
                              condition = value[3]);
      return(object);
   }
)#callFilter<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("unitestFilter", signature(object="UniFilter"),
   function(object) object@unifilter
)#unitestFilter

setReplaceMethod("unitestFilter", signature(object="UniFilter", value="character"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:unitestFilter------")

      if (length(value) == 1) {
         value[2] <- "pval";  #default variable
      } else if (length(value) != 2) {
         stop(paste(sQuote("unifilter"), "must have <cutoff,variable>"));
      }#if

      VAR <- c("stat", "pval", "pcha", "padj");
      if (is.na(match(value[2], VAR))) {
         stop(paste(sQuote("prescall@variable"), "must be <stat,pval,pcha,padj>"));
      }#if

      if (length(object@unitest) == 0) {
         object@numfilters <- object@numfilters + 1;
      }#if
      object@unifilter <- list(cutoff   = as.double(value[1]),
                               variable = value[2]);
      return(object);
   }
)#unitestFilter<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("uniTest", signature(object="UniFilter"),
   function(object) object@unitest
)#uniTest

setReplaceMethod("uniTest", signature(object="UniFilter", value="character"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:uniTest------")

      if (length(value) == 1) {
         value[2] <- "two.sided";
         value[3] <- "none";
         value[4] <- 0;
         value[5] <- 0.0;
         value[6] <- FALSE;
         value[7] <- 0.95;
         value[8] <- FALSE;
      } else if (length(value) != 8) {
         stop(paste(sQuote("unitest"),
              "must have <type,alternative,correction,numperm,mu,paired,conflevel,varequ>"));
      }#if

      ## check type
      TYPE <- c("normal.test", "t.test");
#to do      TYPE <- c("normal.test", "t.test", "wilcox.test", "var.test");
      if (is.na(match(value[1], TYPE))) {
         stop(paste(sQuote("unitest@type"),
              "must be <normal.test,t.test>"));
#to do              "must be <normal.test,t.test,wilcox.test,var.test>"));
      }#if

      ## check alternative
      ALT <- c("two.sided", "less", "greater");
      if (is.na(match(value[2], ALT))) {
         stop(paste(sQuote("unitest@alternative"), "must be <two.sided,less,greater>"));
      }#if

      ## check correction
      COR <- c("none", "bonferroni", "BH", "BY", "fdr", "hochberg", "holm", "wy");
      if (is.na(match(value[3], COR))) {
         stop(paste(sQuote("unitest@correction"),
              "must be <none,bonferroni,BH,BY,fdr,hochberg,holm,wy>"));
      }#if

      ## check paired
      if (is.na(as.logical(value[6]))) {
         stop(paste(sQuote("unitest@paired"), "must be <TRUE,FALSE>"));
      }#if

      ## check conflevel
      if (!(value[7] >= 0 && value[7] <= 1)) {
         stop(paste(sQuote("unitest@conflevel"), "must be in range [0,1]"));
      }#if

      ## check varequ
      if (is.na(as.logical(value[8]))) {
         stop(paste(sQuote("unitest@varequ"), "must be <TRUE,FALSE>"));
      }#if

      object@unitest <- list(type        = as.character(value[1]),
                             alternative = as.character(value[2]),
                             correction  = as.character(value[3]),
                             numperm     = as.integer(value[4]),
                             mu          = as.double(value[5]),
                             paired      = as.logical(value[6]),
                             conflevel   = as.double(value[7]),
                             varequ      = as.logical(value[8]));
      return(object);
   }
)#uniTest<-


#------------------------------------------------------------------------------#
# UniFilter methods:
#------------------------------------------------------------------------------#

#setMethod("show", signature(object="UniFilter"),
#   function(object) {
#      if (debug.xps()) print("------show.UniFilter------")

#   }
#)#show

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


