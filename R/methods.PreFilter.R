#==============================================================================#
# methods.PreFilter.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# madFilter:
# madFilter<-:
# cvFilter:
# cvFilter<-:
# varFilter:
# varFilter<-:
# diffFilter:
# diffFilter<-:
# ratioFilter:
# ratioFilter<-:
# gapFilter:
# gapFilter<-:
# lowFilter:
# lowFilter<-:
# highFilter:
# highFilter<-:
# quantileFilter:
# quantileFilter<-:
# callFilter:
# callFilter<-:
#==============================================================================#


#------------------------------------------------------------------------------#
# PreFilter initialization:
#------------------------------------------------------------------------------#

setMethod("initialize", "PreFilter",
   function(.Object,
            mad         = list(),
            cv          = list(),
            variance    = list(),
            difference  = list(),
            ratio       = list(),
            gap         = list(),
            lothreshold = list(),
            hithreshold = list(),
            quantile    = list(),
            prescall    = list(),
            ...) 
   {
      if (debug.xps()) print("------initialize:PreFilter------")

      .Object@numfilters <- 0;

      if (length(mad))         madFilter(.Object)      <- unlist(mad);
      if (length(cv))          cvFilter(.Object)       <- unlist(cv);
      if (length(variance))    varFilter(.Object)      <- unlist(variance);
      if (length(difference))  diffFilter(.Object)     <- unlist(difference);
      if (length(ratio))       ratioFilter(.Object)    <- unlist(ratio);
      if (length(gap))         gapFilter(.Object)      <- unlist(gap);
      if (length(lothreshold)) lowFilter(.Object)      <- unlist(lothreshold);
      if (length(hithreshold)) highFilter(.Object)     <- unlist(hithreshold);
      if (length(quantile))    quantileFilter(.Object) <- unlist(quantile);
      if (length(prescall))    callFilter(.Object)     <- unlist(prescall);
 
      .Object <- callNextMethod(.Object,
                                mad         = mad,
                                cv          = cv,
                                variance    = variance,
                                difference  = difference,
                                ratio       = ratio,
                                gap         = gap,
                                lothreshold = lothreshold,
                                hithreshold = hithreshold,
                                quantile    = quantile,
                                prescall    = prescall,
                                ...);

      if (length(mad))         madFilter(.Object)      <- unlist(mad);
      if (length(cv))          cvFilter(.Object)       <- unlist(cv);
      if (length(variance))    varFilter(.Object)      <- unlist(variance);
      if (length(difference))  diffFilter(.Object)     <- unlist(difference);
      if (length(ratio))       ratioFilter(.Object)    <- unlist(ratio);
      if (length(gap))         gapFilter(.Object)      <- unlist(gap);
      if (length(lothreshold)) lowFilter(.Object)      <- unlist(lothreshold);
      if (length(hithreshold)) highFilter(.Object)     <- unlist(hithreshold);
      if (length(quantile))    quantileFilter(.Object) <- unlist(quantile);
      if (length(prescall))    callFilter(.Object)     <- unlist(prescall);

      .Object;
   }
)#initialize

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("PreFilter",
   function(object) {
      if (debug.xps()) print("------setValidity:PreFilter------")
      msg <- NULL;

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# PreFilter accessors:
#------------------------------------------------------------------------------#

setMethod("madFilter", signature(object="PreFilter"),
   function(object) object@mad
)#madFilter

setReplaceMethod("madFilter", signature(object="PreFilter", value="numeric"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:madFilter------")

      if (length(value) == 1) {
         value[2] <- 0.01;  #default epsilon
      } else if (length(value) != 2) {
         stop(paste(sQuote("mad"), "must have <cutoff,epsilon>"));
      }#if

      if (length(object@mad) == 0) {
         object@numfilters <- object@numfilters + 1;
      }#if
      object@mad <- list(cutoff  = as.double(value[1]),
                         epsilon = as.double(value[2]));
      return(object);
   }
)#madFilter<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("cvFilter", signature(object="PreFilter"),
   function(object) object@cv
)#cvFilter

setReplaceMethod("cvFilter", signature(object="PreFilter", value="numeric"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:cvFilter------")

      if (length(value) == 1) {
         value[2] <- 0.0;   #default trim
         value[3] <- 0.01;  #default epsilon
      } else if (length(value) == 2) {
         value[3] <- 0.01;  #default epsilon
      } else if (length(value) != 3) {
         stop(paste(sQuote("cv"), "must have <cutoff,trim,epsilon>"));
      }#if

      if (length(object@cv) == 0) {
         object@numfilters <- object@numfilters + 1;
      }#if
      object@cv <- list(cutoff  = as.double(value[1]),
                        trim    = as.double(value[2]),
                        epsilon = as.double(value[3]));
      return(object);
   }
)#cvFilter<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("varFilter", signature(object="PreFilter"),
   function(object) object@variance
)#varFilter

setReplaceMethod("varFilter", signature(object="PreFilter", value="numeric"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:varFilter------")

      if (length(value) == 1) {
         value[2] <- 0.0;   #default trim
         value[3] <- 0.01;  #default epsilon
      } else if (length(value) == 2) {
         value[3] <- 0.01;  #default epsilon
      } else if (length(value) != 3) {
         stop(paste(sQuote("variance"), "must have <cutoff,trim,epsilon>"));
      }#if

      if (length(object@variance) == 0) {
         object@numfilters <- object@numfilters + 1;
      }#if
      object@variance <- list(cutoff  = as.double(value[1]),
                              trim    = as.double(value[2]),
                              epsilon = as.double(value[3]));
      return(object);
   }
)#varFilter<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("diffFilter", signature(object="PreFilter"),
   function(object) object@difference
)#diffFilter

setReplaceMethod("diffFilter", signature(object="PreFilter", value="numeric"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:diffFilter------")

      if (length(value) == 1) {
         value[2] <- 0.0;   #default trim
         value[3] <- 0.01;  #default epsilon
      } else if (length(value) == 2) {
         value[3] <- 0.01;  #default epsilon
      } else if (length(value) != 3) {
         stop(paste(sQuote("difference"), "must have <cutoff,trim,epsilon>"));
      }#if

      if (length(object@difference) == 0) {
         object@numfilters <- object@numfilters + 1;
      }#if
      object@difference <- list(cutoff  = as.double(value[1]),
                                trim    = as.double(value[2]),
                                epsilon = as.double(value[3]));
      return(object);
   }
)#diffFilter<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("ratioFilter", signature(object="PreFilter"),
   function(object) object@ratio
)#ratioFilter

setReplaceMethod("ratioFilter", signature(object="PreFilter", value="numeric"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:ratioFilter------")

      if (length(value) != 1) {
         stop(paste(sQuote("ratio"), "must have <cutoff>"));
      }#if

      if (length(object@ratio) == 0) {
         object@numfilters <- object@numfilters + 1;
      }#if
      object@ratio <- list(cutoff = as.double(value[1]));
      return(object);
   }
)#ratioFilter<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("gapFilter", signature(object="PreFilter"),
   function(object) object@gap
)#gapFilter

setReplaceMethod("gapFilter", signature(object="PreFilter", value="numeric"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:gapFilter------")

      if (length(value) == 1) {
         value[2] <- 0.05;  #default window size
         value[3] <- 0.0;   #default trim value
         value[4] <- 0.01;  #default epsilon
      } else if (length(value) == 2) {
         value[3] <- 0.0;   #default trim value
         value[4] <- 0.01;  #default epsilon
      } else if (length(value) == 3) {
         value[4] <- 0.01;  #default epsilon
      } else if (length(value) != 4) {
         stop(paste(sQuote("gap"), "must have <cutoff,window,trim,epsilon>"));
      }#if

      if (length(object@gap) == 0) {
         object@numfilters <- object@numfilters + 1;
      }#if
      object@gap <- list(cutoff  = as.double(value[1]),
                         window  = as.double(value[2]),
                         trim    = as.double(value[3]),
                         epsilon = as.double(value[4]));
      return(object);
   }
)#gapFilter<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("lowFilter", signature(object="PreFilter"),
   function(object) object@lothreshold
)#lowFilter

setReplaceMethod("lowFilter", signature(object="PreFilter", value="character"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:lowFilter------")

       if (length(value) == 2) {
         value[3] <- "samples";  #default condition
      } else if (length(value) != 3) {
         stop(paste(sQuote("lothreshold"),
              "must have <cutoff,parameter,condition>"));
      }#if

      if (!isFilterCondition(value[length(value)])) {
         stop(paste(sQuote("lothreshold@condition"),
              "must be <percent,samples,mean,percentile>"));
      }#if

      if (length(object@lothreshold) == 0) {
         object@numfilters <- object@numfilters + 1;
      }#if
      object@lothreshold <- list(cutoff    = as.double(value[1]),
                                 parameter = as.double(value[2]),
                                 condition = value[3]);
      return(object);
   }
)#lowFilter<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("highFilter", signature(object="PreFilter"),
   function(object) object@hithreshold
)#highFilter

setReplaceMethod("highFilter", signature(object="PreFilter", value="character"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:highFilter------")

       if (length(value) == 2) {
         value[3] <- "samples";  #default condition
      } else if (length(value) != 3) {
         stop(paste(sQuote("hithreshold"),
              "must have <cutoff,parameter,condition>"));
      }#if

      if (!isFilterCondition(value[length(value)])) {
         stop(paste(sQuote("hithreshold@condition"),
              "must be <percent,samples,mean,percentile>"));
      }#if

      if (length(object@hithreshold) == 0) {
         object@numfilters <- object@numfilters + 1;
      }#if
      object@hithreshold <- list(cutoff    = as.double(value[1]),
                                 parameter = as.double(value[2]),
                                 condition = value[3]);
      return(object);
   }
)#highFilter<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("quantileFilter", signature(object="PreFilter"),
   function(object) object@quantile
)#quantileFilter

setReplaceMethod("quantileFilter", signature(object="PreFilter", value="numeric"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:quantileFilter------")

       if (length(value) == 1) {
         value[2] <- 0.05;  #default low quantile
         value[3] <- 0.95;  #default high quantile
      } else if (length(value) != 3) {
         stop(paste(sQuote("quantile"),
              "must have <cutoff,loquantile,hiquantile>"));
      }#if

      if (length(object@quantile) == 0) {
         object@numfilters <- object@numfilters + 1;
      }#if
      object@quantile <- list(cutoff     = as.double(value[1]),
                              loquantile = as.double(value[2]),
                              hiquantile = as.double(value[3]));
      return(object);
   }
)#quantileFilter<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("callFilter", signature(object="PreFilter"),
   function(object) object@prescall
)#callFilter

setReplaceMethod("callFilter", signature(object="PreFilter", value="character"),
   function(object, value) {
      if (debug.xps()) print("------setReplaceMethod:callFilter------")

       if (length(value) == 2) {
         value[3] <- "samples";  #default condition
      } else if (length(value) != 3) {
         stop(paste(sQuote("prescall"),
              "must have <cutoff,samples,condition>"));
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


#------------------------------------------------------------------------------#
# PreFilter methods:
#------------------------------------------------------------------------------#

#setMethod("show", signature(object="PreFilter"),
#   function(object) {
#      if (debug.xps()) print("------show.PreFilter------")

#   }
#)#show

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


