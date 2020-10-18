#==============================================================================#
# methods.FilterTreeSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# exprTreeset:
# callTreeset:
# validData:
#==============================================================================#


#------------------------------------------------------------------------------#
# FilterTreeSet initialization:
#------------------------------------------------------------------------------#

setMethod("initialize", "FilterTreeSet", 
   function(.Object, ...) {
      if (debug.xps()) print("------initialize:FilterTreeSet------")

      .Object <- callNextMethod(.Object, ...);
      .Object;
   }
)#initialize

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("FilterTreeSet",
   function(object) {
      if (debug.xps()) print("------setValidity:FilterTreeSet------")
      msg <- NULL;

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# FilterTreeSet accessors:
#------------------------------------------------------------------------------#

setMethod("exprTreeset", signature(object="FilterTreeSet"),
   function(object) object@exprset
)#exprTreeset

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("callTreeset", signature(object="FilterTreeSet"),
   function(object) object@callset
)#callTreeset


#------------------------------------------------------------------------------#
# FilterTreeSet methods:
#------------------------------------------------------------------------------#

setMethod("getTreeData", signature(object="FilterTreeSet"),
   function(object,
            treename = character(0),
            treetype = character(0),
            varlist  = "fUnitName:fFlag")
   {
      if (debug.xps()) print("------getTreeData.FilterTreeSet------")

      ds <- export(object,
                   treenames    = treename,
                   treetype     = treetype,
                   varlist      = varlist,
                   as.dataframe = TRUE);
      return(ds);
   }
)#getTreeData

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"dataFilterTreeSet" <-
function(object,
         which = "UnitName") 
{
   if (debug.xps()) print("------dataFilterTreeSet------")

   ## check for presence of data
   data  <- object@data;
   if (min(dim(data)) == 0) {
      stop(paste("slot", sQuote("data"), "has no data"));
   }#if

   ## use names from column "which" as rownames
   if (!is.na(match(which, colnames(data)))) {
      rownames(data) <- data[,which];
      return(data[,-which(colnames(data) == which)]);
   }#if

   return(data);
}#dataFilterTreeSet

setMethod("validData", "FilterTreeSet", dataFilterTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
