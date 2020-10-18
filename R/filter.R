#==============================================================================#
# filter.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# prefilter: prefiltering of expression levels
# unifilter: univariate filtering of expression levels
#==============================================================================#

"prefilter" <-
function(xps.expr,
         filename   = character(0),
         filedir    = getwd(),
         filter     = NULL,
         minfilters = 999,
         logbase    = "log2",
         treename   = "PreFilter",
         xps.call   = NULL,
         verbose    = TRUE)
{
   if (is(xps.expr, "ExprTreeSet")) {
      set <- xpsPreFilter(xps.expr,
                          filename   = filename,
                          filedir    = filedir,
                          filter     = filter,
                          minfilters = minfilters,
                          logbase    = logbase,
                          treename   = treename,
                          xps.call   = xps.call,
                          verbose    = verbose);
   } else {
      stop(paste(sQuote("xps.expr"), "is not a class", sQuote("ExprTreeSet")));
   }#if

   return(set);
}#prefilter

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"unifilter" <-
function(xps.expr,
         filename   = character(0),
         filedir    = getwd(),
         filter     = NULL,
         minfilters = 999,
         logbase    = "log2",
         group      = character(0),
         treename   = "UniTest",
         xps.fltr   = NULL,
         xps.call   = NULL,
         update     = FALSE,
         verbose    = TRUE)
{
   ## check for valid update
   if (!is.logical(update)) {
      stop(paste(sQuote("update"), "must be TRUE or FALSE"));
   }#if

   if (update == TRUE) {
      tmpname  <- filename;
   } else {
      tmpname  <- paste(filename, "_ufr", sep="");
   }#if

   ## new class containing filtered data
   if (is(xps.expr, "ExprTreeSet")) {
      set <- xpsUniFilter(xps.expr,
                          filename   = tmpname,
                          filedir    = filedir,
                          update     = update,
                          filter     = filter,
                          minfilters = minfilters,
                          logbase    = logbase,
                          group      = group,
                          treename   = treename,
                          xps.fltr   = xps.fltr,
                          xps.call   = xps.call,
                          verbose    = verbose);
   } else {
      stop(paste(sQuote("xps.expr"), "is not a class", sQuote("ExprTreeSet")));
   }#if

   return(set);
}#unifilter


