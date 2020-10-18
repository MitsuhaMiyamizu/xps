#==============================================================================#
# dabg.call.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# dabg.call: detection above background call
#==============================================================================#

"dabg.call" <-
function(xps.data,
         filename   = character(0),
         filedir    = getwd(),
         alpha1     = 0.04,
         alpha2     = 0.06,
         option     = "transcript",
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   if (is(xps.data, "DataTreeSet")) {
      set <- xpsDABGCall(xps.data,
                         filename   = filename,
                         filedir    = filedir,
                         alpha1     = alpha1,
                         alpha2     = alpha2,
                         option     = option,
                         exonlevel  = exonlevel,
                         xps.scheme = xps.scheme,
                         add.data   = add.data,
                         verbose    = verbose);
      return(set);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
   }#if
}#dabg.call

