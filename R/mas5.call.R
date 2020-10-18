#==============================================================================#
# mas5.call.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mas5.call: MAS5 detection call
#==============================================================================#

"mas5.call" <-
function(xps.data,
         filename         = character(0),
         filedir          = getwd(),
         tmpdir           = "",
         tau              = 0.015,
         alpha1           = 0.04,
         alpha2           = 0.06,
         ignore.saturated = TRUE,
         bgcorrect.option = "none",
         option           = "transcript",
         exonlevel        = "",
         xps.scheme       = NULL,
         add.data         = TRUE,
         verbose          = TRUE)
{
   if (is(xps.data, "DataTreeSet")) {
      set <- xpsMAS5Call(xps.data,
                         filename         = filename,
                         filedir          = filedir,
                         tmpdir           = tmpdir,
                         tau              = tau,
                         alpha1           = alpha1,
                         alpha2           = alpha2,
                         ignore.saturated = ignore.saturated,
                         bgcorrect.option = bgcorrect.option,
                         option           = option,
                         exonlevel        = exonlevel,
                         xps.scheme       = xps.scheme,
                         add.data         = add.data,
                         verbose          = verbose);
      return(set);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
   }#if
}#mas5.call

