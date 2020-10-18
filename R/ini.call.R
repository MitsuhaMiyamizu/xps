#==============================================================================#
# ini.call.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ini.call: informative/non-informative detection call
#==============================================================================#

"ini.call" <-
function(xps.data,
         filename   = character(0),
         filedir    = getwd(),
         tmpdir     = "",
         weight     = 0.5,
         mu         = 0.0,
         scale      = 1.0,
         tol        = 0.00001,
         cyc        = 100,
         alpha1     = 0.4,
         alpha2     = 0.6,
         version    = "1.3.1",
         option     = "transcript",
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   if (is(xps.data, "DataTreeSet")) {
      set <- xpsINICall(xps.data,
                        filename   = filename,
                        filedir    = filedir,
                        tmpdir     = tmpdir,
                        weight     = weight,
                        mu         = mu,
                        scale      = scale,
                        tol        = tol,
                        cyc        = cyc,
                        alpha1     = alpha1,
                        alpha2     = alpha2,
                        version    = version,
                        option     = option,
                        exonlevel  = exonlevel,
                        xps.scheme = xps.scheme,
                        add.data   = add.data,
                        verbose    = verbose);
      return(set);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
   }#if
}#ini.call

