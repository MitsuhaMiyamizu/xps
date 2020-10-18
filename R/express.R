#==============================================================================#
# express.R: expression functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# express:
#==============================================================================#

"express" <-
function(xps.data,
         ## --
         filename = character(),
         filedir  = getwd(),
         tmpdir   = "",
         update   = FALSE,
         ## --
         bgcorrect.method  = NULL,
         bgcorrect.select  = character(),
         bgcorrect.option  = character(),
         bgcorrect.params  = list(),
         ## --
         normalize.method  = NULL,
         normalize.select  = character(),
         normalize.option  = character(),
         normalize.logbase = character(),
         normalize.params  = list(),
         ## --
         summarize.method  = NULL,
         summarize.select  = character(),
         summarize.option  = character(),
         summarize.logbase = character(),
         summarize.params  = list(),
         ## --
         reference.index   = 0,
         reference.method  = "mean",
         reference.params  = list(0.0),
         ## --
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = TRUE,
         bufsize    = 32000,
         verbose    = TRUE)
{

   if (is(xps.data, "DataTreeSet")) {
      set <- xpsPreprocess(xps.data,
                           filename          = filename,
                           filedir           = filedir,
                           tmpdir            = tmpdir,
                           update            = update,
                           bgcorrect.method  = bgcorrect.method,
                           bgcorrect.select  = bgcorrect.select,
                           bgcorrect.option  = bgcorrect.option,
                           bgcorrect.params  = bgcorrect.params,
                           normalize.method  = normalize.method,
                           normalize.select  = normalize.select,
                           normalize.option  = normalize.option,
                           normalize.logbase = normalize.logbase,
                           normalize.params  = normalize.params,
                           summarize.method  = summarize.method,
                           summarize.select  = summarize.select,
                           summarize.option  = summarize.option,
                           summarize.logbase = summarize.logbase,
                           summarize.params  = summarize.params,
                           reference.index   = reference.index,
                           reference.method  = reference.method,
                           reference.params  = reference.params,
                           exonlevel         = exonlevel,
                           xps.scheme        = xps.scheme,
                           add.data          = add.data,
                           bufsize           = bufsize,
                           verbose           = verbose);
      return(set);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
   }#if
}#express
