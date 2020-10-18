#==============================================================================#
# dfw.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# dfw: DFW (Distribution Free Weighted Fold Change) preprocessing
#==============================================================================#

"dfw" <-
function(xps.data,
         filename   = character(0),
         filedir    = getwd(),
         tmpdir     = "",
         normalize  = TRUE,
         m          = 3,
         n          = 1,
         c          = 0.01,
         option     = "transcript",
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   ## check for valid parameters m, n, c
   if (!(is.numeric(m) && m >= 1)) {
      stop(paste(sQuote("m"), "must be numeric and >= 1"));
   }#if
   if (!(is.numeric(n) && n >= 1)) {
      stop(paste(sQuote("n"), "must be numeric and >= 1"));
   }#if
   if (!is.numeric(c)) {
      stop(paste(sQuote("c"), "must be numeric"));
   }#if

   ## normalize
   if (!is.logical(normalize)) {
      stop(paste(sQuote("normalize"), "must be TRUE or FALSE"));
   } else if (normalize == TRUE) {
      method  <- "quantile";
      normopt <- paste(option, "together", "none", sep=":");
   } else {
      method  <- NULL;
      normopt <- NULL;
   }#if

   if (is(xps.data, "DataTreeSet")) {
      set <- xpsPreprocess(xps.data,
                           filename          = filename,
                           filedir           = filedir,
                           tmpdir            = tmpdir,
                           update            = FALSE,
                           bgcorrect.method  = NULL,
                           bgcorrect.select  = "",
                           bgcorrect.option  = "",
                           bgcorrect.params  = NULL,
                           normalize.method  = method,
                           normalize.select  = "pmonly",
                           normalize.option  = normopt,
                           normalize.logbase = "0",
                           normalize.params  = c(0.0),
                           summarize.method  = "dfw",
                           summarize.select  = "pmonly",
                           summarize.option  = option,
                           summarize.logbase = "log2",
                           summarize.params  = c(m, n, c),
                           exonlevel         = exonlevel,
                           xps.scheme        = xps.scheme,
                           add.data          = add.data,
                           verbose           = verbose);
      return(set);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
   }#if
}#dfw
