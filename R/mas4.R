#==============================================================================#
# mas4.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mas4: MAS4 preprocessing
#==============================================================================#

"mas4" <-
function(xps.data,
         filename   = character(0),
         filedir    = getwd(),
         tmpdir     = "",
         normalize  = FALSE,
         sc         = 500,
         option     = "transcript",
         exonlevel  = "",
         update     = FALSE,
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   ## check for valid normalize
   if (!is.logical(normalize)) {
      stop(paste(sQuote("normalize"), "must be TRUE or FALSE"));
   }#if

   ## check for valid update
   if (!is.logical(update)) {
      stop(paste(sQuote("update"), "must be TRUE or FALSE"));
   }#if

   if (update == TRUE) {
      tmpname  <- filename;
   } else {
      tmpname  <- paste(filename, "_adf", sep="");
   }#if

   ## new class containing MAS5 preprocessed data
   if (is(xps.data, "DataTreeSet")) {
      set <- xpsMAS4(xps.data,
                     filename   = tmpname,
                     filedir    = filedir,
                     tmpdir     = tmpdir,
                     option     = option,
                     exonlevel  = exonlevel,
                     xps.scheme = xps.scheme,
                     add.data   = add.data,
                     verbose    = verbose);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
   }#if

   ## normalize MAS4 data using trimmed mean
   if (normalize) {

      ## check for valid sc
      if (!is.numeric(sc)) {
         stop(paste("parameter", sQuote("sc"), "must be numeric"));
      }#if

      ## new class containing scaled MAS5 data
      set <- xpsNormalize(set,
                          filename  = filename,
                          filedir   = filedir,
                          tmpdir    = tmpdir,
                          update    = update,
                          select    = "separate",
                          method    = "mean",
                          logbase   = "0",
#old                          option    = "all",
                          option    = paste(option, "all", sep=":"),
                          exonlevel = exonlevel,
                          refindex  = 1,
                          refmethod = "mean",
                          params    = c(0.02, sc),
                          add.data  = add.data,
                          verbose   = verbose);
   }#if

   return(set);
}#mas4

