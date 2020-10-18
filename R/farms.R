#==============================================================================#
# farms.R: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# farms: FARMS (Factor Analysis for Robust Microarray Summarization) preprocessing
#==============================================================================#

"farms" <-
function(xps.data,
         filename   = character(0),
         filedir    = getwd(),
         tmpdir     = "",
         normalize  = TRUE,
         weight     = 0.5,
         mu         = 0.0,
         scale      = 1.0,
         tol        = 0.00001,
         cyc        = 100,
         weighted   = TRUE,
         version    = "1.3.1",
         option     = "transcript",
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   ## check for valid version
   version <- as.integer(gsub("\\.","",version));
   if (!(version == 131 || version == 130)) {
      stop(paste("wrong version, currently only", dQuote("1.3.1"), "or", dQuote("1.3.0"),
                 "are allowed"));
   }#if

   ## check for valid parameters
   if (!(is.numeric(weight) && is.numeric(mu) && is.numeric(scale) && is.numeric(tol))) {
      stop("parameters <weight,mu,scale,tol> must be numeric");
   }#if
   if (!is.logical(weighted)) {
      stop(paste(sQuote("weighted"), "must be TRUE or FALSE"));
   }#if
   if (!(is.numeric(cyc) && cyc >= 0)) {
      stop(paste(sQuote("cyc"), "must be integer >= 0"));
   }#if
   cyc <- as.integer(cyc);

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
                           summarize.method  = "farms",
                           summarize.select  = "pmonly",
                           summarize.option  = option,
                           summarize.logbase = "log2",
                           summarize.params  = c(version, weight, mu, scale, tol, cyc, weighted),
                           exonlevel         = exonlevel,
                           xps.scheme        = xps.scheme,
                           add.data          = add.data,
                           verbose           = verbose);
      return(set);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
   }#if
}#farms
