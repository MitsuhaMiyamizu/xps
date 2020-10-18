#==============================================================================#
# qualify.R: quality control functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# qualify:     plm, rlm (medianpolish)
# qualify.plm: PLM-type "probe-level model" quality control
# qualify.rlm: RMA-type "median polish" quality control
#==============================================================================#

"qualify" <-
function(xps.data,
         filename   = character(0),
         filedir    = getwd(),
         tmpdir     = "",
         update     = FALSE,
         select     = "none",
         method     = character(),
         option     = "transcript",
         logbase    = "log2",
         exonlevel  = "",
         params     = list(),
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   if (is(xps.data, "DataTreeSet")) {
      set <- xpsQualify(xps.data,
                        filename   = filename,
                        filedir    = filedir,
                        tmpdir     = tmpdir,
                        update     = update,
                        select     = select,
                        method     = method,
                        option     = paste(option, "huber", sep=":"),
                        logbase    = logbase,
                        exonlevel  = exonlevel,
                        params     = params,
                        xps.scheme = xps.scheme,
                        add.data   = add.data,
                        verbose    = verbose);
      return(set);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
   }#if
}#summarize

#------------------------------------------------------------------------------#

"qualify.plm" <-
function(xps.data,
         filename   = character(0),
         filedir    = getwd(),
         tmpdir     = "",
         update     = FALSE,
         option     = "transcript",
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   set <- qualify(xps.data,
                  filename   = filename,
                  filedir    = filedir,
                  tmpdir     = tmpdir,
                  update     = update,
                  select     = "pmonly",
                  method     = "plm",
                  option     = option,
                  logbase    = "log2",
                  exonlevel  = exonlevel,
                  params     = list(10, 0.01, 1.0), #to do!!!
                  xps.scheme = xps.scheme,
                  add.data   = add.data,
                  verbose    = verbose);
   return(set);
}#qualify.plm

#------------------------------------------------------------------------------#

"qualify.rlm" <-
function(xps.data,
         filename   = character(0),
         filedir    = getwd(),
         tmpdir     = "",
         update     = FALSE,
         option     = "transcript",
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   set <- qualify(xps.data,
                  filename   = filename,
                  filedir    = filedir,
                  tmpdir     = tmpdir,
                  update     = update,
                  select     = "pmonly",
                  method     = "rlm",
                  option     = option,
                  logbase    = "log2",
                  exonlevel  = exonlevel,
                  params     = list(10, 0.01, 1.0),
                  xps.scheme = xps.scheme,
                  add.data   = add.data,
                  verbose    = verbose);
   return(set);
}#qualify.rlm

#------------------------------------------------------------------------------#

