#==============================================================================#
# summarize.R: summarization functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# summarize:      methods: mean, avgdiff, tukeybiweight, medianpolish
# summarize.mas4: MAS4-type "average difference" summarization
# summarize.mas5: MAS5-type "tukey biweight" summarization
# summarize.rma:  RMA-type "median polish" summarization
#==============================================================================#

"summarize" <-
function(xps.data,
         filename   = character(0),
         filedir    = getwd(),
         tmpdir     = "",
         update     = FALSE,
         select     = "none",
         method     = character(),
         option     = "transcript",
         logbase    = "0",
         exonlevel  = "",
         params     = list(),
         xps.scheme = NULL,
         add.data   = TRUE,
         verbose    = TRUE)
{
   if (is(xps.data, "DataTreeSet")) {
      set <- xpsSummarize(xps.data,
                          filename   = filename,
                          filedir    = filedir,
                          tmpdir     = tmpdir,
                          update     = update,
                          select     = select,
                          method     = method,
                          option     = option,
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

"summarize.mas4" <-
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
   set <- summarize(xps.data,
                    filename   = filename,
                    filedir    = filedir,
                    tmpdir     = tmpdir,
                    update     = update,
                    select     = "none",
                    method     = "avgdiff",
                    option     = option,
                    logbase    = "0",
                    exonlevel  = exonlevel,
                    params     = list(3.0),
                    xps.scheme = xps.scheme,
                    add.data   = add.data,
                    verbose    = verbose);
   return(set);
}#summarize.mas4

#------------------------------------------------------------------------------#

"summarize.mas5" <-
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
   set <- summarize(xps.data,
                    filename   = filename,
                    filedir    = filedir,
                    tmpdir     = tmpdir,
                    update     = update,
                    select     = "none",
                    method     = "tukeybiweight",
                    option     = option,
                    logbase    = "log2",
                    exonlevel  = exonlevel,
                    params     = list(0.03, 10.0, 2^(-20), 5.0, 0.0001, 1.0, 0.5),
                    xps.scheme = xps.scheme,
                    add.data   = add.data,
                    verbose    = verbose);
   return(set);
}#summarize.mas5

#------------------------------------------------------------------------------#

"summarize.rma" <-
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
   set <- summarize(xps.data,
                    filename   = filename,
                    filedir    = filedir,
                    tmpdir     = tmpdir,
                    update     = update,
                    select     = "pmonly",
                    method     = "medianpolish",
                    option     = option,
                    logbase    = "log2",
                    exonlevel  = exonlevel,
                    params     = list(10, 0.01, 1.0),
                    xps.scheme = xps.scheme,
                    add.data   = add.data,
                    verbose    = verbose);
   return(set);
}#summarize.rma

