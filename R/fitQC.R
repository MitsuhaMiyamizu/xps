#==============================================================================#
# fitQC.R: quality control functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fitQC:  plm, rlm (medianpolish)
# fitPLM: PLM-type "probe-level model" quality control
# fitRLM: RMA-type "median polish" quality control
# rmaPLM: RMA-type "median polish" quality control
#==============================================================================#

"fitQC" <-
function(xps.data,
         ## --
         filename = character(),
         filedir  = getwd(),
         tmpdir   = "",
         update   = FALSE,
         ## --
         bgcorrect.method  = "rma",
         bgcorrect.select  = "none",
         bgcorrect.option  = "pmonly:epanechnikov",
         bgcorrect.params  = c(16384),
         ## --
         normalize.method  = "quantile",
         normalize.select  = "pmonly",
         normalize.option  = "transcript:together:none",
         normalize.logbase = "0",
         normalize.params  = c(0.0),
         ## --
         qualify.method    = "rlm",
         qualify.select    = "pmonly",
         qualify.qualopt   = "all",
         qualify.option    = "transcript",
         qualify.estimator = "huber",
         qualify.logbase   = "log2",
         qualify.params    = list(10, 0.01, 1.0),
         ## --
         reference.index   = 0,
         reference.method  = "mean",
         reference.params  = list(0.0),
         ## --
         exonlevel  = "",
         xps.scheme = NULL,
         add.data   = FALSE,
         bufsize    = 32000,
         verbose    = TRUE)
{
   if (is(xps.data, "DataTreeSet")) {
      set <- xpsQualityControl(xps.data,
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
                               qualify.method    = qualify.method,
                               qualify.select    = qualify.select,
                               qualify.qualopt   = qualify.qualopt,
                               qualify.option    = qualify.option,
                               qualify.estimator = qualify.estimator,
                               qualify.logbase   = qualify.logbase,
                               qualify.params    = qualify.params,
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
}#fitQC

#------------------------------------------------------------------------------#

#"fitPLM" <-
#function(xps.data,
#         filename   = character(),
#         filedir    = getwd(),
#         tmpdir     = "",
#         update     = FALSE,
#         qualopt    = "all",
#         option     = "transcript",
#         exonlevel  = "",
#         xps.scheme = NULL,
#         add.data   = FALSE,
#         bufsize    = 32000,
#         verbose    = TRUE)
#{
#   if (is(xps.data, "DataTreeSet")) {
#      set <- fitQC(xps.data,
#                  filename          = filename,
#                  filedir           = filedir,
#                  tmpdir            = tmpdir,
#                  update            = update,
#                  bgcorrect.method  = "rma",
#                  bgcorrect.select  = "none",
#                  bgcorrect.option  = "pmonly:epanechnikov",
#                  bgcorrect.params  = c(16384),
#                  normalize.method  = "quantile",
#                  normalize.select  = "pmonly",
#                  normalize.option  = paste(option, "together:none", sep=":"),
#                  normalize.logbase = "0",
#                  normalize.params  = c(0.0),
#                  qualify.method    = "plm",
#                  qualify.select    = "pmonly",
#                  qualify.qualopt   = qualopt,
#                  qualify.option    = option,
#                  qualify.estimator = "huber",
#                  qualify.logbase   = "log2",
#                  qualify.params    = list(10, 0.01, 1.0), #to do!!!
#                  reference.index   = 0,
#                  reference.method  = "mean",
#                  reference.params  = list(0.0),
#                  exonlevel         = exonlevel,
#                  xps.scheme        = xps.scheme,
#                  add.data          = add.data,
#                  bufsize           = bufsize,
#                  verbose           = verbose);
#      return(set);
#   } else {
#      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
#   }#if
#}#fitPLM

#------------------------------------------------------------------------------#

"fitRLM" <-
function(xps.data,
         filename   = character(),
         filedir    = getwd(),
         tmpdir     = "",
         background = "pmonly",
         normalize  = TRUE,
         qualopt    = "all",
         option     = "transcript",
         exonlevel  = "",
         params     = list(16384, 0.0, 1.0, 10, 0.01, 1),
         xps.scheme = NULL,
         add.data   = FALSE,
         bufsize    = 32000,
         verbose    = TRUE)
{
   ## skip normalization/bgrd correction dependent on qualopt
   bg.adjust <- "rma";
   if (qualopt == "raw")    {normalize <- FALSE; bg.adjust <- "none";}
   if (qualopt == "adjust") {normalize <- FALSE;}

   if (is(xps.data, "DataTreeSet")) {
      set <- fitQC(xps.data,
                  filename          = filename,
                  filedir           = filedir,
                  tmpdir            = tmpdir,
                  update            = FALSE,
#                  bgcorrect.method  = "rma",
                  bgcorrect.method  = bg.adjust,
                  bgcorrect.select  = background,
                  bgcorrect.option  = "pmonly:epanechnikov",
                  bgcorrect.params  = unlist(params[1]),
                  normalize.method  = ifelse(normalize, "quantile", "none"),
                  normalize.select  = "pmonly",
                  normalize.option  = paste(option, "together:none", sep=":"),
                  normalize.logbase = "0",
                  normalize.params  = unlist(params[2:3]),
                  qualify.method    = "rlm",
                  qualify.select    = "pmonly",
                  qualify.qualopt   = qualopt,
                  qualify.option    = option,
                  qualify.estimator = "huber",
                  qualify.logbase   = "log2",
                  qualify.params    = unlist(params[4:6]),
                  reference.index   = 0,
                  reference.method  = "mean",
                  reference.params  = list(0.0),
                  exonlevel         = exonlevel,
                  xps.scheme        = xps.scheme,
                  add.data          = add.data,
                  bufsize           = bufsize,
                  verbose           = verbose);
      return(set);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
   }#if
}#fitRLM

#------------------------------------------------------------------------------#

"rmaPLM" <-
function(xps.data,
         filename   = character(),
         filedir    = getwd(),
         tmpdir     = "",
         background = "pmonly",
         normalize  = TRUE,
         qualopt    = "all",
         option     = "transcript",
         exonlevel  = "",
         params     = list(16384, 0.0, 1.0, 10, 0.01, 1),
         xps.scheme = NULL,
         add.data   = FALSE,
         bufsize    = 32000,
         verbose    = TRUE)
{
   if (is(xps.data, "DataTreeSet")) {
      set <- fitRLM(xps.data,
                    filename   = filename,
                    filedir    = filedir,
                    tmpdir     = tmpdir,
                    background = background,
                    normalize  = normalize,
                    qualopt    = qualopt,
                    option     = option,
                    exonlevel  = exonlevel,
                    params     = params,
                    xps.scheme = xps.scheme,
                    add.data   = add.data,
                    bufsize    = bufsize,
                    verbose    = verbose);
      return(set);
   }#if
}#rmaPLM

#------------------------------------------------------------------------------#

