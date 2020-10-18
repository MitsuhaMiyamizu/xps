#==============================================================================#
# normalize.R: normalization functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# normalize:           methods: mean, median, quantile, lowess, supsmu
# normalize.constant:  mean or median normalization
# normalize.quantiles: RMA-tpye quantile normalization
# normxpress.lowess:   lowess normalization
# normxpress.supsmu:   supsmu normalization
#==============================================================================#

"normalize" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         select    = "all",
         method    = "mean",
         option    = "transcript:all",
         logbase   = "0",
         exonlevel = "",
         refindex  = 0,
         refmethod = "mean",
         params    = list(0.02, 0),
         add.data  = TRUE,
         verbose   = TRUE)
{
   if (is(xps.data, "DataTreeSet") || is(xps.data, "ExprTreeSet")) {
      if (tmpdir != "") {
         warning(paste("setting <tmpdir> will result in empty file ",
                 sQuote(paste(filename,".root", sep="")), sep=""));
      }#if 

      set <- xpsNormalize(xps.data,
                          filename  = filename,
                          filedir   = filedir,
                          tmpdir    = tmpdir,
                          update    = update,
                          select    = select,
                          method    = method,
                          option    = option,
                          logbase   = logbase,
                          exonlevel = exonlevel,
                          refindex  = refindex,
                          refmethod = refmethod,
                          params    = params,
                          add.data  = add.data,
                          verbose   = verbose);
      return(set);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class <DataTreeSet,ExprTreeSet>"));
   }#if
}#normalize

#------------------------------------------------------------------------------#

"normalize.constant" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         method    = "mean",
         logbase   = "0",
         exonlevel = "",
         refindex  = 0,
         refmethod = "mean",
         params    = list(0.02, 0),
         add.data  = TRUE,
         verbose   = TRUE)
{
   if (is(xps.data, "DataTreeSet")) {select <- "all"} else {select <- "separate"}

   set <- normalize(xps.data,
                    filename  = filename,
                    filedir   = filedir,
                    tmpdir    = tmpdir,
                    update    = update,
                    select    = select,
                    method    = method,
                    option    = "transcript:all",
                    logbase   = logbase,
                    exonlevel = exonlevel,
                    refindex  = refindex,
                    refmethod = refmethod,
                    params    = params,
                    add.data  = add.data,
                    verbose   = verbose);
   return(set);
}#normalize.constant

#------------------------------------------------------------------------------#

"normalize.quantiles" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         exonlevel = "",
         add.data  = TRUE,
         verbose   = TRUE)
{
   if (is(xps.data, "DataTreeSet")) {select <- "pmonly"} else {select <- "separate"}

   set <- normalize(xps.data,
                    filename  = filename,
                    filedir   = filedir,
                    tmpdir    = tmpdir,
                    update    = update,
                    select    = select,
                    method    = "quantile",
                    option    = "transcript:together:none",
                    logbase   = "0",
                    exonlevel = exonlevel,
                    refindex  = 1,
                    refmethod = "mean",
                    params    = list(0.0),
                    add.data  = add.data,
                    verbose   = verbose);
   return(set);
}#normalize.quantiles

#------------------------------------------------------------------------------#

"normalize.lowess" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         logbase   = "log2",
         exonlevel = "",
         refindex  = 0,
         refmethod = "mean",
         params    = list(0.67, 3, 0.0, 0.0),
         add.data  = TRUE,
         verbose   = TRUE)
{
   if (is(xps.data, "DataTreeSet")) {select <- "all"} else {select <- "separate"}

   set <- normalize(xps.data,
                    filename  = filename,
                    filedir   = filedir,
                    tmpdir    = tmpdir,
                    update    = update,
                    select    = select,
                    method    = "lowess",
                    option    = "transcript:all",
                    logbase   = logbase,
                    exonlevel = exonlevel,
                    refindex  = refindex,
                    refmethod = refmethod,
                    params    = params,
                    add.data  = add.data,
                    verbose   = verbose);
   return(set);
}#normalize.lowess

#------------------------------------------------------------------------------#

"normalize.supsmu" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         logbase   = "log2",
         exonlevel = "",
         refindex  = 0,
         refmethod = "mean",
         params    = list(0.0, 0.0, 0.0, 0.0),
         add.data  = TRUE,
         verbose   = TRUE)
{
   if (is(xps.data, "DataTreeSet")) {select <- "all"} else {select <- "separate"}

   set <- normalize(xps.data,
                    filename  = filename,
                    filedir   = filedir,
                    tmpdir    = tmpdir,
                    update    = update,
                    select    = select,
                    method    = "supsmu",
                    option    = "transcript:all",
                    logbase   = logbase,
                    exonlevel = exonlevel,
                    refindex  = refindex,
                    refmethod = refmethod,
                    params    = params,
                    add.data  = add.data,
                    verbose   = verbose);
   return(set);
}#normalize.supsmu

