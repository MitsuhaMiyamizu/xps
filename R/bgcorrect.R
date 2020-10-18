#==============================================================================#
# bgcorrect.R: background correction functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# bgcorrect:      correction types: sector, weightedsector, rma, gccontent
# bgcorrect.mas4: MAS4-type "sector background" using 4x4 sectors
# bgcorrect.mas5: MAS5-type "weigthed sector background" using 4x4 sectors
# bgcorrect.rma:  RMA-type global background correction
# bgcorrect.gc:   "GC-content background" using anti/genomic probes or MMs 
#==============================================================================#

"bgcorrect" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         select    = "none",
         method    = character(0),
         option    = character(0),
         exonlevel = "",
         params    = list(),
         verbose   = TRUE)
{
   if (is(xps.data, "DataTreeSet")) {
      if (tmpdir != "") {
         stop(paste("setting <tmpdir> will result in empty file ",
              sQuote(paste(filename,".root", sep="")),
              ":  <tmpdir> will be removed in the next release of xps!", sep=""));
      }#if 

      set <- xpsBgCorrect(xps.data,
                          filename  = filename,
                          filedir   = filedir,
                          tmpdir    = tmpdir,
                          update    = update,
                          select    = select,
                          method    = method,
                          option    = option,
                          exonlevel = exonlevel,
                          params    = params,
                          verbose   = verbose);
      return(set);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
   }#if
}#bgcorrect

#------------------------------------------------------------------------------#

"bgcorrect.mas4" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         select    = "all",
         exonlevel = "",
         verbose   = TRUE)
{
   set <- bgcorrect(xps.data,
                    filename  = filename,
                    filedir   = filedir,
                    tmpdir    = tmpdir,
                    update    = update,
                    select    = select,
                    method    = "sector",
                    option    = "subtractbg",
                    exonlevel = exonlevel,
                    params    = c(0.02, 4, 4, 0),
                    verbose   = verbose);
   return(set);
}#bgcorrect.mas4

#------------------------------------------------------------------------------#

"bgcorrect.mas5" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         select    = "both",
         exonlevel = "",
         verbose   = TRUE)
{
   set <- bgcorrect(xps.data,
                    filename  = filename,
                    filedir   = filedir,
                    tmpdir    = tmpdir,
                    update    = update,
                    select    = select,
                    method    = "weightedsector",
                    option    = "correctbg",
                    exonlevel = exonlevel,
                    params    = c(0.02, 4, 4, 0, 100, 0.5),
                    verbose   = verbose);
   return(set);
}#bgcorrect.mas5

#------------------------------------------------------------------------------#

"bgcorrect.rma" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         select    = "none",
         exonlevel = "",
         verbose   = TRUE)
{
   set <- bgcorrect(xps.data, 
                    filename  = filename,
                    filedir   = filedir,
                    tmpdir    = tmpdir,
                    update    = update,
                    select    = select,
                    method    = "rma",
                    option    = "pmonly:epanechnikov",
                    exonlevel = exonlevel,
                    params    = c(16384),
                    verbose   = verbose);
   return(set);
}#bgcorrect.rma

#------------------------------------------------------------------------------#

"bgcorrect.gc" <-
function(xps.data,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         select    = "antigenomic",
         exonlevel = "",
         verbose   = TRUE)
{
   set <- bgcorrect(xps.data,
                    filename  = filename,
                    filedir   = filedir,
                    tmpdir    = tmpdir,
                    update    = update,
                    select    = select,
                    method    = "gccontent",
                    option    = "attenuatebg",
                    exonlevel = exonlevel,
                    params    = c(0.4, 0.005, -1.0),
                    verbose   = verbose);
   return(set);
}#bgcorrect.gc

