#==============================================================================#
# export.R: export functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# export.scheme: 
# export.data: 
# export.expr: 
# export.call: 
# export.filter: 
# export.root: 
#==============================================================================#

"export.scheme" <-
function(xps.scheme,
         treetype     = character(0),
         varlist      = "*",
         outfile      = character(0),
         sep          = "\t",
         as.dataframe = FALSE,
         verbose      = TRUE) 
{
   if (is(xps.scheme, "SchemeTreeSet")) {
      ds <- export(xps.scheme,
                   treetype     = treetype,
                   varlist      = varlist,
                   outfile      = outfile,
                   sep          = sep,
                   as.dataframe = as.dataframe,
                   verbose      = verbose);
      return(ds);
   } else {
      stop(paste(sQuote("xps.scheme"), "is not a class", sQuote("SchemeTreeSet")));
   }#if
}#export.scheme

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"export.data" <-
function(xps.data,
         treename     = "*",
         treetype     = "cel",
         varlist      = "*",
         outfile      = character(0),
         sep          = "\t",
         as.dataframe = FALSE,
         verbose      = TRUE) 
{
   if (is(xps.data, "DataTreeSet")) {
      ds <- export(xps.data,
                   treename     = treename,
                   treetype     = treetype,
                   varlist      = varlist,
                   outfile      = outfile,
                   sep          = sep,
                   as.dataframe = as.dataframe,
                   verbose      = verbose);
      return(ds);
   } else {
      stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
   }#if
}#export.data

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"export.expr" <-
function(xps.expr,
         treename     = "*",
         treetype     = character(0),
         varlist      = "*",
         outfile      = character(0),
         sep          = "\t",
         as.dataframe = FALSE,
         verbose      = TRUE) 
{
   if (is(xps.expr, "ExprTreeSet")) {
      if (length(treetype) == 0) {
         treename1 <- as.character(treeNames(xps.expr)[1]);
         treetype  <- unlist(strsplit(treename1, "\\."))[2];
      }#if

      ds <- export(xps.expr,
                   treename     = treename,
                   treetype     = treetype,
                   varlist      = varlist,
                   outfile      = outfile,
                   sep          = sep,
                   as.dataframe = as.dataframe,
                   verbose      = verbose);
      return(ds);
   } else {
      stop(paste(sQuote("xps.expr"), "is not a class", sQuote("ExprTreeSet")));
   }#if
}#export.expr

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"export.call" <-
function(xps.call,
         treename     = "*",
         treetype     = character(0),
         varlist      = "*",
         outfile      = character(0),
         sep          = "\t",
         as.dataframe = FALSE,
         verbose      = TRUE) 
{
   if (is(xps.call, "CallTreeSet")) {
      if (length(treetype) == 0) {
         treename1 <- as.character(treeNames(xps.call)[1]);
         treetype  <- unlist(strsplit(treename1, "\\."))[2];
      }#if

      ds <- export(xps.call,
                   treename     = treename,
                   treetype     = treetype,
                   varlist      = varlist,
                   outfile      = outfile,
                   sep          = sep,
                   as.dataframe = as.dataframe,
                   verbose      = verbose);
      return(ds);
   } else {
      stop(paste(sQuote("xps.call"), "is not a class", sQuote("CallTreeSet")));
   }#if
}#export.call

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"export.filter" <-
function(xps.fltr,
         treename     = "*",
         treetype     = character(0),
         varlist      = "*",
         outfile      = character(0),
         sep          = "\t",
         as.dataframe = FALSE,
         verbose      = TRUE) 
{
   if (is(xps.fltr, "FilterTreeSet") || is(xps.fltr, "AnalysisTreeSet")) {
      if (length(treetype) == 0) {
         treename1 <- as.character(treeNames(xps.fltr)[1]);
         treetype  <- unlist(strsplit(treename1, "\\."))[2];
      }#if

      ds <- export(xps.fltr,
                   treename     = treename,
                   treetype     = treetype,
                   varlist      = varlist,
                   outfile      = outfile,
                   sep          = sep,
                   as.dataframe = as.dataframe,
                   verbose      = verbose);
      return(ds);
   } else {
      stop(paste(sQuote("xps.fltr"), "must be class", sQuote("FilterTreeSet"),
                 "or", sQuote("AnalysisTreeSet")));
   }#if
}#export.filter

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"export.root" <-
function(datafile     = character(0),
         schemefile   = character(0),
         treeset      = character(0),
         treename     = "*",
         treetype     = character(0),
         varlist      = "*",
         outfile      = character(0),
         sep          = "\t",
         as.dataframe = FALSE,
         verbose      = TRUE) 
{
   ## check for valid root data file
   datafile <- validROOTFile(datafile, "none");

   ## get chiptype (checks also for valid root scheme file)
   chiptype <- getNameType(schemefile)$chiptype;

   ## get datatype
   treetype <- validTreetype(treetype, "all");
   datatype <- getDatatype(treetype);

   ## get tree name "treeset.treename.treetype"
   treename <- paste(treeset, treename, treetype, sep=".");

   ## check for varlist
   if (nchar(varlist) < 1) {
      varlist <- "*";
   }#if

   ## check for presence of outfile=/path/outname
   outname <- "";
   if (treename == "*") {
      outname <- "PivotTable";
   } else {
      outname <- treename;
   }#if
   outfile <- validOutfile(outname, outfile);

   ## check for presence of valid separator
   validSeparator(sep);

   ## export data from root file
   r <- .C("ExportData",
           as.character(datafile),
           as.character(schemefile),
           as.character(chiptype),
           as.character(datatype),
           as.character(treename),
           as.integer(length(treename)),
           as.character(treetype),
           as.character(varlist),
           as.character(outfile),
           as.character(sep),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ExportData")));
      return(NULL);
   }#if

   ## import outfile as dataframe
   ds <- NULL;
   if (as.dataframe) {
      ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep=sep, row.names=NULL, stringsAsFactors=FALSE);
   }#if

   return(ds);
}#export.root

#==============================================================================#

