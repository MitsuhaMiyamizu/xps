#==============================================================================#
# root.files.R: new tree sets for existing root files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# root.scheme:     create new SchemeTreeSet for existing root scheme file
# root.data:       create new DataTreeSet for existing root data file
# root.merge.data: merge ROOT files containing trees from CEL-files
# root.expr:       create new ExprTreeSet for existing root expression file
# root.call:       create new CallTreeSet for existing root present call file
#==============================================================================#

"root.scheme" <-
function(rootfile = character(0),
         add.mask = FALSE)
{
   if (debug.xps()) print("------root.scheme------")

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "none");

   ## check for presence of chip name
   chipinfo <- as.list(getNameType(rootfile));
   chipname <- chipinfo$chipname;
   if (length(chipname) == 0 || nchar(chipname) < 1) {
      stop("missing chip name");
   }#if

   ## check for presence of correct chip type
   chiptype <- chipinfo$chiptype;
   TYPE <- c("GeneChip", "GenomeChip", "ExonChip");
   if (is.na(match(chiptype, TYPE))) {
      stop(paste(sQuote("chiptype"), "must be <GeneChip,GenomeChip,ExonChip>"));
   }#if

   ## get treenames and probe information
   treenames <- as.list(getTreeNames(rootfile));
   probeinfo <- as.list(getProbeInfo(rootfile));

   ## create new class
   set <- new("SchemeTreeSet",
              setname   = chipname,
              settype   = "scheme",
              rootfile  = rootfile,
              filedir   = dirname(rootfile),
              numtrees  = length(treenames),
              treenames = as.list(treenames),
              chipname  = chipname,
              chiptype  = chiptype,
              probeinfo = probeinfo);

   ## get mask for scheme
   if (add.mask) {
      chipMask(set) <- export(set,
                              treetype     = "scm",
                              varlist      = "fMask",
                              as.dataframe = TRUE,
                              verbose      = FALSE);
   }#if

   return(set);
}#root.scheme

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.data" <-
function(xps.scheme,
         rootfile = character(0),
         celnames = "*")
 {
   if (debug.xps()) print("------root.data------")

   ## check for presence of class xps.scheme
   xps.scheme <- validSchemeTreeSet(xps.scheme);

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "none");

   ## check for correct chip name
   chipname  <- getChipName(rootFile(xps.scheme));
   treetitle <- unique(getTreeNames(rootfile, "cel", gettitle=TRUE));
   if (chipname != treetitle) {
      stop(paste("wrong", sQuote("xps.scheme"), "for chip", sQuote(treetitle)));
   }#if

   ## get treenames
   treenames <- getTreeNames(rootfile, "cel");

   if (celnames[1] != "*") {
      celnames  <- paste(namePart(celnames), "cel", sep=".");
      dupnames  <- celnames[duplicated(celnames)];
      if (length(dupnames) > 0) {
         stop(paste("duplicated celnames: ", dupnames, "\n"));
      }#if

      treenames <- treenames[match(celnames, treenames)];
      if (length(treenames[is.na(treenames)]) > 0) {
         stop(paste("some", sQuote("celnames"), " are not found in",
                    sQuote(rootfile)));
      }#if
   }#if

   ## create new class
   set <- new("DataTreeSet",
              setname   = "DataSet",
              settype   = "rawdata",
              rootfile  = rootfile,
              filedir   = dirname(rootfile),
              numtrees  = length(treenames),
              treenames = as.list(treenames),
              scheme    = xps.scheme);
   return(set);
}#root.data

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.merge.data" <- 
function(xps.scheme,
         rootfiles = list(),
         celnames  = "*") 
{
   if (debug.xps()) print("------root.merge.data------")

   ## check for presence of class xps.scheme
   xps.scheme <- validSchemeTreeSet(xps.scheme);

   ## extract chipname from xps.scheme
   chipname <- chipName(xps.scheme);
   setname  <- "DataSet";

   ## check if root files contain trees for chipname
   for (rootfile in rootfiles) {
      treetitle <- getTreeNames(rootfile, treetype = "cel", setname = setname, gettitle = TRUE)[1];
      if (treetitle != chipname) {
         stop(paste("ROOT file", sQuote(rootfile), "has no trees with chipname", sQuote(chipname)));
      }#if
   }#for

   ## get correct treenames from rootfiles
   celnames  <- unique(celnames);
   tmpnames  <- vector();
   fullnames <- vector();

   for (rootfile in rootfiles) {
      treenames <- getTreeNames(rootfile, treetype = "cel", setname = setname, gettitle = FALSE);

      if (length(celnames) == 0) next;
      if (celnames[1] == "*") {
         id <- match(tmpnames, treenames);
         id <- id[!is.na(id)];
         tmpnames <- c(tmpnames, treenames);
         if (length(id) > 0) {
            tmpnames  <- tmpnames[-id];
            treenames <- treenames[-id];
         }#if

         if (length(treenames) > 0) {
            fullnames <- c(fullnames, paste(rootfile, setname, treenames, sep="/"));
         }#if
      } else {
         id <- match(treenames, celnames);
         id <- id[!is.na(id)];

         if (length(id) > 0) {
            treenames <- celnames[id];
            celnames  <- celnames[-id];
            tmpnames  <- c(tmpnames, treenames);
            fullnames <- c(fullnames, paste(rootfile, setname, treenames, sep="/"));
         }#if
      }#if
   }#for

   ## check for presence of celnames
   if (celnames[1] != "*" && length(celnames) > 0) {
      stop(paste(sQuote(celnames), "not found in", sQuote("rootfiles"), "\n"));
   }#if

   ## create new class
   set <- new("DataTreeSet",
              setname   = setname,
              settype   = "rawdata",
              rootfile  = basename(rootfiles[1]),
              filedir   = dirname(rootfiles[1]),
              numtrees  = length(fullnames),
              treenames = as.list(fullnames),
              scheme    = xps.scheme);

   return(set);
}#root.merge.data

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.expr" <-
function(xps.scheme,
         rootfile  = character(0),
         treetype  = character(0),
         treenames = "*")
 {
   if (debug.xps()) print("------root.expr------")

   ## get schemefile and chiptype
   xps.scheme <- validSchemeTreeSet(xps.scheme);
   schemefile <- validROOTFile(rootFile(xps.scheme), "SchemeTreeSet");
   chiptype   <- chipType(xps.scheme);

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "ExprTreeSet");

   ## check for correct treetype
   if (!is.na(match(treetype, NRMTYPE))) {
      setname <- "NormSet";
      settype <- "normation";
   } else if (!is.na(match(treetype, PRETYPE))) {
      setname <- "PreprocesSet";
      settype <- "preprocess";
   } else {
      stop(paste("invalid treetype", sQuote(treetype)));
   }#if

   ## get exprtype
   if (treetype == "adf") {
      exprtype <- "mas4";
   } else if (treetype == "tbw") {
      exprtype <- "mas5";
   } else if (treetype == "mdp") {
      setname  <- "PreprocesSet";
      settype  <- "preprocess";
      exprtype <- "rma";
   } else {
      exprtype <- "none";
   }#if

   ## get normtype
   if (treetype == "tmn") {
      normtype <- "mean";
   } else if (treetype == "med") {
      normtype <- "median";
   } else if (treetype == "low") {
      normtype <- "lowess";
   } else if (treetype == "sup") {
      normtype <- "supsmu";
   } else {
      normtype <- "none";
   }#if

   ## get trees with extension "treetype"
   nametypes <- getTreeNames(rootfile, treetype);
   if (length(nametypes) == 0) {
      stop(paste("root file", sQuote(rootfile), "has no trees of type", sQuote(treetype)));
   }#if

   ## get treenames
   if (treenames[1] == "*") {
      treenames <- nametypes;
   } else {
      treenames <- paste(namePart(treenames), treetype, sep=".");
      dupnames  <- treenames[duplicated(treenames)];
      if (length(dupnames) > 0) {
         stop(paste("duplicated treenames: ", dupnames, "\n"));
      }#if

      nametypes <- nametypes[match(treenames, nametypes)];
      if (length(nametypes[is.na(nametypes)]) > 0) {
         stop(paste("some", sQuote("treenames"), " are not found in",
                    sQuote(rootfile)));
      }#if
   }#if

   ## get tree names "treeset.treename.treetype"
   numtrees <- length(treenames);
   if (numtrees == 1) {
      treeset <- paste(setname, treenames, sep=".");
   } else if (numtrees > 1) {
#      treeset <- paste(rootfile, setname, treenames, sep="/");
      treeset <- paste(setname, treenames, sep="/");
   }#if

   ## export data from root file
   outfile  <- "tmp_expr.txt";

   r <- .C("ExportData",
           as.character(rootfile),
           as.character(schemefile),
           as.character(chiptype),
           as.character(settype),
           as.character(treeset),
           as.integer(numtrees),
           as.character(treetype),
           as.character("fUnitName:fLevel"),
           as.character(outfile),
           as.character("\t"),
           as.integer(FALSE),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ExportData")));
      return(NULL);
   }#if

   ## import outfile as dataframe
   data <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);

   ## create new class
   set <- new("ExprTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = dirname(rootfile),
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = xps.scheme,
              data      = data,
              params    = list(),
              exprtype  = exprtype,
              normtype  = normtype);

   return(set);
}#root.expr

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.call" <-
function(xps.scheme,
         rootfile  = character(0),
         treetype  = character(0),
         treenames = "*")
 {
   if (debug.xps()) print("------root.call------")

   ## get schemefile and chiptype
   xps.scheme <- validSchemeTreeSet(xps.scheme);
   schemefile <- validROOTFile(rootFile(xps.scheme), "SchemeTreeSet");
   chiptype   <- chipType(xps.scheme);

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "CallTreeSet");

   ## check for correct treetype
   if (!is.na(match(treetype, CALTYPE))) {
      setname <- "CallSet";
      settype <- "preprocess";
   } else {
      stop(paste("invalid treetype", sQuote(treetype)));
   }#if

   ## get calltype
   if (treetype == "dc5") {
      calltype <- "mas5";
   } else if (treetype == "dab") {
      calltype <- "dabg";
   } else {
      calltype <- "none";
   }#if

   ## get trees with extension "treetype"
   nametypes <- getTreeNames(rootfile, treetype);
   if (length(nametypes) == 0) {
      stop(paste("root file", sQuote(rootfile), "has no trees of type", sQuote(treetype)));
   }#if

   ## get treenames
   if (treenames[1] == "*") {
      treenames <- nametypes;
   } else {
      treenames <- paste(namePart(treenames), treetype, sep=".");
      dupnames  <- treenames[duplicated(treenames)];
      if (length(dupnames) > 0) {
         stop(paste("duplicated treenames: ", dupnames, "\n"));
      }#if

      nametypes <- nametypes[match(treenames, nametypes)];
      if (length(nametypes[is.na(nametypes)]) > 0) {
         stop(paste("some", sQuote("treenames"), " are not found in",
                    sQuote(rootfile)));
      }#if
   }#if

   ## get tree names "treeset.treename.treetype"
   numtrees <- length(treenames);
   if (numtrees == 1) {
      treeset <- paste(setname, treenames, sep=".");
   } else if (numtrees > 1) {
#      treeset <- paste(rootfile, setname, treenames, sep="/");
      treeset <- paste(setname, treenames, sep="/");
   }#if

   ## export pvalue data from root file
   outfile  <- "tmp_pvalue.txt";

   r <- .C("ExportData",
           as.character(rootfile),
           as.character(schemefile),
           as.character(chiptype),
           as.character(settype),
           as.character(treeset),
           as.integer(numtrees),
           as.character(treetype),
           as.character("fUnitName:fPValue"),
           as.character(outfile),
           as.character("\t"),
           as.integer(FALSE),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ExportData")));
      return(NULL);
   }#if

   ## import outfile as dataframe
   data <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);

   ## export present call from root file
   outfile  <- "tmp_call.txt";

   r <- .C("ExportData",
           as.character(rootfile),
           as.character(schemefile),
           as.character(chiptype),
           as.character(settype),
           as.character(treeset),
           as.integer(numtrees),
           as.character(treetype),
           as.character("fUnitName:fCall"),
           as.character(outfile),
           as.character("\t"),
           as.integer(FALSE),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ExportData")));
      return(NULL);
   }#if

   ## import outfile as dataframe
   dc <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);

   ## create new class
   set <- new("CallTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = dirname(rootfile),
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = xps.scheme,
              data      = data,
              params    = list(),
              calltype  = calltype,
              detcall  = dc);

   return(set);
}#root.call

#==============================================================================#
