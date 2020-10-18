#==============================================================================#
# methods.ExprTreeSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# exprType:
# exprType<-:
# normType:
# normType<-:
# exprs:
# exprs<-:
# se.exprs:
# validExpr:
# validSE:
# attachExpr:
# removeExpr:
# xpsNormalize:
# xpsPreFilter:
# xpsUniFilter:
# corplot:
# madplot:
# mvaplot:
# nuseplot:
# pcaplot:
# rleplot:
#==============================================================================#


#------------------------------------------------------------------------------#
# ExprTreeSet initialization:
#------------------------------------------------------------------------------#

setMethod("initialize", "ExprTreeSet", 
   function(.Object,
            exprtype = "none",
            normtype = "none",
            ...) 
   {
      if (debug.xps()) print("------initialize:ExprTreeSet------")

      ## set default expression/normation type
      if (exprtype == "") exprtype <- "none";
      if (normtype == "") normtype <- "none";

      .Object <- callNextMethod(.Object,
                                exprtype = exprtype,
                                normtype = normtype,
                                ...);
      .Object@exprtype = exprtype;
      .Object@normtype = normtype;
      .Object;
   }
)#initialize

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("ExprTreeSet",
   function(object) {
      if (debug.xps()) print("------setValidity:ExprTreeSet------")

      msg <- NULL;

      ## check for correct settype
      if (!(object@settype == "preprocess" ||
            object@settype == "normation")) {
         msg <- validMsg(msg, paste(sQuote("settype"),
                        "must be of type <preprocess, normation>"));
      }#if

      ## check expression method
      TYPE <- c("none", "rma", "mas4", "mas5", "firma", "custom");
      if (is.na(match(object@exprtype, TYPE))) {
         msg <- validMsg(msg, paste(sQuote("exprtype"), "must be one of",
                         "<none,rma,mas4,mas5,firma,custom>"));
      }#if

      ## check normation method
      TYPE <- c("none", "mean", "median", "lowess", "supsmu");
      if (is.na(match(object@normtype, TYPE))) {
         msg <- validMsg(msg, paste(sQuote("normtype"), "must be one of",
                         "<none,mean,median,lowess,supsmu>"));
      }#if

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# ExprTreeSet accessors:
#------------------------------------------------------------------------------#

setMethod("exprType", signature(object="ExprTreeSet"),
   function(object) object@exprtype
)#exprType

setReplaceMethod("exprType", signature(object="ExprTreeSet", value="character"),
   function(object, value) {
      object@exprtype <- value;
      return(object);
   }
)#exprType<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("normType", signature(object="ExprTreeSet"),
   function(object) object@normtype
)#normType

setReplaceMethod("normType", signature(object="ExprTreeSet", value="character"),
   function(object, value) {
      object@normtype <- value;
      return(object);
   }
)#normType<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("exprs", signature(object="ExprTreeSet"),
   function(object) object@data
)#exprs

setReplaceMethod("exprs", signature(object="ExprTreeSet", value="data.frame"),
   function(object, treenames = NULL, value) {
      ## keep columns "UNIT_ID" and "UnitName"
      ds <- value[, c("UNIT_ID", "UnitName")];

      ## remove columns "UNIT_ID" and "UnitName"
      value <- value[, is.na(match(colnames(value), "UNIT_ID")),  drop=FALSE];
      value <- value[, is.na(match(colnames(value), "UnitName")), drop=FALSE];

      ## result dataframe ds
      if (is.null(treenames)) {
         if (is.null(ds)) {ds <- value} else {ds <-cbind(ds, value)}
      } else if (length(treenames) == 1) {
         value <- value[, "LEVEL", drop=FALSE];
         colnames(value) <- paste(treenames, colnames(value), sep="_");

         if (is.null(ds)) {ds <- value} else {ds <-cbind(ds, value)}
      } else {
         treenames <- namePart(treenames);
         datanames <- namePart(colnames(value));

         pos <- match(datanames, treenames);
         if (length(pos[!is.na(pos)]) != length(treenames)) {
            stop("at least one treename is not present in value");
         }#if

         value <- value[,!is.na(match(datanames, treenames))];

         if (is.null(ds)) {ds <- value} else {ds <-cbind(ds, value)}
      }#if

      ## correct extension
      treenames <- colnames(value);
      treenames <- sub("_LEVEL", "", treenames);

      object@data      <- ds;
      object@treenames <- as.list(treenames);
      object@numtrees  <- length(treenames);
      return(object);
   }
)#exprs<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("se.exprs", signature(object="ExprTreeSet"),
   function(object) {
      if (debug.xps()) print("------se.exprs.ExprTreeSet------")

      se <- export(object,
                   treetype     = extenPart(object@treenames),
                   varlist      = "fUnitName:fStdev",
                   as.dataframe = TRUE,
                   verbose      = FALSE);
      return(se);
   }
)#se.exprs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"exprExprTreeSet" <-
function(object,
         which = "UnitName") 
{
   if (debug.xps()) print("------exprExprTreeSet------")

   return(validData(object, which));
}#exprExprTreeSet

setMethod("validExpr", "ExprTreeSet", exprExprTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"seExprTreeSet" <-
function(object, 
         which = "UnitName") 
{
   if (debug.xps()) print("------seExprTreeSet------")

   ## get se
   se <- se.exprs(object);

   ## get optional subset of columns
   strg <- unlist(strsplit(which,":"));
   if (length(strg) == 2) {
      id   <- grep(strg[2], colnames(se));
      se <- se[,c(strg[1], colnames(se)[id])];
   }#if

   ## use names from column "which" as rownames
   if (!is.na(match(strg[1], colnames(se)))) {
      len <- length(which(duplicated(se[,strg[1]]) == TRUE));
      if (len == 0) {
         rownames(se) <- se[, strg[1]];
      } else {
         warning(paste("cannot use ", sQuote(strg[1]), "as row.names since it has <",
                       len, "> non-unique values.", sep=""));

         ## use names from column "UNIT_ID" as rownames
         if (!is.na(match("UNIT_ID", colnames(se)))) {
            rownames(se) <- se[, "UNIT_ID"];
         }#if
      }#if
   }#if

   ## get name part for matching columns
   treenames <- namePart(object@treenames);
   datanames <- namePart(colnames(se));

   return(se[,!is.na(match(datanames, treenames))]);
}#seExprTreeSet

setMethod("validSE", "ExprTreeSet", seExprTreeSet);


#------------------------------------------------------------------------------#
# ExprTreeSet methods:
#------------------------------------------------------------------------------#

setMethod("attachExpr", signature(object="ExprTreeSet"),
   function(object, treenames="*") {
      if (debug.xps()) print("------attachExpr.ExprTreeSet------")

      oldtrees <- object@treenames;
      treetype <- extenPart(oldtrees);
      if (treenames[1] == "*") treenames <- oldtrees;
      if (length(treenames) > 0) {
         exprs(object, treenames) <- export(object,
                                            treenames    = treenames,
                                            treetype     = treetype,
                                            varlist      = "fUnitName:fLevel",
                                            as.dataframe = TRUE,
                                            verbose      = FALSE);
         ## necessary since "expr()<-" updates slots treenames, numtrees
         object@treenames <- as.list(oldtrees);
         object@numtrees  <- length(oldtrees);
      } else {
         warning("missing data tree names, data will not be added.");
      }#if
      return(object);
   }
)#attachExpr

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("removeExpr", signature(object="ExprTreeSet"),
   function(object) {
      if (debug.xps()) print("------removeExpr.ExprTreeSet------")

      return(removeData(object));
   }
)#removeExpr

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"normalize.ExprTreeSet" <-
function(object,
         filename  = character(0),
         filedir   = getwd(),
         tmpdir    = "",
         update    = FALSE,
         select    = "separate",
         method    = "mean",
         option    = "all",
         logbase   = "0",
         exonlevel = "",
         refindex  = 0,
         refmethod = "mean",
         params    = list(), 
         add.data  = TRUE,
         verbose   = TRUE)
{
   if (debug.xps()) print("------normalize.ExprTreeSet------")

   ## get schemefile and chiptype for object
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chiptype   <- chipType(object);

   ## treenames to normalize as fullnames=/datadir/treenames
   datafile  <- object@rootfile;
   setname   <- object@setname;
   treenames <- as.character(object@treenames);
   numtrees  <- length(treenames);
   fullnames <- paste(setname, treenames, sep="/");

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check for presence of temporary directory
   tmpdir <- validTempDir(tmpdir);

   ## update: need to set filename to rootfile
   if (!is.logical(update)) {
      stop(paste(sQuote("update"), "must be TRUE or FALSE"));
   } else if (update == TRUE) {
      filename <- rootfile;
   } else if (update == FALSE && existsROOTFile(rootfile)) {
      ## check if root file exists (necessary for WinXP to test already here)
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## check for presence of valid selector option
   TYPE <- c("separate", "together");
   if (is.na(match(select, TYPE))) {
      stop(paste(sQuote(select), "is not a valid selector option"));
   }#if

   ## check for valid normalization method
   TYPE <- c("mean", "median", "lowess", "supsmu");
   if (is.na(match(method, TYPE))) {
      stop(paste(sQuote("method"), "must be <mean, median,lowess,supsmu>"));
   }#if

   ## get transcript part and option part
   transcript <- unlist(strsplit(option, ":"))[1];
   transcript <- validTranscriptOption(transcript);
   option     <- unlist(strsplit(option, ":"))[2];

   ## percent units to be selected
   if (option == "all") {
      pc <- 1.0;
   } else if (option == "sel") {
      pc <- 0.3; #default
   } else {
      pc <- as.integer(as.character(option));
   }#if
   if (!(pc > 0 && pc <= 1)) {
      stop(paste(sQuote("option"), "must be <all,sel> or number in range (0,1]"));
   }#if

   ## check parameter list
   if (length(params) == 0) {
      stop(paste("empty parameter list", sQuote("params")));
   }#if

   ## option, tree extension and parameters for normalization method
   logbase <- validLogbase(logbase);
   switch(method,
      mean = {
         if (!(length(params) == 2 &&
               params[[1]] >= 0 && params[[1]] <= 0.5 &&
               params[[2]] >= 0)) {
            warning(paste(sQuote("params"),
                    "are not correct, using default parameters <0.02,0>"));
            params[[1]] <- 0.02;
            params[[2]] <- 0;
         }#if
         params <- list(trim=params[[1]], sc=params[[2]]);
         option <- paste(transcript, ifelse(pc<1, "sel", "all"), logbase, sep=":");
         exten  <- "tmn";
      }, median = {
         if (length(params) == 1) {
            if (params[[1]] < 0) {
               warning(paste(sQuote("params"),
                       "is not correct, using default setting sc=0"));
               params[[1]] <- 0;
            }#if
            params <- list(trim=0.5, sc=params[[1]]);
         } else if (!(length(params) == 2 &&
                    params[[1]] == 0.5 && params[[2]] >= 0)) {
            warning(paste(sQuote("params"),
                    "are not correct, using default parameters <0.5,0>"));
            params[[1]] <- 0.5;
            params[[2]] <- 0;
         }#if
         params <- list(trim=params[[1]], sc=params[[2]]);
         method <- "mean";
         option <- paste(transcript, ifelse(pc<1, "sel", "all"), logbase, sep=":");
         exten  <- "med";
      }, lowess = {
         if (!(length(params) == 4 &&
               params[[1]] >  0 && params[[1]] <= 1 &&
               params[[2]] >= 0 &&
               params[[3]] >= 0 && params[[3]] <= 2  &&
               params[[4]] >= 0 && params[[4]] <= 1)) {
            warning(paste(sQuote("params"),
                          "are not correct, using default parameters <0.67,3,0.,0.>"));
            params[[1]] <- 0.67; # span for lowess
            params[[2]] <- 3;    # iter for lowess
            params[[3]] <- 0.0;  # rule for approx
            params[[4]] <- 0.0;  # f    for approx
         }#if
         params <- list(span=params[[1]], iter=params[[2]], rule=params[[3]], f=params[[4]]);
         option <- paste(transcript, option, logbase, sep=":");
         exten  <- "low";
      }, supsmu = {
         if (!(length(params) == 4 &&
               params[[1]] >= 0 && params[[1]] <= 1  &&
               params[[2]] >= 0 && params[[2]] <= 10 && 
               params[[3]] >= 0 && params[[3]] <= 2  &&
               params[[4]] >= 0 && params[[4]] <= 1)) {
            warning(paste(sQuote("params"),
                          "are not correct, using default parameters <0.,0.,0.,0.>"));
            params[[1]] <- 0.0;  # span for supsmu
            params[[2]] <- 0.0;  # bass for supsmu
            params[[3]] <- 0.0;  # rule for approx
            params[[4]] <- 0.0;  # f    for approx
         }#if
         params <- list(span=params[[1]], bass=params[[2]], rule=params[[3]], f=params[[4]]);
         option <- paste(transcript, option, logbase, sep=":");
         exten  <- "sup";
      }
   )#switch

   ## exon level
   if (exonlevel != "") {
      cat("note:", sQuote("exonlevel"), "has no effect on ExprTreeSet\n");
   }#if
   exonlevel <- 0;

   ## reference tree
   reftree <- "";
   if (!(refindex >= 0 && refindex <= numtrees)) {
      stop(paste(sQuote("refindex"), "must be 0 or less than number of trees"));
   }#if
   if (refindex == 0) {
      xten    <- unlist(strsplit(treenames[[1]], "\\."));
      reftree <- paste(setname, "/*.", xten[2], sep="");
   } else {
      reftree <- paste(setname, treenames[[refindex]], sep="/");
   }#if

   ## reference method
   TYPE <- c("mean", "median");
   if (is.na(match(refmethod, TYPE))) {
      stop(paste(sQuote("refmethod"), "must be <mean, median>"));
   }#if

   ## define setname and settype for new treeset
   setname <- "NormSet";
   settype <- "normation";

   ## normalize
   r <- .C("Normxpress",
           as.character(filename),
           as.character(filedir),
           as.character(chiptype),
           as.character(schemefile),
           as.character(tmpdir),
           as.character(select),
           as.double(pc),
           as.character(method),
           as.character(option),
           as.integer(length(params)),
           as.double(params),
           as.integer(exonlevel),
           as.character(setname),
           as.character(datafile),
           as.character(fullnames),
           as.integer(numtrees),
           as.character(reftree),
           as.character(refmethod),
           as.integer(update),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("Normxpress")));
      return(NULL);
   }#if

   ## export/import result as dataframe
   ds <- data.frame(matrix(nrow=0, ncol=0));
   if (add.data) {
      outfile  <- sub("\\.root", ".txt", rootfile);
      # get treename "treeset.treename.treetype"
      treename <- paste(setname, "*", exten, sep=".");
      numtrees <- 1; # must be one for treename="*"

      r <- .C("ExportData",
              as.character(rootfile),
              as.character(schemefile),
              as.character(chiptype),
              as.character(settype),
              as.character(treename),
              as.integer(numtrees),
              as.character(exten),
              as.character("fUnitName:fLevel"),
              as.character(outfile),
              as.character("\t"),
              as.integer(verbose),
              err=integer(1),
              PACKAGE="xps")$err;

      if (r != 0) {
         stop(paste("error in function", sQuote("ExportData")));
         return(NULL);
      }#if

      if (file.exists(outfile)) {
         ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
      } else {
         warning(paste("could not export results as", sQuote(outfile)));
      }#if
   }#if

   ## get treenames after normalization
   treenames <- getTreeNames(rootfile, exten);
   numtrees  <- length(treenames);

   ## get exprtype of object
   exprtype = "none";
   if (class(object) == "ExprTreeSet") {
      exprtype = object@exprtype;
   }#if

   ## create new class ExprTreeSet
   set <- new("ExprTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = scheme,
              data      = ds,
              params    = params,
              exprtype  = exprtype,
              normtype  = method);

   return(set);
}#normalize.ExprTreeSet

setMethod("xpsNormalize", "ExprTreeSet", normalize.ExprTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"prefilter.ExprTreeSet" <-
function(object,
         filename   = character(0),
         filedir    = getwd(),
         filter     = NULL,
         minfilters = 999,
         logbase    = "log2",
         treename   = "PreFilter",
         xps.call   = NULL,
         verbose    = TRUE)
{
   if (debug.xps()) print("------prefilter.ExprTreeSet------")

   ## get schemefile and chiptype for object
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chiptype   <- chipType(object);
   chipname   <- chipName(scheme);

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## check if root file exists (necessary for WinXP to test already here)
   if (existsROOTFile(rootfile)) {
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## check for prefilter
   if (is.null(filter) || !is(filter, "PreFilter")) {
      stop(paste("missing filter of class", sQuote("PreFilter")));
   }#if

   ## check for number of filters
   if (filter@numfilters <= 0) {
      stop(paste(sQuote("PreFilter"), "must have at least one initialized filter"));
   }#if

   ## check for minfilters
   if (minfilters <= 0) {
      stop(paste(sQuote("minfilters"), "must be at least one"));
   }#if

   ## check for valid logbase
   logbase <- validLogbase(logbase);

   ## get exprtree names to add to filter set
   datafile  <- object@rootfile;
   setname   <- object@setname;
   treenames <- as.character(object@treenames);
   numexprs  <- length(treenames);
   exprnames <- paste(datafile, setname, treenames, sep="/");

   ## check for presence of optional root detection-call file
   numcalls  <- 0;
   callnames <- "";
   if (!is.null(xps.call)) {
      if (!(is(xps.call, "CallTreeSet") &&
            (schemeFile(xps.call) == schemefile) &&
            (chipType(xps.call)   == chiptype) &&
            (file.exists(rootFile(xps.call))))) {
         stop(paste("wrong detection call file", sQuote("xps.call")));
      }#if

      ## get calltree names to add to filter set
      datafile  <- rootFile(xps.call);
      setname   <- setName(xps.call);
      treenames <- treeNames(xps.call);
      numcalls  <- length(treenames);
      callnames <- paste(datafile, setname, treenames, sep="/");
   }#if

   ## prefilter data
   nvarparams  <- 0; varparams  <- 0; varoption  = "";
   nlowparams  <- 0; lowparams  <- 0; lowoption  = "";
   nhiparams   <- 0; hiparams   <- 0; hioption   = "";
   nquanparams <- 0; quanparams <- 0; quanoption = "";
   ncallparams <- 0; callparams <- 0; calloption = "";
   trim    <- NULL;
   epsilon <- NULL;

   ## variance filter
   madev      <- madFilter(filter);
   cv         <- cvFilter(filter);
   variance   <- varFilter(filter);
   difference <- diffFilter(filter);
   ratio      <- ratioFilter(filter);
   gap        <- gapFilter(filter);
   if (length(madev) == 2) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- madev[1];
      epsilon    <- madev[2];
      varoption  <- paste(varoption, "mad", sep=":", collapse=":");
   }#if
   if (length(cv) == 3) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- cv[1];
      trim       <- cv[2];
      epsilon    <- cv[3];
      varoption  <- paste(varoption, "cov2mn", sep=":", collapse=":");
   }#if
   if (length(variance) == 3) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- variance[1];
      trim       <- variance[2];
      epsilon    <- variance[3];
      varoption  <- paste(varoption, "var2mn", sep=":", collapse=":");
   }#if
   if (length(difference) == 3) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- difference[1];
      trim       <- difference[2];
      epsilon    <- difference[3];
      varoption  <- paste(varoption, "dif2mn", sep=":", collapse=":");
   }#if
   if (length(ratio) == 1) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- ratio[1];
      varoption  <- paste(varoption, "max2min", sep=":", collapse=":");
   }#if
   if (length(gap) == 4) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- gap[1];
      nvarparams <- nvarparams + 1;
      varparams[nvarparams]  <- gap[2];
      trim       <- gap[3];
      epsilon    <- gap[4];
      varoption  <- paste(varoption, "gap2mn", sep=":", collapse=":");
   }#if
   if (!is.null(trim)) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams] <- trim;
   }#if
   if (!is.null(epsilon)) {
      nvarparams <- nvarparams + 1;
      varparams[nvarparams] <- epsilon;
   }#if
   varoption <- substring(varoption, 2);

   ## lowFilter
   lothreshold <- lowFilter(filter);
   if (length(lothreshold) == 3) {
      nlowparams <- 2;
      lowparams  <- lothreshold[1:2];
      lowoption  <- lothreshold[3];
   }#if

   ## highFilter
   hithreshold <- highFilter(filter);
   if (length(hithreshold) == 3) {
      nhiparams <- 2;
      hiparams  <- hithreshold[1:2];
      hioption  <- hithreshold[3];
   }#if

   ## quantileFilter
   quantil <- quantileFilter(filter);
   if (length(quantil) == 3) {
      nquanparams <- 3;
      quanparams  <- quantil;
   }#if

   ## callFilter
   prescall <- callFilter(filter);
   if (length(prescall) == 3) {
      ncallparams <- 2;
      callparams  <- prescall[1:2];
      calloption  <- prescall[3];
   }#if

   ## define setname and settype for new treeset
   setname <- "PreFilterSet";
   settype <- "prefilter";
   exten   <- "pfr";

   ## prefilter
   r <- .C("PreFilter",
           as.character(filename),
           as.character(filedir),
           as.character(chiptype),
           as.character(chipname),
           as.character(schemefile),
           as.character(setname),
           as.character(treename),
           as.integer(minfilters),
           as.character(logbase),
           as.integer(nvarparams),
           as.double(varparams),
           as.character(varoption),
           as.integer(nlowparams),
           as.double(lowparams),
           as.character(lowoption),
           as.integer(nhiparams),
           as.double(hiparams),
           as.character(hioption),
           as.integer(nquanparams),
           as.double(quanparams),
           as.character(quanoption),
           as.integer(ncallparams),
           as.double(callparams),
           as.character(calloption),
           as.character(exprnames),
           as.integer(numexprs),
           as.character(callnames),
           as.integer(numcalls),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("PreFilter")));
      return(NULL);
   }#if

   ## export/import filter mask as dataframe
   outfile  <- sub("\\.root", ".txt", rootfile);
   # get treename "treeset.treename.treetype"
   treename <- paste(setname, treename, exten, sep=".");
   numtrees <- 1; 

   r <- .C("ExportData",
           as.character(rootfile),
           as.character(schemefile),
           as.character(chiptype),
           as.character(settype),
           as.character(treename),
           as.integer(numtrees),
           as.character(exten),
           as.character("fUnitName:fFlag"),
           as.character(outfile),
           as.character("\t"),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ExportData")));
      return(NULL);
   }#if

   ## get treenames after normalization
   treenames <- getTreeNames(rootfile, exten);
   numtrees  <- length(treenames);

   ds <- data.frame(matrix(nrow=0, ncol=0));
   if (file.exists(outfile)) {
      ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
   } else {
      warning(paste("could not export results as", sQuote(outfile)));
   }#if

   ## create new class FilterTreeSet
   set <- new("FilterTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = scheme,
              data      = ds,
#?              params    = params,
              filter    = filter,
              exprset   = object);

   ## add callset
   if (!is.null(xps.call) && is(xps.call, "CallTreeSet")) {
      set@callset <- xps.call;
   }#if

   return(set);
}#prefilter.ExprTreeSet

setMethod("xpsPreFilter", "ExprTreeSet", prefilter.ExprTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"unifilter.ExprTreeSet" <-
function(object,
         filename   = character(0),
         filedir    = getwd(),
         update     = FALSE,
         filter     = NULL,
         minfilters = 999,
         logbase    = "log2",
         group      = character(0),
         treename   = "UniTest",
         xps.fltr   = NULL,
         xps.call   = NULL,
         verbose    = TRUE)
{
   if (debug.xps()) print("------unifilter.ExprTreeSet------")

   ## get schemefile and chiptype for object
   scheme     <- object@scheme;
   schemefile <- schemeFile(object);
   chiptype   <- chipType(object);
   chipname   <- chipName(scheme);

   ## root file /filedir/filename.root
   rootfile <- rootDirFile(filename, filedir);

   ## update: need to set filename to rootfile
   if (!is.logical(update)) {
      stop(paste(sQuote("update"), "must be TRUE or FALSE"));
   } else if (update == TRUE) {
      filename <- rootfile;
   } else if (update == FALSE && existsROOTFile(rootfile)) {
      ## check if root file exists (necessary for WinXP to test already here)
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## check for unifilter
   if (is.null(filter) || !is(filter, "UniFilter")) {
      stop(paste("missing filter of class", sQuote("UniFilter")));
   }#if

   ## check for minfilters
   if (minfilters <= 0) {
      stop(paste(sQuote("minfilters"), "must be at least one"));
   }#if

   ## check for valid logbase
   logbase <- validLogbase(logbase);

   ## check for presence of optional root prefilter file
   fltrname <- "";
   numfltr  <- 0;
   if (!is.null(xps.fltr)) {
      if (setType(xps.fltr) != "prefilter") {
         stop(paste(sQuote("xps.fltr"), "must be of type", sQuote("PreFilter")));
      }#if
      if (chipName(xps.fltr@scheme) != chipname) {
         stop(paste(sQuote("xps.fltr"), "was not created for chip", sQuote(chipname)));
      }#if

      # get treename "rootfile.treeset.treename"
      rootname <- rootFile(xps.fltr);
      setname  <- setName(xps.fltr);
      fltrname <- treeNames(xps.fltr);
      fltrname <- paste(rootname, setname, fltrname, sep="/");
      numfltr  <- 1;
   }#if

   ## get exprtree names to add to filter set
   datafile  <- object@rootfile;
   setname   <- object@setname;
   treenames <- as.character(object@treenames);
   numexprs  <- length(treenames);
   exprnames <- paste(datafile, setname, treenames, sep="/");

   ## check for presence of optional root detection call file
   numcalls  <- 0;
   callnames <- "";
   if (!is.null(xps.call)) {
      if (!(is(xps.call, "CallTreeSet") &&
            (schemeFile(xps.call) == schemefile) &&
            (chipType(xps.call)   == chiptype) &&
            (file.exists(rootFile(xps.call))))) {
         stop(paste("wrong detection call file", sQuote("xps.call")));
      }#if

      ## get calltree names to add to filter set
      datafile  <- rootFile(xps.call);
      setname   <- setName(xps.call);
      treenames <- treeNames(xps.call);
      numcalls  <- length(treenames);
      if (numcalls != numexprs) {
         stop("number of call trees must be equal to number of expression trees");
      }#if
      callnames <- paste(datafile, setname, treenames, sep="/");
   }#if

   ## unifilter data
   nuniparams  <- 0; uniparams  <- 0; unioption  = "";
   nfcparams   <- 0; fcparams   <- 0; fcoption   = "";
   nufparams   <- 0; ufparams   <- 0; ufoption   = "";
   ncallparams <- 0; callparams <- 0; calloption = "";

   ## uniTest
   unitest <- uniTest(filter);
   if (length(unitest) == 8) {
      nuniparams <- 5;
      uniparams <- as.double(unitest[4:8]);
      unioption <- sub("\\.","",paste(unitest[2], unitest[3], sep=":"));
      unitest   <- sub("\\.","",unitest[1]);
   } else {
         stop(paste("missing or incorrect parameter", sQuote("UniFilter@unitest")));
   }#if

   ## fcFilter
   foldchange <- fcFilter(filter);
   if (length(foldchange) == 2) {
      nfcparams <- 2;
      fcparams[1] <- foldchange[1];
      fcparams[2] <- ifelse(foldchange[2]=="up",1, ifelse(foldchange[2]=="down",-1,0));
   }#if

   ## unitestFilter
   unifltr <- unitestFilter(filter);
   if (length(unifltr) > 0 & length(foldchange) == 2) {
      nufparams <- 1;
      ufparams <- unifltr[1];
      ufoption <- unifltr[2];
   }#if

   ## callFilter
   prescall <- callFilter(filter);
   if (length(prescall) == 3) {
      ncallparams <- 3;
      callparams[1] <- prescall[1];
      callparams[2] <- prescall[2];
      callparams[3] <- prescall[2];
      calloption    <- paste(prescall[3], prescall[3], sep=":");
   }#if

   ## check group and get group index
   if (length(unique(group)) != 2) {
         stop("number of groups must be two");
   }#if
   grpidx <- as.integer(as.factor(group));
   
   ## define setname and settype for unifilter
   setname <- "UniFilterSet";
   settype <- "unifilter";
   exten   <- "ufr";

   ## unifilter
   r <- .C("UniFilter",
           as.character(filename),
           as.character(filedir),
           as.character(chiptype),
           as.character(chipname),
           as.character(schemefile),
           as.character(setname),
           as.character(treename),
           as.integer(minfilters),
           as.character(logbase),
           as.character(unitest),
           as.integer(nuniparams),
           as.double(uniparams),
           as.character(unioption),
           as.integer(nfcparams),
           as.double(fcparams),
           as.character(fcoption),
           as.integer(nufparams),
           as.double(ufparams),
           as.character(ufoption),
           as.integer(ncallparams),
           as.double(callparams),
           as.character(calloption),
           as.character(exprnames),
           as.integer(numexprs),
           as.character(callnames),
           as.integer(numcalls),
           as.character(group),
           as.integer(grpidx),
           as.character(fltrname),
           as.integer(numfltr),
           as.character("*"),
           as.integer(update),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("UniFilter")));
      return(NULL);
   }#if

   ## export optional unifilter
   ds <- data.frame(matrix(nrow=0, ncol=0));
   if (numberFilters(filter) > 0) {
      ## export/import unifilter as dataframe
      outfile  <- sub("\\.root", ".txt", rootfile);
      # get treename "treeset.treename.treetype"
      ufltrname <- paste(setname, treename, exten, sep=".");
      numtrees  <- 1;

      r <- .C("ExportData",
              as.character(rootfile),
              as.character(schemefile),
              as.character(chiptype),
              as.character(settype),
              as.character(ufltrname),
              as.integer(numtrees),
              as.character(exten),
              as.character("fUnitName:fFlag"),
              as.character(outfile),
              as.character("\t"),
              as.integer(verbose),
              err=integer(1),
              PACKAGE="xps")$err;

      if (r != 0) {
         stop(paste("error in function", sQuote("ExportData")));
         return(NULL);
      }#if

      if (file.exists(outfile)) {
         ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
      } else {
         warning(paste("could not export results as", sQuote(outfile)));
      }#if
   }#if

   ## get treenames after normalization
   treenames <- getTreeNames(rootfile, exten);
   numtrees  <- length(treenames);

   ## create new class FilterTreeSet
   fltrset <- new("FilterTreeSet",
                  setname   = setname,
                  settype   = settype,
                  rootfile  = rootfile,
                  filedir   = filedir,
                  numtrees  = numtrees,
                  treenames = as.list(treenames),
                  scheme    = scheme,
                  data      = ds,
#?                  params    = params,
                  filter    = filter,
                  exprset   = object);

   ## add callset to class FilterTreeSet
   if (!is.null(xps.call) && is(xps.call, "CallTreeSet")) {
      fltrset@callset <- xps.call;
   }#if

   ## varlist for export
#ev depending on uniparams, unioption
   varlist <- "fUnitName:df:mn1:mn2:fc:pval:padj:pcha:se:stat";

   ## define setname and settype for new treeset
   setname <- "UniFilterSet";
   settype <- "UnivariateAnalysis";
#??   settype <- "unifilter";
   exten   <- "stt";

   ## export/import unitest result as dataframe
   outfile  <- sub("\\.root", ".txt", rootfile);
   # get treename "treeset.treename.treetype"
   treename <- paste(setname, treename, exten, sep=".");
   numtrees <- 1;

   r <- .C("ExportData",
           as.character(rootfile),
           as.character(schemefile),
           as.character(chiptype),
           as.character(settype),
           as.character(treename),
           as.integer(numtrees),
           as.character(exten),
           as.character(varlist),
           as.character(outfile),
           as.character("\t"),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ExportData")));
      return(NULL);
   }#if

   ## get treenames
   treenames <- getTreeNames(rootfile, exten);
   numtrees  <- length(treenames);

   ds <- data.frame(matrix(nrow=0, ncol=0));
   if (file.exists(outfile)) {
      ds <- read.table(outfile, header=TRUE, check.names=FALSE, sep="\t", row.names=NULL);
   } else {
      warning(paste("could not export results as", sQuote(outfile)));
   }#if

   ## create new class AnalysisTreeSet
   set <- new("AnalysisTreeSet",
              setname   = setname,
              settype   = settype,
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = numtrees,
              treenames = as.list(treenames),
              scheme    = scheme,
              data      = ds,
#?              params    = params,
              fltrset   = fltrset);

   return(set);
}#unifilter.ExprTreeSet

setMethod("xpsUniFilter", "ExprTreeSet", unifilter.ExprTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("corplot", signature(x="ExprTreeSet"),
   function(x,
            which      = "UnitName",
            transfo    = log2,
            method     = "spearman",
            col        = NULL,
            names      = "namepart",
            sort       = FALSE,
            reverse    = TRUE,
            bmar       = NULL,
            add.legend = FALSE,
            ...) 
   {
      if (debug.xps()) print("------corplot.ExprTreeSet------")

      m <- validData(x, which=which);
      if (is.function(transfo)) m <- transfo(m);

      if (is.null(names))              names <- colnames(m)
      else if (names[1] == "namepart") names <- namePart(colnames(m))
      else                             m     <- m[, names, drop=F];

      m <- cor(m, method = method);
      if (reverse) 1 - m;
      if (sort)    m[order(m[,1], decreasing = TRUE),];

      if (is.null(col)) col <- heat.colors(12);

      if (is.null(bmar)) bmar <- adjustXLabMargins(names, bottom=6, cex=1.0);

      oldpar <- par(no.readonly=TRUE, mar=c(bmar$b, bmar$b, 3, 1));
      if (add.legend) {
         layout(matrix(c(1, 2), 1, 2, byrow=TRUE), widths=c(7,1), heights=8, TRUE);
         par(mar = c(bmar$b, bmar$b, 3, 0));
      }#if

      ## plot image
      graphics::image(x    = 1:ncol(m),
                      y    = 1:ncol(m),
                      z    = m,
                      xlab = "",
                      ylab = "",
                      zlim = range(m, na.rm=TRUE),
                      main = "Array-Array Expression Level Correlation",
                      col  = col,
                      axes = FALSE);
      box();
      axis(1, at=1:ncol(m), labels=names, las=3);
      axis(2, at=1:ncol(m), labels=names, las=2);

      if (add.legend) {
         par(mar = c(bmar$b, 1, 3, 3));
         y <- pretty(range(m, na.rm=TRUE), 10);
         m <- matrix(y, nrow=1, ncol=length(y));
         graphics::image(m, xaxt="n", yaxt="n", col=col);
         axis(4, labels=y, at=seq(0, 1, by=(1/(length(y)-1))), las=2, cex.axis=0.8);
         layout(1);
         par(mar = c(5, 4, 4, 2) + 0.1);
      }#if

      par(oldpar);
   }
)#corplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("madplot", signature(x="ExprTreeSet"),
   function(x,
            which      = "UnitName",
            transfo    = log2,
            col        = NULL,
            names      = "namepart",
            sort       = FALSE,
            bmar       = NULL,
            add.legend = FALSE,
            ...) 
   {
      if (debug.xps()) print("------madplot.ExprTreeSet------")

      m <- validData(x, which=which);
      if (is.function(transfo)) m <- transfo(m);

      if (is.null(names))              names <- colnames(m)
      else if (names[1] == "namepart") names <- namePart(colnames(m))
      else                             m     <- m[, names, drop=F];

      m <- distMAD(m);
      if (sort) m[order(m[,1], decreasing = TRUE),];

      ## color adapted from brewer.pal(9, "RdPu")
      if (is.null(col)) col <- colorRampPalette(rgb(c(255,253,252,250,247,221,174,122,73),
                                                    c(247,224,197,159,104,52,1,1,0),
                                                    c(243,221,192,181,161,151,126,119,106),
                                                    maxColorValue=255)
                                               )(256);

      if (is.null(bmar)) bmar <- adjustXLabMargins(names, bottom=6, cex=1.0);

      oldpar <- par(no.readonly=TRUE, mar=c(bmar$b, bmar$b, 3, 1));
      if (add.legend) {
         layout(matrix(c(1, 2), 1, 2, byrow=TRUE), widths=c(7,1), heights=8, TRUE);
         par(mar = c(bmar$b, bmar$b, 3, 0));
      }#if

      ## plot image
      graphics::image(x    = 1:ncol(m),
                      y    = 1:ncol(m),
                      z    = m,
                      xlab = "",
                      ylab = "",
                      zlim = range(m, na.rm=TRUE),
                      main = "",
                      col  = col,
                      axes = FALSE);
      box();
      axis(1, at=1:ncol(m), labels=names, las=3);
      axis(2, at=1:ncol(m), labels=names, las=2);

      if (add.legend) {
         par(mar = c(bmar$b, 1, 3, 3));
         y <- pretty(range(m, na.rm=TRUE), 10);
         m <- matrix(y, nrow=1, ncol=length(y));
         graphics::image(m, xaxt="n", yaxt="n", col=col);
         axis(4, labels=y, at=seq(0, 1, by=(1/(length(y)-1))), las=2, cex.axis=0.8);
         layout(1);
         par(mar = c(5, 4, 4, 2) + 0.1);
      }#if

      par(oldpar);
   }
)#madplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("mvaplot", signature(x="ExprTreeSet"),
   function(x,
            which   = "UnitName",
            transfo = log2,
            method  = "median",
            names   = "namepart",
            ylim    = c(-6,6),
            xlab    = "A",
            ylab    = "M",
            pch     = '.',
            las     = 2,
            ...) 
   {
      if (debug.xps()) print("------mvaplot.ExprTreeSet------")

      m <- validData(x, which=which);
      if (is.function(transfo)) m <- transfo(m);

      if (method == "median") mn  <- apply(m, 1, median)
      else                    mn  <- rowMeans(m);

      if (is.null(names))              names <- colnames(m)
      else if (names[1] == "namepart") names <- namePart(colnames(m))
      else                             m     <- m[, names, drop=F];

      ## plot images
      if (ncol(m) > 1 && interactive()) par(ask=TRUE) else par(ask=FALSE);
      for (i in 1:ncol(m)) {
         plot(x    = mn,
              y    = (m[,i] - mn),
              main = names[i],
              ylim = ylim,
              xlab = xlab,
              ylab = ylab,
              pch  = pch,
              las  = las,
              ...);
         abline(h=0);
      }#for
      par(ask=FALSE);
   }
)#mvaplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("nuseplot", signature(x="ExprTreeSet"),
   function(x,
            which    = "UnitName",
            size     = 0,
            range    = 0,
            names    = "namepart",
            main     = "NUSE Plot",
            ylim     = c(0.8, 1.2),
            las      = 2,
            add.line = TRUE,
            outline  = FALSE,
            ...) 
   {
      if (debug.xps()) print("------nuseplot.ExprTreeSet------")

      m <- validSE(x, which=which);
      if (size > 1) m <- m[seq(1,nrow(m),len=size),];
      m <- log2(m);

      if (is.null(names))              names <- colnames(m)
      else if (names[1] == "namepart") names <- namePart(colnames(m))
      else                             m     <- m[, names, drop=F];

      med  <- apply(m, 1, median, na.rm=TRUE);
      nuse <- data.frame(sweep(m, MARGIN=1, STATS=med, FUN='/'));

      boxplot(nuse,
              range   = range,
              xlab    = '',
              xaxt    = 'n',
              main    = main,
              ylim    = ylim,
              outline = outline);
      if (add.line) {
         abline(1, 0, lty=2, col="gray70");
      }#if
      axis(side=1, outer=F, at=1:length(names), labels=names, las=las, ...);
   }
)#nuseplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("pcaplot", signature(x="ExprTreeSet"),
   function(x,
            which      = "UnitName",
            transfo    = log2,
            method     = "none",
            groups     = NULL,
            screeplot  = FALSE,
            squarepca  = FALSE,
            pcs        = c(1,2),
            add.labels = FALSE,
            add.legend = FALSE,
            col        = NULL,
            names      = "namepart",
            as.list    = FALSE,
            ...) 
   {
      if (debug.xps()) print("------pcaplot.ExprTreeSet------")

      if (length(pcs) != 2) stop("only two principal components can be plotted.\n");

      m <- validData(x, which=which);
      if (is.function(transfo)) m <- transfo(m);

      if (is.null(names))              names <- colnames(m)
      else if (names[1] == "namepart") names <- namePart(colnames(m))
      else                             m     <- m[, names, drop=F];

      ## PCA
      if (method == "none") pca <- prcomp(t(m))
      else                  pca <- prcomp((1 - cor(m, method = method)));

      ## plot PCA
      if (screeplot) {
         screeplot(pca,
                   xaxt = "n",
                   main = "Screeplot for PCA",
                   ...);
         axis(1, at=1:ncol(m), labels=colnames(pca[[2]]), las=3);
      }else{
         if (squarepca) {
            ylim <- max(abs(range(pca$x[,pcs[1]])));
            ylim <- c(-ylim, ylim);
         } else {
            ylim <- NULL;
         }#if


         if (!is.null(groups)) {
            if (length(groups) != length(names)) {
               stop("length of groups must be equal to number of samples.\n");
            }#if

            groupnames <- as.character(groups);
            groups     <- unique(groups);
            plotsymbol <- 0:25;
            plotsymbol <- split(plotsymbol[1:length(groups)], groups);
            plotsymbol <- unlist(sapply(as.factor(groupnames),function(x)plotsymbol[[x]]));
         }else{
            groupnames <- names;
            plotsymbol <- 1;
         }#if

         if (is.null(col)) {
            col <- c("blue3", "blue2", "blue1", "steelblue3", "steelblue2", "steelblue1",
                     "lightblue3", "lightblue2", "lightblue1", "gray60",
                     "red3", "red2", "red1", "orange3", "orange2", "orange1",
                     "yellow3", "yellow2", "yellow1", "black");
            col <-rep(col, (floor(length(names)/length(col)) + 1));
         }#if

         plot(pca$x[,pcs],
              ylim = ylim,
              col  = col,
              pch  = plotsymbol,
              main = "Principal Components Plot",
              ...);

         if (add.labels) {
            offset <-  (par("usr")[4] - par("usr")[3])/50;
            text(pca$x[,pcs[1]], pca$x[,pcs[2]] + offset, label=names, cex=0.7);
         }#if

         pos.legend <- "topleft";
         if (is.character(add.legend)) {
            pos.legend <- add.legend;
            add.legend <- TRUE;
         }#if
         if (!is.null(groups) & add.legend) {
            legend(x      = pos.legend,
                   legend = unique(groupnames),
                   pch    = unique(plotsymbol),
                   pt.bg  = "white",
                   cex    = 0.6
                  );
         }#if
      }#if

      if (as.list) return(pca);
   }
)#pcaplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("rleplot", signature(x="ExprTreeSet"),
   function(x,
            which    = "UnitName",
            size     = 0,
            range    = 0,
            names    = "namepart",
            main     = "RLE Plot",
            ylim     = c(-1.0, 1.0),
            las      = 2,
            add.line = TRUE,
            outline  = FALSE,
            ...) 
   {
      if (debug.xps()) print("------rleplot.ExprTreeSet------")

      mboxplot(x,
               which    = which,
               size     = size,
               transfo  = log2,
               method   = "median",
               range    = range,
               names    = names,
               main     = main,
               ylim     = ylim,
               las      = las,
               add.line = add.line,
               outline  = outline,
               ...);
   }
)#rleplot

#------------------------------------------------------------------------------#
