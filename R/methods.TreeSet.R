#==============================================================================#
# methods.TreeSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# rootFile:
# rootFile<-:
# fileDir:
# fileDir<-:
# setName:
# setName<-:
# setType:
# setType<-:
# treeNames:
# treeInfo:
# export:
# root.browser: 
#==============================================================================#


#------------------------------------------------------------------------------#
# TreeSet initialization:
#------------------------------------------------------------------------------#

setMethod("initialize", "TreeSet", 
   function(.Object,
            setname   = character(),
            settype   = character(),
            rootfile  = NULL,
            filedir   = getwd(),
            numtrees  = 0,
            treenames = list(),
            ...) 
   {
      if (debug.xps()) print("------initialize:TreeSet------")

      ## prevent effects of multiple initialization
      if (is.null(rootfile)) return(.Object);
      .Object@setname  <- setname;
      .Object@settype  <- settype;
      .Object@rootfile <- rootfile;
      .Object@filedir  <- filedir;

      rootfile <- rootDirFile(rootfile, filedir);

      .Object <- callNextMethod(.Object,
                                setname   = setname,
                                settype   = settype,
                                rootfile  = rootfile,
                                filedir   = filedir,
                                numtrees  = numtrees,
                                treenames = treenames,
                                ...);
      .Object@setname   <- setname;
      .Object@settype   <- settype;
      .Object@rootfile  <- rootfile;
      .Object@filedir   <- filedir;
      .Object@numtrees  <- numtrees;
      .Object@treenames <- treenames;
      .Object;
   }
)#initialize

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("TreeSet",
   function(object) {
      if (debug.xps()) print("------setValidity:TreeSet------")

      msg <- NULL;

      ## check rootfile
      if (!isROOTFile(object@rootfile)) {
         msg <- validMsg(msg, paste(sQuote("rootfile"),
                                   "is missing or is not a ROOT file"));
      }#if

      ## check for presence of directory for file
      if (!file.exists(object@filedir)) {
         msg <- validMsg(msg, paste(sQuote("filedir"),
                                   "is not a system directory"));
      }#if

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# TreeSet accessors:
#------------------------------------------------------------------------------#

setMethod("rootFile", signature(object="TreeSet"),
   function(object) object@rootfile
)#rootFile

setReplaceMethod("rootFile", signature(object="TreeSet", value="character"),
   function(object, value) {
      object@rootfile <- value;
      return(object);
   }
)#rootFile<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("fileDir", signature(object="TreeSet"),
   function(object) object@filedir
)#fileDir

setReplaceMethod("fileDir", signature(object="TreeSet", value="character"),
   function(object, value) {
      object@filedir <- value;
      return(object);
   }
)#fileDir<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("setName", signature(object="TreeSet"),
   function(object) object@setname
)#setName

setReplaceMethod("setName", signature(object="TreeSet", value="character"),
   function(object, value) {
      object@setname <- value;
      return(object);
   }
)#setName<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("setType", signature(object="TreeSet"),
   function(object) object@settype
)#setType

setReplaceMethod("setType", signature(object="TreeSet", value="character"),
   function(object, value) {
      object@settype <- value;
      return(object);
   }
)#setType<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("treeNames", signature(object="TreeSet"),
   function(object) object@treenames
)#treeNames


#------------------------------------------------------------------------------#
# TreeSet methods:
#------------------------------------------------------------------------------#

"exportUserInfo" <-
function(object,
         treename = "*",
         treetype = character(0),
         varlist  = "*",
         qualopt  = NULL,
         verbose  = FALSE,
         ...)
{
   if (debug.xps()) print("------exportUserInfo.TreeSet------")

   ds <- export(object,
                treename     = treename,
                treetype     = treetype,
                varlist      = paste("userinfo", varlist, sep=":"), 
                as.dataframe = TRUE,
                verbose      = verbose);

   if (nrow(ds) == 0) {
      stop(paste(sQuote("varlist"), "is not available in userinfo of current tree(s)"));
   }#if

   rownames(ds) <- ds[,"Parameter"]; ds <- ds[,-1, drop=FALSE];

   if (!is.null(qualopt)) {
      qualopt <-validQualityOption(qualopt);
      ds <- ds[, grep(qualopt, colnames(ds)), drop=FALSE];
   }#if

   return(ds);
}#exportUserInfo

setMethod("treeInfo", "TreeSet", exportUserInfo);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"exportTreeSet" <-
function(object, ...) 
{
   if (debug.xps()) print("------export.TreeSet------")

   ## check for presence of root file
   rootfile <- validROOTFile(rootFile(object), class(object));

   ## check for number of trees
   if (nchar(object@numtrees) < 1) {
      stop("number of trees is null");
   }#if
}#exportTreeSet

setMethod("export", "TreeSet", exportTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"rootBrowser" <-
function(object) 
{
   if (debug.xps()) print("------rootBrowser.TreeSet------")

   savedir <- getwd();
   setwd(object@filedir);

   ## unix vs windows settings
   is.win <- (.Platform$OS.type == "windows");
   if (is.win) {
      sq  <- as.character("");  #
      dq  <- as.character("\\\""); #\"
   } else {
      sq  <- as.character("'");  #'
      dq  <- as.character("\""); #"
   }#if

   macro <- paste(system.file(package="xps"), "rootsrc/macroOpenBrowser.C", sep="/");
   xpsso <- system.file("libs", .Platform$r_arch,
                  paste("xps", .Platform$dynlib.ext, sep=''), package="xps");
   xpsso <- gsub("//", "/", xpsso);
   xpsso <- paste(dq, xpsso, dq, sep="");
   macro <- paste(sq, macro, "(", xpsso, ")", sq, sep="");
   system(paste("root -l", macro));

   setwd(savedir);
}#rootBrowser

setMethod("root.browser", "TreeSet", rootBrowser);


#------------------------------------------------------------------------------#
