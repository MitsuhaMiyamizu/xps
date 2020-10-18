#==============================================================================#
# utils.R: contains utilty functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# xpsOptions:
# isChipType:
# isFilterCondition:
# isROOTFile:
# existsROOTFile:
# rootDirFile:
# rootDrawName:
# validEstimatorOption: 
# validLogbase:
# validMsg:
# validOption:
# validOutfile:
# validQualityOption: 
# validROOTFile:
# validSchemeTreeSet:
# validSeparator:
# validTempDir: 
# validTranscriptOption:
# validTreenames:
# validTreetype:
# getDatatype:
# getDataXY:
# CELNames:
# CELHeader:
# exonLevel:
# exonLevelIDs:
# extenPart:
# namePart:
# type2Exten:
# listTreeNames:
# circle:
# distMAD:
# plotDensities:
# pseudoPalette: 
# adjustXLabMargins:
# getNameType:
# getChipName:
# getChipType:
# getProbeInfo:
# getNumberTrees:
# getTreeNames:
# metaProbesets:
#==============================================================================#


#------------------------------------------------------------------------------#
# global constants: 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   eINITWEIGHT <- -16384;  # see XPSProcessing.h


#------------------------------------------------------------------------------#
# treetype: parameter used to get tree with corresponding tree extension
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ## scheme tree extensions are defined in src/XPSSchemes.cxx
   SCMTYPE <- c("scm", "idx", "prb", "ann");
   EXNTYPE <- c("cxy", "exn", "pbs", "anx", "anp");

   ## data tree extensions are defined in src/XPSProcessing.cxx
   RAWTYPE <- c("cel");
   ADJTYPE <- c("int");
   BGDTYPE <- c("sbg", "wbg", "rbg", "gbg");
   CNRTYPE <- c("cmn", "cmd", "clw", "css", "cqu");
   CALTYPE <- c("dc5", "dab", "ini");
   EXPTYPE <- c("amn", "gmn", "wmn", "wdf", "adf", "tbw", "mdp", "frm", "dfw", "fir");
   QUATYPE <- c("plm", "rlm", "brd", "res");
   PRETYPE <- c(ADJTYPE, BGDTYPE, CNRTYPE, CALTYPE, EXPTYPE, QUATYPE);
   NRMTYPE <- c("tmn", "med", "ksm", "low", "sup", "qua", "mdp");
   FLRTYPE <- c("pfr", "ufr", "mfr");
   UNITYPE <- c("uvt", "stt", "wil", "var");


#------------------------------------------------------------------------------#
# options: parameters for images
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   QUALOPT <- c("raw", "adjusted", "normalized");
   RESIOPT <- c("resids", "pos.resids", "neg.resids", "sign.resids", "weights");


################################################################################
# general utility functions
################################################################################

#------------------------------------------------------------------------------#
# xpsOptions: set options for xps
# currently only used to enable/disable debug messages
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
OPTIONS <- new.env(hash=TRUE, parent=emptyenv());
OPTIONS$debug <- FALSE;

xpsOptions <- function(debug=FALSE) {
   if (debug == FALSE) {
      OPTIONS$debug <- FALSE;
   } else {
      OPTIONS$debug <- TRUE;
   }#if
}#xpsOptions

debug.xps <- function(x) {OPTIONS$debug};

#------------------------------------------------------------------------------#
# isChipType: check if chiptype is "GeneChip", "GenomeChip" or "ExonChip"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
isChipType <- function(chiptype) {
   if (debug.xps()) print("------isChipType------")

   TYPE <- c("GeneChip", "GenomeChip", "ExonChip");
   if (is.na(match(chiptype, TYPE))) return(FALSE);
   return(TRUE);
}#isChipType

#------------------------------------------------------------------------------#
# isFilterCondition: check for correct filter condition
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
isFilterCondition <- function(condition) {
   if (debug.xps()) print("------isFilterCondition------")

   CONDITION <- c("percent", "samples", "mean", "percentile");
   if (is.na(match(condition, CONDITION))) return(FALSE);
   return(TRUE);
}#isFilterCondition

#------------------------------------------------------------------------------#
# isROOTFile: check if filename is a ROOT file, i.e. "filename.root"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
isROOTFile <- function(filename) {
   if (debug.xps()) print("------isROOTFile------")

   if (!file.exists(filename)) {
      return(FALSE);
   }#if

   rtf <- unlist(strsplit(filename, "/"));
   rtf <- rtf[length(rtf)];
   xtn <- substr(rtf, nchar(rtf)-3, nchar(rtf));
   if (length(xtn) == 0 || xtn != "root") {
      return(FALSE);
   }#if

   return(TRUE);
}#isROOTFile

#------------------------------------------------------------------------------#
# existsROOTFile: check if ROOT file exists already
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
existsROOTFile <- function(filename, tmp.rm = TRUE) {
   if (debug.xps()) print("------existsROOTFile------")

   if (!isROOTFile(filename)) {
      return(FALSE);
   }#if

   rtf <- unlist(strsplit(filename, "/"));
   rtf <- rtf[length(rtf)];
   rtf <- unlist(strsplit(rtf, "_"));
   tmp <- substr(rtf[1], 1, 3);

   if (length(rtf) == 1 && tmp != "tmp") {
      return(TRUE);
   }#if

   if (tmp.rm && tmp == "tmp") {
      return(FALSE);
   }#if

   return(TRUE);
}#existsROOTFile

#------------------------------------------------------------------------------#
# rootDirFile: return root file as /filedir/filename.root
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootDirFile <- function(filename=character(0), filedir=character(0)) {
   if (debug.xps()) print("------rootDirFile------")

   ## for case that filename=/path/filename
   if (length(filename) == 0) {
      filename <- "ERROR_MISSING_FILENAME";
   } else {
      filename <- basename(filename);
   }#if

   ## set default directory
   if (dirname(filename) != ".") {
      filedir <- dirname(filename);
   }#if
   if (filedir == "" || filedir == ".") {
      filedir <- as.character(getwd());
   }#if
   ## remove trailing "/"
   if (substr(filedir, nchar(filedir), nchar(filedir)) == "/") {
      filedir <- substr(filedir, 0, nchar(filedir)-1);
   }#if

   ## set rootfile to /fullpath/filename.root
   exten <- unlist(strsplit(filename, "\\."));
   if (exten[length(exten)] == "root") {
      rootfile <- paste(filedir, "/", filename, sep="");
   } else {
      rootfile <- paste(filedir, "/", filename, ".root", sep="");
   }#if

   return(rootfile);
}#rootDirFile

#------------------------------------------------------------------------------#
# rootDrawName: check for valid graphics type and return name to save
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootDrawName <- function(canvasname, graphtype) {
   if (debug.xps()) print("------rootDrawName------")

   if (graphtype == "") return(as.character(""));

   TYPE <- c("ps", "eps", "pdf", "jpg", "gif", "png", "tiff");
   if (is.na(match(graphtype, TYPE))) {
      stop(paste(sQuote(graphtype), "is not a valid graphics type"));
   }#if

   return(paste(canvasname, graphtype, sep="."));
}#rootDrawName

#------------------------------------------------------------------------------#
# validEstimatorOption: check for presence of valid M-estimator option
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validEstimatorOption <- function(estimator) {
   if (debug.xps()) print("------validEstimatorOption------")

   TYPE <- c("huber", "fair", "cauchy", "gemanmcclure", "welsch", "tukey", "andrew");
   if (is.na(match(estimator, TYPE))) {
      stop(paste(sQuote(estimator), "is not a valid M-estimator type"));
   }#if

   return(as.character(estimator));
}#validEstimatorOption

#------------------------------------------------------------------------------#
# validLogbase: check for presence of valid logbase
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validLogbase <- function(logbase) {
   if (debug.xps()) print("------validLogbase------")

   if (!(identical(logbase, "0")    || identical(logbase, "log") ||
         identical(logbase, "log2") || identical(logbase, "log10"))) {
      stop(paste(sQuote(logbase), "is not a valid logbase option"));
   }#if

   return(as.character(logbase));
}#validLogbase

#------------------------------------------------------------------------------#
# validMsg: utility function for setValidity
# taken from package BioBase: tools.R
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validMsg <- function(msg, result) {
   if (is.character(result)) {
      append(msg, result);
   } else {
      msg;
   }#if
}#validMsg

#------------------------------------------------------------------------------#
# validOption: check for presence of valid option
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validOption <- function(option) {
   if (debug.xps()) print("------validOption------")

   TYPE <- c("all", "together:none", "separate:none");
   if (is.na(match(option, TYPE))) {
      transcript <- unlist(strsplit(option, ":"))[1];
      transcript <- validTranscriptOption(transcript);

      opt <- substring(option, nchar(transcript)+2);
      if (is.na(match(opt, TYPE))) {
         stop(paste(sQuote("option"), "is not a valid method option"));
      }#if
   }#if

   return(as.character(option));
}#validOption

#------------------------------------------------------------------------------#
# validOutfile: check for presence of valid outfile = /filepath/filename
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validOutfile <- function(outname, outfile=character(0)) {
   if (debug.xps()) print("------validOutfile------")

   filename <- "";
   filepath <- paste(getwd(), "/", sep="");

   if (!(length(outfile) <= 0 || outfile == "")) {
      nsep <- gregexpr("/", outfile)[[1]];

      if (!(length(nsep) ==1 && nsep < 0)) {
         nsep <- nsep[length(nsep)];
         filepath <- substr(outfile, 1, nsep);
      }#if

      if (debug.xps()) print(paste("filepath = ",filepath));

#      if (!(filepath != "" && file.exists(filepath))) {
      if (!(filepath != "" && file.access(filepath) == 0)) {
         stop(paste("path of", sQuote("outfile"), "is not a valid path"));
      }#if

      filename <- substring(outfile, nsep+1);
   }#if

   if (debug.xps()) print(paste("filename = ",filename));

   if (filename == "") {
      filename <- paste(outname, "txt", sep=".");
      outfile  <- paste(filepath, filename, sep="");
   }#if

   if (debug.xps()) print(paste("outfile = ",outfile));

   return(outfile);
}#validOutfile

#------------------------------------------------------------------------------#
# validQualityOption: check for presence of valid quality option
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validQualityOption <- function(option, as.logical=FALSE) {
   if (debug.xps()) print("------validQualityOption------")

   TYPE <- c("raw", "adjusted", "normalized", "all");
#   if (is.na(match(option, TYPE))) {
   if (is.na(match(option, TYPE)[1])) {
      if (as.logical) return(FALSE);
      stop(paste(sQuote("option"), "is not a valid quality option"));
   }#if

   if (as.logical) return(TRUE);
   return(as.character(option));
}#validQualityOption

#------------------------------------------------------------------------------#
# validROOTFile: check if rootfile is a valid ROOT TFile, i.e. "myfile.root"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validROOTFile <- function(rootfile, treeset) {
   if (debug.xps()) print("------validROOTFile------")

   if (!isROOTFile(rootfile)) {
      stop(paste(sQuote("rootfile"), "of class", sQuote(treeset),
          "is missing or not a ROOT file *.root"));
   }#if

   filedir <- dirname(rootfile);
   if (!(is(filedir, "character") && file.exists(filedir))) {
      stop(paste(sQuote("filedir"), "is not a system directory"));
   }#if

   return(rootfile);
}#validROOTFile

#------------------------------------------------------------------------------#
# validSchemeTreeSet: check if "object" is a valid class with correct parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validSchemeTreeSet <- function(object) {
   if (debug.xps()) print("------validSchemeTreeSet------")

   if (!is (object, "SchemeTreeSet")) {
      stop(paste(sQuote("object"), "is not class", sQuote("SchemeTreeSet")));
   }#if

   ## check for correct settype
   settype  <- as.character(object@settype);
   if (settype != "scheme") {
      stop(paste(sQuote("object"), "is not of settype", sQuote("scheme")));
   }#if

   ## check for presence of root scheme file
   scheme   <- as.character(object@rootfile);
   scheme <- validROOTFile(scheme);

   ## check for presence of chip name
   chipname <- as.character(object@chipname);
   if (length(chipname) == 0 || nchar(chipname) < 1) {
      stop("missing chip name");
   }#if

   ## check for presence of correct chip type
   if(!isChipType(object@chiptype)) {
      stop(paste(sQuote("chiptype"), "of class", sQuote("SchemeTreeSet"),
                 "must be", sQuote("<GeneChip,GenomeChip,ExonChip>")));
   }#if

   return(object);
}#validSchemeTreeSet

#------------------------------------------------------------------------------#
# validSeparator: check for presence of valid separator and return extension
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validSeparator <- function(sep) {
   if (debug.xps()) print("------validSeparator------")

   if (identical(sep, " ") || identical(sep, "\t")) {
      return("txt");
   } else if (identical(sep, ",") || identical(sep, ";")) {
      return("csv");
   } else {
      stop(paste(sQuote(sep), "is not a valid separator"));
   }#if
}#validSeparator

#------------------------------------------------------------------------------#
# validTempDir: check for presence of valid temporary directory
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validTempDir <- function(tmpdir) {
   if (debug.xps()) print("------validTempDir------")

   if (tmpdir != "") {
      if (!file.exists(tmpdir)) {
         stop(paste(sQuote("tmpdir"), "is not a system directory"));
      }#if
      if (substr(tmpdir,nchar(tmpdir), nchar(tmpdir)) == "/") {
         tmpdir <- substr(tmpdir, 0, nchar(tmpdir)-1);
      }#if
   }#if

   return(tmpdir);
}#validTempDir

#------------------------------------------------------------------------------#
# validTranscriptOption: check for presence of valid transcript option
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validTranscriptOption <- function(transcript) {
   if (debug.xps()) print("------validTranscriptOption------")

   if (!(identical(transcript, "transcript") ||
         identical(transcript, "exon") ||
         identical(transcript, "probeset"))) {
      stop(paste(sQuote(transcript), "is not a valid transcript option"));
   }#if

   return(as.character(transcript));
}#validTranscriptOption

#------------------------------------------------------------------------------#
# validTreenames: remove potential reference tree
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validTreenames <- function(treenames, refname) {
   if (debug.xps()) print("------validTreenames------")

   treexten  <- extenPart(treenames);
   treenames <- namePart(treenames);
   reference <- match(refname, treenames);

   if (!is.na(reference)) treenames <- treenames[-reference];

   return(paste(treenames, treexten, sep="."));
}#validTreenames

#------------------------------------------------------------------------------#
# validTreetype: check for presence of valid tree extension
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
validTreetype <- function(treetype, datatype) {
   if (debug.xps()) print("------validTreetype------")

   if (datatype == "scheme") {
      TYPE <- c(SCMTYPE, EXNTYPE);
   } else if (datatype == "rawdata") {
      TYPE <- RAWTYPE;
   } else if (datatype == "preprocess") {
      TYPE <- PRETYPE;
   } else if (datatype == "normation") {
      TYPE <- NRMTYPE;
   } else if (datatype == "prefilter") {
      TYPE <- FLRTYPE[1];
   } else if (datatype == "unifilter") {
      TYPE <- FLRTYPE[2];
   } else if (datatype == "UnivariateAnalysis") {
      TYPE <- c(UNITYPE, FLRTYPE[2]);
   } else if (datatype == "all") {
      TYPE <- c(RAWTYPE, PRETYPE, NRMTYPE, FLRTYPE, UNITYPE);
   } else {
      stop(paste("invalid parameter", sQuote("datatype")));
   }#if

   if (is.na(match(treetype, TYPE))) {
      stop(paste("invalid treetype", sQuote(treetype),
                 "for datatype", sQuote(datatype)));
   }#if

   return(treetype);
}#validTreetype

#------------------------------------------------------------------------------#
# getDatatype: get data type for tree extension treetype
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getDatatype <- function(treetype) {
   if (debug.xps()) print("------getDatatype------")

   if (!is.na(match(treetype, RAWTYPE))) {
      datatype <- "rawdata";
   } else if (!is.na(match(treetype, PRETYPE))) {
      datatype <- "preprocess";
   } else if (!is.na(match(treetype, NRMTYPE))) {
      datatype <- "normation";
   } else {
      stop(paste("invalid parameter", sQuote("treetype")));
   }#if

   return(datatype);
}#getDatatype

#------------------------------------------------------------------------------#
# getDataXY: get (X,Y) from data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getDataXY <- function(object) {
   if (debug.xps()) print("------getDataXY------")

   ## get (x,y) from data
   if (min(dim(object@data)) == 0) {
      treename <- treeNames(object)[[1]];
      treetype <- extenPart(treename);
      data <- export(object,
                     treenames    = treename,
                     treetype     = treetype,
                     varlist      = "fX:fY",
                     outfile      = "dataXY.txt",
                     as.dataframe = TRUE,
                     verbose      = FALSE);
   } else {
      data  <- object@data[,c("X","Y")];
   }#if

   return(data);
}#getDataXY

#------------------------------------------------------------------------------#
# CELNames: utility function to get name part of CEL-file only
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CELNames <- function(celnames) {
   if (debug.xps()) print("------CELNames------")

   cel <- sapply(celnames,function(x){y<-unlist(strsplit(x, "/")); y[length(y)]});
   cel <- unlist(strsplit(cel, "\\.[cC][eE][lL]"));
   return(cel);
}#CELNames

#------------------------------------------------------------------------------#
# CELHeader: utility function to create header part for CEL-file (Version 3)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CELHeader <- function(celname, scheme) {
   if (debug.xps()) print("------CELHeader------")

   dathdr <- paste("DatHeader=[0..65534]  ",
                   celname,
                   ":CLS=1167 RWS=1167 XIN=3  YIN=3  VE=17        ",
                   "2.0 10/13/02 10:28:08       ",
                   chipName(scheme), ".1sq", 
                   "                  6",
                   sep="");

   header <- data.frame(matrix(nrow=23,ncol=1));

   header[1, 1] <- "[CEL]";
   header[2, 1] <- "Version=3";
   header[3, 1] <- "";
   header[4, 1] <- "[HEADER]";
   header[5, 1] <- paste("Cols=", ncols(scheme));
   header[6, 1] <- paste("Rows=", nrows(scheme));
   header[7, 1] <- paste("TotalX=", ncols(scheme));
   header[8, 1] <- paste("TotalY=", nrows(scheme));
   header[9, 1] <- "OffsetX=0";
   header[10,1] <- "OffsetY=0";
   header[11,1] <- "GridCornerUL=100 100";
   header[12,1] <- "GridCornerUR=1000 100";
   header[13,1] <- "GridCornerLR=1000 1000";
   header[14,1] <- "GridCornerLL=100 1000";
   header[15,1] <- "Axis-invertX=0";
   header[16,1] <- "AxisInvertY=0";
   header[17,1] <- "swapXY=0";
   header[18,1] <- dathdr;
   header[19,1] <- "Algorithm=Percentile";
   header[20,1] <- "AlgorithmParameters=Percentile:75;CellMargin:2;OutlierHigh:1.500;OutlierLow:1.004";
   header[21,1] <- "";
   header[22,1] <- "[INTENSITY]";
   header[23,1] <- paste("NumberCells=", nrows(scheme)*ncols(scheme));

   return(header);
}#CELHeader

#------------------------------------------------------------------------------#
# exonLevel: utility function to convert exonlevel to mask
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exonLevel <- function(exonlevel="", chiptype="GeneChip", as.sum=TRUE) {
   if (debug.xps()) print("------exonLevel------")

   exlevel <- 0; if (as.sum) exlevel <- c(0, 0, 0);

   if (chiptype == "GenomeChip") {
      ## separate levels for bgrd, norm, expr must be given as integers
      if (is.numeric(exonlevel)) {
         if (length(exonlevel) == 3 && identical(c(TRUE,TRUE,TRUE), exonlevel > 0)) {
            return(exonlevel);
         } else {
            stop(paste("numeric argument", sQuote("exonlevel"), "must be of length 3 and positive."));
         }#if
      } else if (exonlevel == "") {
         stop(paste("invalid argument", sQuote("exonlevel")));
      }#if

      ## levels are defined in src/XPSSchemes.h
      level <- unlist(strsplit(exonlevel, "\\+"));
      LEVEL <-    c("core", "metacore", "affx",           "all");
      CODE  <- list( 1024,   8192,       c(4,8,16,32),     9276);
      level <- match(level, LEVEL);

      if (is.na(all(level))) {
         stop(paste("invalid argument", sQuote("exonlevel")));
      }#if

      ## replace "core" with "core+metacore"
      if (any(level == 1)) level <- c(level, 1, 2);
      if (any(level == 4)) level <- c(1:3); 
      level <- unique(level);

      if (as.sum) {
         exlevel <- sum(unlist(CODE[level]));
         exlevel <- c(exlevel, exlevel, exlevel);
      } else {
         exlevel <- unlist(CODE[level]);
      }#if
   } else if (chiptype == "ExonChip") {
      ## separate levels for bgrd, norm, expr must be given as integers
      if (is.numeric(exonlevel)) {
         if (length(exonlevel) == 3 && identical(c(TRUE,TRUE,TRUE), exonlevel > 0)) {
            return(exonlevel);
         } else {
            stop(paste("numeric argument", sQuote("exonlevel"), "must be of length 3 and positive."));
         }#if
      } else if (exonlevel == "") {
         stop(paste("invalid argument", sQuote("exonlevel")));
      }#if

      ## levels are defined in src/XPSSchemes.h
      level <- unlist(strsplit(exonlevel, "\\+"));
      LEVEL <-    c("core", "metacore", "extended", "metaextended", "full", "metafull", "ambiguous", "affx",           "all");
      CODE  <- list( 1024,   8192,       512,        4096,           256,    2048,       128,         c(4,8,16,32),    16316);
      level <- match(level, LEVEL);

      if (is.na(all(level))) {
         stop(paste("invalid argument", sQuote("exonlevel")));
      }#if

      ## replace e.g. "core" with "core+metacore"
      if (any(level == 1)) level <- c(level, 1, 2);
      if (any(level == 3)) level <- c(level, 3, 4);
      if (any(level == 5)) level <- c(level, 5, 6);
      if (any(level == 9)) level <- c(1:8);
      level <- unique(level);

      if (as.sum) {
         exlevel <- sum(unlist(CODE[level]));
         exlevel <- c(exlevel, exlevel, exlevel);
      } else {
         exlevel <- unlist(CODE[level]);
      }#if
   }#if

   return(exlevel);
}#exonLevel

#------------------------------------------------------------------------------#
# exonLevelIDs: utility function to get row numbers for exonlevel
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exonLevelIDs <- function(exonlevel, data, mask, ncol) {
   if (debug.xps()) print("------exonLevelIDs------")

   mask <- cbind((mask[, "X"] + 1 + ncol*mask[,"Y"]), mask);
   data <- cbind((data[, "X"] + 1 + ncol*data[,"Y"]), data);
   colnames(mask)[1] <- "XY";
   colnames(data)[1] <- "XY";

   id <- sapply(unique(exonlevel),
                function(x) {
                   xy <- mask[mask[, "Mask"] == x,];
                   id <- match(xy[,"XY"], data[,"XY"]);
                }
         );
   id <- unique(unlist(id));

   return(id[order(id)]);
}#exonLevelIDs

#exonLevelIDs <- function(exonlevel, data, mask, ncol) {
#   if (debug.xps()) print("------exonLevelIDs------")
#
#   id <- sapply(exonlevel,
#                function(x) {
#                   xy <- mask[mask[, "Mask"] == x, c("X","Y")];
#                   id <- match(xy[,"X"]   + ncol*xy[,"Y"],
#                               data[,"X"] + ncol*data[,"Y"]);
#                }
#         );
#   id <- unique(unlist(id));
#
#   return(id[order(id)]);
#}#exonLevelIDs

#------------------------------------------------------------------------------#
# extenPart: utility function to extract extension part from "name.exten"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
extenPart <- function(names, as.unique=TRUE) {
   if (debug.xps()) print("------extenPart------")

   subPart <- function(strg, pat) {
      pos <- which(unlist(strsplit(strg,"")) == pat);
      pos <- ifelse(length(pos), pos[length(pos)] + 1, 1);
      return(substr(strg, pos, 99999));
   }

   strg <- unlist(names);
   strg <- sapply(strg, function(x) subPart(x, "."));
   if (as.unique == TRUE) {
      strg <- unique(strg);
   }#if
 
   return(strg[!is.na(strg)]);
}#extenPart

#------------------------------------------------------------------------------#
# namePart: utility function to extract name part from "name.exten"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
namePart <- function(names) {
   if (debug.xps()) print("------namePart------")

   subPart <- function(strg, pat) {
      pos <- which(unlist(strsplit(strg,"")) == pat);
      pos <- ifelse(length(pos), pos[length(pos)] - 1, 99999);
      return(substr(strg, 1, pos));
   }

   strg <- unlist(names);
   strg <- sapply(strg, function(x) subPart(x, "."));
 
   return(strg[!is.na(strg)]);
}#namePart

#------------------------------------------------------------------------------#
# type2Exten: utility function to get extension for method type
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type2Exten <- function(type, datatype) {
   if (debug.xps()) print("------type2Exten------")

   ## type and extensions are defined in src/XPSProcessing.cxx
   if (datatype == "preprocess") {
      TYPE <- c("mean", "median", "lowess", "supsmu", "quantile",
                "avgdiff", "tukeybiweight", "medianpolish",
                "farms", "dfw", "firma", "plm", "rlm");
      XTEN <- c("cmn", "cmd", "clw", "css", "cqu", "adf", "tbw",
                "mdp", "frm", "dfw", "fir", "plm", "rlm");
   } else if (datatype == "normation") {
      TYPE <- c("mean", "median", "ksmooth", "lowess", "supsmu",
                "quantile", "medianpolish");
      XTEN <- c("tmn", "med", "ksm", "low", "sup", "qua", "mdp");
   } else {
      stop(paste("invalid parameter", sQuote("datatype")));
   }#if
 
   return(XTEN[match(type, TYPE)]);
}#type2Exten

#------------------------------------------------------------------------------#
# listTreeNames: utility function to get treenames and /path/treenames
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"listTreeNames" <- function(object) {
   if (debug.xps()) print("------listTreeNames------")

   datafile  <- object@rootfile;
   setname   <- object@setname;
   treenames <- as.character(object@treenames);

   ## get treenames as fullnames=/datafile/setname/treenames
   if (length(unlist(strsplit(treenames[1],"/"))) == 1) {
      fullnames <- paste(datafile, setname, treenames, sep="/");
   } else {
      fullnames <- treenames;
      treenames <- as.character(sapply(treenames,
                                       function(x){tree <- unlist(strsplit(x, "/"));
                                                   tree[length(tree)];}));
   }#if

   return(list(treenames = treenames, fullnames = fullnames));
}#listTreeNames


################################################################################
# utility functions for graphics
################################################################################

#------------------------------------------------------------------------------#
# circle: plot circle
# adapted from package calibrate: circle.R
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"circle" <- function(radius = 1.0, col = "black", lty = 1, ...) {
   x <- seq(-radius, radius, by=0.01);
   y <- sqrt(radius^2 - x^2);
   lines(x,  y, col=col, lty=lty, ...);
   lines(x, -y, col=col, lty=lty, ...);
}#circle

#------------------------------------------------------------------------------#
# distMAD: calculate matrix of the mad of all pairwise differences of columns of x
# adapted from package genefilter: dist2.R
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"distMAD" <- function(x) {
   m <- matrix(0, ncol=ncol(x), nrow=ncol(x));
   rownames(m) <- colnames(x);
   colnames(m) <- colnames(x);

   if (ncol(x) == 1) return(m);

   for(j in 2:ncol(x)) {
      for(i in 1:(j-1)) {
         m[i, j] = m[j, i] = mad(x[,i] - x[,j]);
      }#for_i
   }#for_j

   return(m);
}#distMAD

#------------------------------------------------------------------------------#
# plotDensities: function to plot density
# taken from package affy: plot.density.R
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plotDensities <- function(mat,
                          ylab="density", xlab="x", type="l", lty=1:5, col=1:6,
                          ...) 
{
   x.density <- apply(mat, 2, density);

   all.x <- do.call("cbind", lapply(x.density, function(x) x$x));
   all.y <- do.call("cbind", lapply(x.density, function(x) x$y));
  
   matplot(all.x, all.y, ylab=ylab, xlab=xlab, type=type, lty=lty, col=col, ...);

   invisible(list(all.x=all.x, all.y=all.y));
}#plotDensities

#------------------------------------------------------------------------------#
# pseudoPalette: function for coloring pseudo chip images.
# taken from package affyPLM: PLMset.R
# Copyright (C) 2003-2008     Ben Bolstad
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pseudoPalette <-function (low="white", high=c("green", "red"), mid=NULL, k=50) {
   if (debug.xps()) print("------pseudoPalette------")

   low  <- col2rgb(low)/255;
   high <- col2rgb(high)/255;

   if (is.null(mid)) {
      r <- seq(low[1], high[1], len = k);
      g <- seq(low[2], high[2], len = k);
      b <- seq(low[3], high[3], len = k);
   }#if

   if (!is.null(mid)) {
      k2  <- round(k/2);
      mid <- col2rgb(mid)/255;

      r <- c(seq(low[1], mid[1], len = k2), seq(mid[1], high[1], len = k2));
      g <- c(seq(low[2], mid[2], len = k2), seq(mid[2], high[2], len = k2));
      b <- c(seq(low[3], mid[3], len = k2), seq(mid[3], high[3], len = k2));
   }#if

   rgb(r, g, b);
}#pseudoPalette

#------------------------------------------------------------------------------#
# adjustXLabMargins: adjust bottom margin, pointsize, plot width
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
adjustXLabMargins <- function(xlab, bottom=6, cex=1.0, width=800) {
   if (debug.xps()) print("------adjustXLabMargins------")

   b <- bottom;
   x <- cex;
   w <- width;

   ## get bottom margin
   maxchar <- max(nchar(xlab));
   if (maxchar > 8)  {
      bottom <- max(as.integer(maxchar/2.0), b);
      cex    <- x;
   }#if
   if (maxchar > 16) {
      bottom <- max(as.integer(maxchar/1.8), b);
      cex    <- x - 0.2;
   }#if

   ## get plot width
   numlab <- length(xlab);
   if (numlab > 60)  {
      width <- numlab * 16 + 120;
      cex   <- x - 0.2;
   }#if
   if (numlab > 120) {
      width <- numlab * 14 + 120;
      cex   <- x - 0.3;
   }#if
   if (width < w) width <- w;

   return(list(b=bottom, cex=cex, w=width));
}#adjustXLabMargins


################################################################################
# utility functions based on rwrapper.h C-functions
################################################################################

#------------------------------------------------------------------------------#
# getNameType: get chip name and type from root scheme file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getNameType <- function(rootfile) {
   if (debug.xps()) print("------getNameType------")

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "none");

   ## get chip name and type
   info <- .C("ChipNameType",
              as.character(rootfile),
              nametype=character(2),
              PACKAGE="xps")$nametype;

   return(list(chipname=info[1], chiptype=info[2]));
}#getNameType

#------------------------------------------------------------------------------#
# getChipName: get chip name from root scheme file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getChipName <- function(rootfile) {
   return(getNameType(rootfile)$chipname);
}#getChipName

#------------------------------------------------------------------------------#
# getChipType: get chip type from root scheme file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getChipType <- function(rootfile) {
   return(getNameType(rootfile)$chiptype);
}#getChipType

#------------------------------------------------------------------------------#
# getProbeInfo: get GeneChip probe information from root scheme file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getProbeInfo <- function(rootfile) {
   if (debug.xps()) print("------probeInfo------")

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "none");

   ## get probe info
   info <- .C("GeneChipProbeInfo",
              as.character(rootfile),
              value=integer(8),
              PACKAGE="xps")$value;

   return(list(nrows=info[1], ncols=info[2], nprobes=info[3], ncontrols=info[4],
            ngenes=info[5], nunits=info[6], nprobesets=info[7], naffx=info[8]));
}#getProbeInfo

#------------------------------------------------------------------------------#
# getNumberTrees: get number of trees saved in root file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getNumberTrees <- function(rootfile, treetype="*", setname=NULL) {
   if (debug.xps()) print("------getNumberTrees------")

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "none");

   ## get number of trees
   if (is.null(setname)) {
      ntrees <- .C("GetNumberOfTrees4Exten",
                   as.character(rootfile),
                   as.character(treetype),
                   numtrees=integer(1),
                   PACKAGE="xps")$numtrees;
   } else {
      ntrees <- .C("GetNumberOfTrees",
                   as.character(rootfile),
                   as.character(setname),
                   as.character(treetype),
                   numtrees=integer(1),
                   PACKAGE="xps")$numtrees;
   }#if
   return(ntrees);
}#getNumberTrees

#------------------------------------------------------------------------------#
# getTreeNames: get names of trees saved in root file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getTreeNames <- function(rootfile, treetype="*", setname=NULL, gettitle=FALSE) {
   if (debug.xps()) print("------getTreeNames------")

   ## check for presence of root file
   rootfile <- validROOTFile(rootfile, "none");

   ## get tree names
   if (is.null(setname)) {
      ntrees <- .C("GetNumberOfTrees4Exten",
                   as.character(rootfile),
                   as.character(treetype),
                   numtrees=integer(1),
                   PACKAGE="xps")$numtrees;

      tnames <- .C("GetTreeNames4Exten",
                   as.character(rootfile),
                   as.character(treetype),
                   as.integer(gettitle),
                   treenames=character(ntrees),
                   PACKAGE="xps")$treenames;
   } else {
      ntrees <- .C("GetNumberOfTrees",
                   as.character(rootfile),
                   as.character(setname),
                   as.character(treetype),
                   numtrees=integer(1),
                   PACKAGE="xps")$numtrees;

      tnames <- .C("GetTreeNames",
                   as.character(rootfile),
                   as.character(setname),
                   as.character(treetype),
                   as.integer(gettitle),
                   treenames=character(ntrees),
                   PACKAGE="xps")$treenames;
   }#if

#   return(tnames);
   return(as.vector(tnames, mode="character"));
}#getTreeNames


################################################################################
# other utility functions
################################################################################

#------------------------------------------------------------------------------#
# metaProbesets: create file containing meta probeset definitions for "apt-probeset-summarize"
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
metaProbesets <- function(xps.scheme, infile = character(0),
                          outfile = character(0), exonlevel="metacore") {
   if (debug.xps()) print("------metaProbesets------")

   ## check for presence of exon scheme file
   if (!(is(xps.scheme, "SchemeTreeSet") &&
       (chipType(xps.scheme) == "ExonChip") &&
       (file.exists(rootFile(xps.scheme))))) {
      stop(paste("File", sQuote("xps.scheme"), "is not a SchemeTreeSet of type ExonChip."));
   }#if
   schemefile <- rootFile(xps.scheme);

   ## check for presence of infile
   if (!file.exists(infile)) {
      stop(paste("File", sQuote(infile), "does not exist."));
   }#if

    ## levels are defined in src/XPSSchemes.h
   level <- unlist(strsplit(exonlevel, "\\+"));
   LEVEL <- c("core", "metacore", "extended", "metaextended", "full", "metafull", "ambiguous", "affx", "all");
#   CODE  <- c( 1024,   8192,       512,        4096,           256,    2048,       128,         60,    16316);
   CODE  <- c( 1024,   1024,       512,        512,            256,    256,        128,         0,     1920);
   level <- match(level, LEVEL);

   if (is.na(all(level))) {
      stop(paste("invalid argument", sQuote("exonlevel")));
   }#if

   ## replace e.g. "all" with "core+extended+full+ambiguous"
   if (any(level == 9)) level <- c(1, 3, 5, 7);
   level <- unique(level);
   level <- sum(CODE[level]);

   ## check for "meta"-level
   meta <- 0;
   if (substr(exonlevel, 1, 4) == "meta") {
      meta <- 1;
   }#if

   ## create metaprobeset file
   r <- .C("MetaProbesets",
           as.character(schemefile),
           as.character(infile),
           as.character(outfile),
           as.integer(level),
           as.integer(meta),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in utility function", sQuote("MetaProbesets")));
   }#if
}#metaProbesets

#------------------------------------------------------------------------------#

