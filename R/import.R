#==============================================================================#
# import.R: import scheme file and data file functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# import.expr.scheme:   import CDF-file and annotation files for exression array
# import.genome.scheme: import CLF, PGF and annotation files for whole genome array
# import.exon.scheme:   import CLF, PGF and annotation files for exon array
# import.data:          import CEL-files
#==============================================================================#

"import.expr.scheme" <-
function(filename   = character(0),
         filedir    = getwd(),
         schemefile = character(0),
         probefile  = character(0),
         annotfile  = character(0),
         chipname   = NULL,
         add.mask   = FALSE,
         verbose    = TRUE)
{
   if (debug.xps()) print("------import.expr.scheme------")

   ## check for presence of filename
   if (filename == "") {
      stop(paste(sQuote("filename"), "is missing"))
   }#if

   ## root scheme file
   rootfile <- rootDirFile(filename, filedir);

   ## check for presence of parameters
   if (!(is.character(schemefile) &&
         is.character(probefile) &&
         is.character(annotfile))) {
      stop("all arguments must be of type character");
   }#if

   ## check for presence of scheme file
   if (!file.exists(schemefile)) {
      stop("missing chip definition file *.CDF");
   }#if

   ## get chipname and extension
   scheme <- unlist(strsplit(schemefile, "/"));
   scheme <- scheme[length(scheme)];
   scheme <- unlist(strsplit(scheme, "\\."));

   ## get chipname from scheme if not given
   if (is.null(chipname)) {
      chipname <- scheme[1];
   }#if
   exten <- scheme[length(scheme)];

   if (!((exten == "CDF") || (exten == "cdf"))) {
      stop(paste(sQuote("schemefile"), "is not a CDF-file"))
   }#if

   ## check for presence of probe file
   if (!file.exists(probefile)) {
      stop("missing probe file");
   }#if

   probe <- unlist(strsplit(probefile, "/"));
   probe <- probe[length(probe)];
   exten <- substr(probe, nchar(probe)-8, nchar(probe));
   if (exten != "probe.tab") {
      stop(paste(sQuote("probefile"),
                 "does not end with", sQuote("probe.tab")));
   }#if

   ## check for presence of annotation file
   if (length(annotfile) == 0 || annotfile == "") {
      warning("no annotation file is given");
      annotfile <- "";
   } else if (!file.exists(annotfile)) {
      stop("missing chip annotation file");
   } else {
      ann   <- unlist(strsplit(annotfile, "/"));
      ann   <- ann[length(ann)];
      exten <- substr(ann, nchar(ann)-8, nchar(ann));

      if (exten != "annot.csv") {
         stop(paste(sQuote("annotfile"),
                    "does not end with", sQuote("annot.csv")));
      }#if
   }#if

   ## create root scheme file for expression array
   r <- .C("ImportExprSchemes",
           as.character(filename),
           as.character(filedir),
           as.character(chipname),
           as.character(schemefile),
           as.character(probefile),
           as.character(annotfile),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ImportExprSchemes")));
      return(NULL);
   }#if

   ## get treenames and probe information
   treenames <- as.list(getTreeNames(rootfile));
   probeinfo <- as.list(getProbeInfo(rootfile));

   ## create new class
   set <- new("SchemeTreeSet",
              setname   = chipname,
              settype   = "scheme",
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = length(treenames),
              treenames = as.list(treenames),
              chipname  = chipname,
              chiptype  = "GeneChip",
              probeinfo = probeinfo);

   ## get mask for scheme
   if (add.mask) {
      chipMask(set) <- export(set,
                              treetype     = "scm",
                              varlist      = "fMask",
                              as.dataframe = TRUE,
                              verbose      = TRUE);
   }#if

   return(set);
}#import.expr.scheme

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"import.genome.scheme" <-
function(filename   = character(0),
         filedir    = getwd(),
         layoutfile = character(0),
         schemefile = character(0),
         transcript = character(0),
         add.mask   = FALSE,
         verbose    = TRUE)
{
   if (debug.xps()) print("------import.genome.scheme------")

   ## check for presence of filename
   if (filename == "") {
      stop(paste(sQuote("filename"), "is missing"))
   }#if

   ## root scheme file
   rootfile <- rootDirFile(filename, filedir);

   ## check for presence of parameters
   if (!(is.character(layoutfile) && is.character(schemefile) &&
         is.character(transcript))) {
      stop("all arguments must be of type character");
   }#if

   ## check for presence of layout file
   if (!file.exists(layoutfile)) {
      stop("missing CEL Layout File *.CLF");
   }#if

   scheme <- unlist(strsplit(layoutfile, "/"));
   scheme <- scheme[length(scheme)];
   scheme <- unlist(strsplit(scheme, "\\."));
   exten  <- scheme[length(scheme)];

   if (!((exten == "CLF") || (exten == "clf"))) {
      stop(paste("argument", sQuote("layoutfile"), "is not a CLF-file"))
   }#if

   ## check for presence of scheme file
   if (!file.exists(schemefile)) {
      stop("missing Probe Group File *.PGF");
   }#if

   ## get chipname and extension
   scheme <- unlist(strsplit(schemefile, "/"));
   scheme <- scheme[length(scheme)];
   scheme <- unlist(strsplit(scheme, "\\."));

   chipname <- scheme[1];
   exten    <- scheme[length(scheme)];

   if (!((exten == "PGF") || (exten == "pgf"))) {
      stop(paste(sQuote("schemefile"), "is not a PGF-file"))
   }#if

   ## check for presence of transcript annotation file
   if (!file.exists(transcript)) {
      stop("missing transcript annotation file");
   }#if

   annot <- unlist(strsplit(transcript, "/"));
   annot <- annot[length(annot)];
   exten <- substr(annot, nchar(annot)-2, nchar(annot));
   if (!((length(grep("transcript",annot)) > 0) && (exten == "csv"))) {
      stop(paste(sQuote("transcript"), "is not a transcript annotation file"));
   }#if

   ## create root scheme file for whole genome array
   r <- .C("ImportGenomeSchemes",
           as.character(filename),
           as.character(filedir),
           as.character(chipname),
           as.character(layoutfile),
           as.character(schemefile),
           as.character(transcript),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ImportGenomeSchemes")));
      return(NULL);
   }#if

   ## get treenames and probe information
   treenames <- as.list(getTreeNames(rootfile));
   probeinfo <- as.list(getProbeInfo(rootfile));

   ## create new class
   set <- new("SchemeTreeSet",
              setname   = chipname,
              settype   = "scheme",
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = length(treenames),
              treenames = as.list(treenames),
              chipname  = chipname,
              chiptype  = "GenomeChip",
              probeinfo = probeinfo);

   ## get mask for scheme
   if (add.mask) {
      chipMask(set) <- export(set,
                              treetype     = "scm",
                              varlist      = "fMask",
                              as.dataframe = TRUE,
                              verbose      = TRUE);
   }#if

   return(set);
}#import.genome.scheme

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"import.exon.scheme" <-
function(filename   = character(0),
         filedir    = getwd(),
         layoutfile = character(0),
         schemefile = character(0),
         probeset   = character(0),
         transcript = character(0),
         control    = "",
         add.mask   = FALSE,
         verbose    = TRUE)
{
   if (debug.xps()) print("------import.exon.scheme------")

   ## check for presence of filename
   if (filename == "") {
      stop(paste(sQuote("filename"), "is missing"))
   }#if

   ## root scheme file
   rootfile <- rootDirFile(filename, filedir);

   ## check for presence of parameters
   if (!(is.character(layoutfile) && is.character(schemefile) &&
         is.character(probeset) && is.character(transcript))) {
      stop("all arguments must be of type character");
   }#if

   ## check for presence of layout file
   if (!file.exists(layoutfile)) {
      stop("missing CEL Layout File *.CLF");
   }#if

   scheme <- unlist(strsplit(layoutfile, "/"));
   scheme <- scheme[length(scheme)];
   scheme <- unlist(strsplit(scheme, "\\."));
   exten  <- scheme[length(scheme)];

   if (!((exten == "CLF") || (exten == "clf"))) {
      stop(paste("argument", sQuote("layoutfile"), "is not a CLF-file"))
   }#if

   ## check for presence of scheme file
   if (!file.exists(schemefile)) {
      stop("missing Probe Group File *.PGF");
   }#if

   ## get chipname and extension
   scheme <- unlist(strsplit(schemefile, "/"));
   scheme <- scheme[length(scheme)];
   scheme <- unlist(strsplit(scheme, "\\."));

   chipname <- scheme[1];
   exten    <- scheme[length(scheme)];

   if (!((exten == "PGF") || (exten == "pgf"))) {
      stop(paste(sQuote("schemefile"), "is not a PGF-file"))
   }#if

   ## check for presence of probeset annotation file
   if (!file.exists(probeset)) {
      stop("missing probeset annotation file");
   }#if

   annot <- unlist(strsplit(probeset, "/"));
   annot <- annot[length(annot)];
   exten <- substr(annot, nchar(annot)-2, nchar(annot));
   if (!((length(grep("probeset",annot)) > 0) && (exten == "csv"))) {
      stop(paste(sQuote("probeset"), "is not a probeset annotation file"));
   }#if

   ## check for presence of transcript annotation file
   if (!file.exists(transcript)) {
      stop("missing transcript annotation file");
   }#if

   annot <- unlist(strsplit(transcript, "/"));
   annot <- annot[length(annot)];
   exten <- substr(annot, nchar(annot)-2, nchar(annot));
   if (!((length(grep("transcript",annot)) > 0) && (exten == "csv"))) {
      stop(paste(sQuote("transcript"), "is not a transcript annotation file"));
   }#if

   ## check for presence of control annotation file
   if (control != "" && !file.exists(control)) {
      stop("control file does not exist");
   }#if

   if (control != "") {
      warning("control file must only be given for old annotation files!!!");
   }#if

   ## create root scheme file for exon array
   r <- .C("ImportExonSchemes",
           as.character(filename),
           as.character(filedir),
           as.character(chipname),
           as.character(layoutfile),
           as.character(schemefile),
           as.character(probeset),
           as.character(transcript),
           as.character(control),
           as.integer(verbose),
           err=integer(1),
           PACKAGE="xps")$err;

   if (r != 0) {
      stop(paste("error in function", sQuote("ImportExonSchemes")));
      return(NULL);
   }#if

   ## get treenames and probe information
   treenames <- as.list(getTreeNames(rootfile));
   probeinfo <- as.list(getProbeInfo(rootfile));

   ## create new class
   set <- new("SchemeTreeSet",
              setname   = chipname,
              settype   = "scheme",
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = length(treenames),
              treenames = as.list(treenames),
              chipname  = chipname,
              chiptype  = "ExonChip",
              probeinfo = probeinfo);

   ## get mask for scheme
   if (add.mask) {
      chipMask(set) <- export(set,
                              treetype     = "scm",
                              varlist      = "fMask",
                              as.dataframe = TRUE,
                              verbose      = TRUE);
   }#if

   return(set);
}#import.exon.scheme

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"import.data" <-
function(xps.scheme,
         filename = character(0),
         filedir  = getwd(),
         celdir   = NULL,
         celfiles = "*",
         celnames = NULL,
         project  = NULL,
         verbose  = TRUE)
{
   if (debug.xps()) print("------import.data------")

   ## check for presence of class xps.scheme
   xps.scheme <- validSchemeTreeSet(xps.scheme);

   ## extract scheme, chipname and chiptype from xps.scheme
   settype  <- setType(xps.scheme);
   scheme   <- rootFile(xps.scheme);
   chipname <- chipName(xps.scheme);
   chiptype <- chipType(xps.scheme);

   ## check for presence of parameters
   if (!(is.character(filename) &&
         is.character(filedir))) {
      stop("arguments <filename,filedir> must be of type character");
   }#if
   if (!file.exists(filedir)) {
      stop(paste("directory", sQuote(filedir), "does not exist."));
   }#if

   ## root data file
   rootdata <- paste(filename, "_cel.root", sep="");
   rootfile <- rootDirFile(rootdata, filedir);

   ## check if root file exists (necessary for WinXP to test already here)
   if (existsROOTFile(rootfile)) {
      stop(paste("ROOT file", sQuote(rootfile), "does already exist."));
   }#if

   ## check for presence of cel-file directory
   if (!is.null(celdir)) {
      if (celdir == "") {
         celdir <- getwd();
      }#if
      if (!file.exists(celdir)) {
         stop(paste("argument", sQuote("celdir"), "is not a system directory"));
      }#if
      if (substr(celdir,nchar(celdir),nchar(celdir)) != "/") {
         celdir <- paste(celdir, "/", sep="");
      }#if
   }#if

   ## get cel-files
   if (celfiles[1] == "*") {
      celfiles <- list.files(celdir, pattern = "\\.[cC][eE][lL]");
   } else if (is.null(celdir)) {
      existing <- file.exists(celfiles);
      if (length(existing[existing==TRUE]) != length(celfiles)) {
         stop(paste("some", sQuote("celfiles"), " are not found."));
      }#if
   } else {
      tmpfiles <- list.files(celdir, pattern = "\\.[cC][eE][lL]");
      tmpfiles <- match(celfiles, tmpfiles);
      if (length(tmpfiles[is.na(tmpfiles)]) > 0) {
         stop(paste("some", sQuote("celfiles"), " are not found in directory",
                    sQuote(celdir)));
      }#if
   }#if

   ## get number of cel-files
   ncel <- length(celfiles);
   if (ncel == 0) {
      stop(paste("directory", sQuote(celdir), "does not contain CEL-files"));
   }#if

   ## check for presence of potential celnames
   if (!is.null(celnames) && (length(unique(celnames)) != ncel)) {
      stop(paste(sQuote("celnames"), "must have same number of entries as",
                 sQuote("celfiles")));
   }#if
   if (is.null(celnames)) {
      celnames <- CELNames(celfiles);
   } else {
      celnames <- unlist(strsplit(celnames, "\\.[cC][eE][lL]"));
   }#if

   ## concatenate celfiles with celdir
   tmpfiles <- unlist(celfiles);
   celfiles <- paste(celdir, celfiles, sep="");

   ## tree names from celfiles
   if (length(celnames) == 1 && celnames == "") {
      celnames <- unique(unlist(strsplit(tmpfiles, "\\.[cC][eE][lL]")));
   }#if

   ## check celnames for presence of characters [](){}.:# etc and replace with "_"
   tmpnames <- make.names(celnames);
   if (length(unique(tmpnames %in% celnames)) > 1 || unique(tmpnames %in% celnames) == FALSE) {
      if (verbose) warning("characters [](){}.:# etc in 'celnames' will be replaced  with '_'");
      celnames <- gsub("\\.", "_", tmpnames);
   }#if

   ## project information
   projct  <- ""; nproject <- 0;
   author  <- ""; nauthor  <- 0;
   dataset <- ""; ndataset <- 0;
   source  <- ""; nsource  <- 0;
   sample  <- ""; nsample  <- 0;
   cell    <- ""; ncell    <- 0;
   pcell   <- ""; npcell   <- 0;
   tissue  <- ""; ntissue  <- 0;
   biopsy  <- ""; nbiopsy  <- 0;
   array   <- ""; narray   <- 0;
   hyb     <- ""; nhyb     <- 0;
   treat   <- ""; ntreat   <- 0;

   if (!is.null(project) && is(project, "ProjectInfo")) {
      nproject <- length(projectInfo(project));
      nauthor  <- length(authorInfo(project));
      ndataset <- length(datasetInfo(project));
      nsource  <- length(sourceInfo(project));
      nsample  <- length(sampleInfo(project));
      ncell    <- length(cellineInfo(project));
      npcell   <- length(primcellInfo(project));
      ntissue  <- length(tissueInfo(project));
      nbiopsy  <- length(biopsyInfo(project));
      narray   <- length(arrayInfo(project));
      nhyb     <- length(unlist(hybridizInfo(project)));
      ntreat   <- length(unlist(treatmentInfo(project)));

      if (nproject == 5)  projct  <- projectInfo(project)  else projct  <- "";
      if (nauthor  == 8)  author  <- authorInfo(project)   else author  <- "";
      if (ndataset == 7)  dataset <- datasetInfo(project)  else dataset <- "";
      if (nsource  == 6)  source  <- sourceInfo(project)   else source  <- "";
      if (nsample  == 12) sample  <- sampleInfo(project)   else sample  <- "";
      if (ncell    == 15) cell    <- cellineInfo(project)  else cell    <- "";
      if (npcell   == 14) pcell   <- primcellInfo(project) else pcell   <- "";
      if (ntissue  == 19) tissue  <- tissueInfo(project)   else tissue  <- "";
      if (nbiopsy  == 19) biopsy  <- biopsyInfo(project)   else biopsy  <- "";
      if (narray   == 4)  array   <- arrayInfo(project)    else array   <- "";

      if ((narray == 4) & (array$chipname != chipname)) {
         stop(paste(sQuote("arrayInfo$chipname"), " must be", sQuote(chipname)));
      }#if
      if ((narray == 4) & (array$chiptype != chiptype)) {
         stop(paste(sQuote("arrayInfo$chiptype"), " must be", sQuote(chiptype)));
      }#if

###########
# to do: check if hybnames and treatnames are identical to celnames
###########
      if (nhyb >= 9) {
         hyb <- as.vector(unlist(t(hybridizInfo(project))));
      } else {
         hyb <- "";
      }#if
      if (ntreat >= 8) {
         treat <- as.vector(unlist(t(treatmentInfo(project))));
      } else {
         treat <- "";
      }#if
   }#if

   ## define treeset
   setname <- "DataSet";

   ## create root data file
   r <- .C("ImportData",
           as.character(filename),
           as.character(filedir),
           as.character(chiptype),
           as.character(chipname),
           as.character(scheme),
           as.character(setname),
           as.character(celfiles),
           as.character(celnames),
           as.integer(ncel),
           as.character(projct),
           as.integer(nproject),
           as.character(author),
           as.integer(nauthor),
           as.character(dataset),
           as.integer(ndataset),
           as.character(source),
           as.integer(nsource),
           as.character(sample),
           as.integer(nsample),
           as.character(cell),
           as.integer(ncell),
           as.character(pcell),
           as.integer(npcell),
           as.character(tissue),
           as.integer(ntissue),
           as.character(biopsy),
           as.integer(nbiopsy),
           as.character(array),
           as.integer(narray),
           as.character(hyb),
           as.integer(nhyb),
           as.character(treat),
           as.integer(ntreat),
           as.integer(0),
           as.integer(0),
           as.integer(verbose),
           result=character(2),
           PACKAGE="xps")$result;

   ## returned result: saved rootfile and error
   rootfile <- r[1];
   error    <- as.integer(r[2]);

   if (error != 0) {
      stop(paste("error in function", sQuote("ImportData")));
      return(NULL);
   }#if

   celnames <- paste(celnames, ".cel", sep="");

   ## create new class
   set <- new("DataTreeSet",
              setname   = setname,
              settype   = "rawdata",
              rootfile  = rootfile,
              filedir   = filedir,
              numtrees  = as.numeric(ncel),
              treenames = as.list(celnames),
              scheme    = xps.scheme);

   ## add project information
   if (!is.null(project) && is(project, "ProjectInfo")) {
      projectInfo(set) <- project;
   }#if

   return(set);
}#import.data

#==============================================================================#
