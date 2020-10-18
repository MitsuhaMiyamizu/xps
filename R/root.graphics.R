#==============================================================================#
# root.graphics.R: root functions for drawing root graphics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# root.density: 
# root.image:
# root.profile: 
# root.draw: 
# root.draw.leaves:
# root.graph1D: 
# root.graph2D: 
# root.mvaplot: 
# root.hist1D: 
# root.hist2D: 
# root.hist3D: 
#==============================================================================#

"root.density" <-
function(x,
         treename   = "*",
         logbase    = "log2",
         canvasname = "DensityPlot",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.density------")

   rootfile  <- rootFile(x);
   treenames <- treeNames(x);
   savename  <- rootDrawName(canvasname, save.as);

   ## get tree extension
   exten <- unlist(treenames)[1];
   exten <- unlist(strsplit(exten, "\\."))[2];

   ## test for valid treename
   namepart <- unlist(strsplit(treename, "\\."))[1];
   if (namepart == "*") {
      treename <- paste(namepart, exten, sep=".");
   } else if (treename%in%treenames == FALSE) {
      stop(paste("treename", sQuote(treename), "is not present in", sQuote(rootfile)));
   }#if

   ## add ROOT directory in rootfile and select leafname
   treename <- paste(setName(x), treename, sep="/");
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else {
      stop(paste(sQuote("x"), "must be  class <DataTreeSet,ExprTreeSet>"));
   }#if

   varlist <- leafname;
   option  <- "L";

   ## unix vs windows settings
   is.win <- (.Platform$OS.type == "windows");
   if (is.win) {
      sq  <- as.character("");  #
      dq  <- as.character("\\\""); #\"
   } else {
      sq  <- as.character("'");  #'
      dq  <- as.character("\""); #"
   }#if

   ## get path to library xps.so
   xpsso <- system.file("libs", .Platform$r_arch,
                  paste("xps", .Platform$dynlib.ext, sep=''), package="xps");
   xpsso <- gsub("//", "/", xpsso);

   ## need to put character variables in parenthesis
   xpsso      <- paste(dq, xpsso, dq, sep="");
   rootfile   <- paste(dq, rootfile, dq, sep="");
   canvasname <- paste(dq, canvasname, dq, sep="");
   treename   <- paste(dq, treename, dq, sep="");
   varlist    <- paste(dq, varlist, dq, sep="");
   logbase    <- paste(dq, logbase, dq, sep="");
   option     <- paste(dq, option, dq, sep="");
   savename   <- paste(dq, savename, dq, sep="");
   w          <- as.integer(w);
   h          <- as.integer(h);

   ## call ROOT macro
   macro <- paste(system.file(package="xps"), "rootsrc/macroDrawDensity.C", sep="/");
   macro <- paste(sq, macro, "(", xpsso,",", rootfile,",", canvasname,",",
                  treename,",", varlist,",", logbase,",", option,",", savename,",", 
                  w,",", h, ")", sq, sep="");

   if (save.as == "") {
      system(paste("root -l", macro));
   } else {
      system(paste("root -l -q", macro));
   }#if
}#root.density

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.image" <-
function(x,
         treename   = character(0),
         leafname   = "fInten",
         logbase    = "log2",
         option     = "COLZ",
         zlim       = NULL,
         canvasname = "Image",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.image------")

   ## check for correct class
   if (!(is(x, "DataTreeSet") | is(x, "QualTreeSet"))) {
      stop(paste(sQuote("x"), "must be one of class <DataTreeSet, QualTreeSet>"));
   }#if

   rootfile <- rootFile(x);
   varlist  <- paste("fX", "fY", leafname, sep=":");
   savename <- rootDrawName(canvasname, save.as);

   ## test for valid treename
   treenames <- getTreeNames(rootfile);
   if (treename%in%treenames == FALSE) {
      stop(paste("treename", sQuote(treename), "is not present in", sQuote(rootfile)));
   }#if

### to do ###
# test for extenPart suitable for images, ev validImageExtension()
#############

   ## need to add ROOT directory in rootfile
   treename <- paste(setName(x), treename, sep="/");

   ## test for zlim
   if (is.null(zlim)) {
      minz <- 0.0;
      maxz <- 0.0;
   } else {
      minz <- zlim[1];
      maxz <- zlim[2];
   }#if

   ## unix vs windows settings
   is.win <- (.Platform$OS.type == "windows");
   if (is.win) {
      sq  <- as.character("");  #
      dq  <- as.character("\\\""); #\"
   } else {
      sq  <- as.character("'");  #'
      dq  <- as.character("\""); #"
   }#if

   ## get path to library xps.so
   xpsso <- system.file("libs", .Platform$r_arch,
                  paste("xps", .Platform$dynlib.ext, sep=''), package="xps");
   xpsso <- gsub("//", "/", xpsso);

   ## need to put character variables in parenthesis
   xpsso      <- paste(dq, xpsso, dq, sep="");
   rootfile   <- paste(dq, rootfile, dq, sep="");
   canvasname <- paste(dq, canvasname, dq, sep="");
   treename   <- paste(dq, treename, dq, sep="");
   varlist    <- paste(dq, varlist, dq, sep="");
   logbase    <- paste(dq, logbase, dq, sep="");
   option     <- paste(dq, option, dq, sep="");
   savename   <- paste(dq, savename, dq, sep="");
   minz       <- as.double(minz);
   maxz       <- as.double(maxz);
   w          <- as.integer(w);
   h          <- as.integer(h);

   ## call ROOT macro
   macro <- paste(system.file(package="xps"), "rootsrc/macroDrawImage.C", sep="/");
   macro <- paste(sq, macro, "(", xpsso,",", rootfile,",", canvasname,",",
                  treename,",", varlist,",", logbase,",", option,",",
                  savename,",", minz,",", maxz,",", w,",", h, ")", sq, sep="");

   if (save.as == "") {
      system(paste("root -l", macro));
   } else {
      system(paste("root -l -q", macro));
   }#if
}#root.image

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.profile" <-
function(x,
         treename    = "*",
         varlist     = NULL,
         as.log      = TRUE,
         globalscale = TRUE,
         boxes       = TRUE,
         ylim        = NULL,
         canvasname  = "ProfilePlot",
         save.as     = "",
         w           = 800,
         h           = 600) 
{
   if (debug.xps()) print("------root.profile------")

   rootfile  <- rootFile(x);
   treenames <- treeNames(x);
   savename  <- rootDrawName(canvasname, save.as);

   ## get tree extension
   exten <- unlist(treenames)[1];
   exten <- unlist(strsplit(exten, "\\."))[2];

   ## test for valid treename
   namepart <- unlist(sapply(treename, function(x)strsplit(x, "\\.")[[1]][1]));
   hastrees <- unique(treename%in%treenames);
   if (namepart[1] == "*") {
      namepart <- namepart[1];
   } else if (length(hastrees) == 1 & hastrees == TRUE) {
      namepart <- paste(namepart, collapse=":");
   } else {
      stop(paste("at least one treename is not present in", sQuote(rootfile)));
   }#if

   ## select leafname as varlist
   if (is.null(varlist)) {
      if (is (x, "DataTreeSet")) {
         varlist <- "fInten";
      } else if (is (x, "ExprTreeSet")) {
         varlist <- "fLevel";
      } else {
         stop(paste("for <varlist=NULL>", sQuote("x"), 
                    "must be  class <DataTreeSet,ExprTreeSet>"));
      }#if
   }#if

   ## test for ylim
   if (is.null(ylim)) {
      miny <- 0.0;
      maxy <- 0.0;
   } else {
      miny <- ylim[1];
      maxy <- ylim[2];
   }#if

   ## unix vs windows settings
   is.win <- (.Platform$OS.type == "windows");
   if (is.win) {
      sq  <- as.character("");  #
      dq  <- as.character("\\\""); #\"
   } else {
      sq  <- as.character("'");  #'
      dq  <- as.character("\""); #"
   }#if

   ## get path to library xps.so
   xpsso <- system.file("libs", .Platform$r_arch,
                  paste("xps", .Platform$dynlib.ext, sep=''), package="xps");
   xpsso <- gsub("//", "/", xpsso);

   ## need to put character variables in parenthesis
   xpsso       <- paste(dq, xpsso, dq, sep="");
   rootfile    <- paste(dq, rootfile, dq, sep="");
   canvasname  <- paste(dq, canvasname, dq, sep="");
   setname     <- paste(dq, setName(x), dq, sep="");
   treename    <- paste(dq, namepart, dq, sep="");
   exten       <- paste(dq, exten, dq, sep="");
   varlist     <- paste(dq, varlist, dq, sep="");
   savename    <- paste(dq, savename, dq, sep="");
   miny        <- as.double(miny);
   maxy        <- as.double(maxy);
   aslog       <- as.integer(as.log);
   globalscale <- as.integer(globalscale);
   boxes       <- as.integer(boxes);
   w           <- as.integer(w);
   h           <- as.integer(h);

   ## call ROOT macro
   macro <- paste(system.file(package="xps"), "rootsrc/macroDrawProfilePlot.C", sep="/");
   macro <- paste(sq, macro, "(", xpsso,",", rootfile,",", canvasname,",", setname,",",
                  treename,",", exten,",", varlist,",", savename,",", miny,",", maxy,",",
                  aslog,",", globalscale,",", boxes,",", w,",", h, ")", sq, sep="");

   if (save.as == "") {
      system(paste("root -l", macro));
   } else {
      system(paste("root -l -q", macro));
   }#if
}#root.profile

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.draw" <-
function(x,
         treename   = "",
         leafname   = "fInten",
         logbase    = "log2",
         type       = "graph",
         option     = "L",
         canvasname = "DrawCanvas",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.draw------")

   ## check for correct class
   if (!extends(class(x), "ProcesSet")) {
      stop(paste(sQuote("x"), "is not derived from class", sQuote("ProcesSet")));
   }#if

   rootfile <- rootFile(x);
   varlist  <- leafname;
   savename <- rootDrawName(canvasname, save.as);

   ## test for valid treenames
   ntree <- length(treename);
   if (ntree > 3) {
      stop(paste("function", dQuote("root.draw"), "can only be applied to at most 3 trees."));
   }#if
   treenames <- getTreeNames(rootfile);
   validtree <- treename%in%treenames;
   if (length(validtree[validtree==TRUE]) < ntree) {
      stop(paste("at least one tree", dQuote("treename"),
                 "is missing or not present in", sQuote(rootfile)));
   }#if

   ## need to add ROOT directory in rootfile
   if (is (x, "DataTreeSet")) {
      treename <- paste("DataSet", treename, sep="/");
   } else if (is (x, "ExprTreeSet")) {
      treename <- paste("PreprocesSet", treename, sep="/");
   } else if (is (x, "CallTreeSet")) {
      treename <- paste("CallSet", treename, sep="/");
   }#if

   ## unix vs windows settings
   is.win <- (.Platform$OS.type == "windows");
   if (is.win) {
      sq  <- as.character("");  #
      dq  <- as.character("\\\""); #\"
   } else {
      sq  <- as.character("'");  #'
      dq  <- as.character("\""); #"
   }#if

   ## get path to library xps.so
   xpsso <- system.file("libs", .Platform$r_arch,
                  paste("xps", .Platform$dynlib.ext, sep=''), package="xps");
   xpsso <- gsub("//", "/", xpsso);

   ## need to put character variables in parenthesis
   xpsso      <- paste(dq, xpsso, dq, sep="");
   rootfile   <- paste(dq, rootfile, dq, sep="");
   canvasname <- paste(dq, canvasname, dq, sep="");
   ntree      <- as.integer(ntree);
   treename1  <- paste(dq, treename[1], dq, sep="");
   treename2  <- paste(dq, treename[2], dq, sep="");
   treename3  <- paste(dq, treename[3], dq, sep="");
   varlist    <- paste(dq, varlist, dq, sep="");
   logbase    <- paste(dq, logbase, dq, sep="");
   type       <- paste(dq, type, dq, sep="");
   option     <- paste(dq, option, dq, sep="");
   savename   <- paste(dq, savename, dq, sep="");
   w          <- as.integer(w);
   h          <- as.integer(h);

   ## call ROOT macro
   macro <- paste(system.file(package="xps"), "rootsrc/macroDraw.C", sep="/");
   macro <- paste(sq, macro, "(", xpsso,",", rootfile,",", canvasname,",",
                  ntree,",", treename1,",", treename2,",", treename3,",",
                  varlist,",", logbase,",", type,",", option,",", savename,",", 
                  w,",", h, ")", sq, sep="");

   if (save.as == "") {
      system(paste("root -l", macro));
   } else {
      system(paste("root -l -q", macro));
   }#if
}#root.draw

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.draw.leaves" <-
function(x,
         treename   = "*",
         leafname   = "fInten",
         logbase    = "log2",
         type       = "graph",
         option     = "L",
         sortopt    = 0,
         canvasname = "DrawCanvas",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.draw.leaves------")

   ## check for correct class
   if (!extends(class(x), "ProcesSet")) {
      stop(paste(sQuote("x"), "is not derived from class", sQuote("ProcesSet")));
   }#if

   rootfile  <- rootFile(x);
   treenames <- treeNames(x);
   varlist   <- leafname;
   savename  <- rootDrawName(canvasname, save.as);

   ## get tree extension
   exten <- unlist(treenames)[1];
   exten <- unlist(strsplit(exten, "\\."))[2];

   ## test for valid treename
   namepart <- unlist(strsplit(treename, "\\."))[1];
   if (namepart == "*") {
      treename <- paste(namepart, exten, sep=".");
   } else if (treename%in%treenames == FALSE) {
      stop(paste("treename", sQuote(treename), "is not present in", sQuote(rootfile)));
   }#if

   ## need to add ROOT directory in rootfile
   if (is (x, "DataTreeSet")) {
      treename <- paste("DataSet", treename, sep="/");
   } else if (is (x, "ExprTreeSet")) {
      treename <- paste("PreprocesSet", treename, sep="/");
   } else if (is (x, "CallTreeSet")) {
      treename <- paste("CallSet", treename, sep="/");
   }#if

   ## unix vs windows settings
   is.win <- (.Platform$OS.type == "windows");
   if (is.win) {
      sq  <- as.character("");  #
      dq  <- as.character("\\\""); #\"
   } else {
      sq  <- as.character("'");  #'
      dq  <- as.character("\""); #"
   }#if

   ## get path to library xps.so
   xpsso <- system.file("libs", .Platform$r_arch,
                  paste("xps", .Platform$dynlib.ext, sep=''), package="xps");
   xpsso <- gsub("//", "/", xpsso);

   ## need to put character variables in parenthesis
   xpsso      <- paste(dq, xpsso, dq, sep="");
   rootfile   <- paste(dq, rootfile, dq, sep="");
   canvasname <- paste(dq, canvasname, dq, sep="");
   treename   <- paste(dq, treename, dq, sep="");
   varlist    <- paste(dq, varlist, dq, sep="");
   logbase    <- paste(dq, logbase, dq, sep="");
   type       <- paste(dq, type, dq, sep="");
   option     <- paste(dq, option, dq, sep="");
   savename   <- paste(dq, savename, dq, sep="");
   sortopt    <- as.integer(sortopt);
   w          <- as.integer(w);
   h          <- as.integer(h);

   ## call ROOT macro
   macro <- paste(system.file(package="xps"), "rootsrc/macroDrawLeaves.C", sep="/");
   macro <- paste(sq, macro, "(", xpsso,",", rootfile,",", canvasname,",",
                  treename,",", varlist,",", logbase,",", type,",", option,",",
                  savename,",", sortopt,",", w,",", h, ")", sq, sep="");

   if (save.as == "") {
      system(paste("root -l", macro));
   } else {
      system(paste("root -l -q", macro));
   }#if
}#root.draw.leaves

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.graph1D" <-
function(x,
         treename   = character(0),
         logbase    = "log2",
         option     = "P",
         canvasname = "Graph1D",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.graph1D------")

   ## select leafname for class
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else if (is (x, "CallTreeSet")) {
      leafname <- "fPValue";
   }#if

   root.draw(x,
             treename   = treename,
             leafname   = leafname,
             logbase    = logbase,
             type       = "graph",
             option     = option,
             canvasname = canvasname,
             save.as    = save.as,
             w          = w,
             h          = h);
}#root.graph1D

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.graph2D" <-
function(x,
         treename1  = character(0),
         treename2  = character(0),
         logbase    = "log2",
         option     = "P",
         canvasname = "Graph2D",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.graph2D------")

   ## select leafname for class
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else if (is (x, "CallTreeSet")) {
      leafname <- "fPValue";
   }#if

   root.draw(x,
             treename   = c(treename1, treename2),
             leafname   = leafname,
             logbase    = logbase,
             type       = "graph",
             option     = option,
             canvasname = canvasname,
             save.as    = save.as,
             w          = w,
             h          = h);
}#root.graph2D

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.mvaplot" <-
function(x,
         treename1  = character(0),
         treename2  = character(0),
         logbase    = "log2",
         option     = "P",
         canvasname = "MvAPlot",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.mvaplot------")

   ## select leafname for class
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else if (is (x, "CallTreeSet")) {
#      leafname <- "fPValue";
      stop(paste(sQuote("mvaplot"), "does not make sense for class", sQuote("CallTreeSet")));
   }#if

   root.draw(x,
             treename   = c(treename1, treename2),
             leafname   = leafname,
             logbase    = logbase,
             type       = "mvaplot",
             option     = option,
             canvasname = canvasname,
             save.as    = save.as,
             w          = w,
             h          = h);
}#root.mvaplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.hist1D" <-
function(x,
         treename   = character(0),
         logbase    = "log2",
         type       = "hist",
         option     = "HIST",
         canvasname = "Histogram1D",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.hist1D------")

   if (!(identical(type, "hist") || identical(type, "density"))) {
      stop(paste(sQuote(type), "is not a valid type option"));
   }#if

   ## density requires special option
   if (identical(type, "density")) {
      option = "AL";
   }#if

   ## select leafname for class
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else if (is (x, "CallTreeSet")) {
      leafname <- "fPValue";
   }#if

   root.draw(x,
             treename   = treename,
             leafname   = leafname,
             logbase    = logbase,
             type       = type,
             option     = option,
             canvasname = canvasname,
             save.as    = save.as,
             w          = w,
             h          = h);
}#root.hist1D

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.hist2D" <-
function(x,
         treename1  = character(0),
         treename2  = character(0),
         logbase    = "log2",
         option     = "COLZ",
         canvasname = "Histogram2D",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.hist2D------")

   ## select leafname for class
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else if (is (x, "CallTreeSet")) {
      leafname <- "fPValue";
   }#if

   root.draw(x,
             treename   = c(treename1, treename2),
             leafname   = leafname,
             logbase    = logbase,
             type       = "hist",
             option     = option,
             canvasname = canvasname,
             save.as    = save.as,
             w          = w,
             h          = h);
}#root.hist2D

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"root.hist3D" <-
function(x,
         treename1  = character(0),
         treename2  = character(0),
         treename3  = character(0),
         logbase    = "log2",
         option     = "HIST",
         canvasname = "Histogram3D",
         save.as    = "",
         w          = 540,
         h          = 540) 
{
   if (debug.xps()) print("------root.hist3D------")

   ## select leafname for class
   if (is (x, "DataTreeSet")) {
      leafname <- "fInten";
   } else if (is (x, "ExprTreeSet")) {
      leafname <- "fLevel";
   } else if (is (x, "CallTreeSet")) {
      leafname <- "fPValue";
   }#if

   root.draw(x,
             treename   = c(treename1, treename2, treename3),
             leafname   = leafname,
             logbase    = logbase,
             type       = "hist",
             option     = option,
             canvasname = canvasname,
             save.as    = save.as,
             w          = w,
             h          = h);
}#root.hist3D

#==============================================================================#
