#==============================================================================#
# methods.ProcesSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# schemeFile:
# schemeFile<-:
# schemeSet:
# schemeSet<-:
# chipName:
# chipType:
# attachData: 
# removeData: 
# getTreeData:
# treeData:
# validData:
# export:
# boxplot:
# mboxplot:
# hist:
# image:
#==============================================================================#


#------------------------------------------------------------------------------#
# ProcesSet initialization:
#------------------------------------------------------------------------------#

setMethod("initialize", "ProcesSet", 
   function(.Object,
            scheme = NULL,
            data   = data.frame(),
            params = list(),
            ...) 
   {
      if (debug.xps()) print("------initialize:ProcesSet------")

      ## prevent effects of multiple initialization
      if (is.null(scheme))  return(.Object);
      if (is.null(data))    data   <- data.frame(matrix(nrow=0,ncol=0));
      if (!is.list(params)) params <- as.list(params);

      .Object@scheme = scheme;
      .Object@data   = data;
      .Object@params = params;

      .Object <- callNextMethod(.Object,
                                scheme = scheme,
                                data   = data,
                                params = params,
                                ...);
      .Object@scheme = scheme;
      .Object@data   = data;
      .Object@params = params;
      .Object;
   }
)#initialize

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("ProcesSet",
   function(object) 
   {
      if (debug.xps()) print("------setValidity:ProcesSet------")

      msg <- NULL;

      ## check setname
      if (!(is(object@setname, "character") && nchar(object@setname) > 0)) {
         msg <- validMsg(msg, paste(sQuote("setname"), "is missing"));
      }#if

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# ProcesSet accessors:
#------------------------------------------------------------------------------#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("schemeFile", signature(object="ProcesSet"),
   function(object) {
      scheme <- object@scheme;
      return(scheme@rootfile);
   }
)#schemeFile

setReplaceMethod("schemeFile", signature(object="ProcesSet", value="character"),
   function(object, value) {
      scheme <- object@scheme;
      scheme@rootfile <- value;
      object@scheme <- scheme;
      return(object);
   }
)#schemeFile<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("schemeSet", signature(object="ProcesSet"),
   function(object) object@scheme
)#schemeSet

setReplaceMethod("schemeSet", signature(object="ProcesSet", value="SchemeTreeSet"),
   function(object, value) {
      object@scheme <- value;
      return(object);
   }
)#schemeSet<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("chipName", signature(object="ProcesSet"),
   function(object) {
      scheme <- object@scheme;
      return(scheme@chipname);
   }
)#chipName

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("chipType", signature(object="ProcesSet"),
   function(object) {
      scheme <- object@scheme;
      return(scheme@chiptype);
   }
)#chipType


#------------------------------------------------------------------------------#
# ProcesSet methods:
#------------------------------------------------------------------------------#

setMethod("attachData", signature(object="ProcesSet"),
   function(object,
            treenames = character(0),
            varlist   = character(0),
            outfile   = "data.txt")
   {
      if (debug.xps()) print("------attachData.ProcesSet------")

      treetype <- extenPart(object@treenames);
      if (treenames[1] == "*") treenames <- object@treenames;
      if (length(treenames) > 0) {
         object@data <- export(object,
                               treenames    = treenames,
                               treetype     = treetype,
                               varlist      = varlist,
                               outfile      = outfile,
                               as.dataframe = TRUE,
                               verbose      = FALSE);
      } else {
         warning("missing tree names, data will not be added.");
      }#if
      return(object);
   }
)#attachData

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("removeData", signature(object="ProcesSet"),
   function(object) {
      if (debug.xps()) print("------removeData.ProcesSet------")

      object@data <- data.frame(matrix(nrow=0,ncol=0));
      gc(); #????
      return(object);
   }
)#removeData

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("getTreeData", signature(object="ProcesSet"),
   function(object,
            treetype = "cel",
            varlist  = "fInten")
   {
      if (debug.xps()) print("------getTreeData.ProcesSet------")

      ds <- export(object,
                   treetype     = treetype,
                   varlist      = varlist,
                   as.dataframe = TRUE);
      return(ds);
   }
)#getTreeData

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("treeData", signature(object="ProcesSet"),
   function(object) object@data
)#treeData

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"dataProcesSet" <-
function(object,
         which = "UnitName") 
{
   if (debug.xps()) print("------dataProcesSet------")

   ## check for presence of data
   data  <- object@data;
   if (min(dim(data)) == 0) {
      stop(paste("slot", sQuote("data"), "has no data"));
   }#if

   ## get optional subset of columns
   strg <- unlist(strsplit(which,":"));
   if (length(strg) == 2) {
      id   <- grep(strg[2], colnames(data));
      data <- data[,c(strg[1], colnames(data)[id])];
   }#if

   ## use names from column "which" as rownames
   if (!is.na(match(which, colnames(data)))) {
      len <- length(which(duplicated(data[,which])==TRUE));
      if (len == 0) {
         rownames(data) <- data[, which];
      } else {
         warning(paste("cannot use ", sQuote(which), "as row.names since it has <",
                       len, "> non-unique values.", sep=""));

         ## use names from column "UNIT_ID" as rownames
         if (!is.na(match("UNIT_ID", colnames(data)))) {
            rownames(data) <- data[, "UNIT_ID"];
         }#if
      }#if
   }#if

   ## get name part for matching columns
   treenames <- namePart(object@treenames);
   datanames <- namePart(colnames(data));

   return(data[,!is.na(match(datanames, treenames))]);
}#dataProcesSet

setMethod("validData", "ProcesSet", dataProcesSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"exportProcesSet" <-
function(object,
         treenames    = "*",
         treetype     = character(0),
         varlist      = "*",
         outfile      = character(0),
         sep          = "\t",
         as.dataframe = FALSE,
         verbose      = TRUE,
         ...) 
{
   if (debug.xps()) print("------export.ProcesSet------")

   callNextMethod();

   ## get name of root file
   filename <- rootFile(object);

   ## get name of root scheme file
   schemefile <- validROOTFile(schemeFile(object), class(object));

   ## get chiptype
   chiptype <- chipType(object);

   ## get treeset
   treeset <- setName(object);

   ## check for treenames
   numtrees <- length(treenames);
   if (treenames[1] != "*") {
      treenames <- namePart(treenames);
   }#if

   ## check for presence of correct tree type
   datatype <- setType(object);
   treetype <- validTreetype(treetype, datatype);

   ## check for varlist
   if (nchar(varlist) < 1) {
      varlist <- "*";
   }#if

   ## check for presence of valid separator and outfile
   exten <- validSeparator(sep);
   if (length(outfile) && outfile != "") {
      if (regexpr(".txt",outfile) > 0 || regexpr(".csv",outfile) > 0) {
         outfile <- substr(outfile, 1, nchar(outfile)-4);
      }#if
      outfile <- paste(outfile, exten, sep=".");
   }#if

   ## check for presence of outfile=/path/outname
   outname <- "";
   if (numtrees > 1 || treenames[1] == "*") {
      outname <- "PivotTable";
   } else {
      outname <- treenames;
   }#if
   outfile <- validOutfile(outname, outfile);

   ## get tree names "treeset.treename.treetype"
   if (numtrees == 1) {
      treenames <- paste(treeset, treenames, treetype, sep=".");
   } else if (numtrees > 1) {
      treenames <- paste(treenames, treetype, sep=".");
#      treenames <- paste(filename, treeset, treenames, sep="/");
      treenames <- paste(treeset, treenames, sep="/");
   }#if

   ## export data from root file
   r <- .C("ExportData",
           as.character(filename),
           as.character(schemefile),
           as.character(chiptype),
           as.character(datatype),
           as.character(treenames),
           as.integer(numtrees),
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
      ds <- read.table(outfile, header=TRUE, row.names=NULL, sep=sep, 
                       check.names=FALSE, stringsAsFactors=FALSE, comment.char="", ...);
   }#if

   return(ds);
}#exportProcesSet

setMethod("export", "ProcesSet", exportProcesSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("boxplot", signature(x="ProcesSet"),
   function(x,
            which   = "",
            size    = 0,
            transfo = log2,
            range   = 0,
            names   = "namepart",
            bmar    = NULL,
            las     = 2,
            ...) 
   {
      if (debug.xps()) print("------boxplot.ProcesSet------")

      has.userinfo <- (unlist(strsplit(which, ":"))[1] == "userinfo");
      has.userinfo <- ifelse(is.na(has.userinfo), FALSE, has.userinfo);

      if (has.userinfo) {
         ds <- treeInfo(x,
                        treetype = extenPart(treeNames(x)),
                        varlist  = sub("userinfo:", "", which),
                        verbose  = FALSE
                       );
      } else {
         ds <- validData(x, which=which);
         if (size > 1) ds <- ds[seq(1,nrow(ds),len=size),];
      }#if

      if (is.function(transfo)) ds <- transfo(ds);

      if (is.null(names)) {
         names <- colnames(ds);
      } else if (validQualityOption(names, as.logical = TRUE)) {
         ds <- ds[, grep(names, colnames(ds))];
         names <- colnames(ds);
      } else if (unlist(strsplit(names,":"))[1] == "namepart") {
         qualopt <- unlist(strsplit(names,":"))[2];
         if (validQualityOption(qualopt, as.logical = TRUE)) {
            ds <- ds[, grep(qualopt, colnames(ds))];
         }#if
         names <- namePart(colnames(ds));
      } else {
         ds    <- ds[, names, drop=FALSE];
      }#if

      if (is.null(bmar)) bmar <- adjustXLabMargins(names, bottom=6, cex=1.0);

      oldpar <- par(mar=c(bmar$b,5,2,1), pty="m", cex.axis=bmar$cex);

      if (has.userinfo & range == 0) {
         bx.p <- list(stats = as.matrix(ds[c(1,3:5,7),]),
                      names = names,
                      out   = numeric(0),
                      group = numeric(0)
                     );
      } else if (has.userinfo & range > 0) {
         bx.p <- list(stats = as.matrix(ds[2:6,]),
                      names = names,
                      out   = as.vector(as.matrix(ds[c(1,7),])),
                      group = rep(1:ncol(ds), each = 2)
                     );
      } else {
         bx.p <- boxplot(ds, plot = FALSE, names = names, range = range);
      }#if

      bxp(bx.p,
          names = names,
          las    = las,
          ...);

      par(oldpar);
   }
)#boxplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("mboxplot", signature(x="ProcesSet"),
   function(x,
            which    = "",
            size     = 0,
            transfo  = log2,
            method   = "mean",
            range    = 0,
            names    = "namepart",
            main     = "RLE Plot",
            ylim     = c(-1,1),
            las      = 2,
            add.line = TRUE,
            outline  = FALSE,
            ...) 
   {
      if (debug.xps()) print("------mboxplot.ProcesSet------")

      m <- validData(x, which=which);
      if (size > 1)             m <- m[seq(1,nrow(m),len=size),];
      if (is.function(transfo)) m <- transfo(m);

      if (is.null(names))              names <- colnames(m)
      else if (names[1] == "namepart") names <- namePart(colnames(m))
      else                             m     <- m[, names, drop=FALSE];

      if (method == "median") {
         mn  <- apply(m, 1, median);
      } else {
         mn  <- rowMeans(m);
      }#if
      rle <- data.frame(sweep(m, MARGIN=1, STATS=mn, FUN='-'));

      boxplot(rle,
              range   = range,
              main    = main,
              xlab    = '',
              xaxt    = 'n',
              ylim    = ylim,
              outline = outline);

      if (add.line) {
         abline(0, 0, lty=2, col="gray70");
      }#if
      axis(side=1, outer=F, at=1:length(names), labels=names, las=las, ...);
   }
)#mboxplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("hist", signature(x="ProcesSet"),
   function(x,
            which      = "",
            size       = 0,
            transfo    = log2,
            xlab       = "log intensity",
            ylab       = "density",
            names      = "namepart",
            type       = "l",
            col        = 1:6,
            lty        = 1:5,
            add.legend = FALSE,
            verbose    = TRUE,
            ...) 
   {
      if (debug.xps()) print("------hist.ProcesSet------")

      if (verbose) {verbose=FALSE;} ## for compatibility to hist.DataTreeSet

      m <- validData(x, which=which);
      if (size > 1) m <- m[seq(1,nrow(m),len=size),];
      if (is.function(transfo)) m <- transfo(m);

      treenames <- colnames(m);
      if (is.null(names)) {
         names <- treenames;
      } else if (names[1] == "namepart") {
         names <- namePart(treenames);
      } else {
         treenames <- namePart(treeNames(x));

         ok <- match(namePart(names), treenames);
         ok <- ok[!is.na(ok)];
         if (length(ok) == 0) {
            stop(paste(sQuote("names"), "is not a valid tree name"));
         }#if

         m     <- m[,ok];
         names <- treenames[ok];
      }#if

      plotDensities(m,
                    ylab = ylab,
                    xlab = xlab,
                    type = type,
                    lty  = lty,
                    col  = col,
                    ...);

      if (add.legend) {
         legend(x      = "topright",
                legend = names[1:ncol(m)],
                lty    = lty,
                pt.bg  = "white",
                col    = col,
                cex    = 0.6
               );
      }#if
   }
)#hist

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("image", signature(x="ProcesSet"),
   function(x,
            bg         = FALSE,
            transfo    = log,
            col        = NULL,
            names      = "namepart",
            xlab       = "",
            ylab       = "",
            add.legend = FALSE,
            ...) 
   {
      if (debug.xps()) print("------image.ProcesSet------")

      treenames <- getTreeNames(rootFile(x));
      treetype  <- extenPart(treenames);

      if (bg == FALSE) {
         varlist <- "fInten";
         colname <- "MEAN";
         treetype <- treetype[!is.na(sapply(treetype, function(x) match(x, c(RAWTYPE, ADJTYPE, CNRTYPE))))];
         treenames <- treenames[grep(unique(treetype), treenames)];
      } else {
         varlist  <- "fBg";
         colname  <- "BGRD";
         treetype <- treetype[!is.na(sapply(treetype, function(x) match(x, BGDTYPE)))];
         treenames <- treenames[grep(unique(treetype), treenames)];
      }#if

      if (is.null(names)) {
         names <- treenames;
      } else if (names[1] == "namepart") {
         names <- namePart(treenames);
      } else {
         ok <- match(names, treenames);
         ok <- ok[!is.na(ok)];

         if (length(ok) == 0) {
            stop(paste(sQuote("names"), "is not a valid data tree name"));
         }#if

         names <- unlist(treenames[ok]);
      }#if

      ## colors
      if (is.null(col)) col = gray((0:64)/64);

      ntrees <- length(names);
      nrows  <- nrows(schemeSet(x));
      ncols  <- ncols(schemeSet(x));

      ## plot images
      if (ntrees > 1 && interactive()) par(ask=TRUE) else par(ask=FALSE);

      for (i in 1:ntrees) {
         ds <- export(x,
                      treenames    = names[i],
                      treetype     = unique(treetype),
                      varlist      = varlist,
                      as.dataframe = TRUE,
                      verbose      = FALSE);
         m  <- matrix(ds[,colname], nrow=ncols, ncol=nrows, byrow=FALSE);
         m  <- m[,nrows:1];
         if (is.function(transfo)) m <- transfo(m);

         par(mar = c(1, 1, 2, 1));
         if (add.legend) {
            layout(matrix(c(1, 2), 1, 2, byrow=TRUE), widths=c(7,1), heights=8, TRUE);
            par(mar = c(1, 1, 2, 0));
         }#if

         graphics::image(as.matrix(m),
                         col  = col,
                         main = names[i],
                         xlab = xlab,
                         ylab = ylab,
                         xaxt = 'n',
                         yaxt = 'n',
                         ...);
         box();

         if (add.legend) {
            if (is.function(transfo)) {
               lo <- 0.1;
               hi <- transfo(max(ds[,colname], na.rm=TRUE));
            } else {
               lo <- 0.0;
               hi <- max(ds[,colname], na.rm=TRUE);
            }#if
            y <- pretty(c(lo, hi), 10);
            m <- matrix(y, nrow=1, ncol=length(y));

            par(mar = c(1, 1, 2, 3));
            graphics::image(m, xaxt="n", yaxt="n", col=col);
            axis(4, labels=y, at=seq(0, 1, by=(1/(length(y)-1))), las=2, cex.axis=0.8);
            layout(1);
            par(mar = c(1, 1, 2, 1));
         }#if
      }#for

      par(ask=FALSE);
   }
)#image

#------------------------------------------------------------------------------#
