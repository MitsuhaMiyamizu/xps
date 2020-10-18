#==============================================================================#
# methods.QualTreeSet.R: initialization, accessors, methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# initialize:
# setValidity:
# qualType:
# qualType<-:
# qualOption:
# qualOption<-:
# borders:
# residuals:
# weights:
# xpsRNAdeg:
# image:
# borderplot:
# coiplot:
# nuseplot:
# rleplot:
# NUSE:
# RLE:
#==============================================================================#


#------------------------------------------------------------------------------#
# QualTreeSet initialization:
#------------------------------------------------------------------------------#

setMethod("initialize", "QualTreeSet", 
   function(.Object,
            qualtype = "rlm",
            qualopt  = "raw",
            ...) 
   {
      if (debug.xps()) print("------initialize:QualTreeSet------")

      ## set default quality type and option
      if (qualtype == "") qualtype <- "rlm";
      if (qualopt  == "") qualopt  <- "raw";

      .Object <- callNextMethod(.Object,
                                qualtype = qualtype,
                                qualopt  = qualopt,
                                ...);
      .Object@qualtype = qualtype;
      .Object@qualopt  = qualopt;
      .Object;
   }
)#initialize

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setValidity("QualTreeSet",
   function(object) {
      if (debug.xps()) print("------setValidity:QualTreeSet------")

      msg <- NULL;

      ## check for correct settype
      if (!(object@settype == "preprocess")) {
         msg <- validMsg(msg, paste(sQuote("settype"),
                        "must be of type <preprocess>"));
      }#if

      ## check quality method
      TYPE <- c("rlm", "plm");
      if (is.na(match(object@qualtype, TYPE))) {
         msg <- validMsg(msg, paste(sQuote("qualtype"), "must be one of",
                         "<rlm,plm>"));
      }#if

      ## check quality option
      TYPE <- c("raw", "adjusted", "normalized", "all");
      if (is.na(match(object@qualopt, TYPE))) {
         msg <- validMsg(msg, paste(sQuote("qualopt"), "must be one of",
                         "<raw,adjusted,normalized,all>"));
      }#if

      if (is.null(msg)) TRUE else msg;
   }
)#setValidity


#------------------------------------------------------------------------------#
# QualTreeSet accessors:
#------------------------------------------------------------------------------#

setMethod("qualType", signature(object="QualTreeSet"),
   function(object) object@qualtype
)#qualType

setReplaceMethod("qualType", signature(object="QualTreeSet", value="character"),
   function(object, value) {
      object@qualtype <- value;
      return(object);
   }
)#qualType<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("qualOption", signature(object="QualTreeSet"),
   function(object) object@qualopt
)#qualOption

setReplaceMethod("qualOption", signature(object="QualTreeSet", value="character"),
   function(object, value) {
      object@qualopt <- value;
      return(object);
   }
)#qualOption<-

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("borders", signature(object="QualTreeSet"),
   function(object, treename = "*") {
      if (debug.xps()) print("------borders.QualTreeSet------")

      brd <- export(object,
                    treename     = treename,
                    treetype     = QUATYPE[3],  # "brd"
                    varlist      = "fFlag:fInten",
                    as.dataframe = TRUE,
                    verbose      = FALSE);
      return(brd);
   }
)#borders

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("residuals", signature(object="QualTreeSet"),
   function(object, treename = "*") {
      if (debug.xps()) print("------residuals.QualTreeSet------")

      res <- export(object,
                    treename     = treename,
                    treetype     = QUATYPE[4],  # "res"
                    varlist      = "fResidual",
                    as.dataframe = TRUE,
                    verbose      = FALSE);
      return(res[order(res[,"Y"], res[,"X"]),]);
   }
)#residuals

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("weights", signature(object="QualTreeSet"),
   function(object, treename = "*") {
      if (debug.xps()) print("------weights.QualTreeSet------")

      w <- export(object,
                  treename     = treename,
                  treetype     = QUATYPE[4],  # "res"
                  varlist      = "fWeight",
                  as.dataframe = TRUE,
                  verbose      = FALSE);

      w[w == eINITWEIGHT] <- NA;
      w <- w[order(w[,"Y"], w[,"X"]),];

      return(w);
   }
)#weights


#------------------------------------------------------------------------------#
# QualTreeSet methods:
#------------------------------------------------------------------------------#

"RNAdeg.QualTreeSet" <-
function(object,
         treename = "*",
         qualopt  = NULL)
{
   if (debug.xps()) print("------RNAdeg.QualTreeSet------")

   ds <- treeInfo(object,
                  treename = treename,
                  treetype = QUATYPE[2],  #rlm
                  varlist  = "userinfo:fNDegUnits:fNCells:fMNS:fSES",
                  qualopt  = qualopt
                 );

   nunits <- as.integer(ds[1,1]);
   ncells <- as.integer(ds[2,1]);

   res <- list(N            = nunits,
               sample.names = colnames(ds),
               mns          = t(as.matrix(ds[3:(ncells+2),])),
               ses          = t(as.matrix(ds[(ncells+3):nrow(ds),]))
              );

   return(res);
}

setMethod("xpsRNAdeg", "QualTreeSet", RNAdeg.QualTreeSet);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("image", signature(x="QualTreeSet"),
   function(x,
            type       = c("resids", "pos.resids", "neg.resids", "sign.resids", "weights"),
            qualopt    = c("raw", "adjusted", "normalized"),
            transfo    = log2,
            col        = NULL,
            names      = "namepart",
            xlab       = "",
            ylab       = "",
            add.legend = FALSE,
            ...) 
   {
      if (debug.xps()) print("------image.QualTreeSet------")

      type    <- match.arg(type);
      qualopt <- match.arg(qualopt);
      resname <- getTreeNames(rootFile(x), treetype=QUATYPE[4]);
      if (qualopt != "all") {
         resname <- resname[grep(qualopt, resname)];
      }#if

      if (is.null(names)) {
         names <- resname;
      } else if (names[1] == "namepart") {
         names <- namePart(resname);
      } else {
         ok <- match(names, resname);
         ok <- ok[!is.na(ok)];

         if (length(ok) == 0) {
            stop(paste(sQuote("names"), "is not a valid residual tree name"));
         }#if

         names <- resname[ok];
      }#if

      nres  <- length(names);
      nrows <- nrows(schemeSet(x));
      ncols <- ncols(schemeSet(x));

      ## get range across all arrays
      if (type != "weights") {
         rg <- treeInfo(x, treetype="res", varlist="fResiduQuant", qualopt=qualopt);
         rg <- range(rg);
      }#if

      ## colors
      if (is.null(col)) {
         if (type == "weights") {
            col <- terrain.colors(25);
         } else if (type == "resids") {
            col <- pseudoPalette(low="blue",  high="red", mid="white");
         } else if (type == "pos.resids") {
            col <- pseudoPalette(low="white", high="red");
         } else if (type == "neg.resids") {
            col <- pseudoPalette(low="blue",  high="white");
         } else if (type == "sign.resids") {
            col <- pseudoPalette(low="blue",  high="red", mid="white");
         }#if
      }#if

      ## plot images
      if (nres > 1 && interactive()) par(ask=TRUE) else par(ask=FALSE);

      for (i in 1:nres) {
         if (type == "weights") {
            ds <- weights(x, treename=names[i]);
            m  <- matrix(ds[,"WEIGHT"], nrow=ncols, ncol=nrows, byrow=FALSE);
            m  <- m[,nrows:1];
            hi <- 1.0;
            lo <- 0.0;
         } else {
            ds <- residuals(x, treename=names[i]);
            m  <- matrix(ds[,"RESIDUAL"], nrow=ncols, ncol=nrows, byrow=FALSE);
            m  <- m[,nrows:1];
         }#if

         if (type == "resids") {
            hi <- max(abs(rg));
            lo <- -max(abs(rg));
            if (is.function(transfo)) {
               m  <- sign(m) * transfo(abs(m) + 1);
               hi <- max(transfo(abs(rg) + 1));
               lo <- -max(transfo(abs(rg) + 1));
            }#if
         } else if (type == "pos.resids") {
            m  <- pmax(m, 0);
            hi <- max(rg);
            lo <- 0.0;
            if (is.function(transfo)) {
               m  <- sign(m) * transfo(abs(m) + 1);
               hi <- max(transfo(pmax(rg, 0) + 1));
               lo <- 0.0;
            }#if
         } else if (type == "neg.resids") {
            m  <- pmin(m, 0);
            hi <- 0.0;
            lo <- -abs(min(rg));
            if (is.function(transfo)) {
               m  <- sign(m) * transfo(abs(m) + 1);
               hi <- 0.0;
               lo <- -transfo(abs(min(rg)) + 1);
            }#if
         } else if (type == "sign.resids") {
            m  <- sign(m);
            hi <- 1.0;
            lo <- -1.0;
         }#if

         par(mar = c(1, 1, 2, 1));
         if (add.legend) {
            layout(matrix(c(1, 2), 1, 2, byrow=TRUE), widths=c(7,1), heights=8, TRUE);
            par(mar = c(1, 1, 2, 0));
         }#if

         graphics::image(m,
                         zlim = c(lo, hi),
                         col  = col,
                         main = names[i],
                         xlab = xlab,
                         ylab = ylab,
                         xaxt = 'n',
                         yaxt = 'n',
                         ...);
         box();

         if (add.legend) {
            par(mar = c(1, 1, 2, 3));
            y <- pretty(c(lo, hi), 10);
            m <- matrix(y, nrow=1, ncol=length(y));
            graphics::image(m, xaxt="n", yaxt="n", col=col);
            axis(4, labels=y, at=seq(0, 1, by=(1/(length(y)-1))), las=2, cex.axis=0.8);
            layout(1);
            par(mar = c(1, 1, 2, 1));
         }#if
      }#for

      par(ask=FALSE);
   }
)#image

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("borderplot", signature(x="QualTreeSet"),
   function(x,
            type    = c("pos", "neg"),
            qualopt = "raw",
            transfo = log2,
            range   = 0,
            names   = "namepart",
            ylim    = NULL,
            bmar    = NULL,
            las     = 2,
            ...) 
   {
      if (debug.xps()) print("------borderplot.QualTreeSet------")

      brd <- borders(x);
      pos <- brd[brd[,"FLAG"] ==  1, grep(qualopt, colnames(brd))];
      neg <- brd[brd[,"FLAG"] == -1, grep(qualopt, colnames(brd))];

      if (is.function(transfo)) pos <- transfo(pos);
      if (is.function(transfo)) neg <- transfo(neg);

      if (is.null(names))              names <- colnames(pos)
      else if (names[1] == "namepart") names <- namePart(colnames(pos))
      else {
         pos <- pos[, names, drop=F];
         neg <- neg[, names, drop=F];
      }#if

      if (is.null(ylim)) ylim <- c(floor(min(neg)), ceiling(max(pos)));
      if (is.null(bmar)) bmar <- adjustXLabMargins(names, bottom=6, cex=1.0);

      oldpar <- par(no.readonly = TRUE);

      if (length(type) == 2) {
         layout(matrix(c(1,2), 1, 2, byrow=FALSE), c(1.05,0.95), 2, FALSE);
         par(mar=c(bmar$b,5,2,0), pty="m", cex.axis=bmar$cex);
      } else {
         par(mar=c(bmar$b,5,2,1), pty="m", cex.axis=bmar$cex);
      }#if

      if (any(type == "pos")) {
         graphics::boxplot(pos,
                           range = range,
                           names = names,
                           ylim  = ylim,
                           las   = las,
                           ylab  = "Intensity",
                           main  = "Positive Border Elements",
                           ...);
      }#if

      if (length(type) == 2) {
         par(mar=c(bmar$b,0,2,1), yaxt="n", pty="m", cex.axis=bmar$cex);
      }#if

      if (any(type == "neg")) {
         graphics::boxplot(neg,
                           range = range,
                           names = names,
                           ylim  = ylim,
                           las   = las,
                           ylab  = "Intensity",
                           main  = "Negative Border Elements",
                           ...);
      }#if

      par(oldpar);
   }
)#borderplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("coiplot", signature(x="QualTreeSet"),
   function(x,
            type    = c("pos", "neg"),
            qualopt = "raw",
            radius  = 0.5,
            linecol = "gray70",
            visible = TRUE,
            ...) 
   {
      if (debug.xps()) print("------coiplot.QualTreeSet------")

      ds <- export(x,
                   treetype     = "brd",
                   varlist      = "userinfo:fCOIXhi:fCOIYhi:fCOIXlo:fCOIYlo",
                   as.dataframe = TRUE,
                   verbose      = FALSE);

      pos   <- t(ds[1:2,grep(qualopt, colnames(ds))]);
      neg   <- t(ds[3:4,grep(qualopt, colnames(ds))]);
      names <- namePart(rownames(pos));
      colnames(pos) <- c("X", "Y");
      colnames(neg) <- c("X", "Y");

      ds <- NULL;

      if (length(type) == 2) {
         oldpar <- par(no.readonly = TRUE);
         layout(matrix(c(1,2), 1, 2, byrow=FALSE), c(1.05,0.95), 2, FALSE);
         par(mar=c(5,5,2,0));
      }#if

      if (any(type == "pos")) {
         plot(pos,
              xlim = c(-1,1),
              ylim = c(-1,1),
              xlab = "X Center of Intensity position",
              ylab = "Y Center of Intensity position",
              main = "Positive Elements",
              ...);

         abline(h=0,    col=linecol, lty=3);
         abline(v=0,    col=linecol, lty=3);
         circle(radius, col=linecol, lty=3);

         id <- (sqrt(pos[,1]*pos[,1] + pos[,2]*pos[,2]) >= radius);
         if (visible & any(id == TRUE)) {
            text(pos[id,,drop=F], labels=names[id], cex=0.8, adj=c(1,1));

            ds <- rbind(ds, pos[id,,drop=F]);
         }#if
      }#if

      if (length(type) == 2) {
         par(mar=c(5,0,2,1), yaxt="n");
      }#if

      if (any(type == "neg")) {
         plot(neg,
              xlim = c(-1,1),
              ylim = c(-1,1),
              xlab = "X Center of Intensity position",
              ylab = "Y Center of Intensity position",
              main = "Negative Elements",
              ...);

         abline(h=0,    col=linecol, lty=3);
         abline(v=0,    col=linecol, lty=3);
         circle(radius, col=linecol, lty=3);

         id <- (sqrt(neg[,1]*neg[,1] + neg[,2]*neg[,2]) >= radius);
         if (visible & any(id == TRUE)) {
            text(neg[id,,drop=F], labels=names[id], cex=0.8, adj=c(1,1));

            ds <- rbind(ds, neg[id,,drop=F]);
         }#if
      }#if

      if (length(type) == 2) {
         par(oldpar);
      }#if

      if (visible) return(ds);
   }
)#coiplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("nuseplot", signature(x="QualTreeSet"),
   function(x,
            range    = 0,
            names    = "namepart",
            main     = "NUSE Plot",
            ylim     = c(0.8, 1.2),
            las      = 2,
            add.line = TRUE,
            ...) 
   {
      if (debug.xps()) print("------nuseplot.QualTreeSet------")

      boxplot(x,
              which   = "userinfo:fNUSEQuant",
              transfo = NULL,
              range   = range,
              names   = names,
              main    = main,
              ylim    = ylim,
              las     = las,
              ...);

      if (add.line) {
         abline(1, 0, lty=2, col="gray70");
      }#if
   }
)#nuseplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("rleplot", signature(x="QualTreeSet"),
   function(x,
            range    = 0,
            names    = "namepart",
            main     = "RLE Plot",
            ylim     = c(-1.0,1.0),
            las      = 2,
            add.line = TRUE,
            ...) 
   {
      if (debug.xps()) print("------rleplot.QualTreeSet------")

      boxplot(x,
              which   = "userinfo:fRLEQuant",
              transfo = NULL,
              range   = range,
              names   = names,
              main    = main,
              ylim    = ylim,
              las     = las,
              ...);

      if (add.line) {
         abline(0, 0, lty=2, col="gray70");
      }#if
   }
)#rleplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("NUSE", signature(x="QualTreeSet"),
   function(x,
            treename = "*",
            type     = c("plot", "stats", "values"),
            qualopt  = NULL,
            ...)
   {
      if (debug.xps()) print("------NUSE.QualTreeSet------")

      if (treename == "*") treename <- unlist(treeNames(x));

      exten <- extenPart(treename);
      if (is.na(match(exten, QUATYPE[1:2]))) {
         stop(paste(sQuote(exten), "is not a quality tree extension"));
      }#if

      type <- match.arg(type);
      if (type == "stats") {
         stats <- treeInfo(x,
                           treename = treename,
                           treetype = exten,
                           varlist  = "fNQuantiles:fNUSEQuant",
                           qualopt  = qualopt,
                           ...);
         return(stats);
      } else if (type == "values") {
         nuse <- export(x,
                        treename     = treename,
                        treetype     = exten,  # "rlm" or "plm"
                        varlist      = "fUnitName:fNUSE",
                        as.dataframe = TRUE,
                        verbose      = FALSE);
         return(nuse);
      }#if

      if (type == "plot") {
         nuseplot(x, 
                  names = qualopt,
                  ...);
      }#if
   }
)#NUSE

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("RLE", signature(x="QualTreeSet"),
   function(x,
            treename = "*",
            type     = c("plot", "stats", "values"),
            qualopt  = NULL,
            ...) 
   {
      if (debug.xps()) print("------RLE.QualTreeSet------")

      if (treename == "*") treename <- unlist(treeNames(x));

      exten <- extenPart(treename);
      if (is.na(match(exten, QUATYPE[1:2]))) {
         stop(paste(sQuote(exten), "is not a quality tree extension"));
      }#if

      type <- match.arg(type);
      if (type == "stats") {
         stats <- treeInfo(x,
                           treename = treename,
                           treetype = exten,
                           varlist  = "fNQuantiles:fRLEQuant",
                           qualopt  = qualopt,
                           ...);
         return(stats);
      } else if (type == "values") {
         rle <- export(x,
                       treename     = treename,
                       treetype     = exten,  # "rlm" or "plm"
                       varlist      = "fUnitName:fRLE",
                       as.dataframe = TRUE,
                       verbose      = FALSE);
         return(rle);
      }#if

      if (type == "plot") {
         rleplot(x,
                 names = qualopt,
                 ...);
      }#if
   }
)#RLE

#------------------------------------------------------------------------------#
