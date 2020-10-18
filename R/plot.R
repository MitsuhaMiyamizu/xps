#==============================================================================#
# plot.R: plot to device functions:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# plotBorder: 
# plotBoxplot: 
# plotCall: 
# plotCOI: 
# plotCorr: 
# plotDensity: 
# plotImage: 
# plotIntensity2GC: 
# plotMA: 
# plotMAD: 
# plotNUSE: 
# plotPCA: 
# plotPM: 
# plotProbeset: 
# plotRLE: 
# plotVolcano: 
#==============================================================================#

"plotBorder" <-
function(x,
         type    = c("pos", "neg"),
         qualopt = "raw",
         transfo = log2,
         range   = 0,
         names   = "namepart",
         ylim    = NULL,
         bmar    = NULL,
         las     = 2,
         dev     = "screen",
         outfile = "BorderPlot",
         w       = 800,
         h       = 540,
         ...) 
{
   if (debug.xps()) print("------plotBorder------")

   ## check for correct class
   if (!is(x, "QualTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("QualTreeSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot data
   borderplot(x,
              type    = type,
              qualopt = qualopt,
              transfo = transfo,
              range   = range,
              names   = names,
              ylim    = ylim,
              bmar    = bmar,
              las     = las,
              ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotBorder

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotBoxplot" <-
function(x,
         which   = "",
         size    = 0,
         transfo = log2,
         range   = 0,
         names   = "namepart",
         mar     = NULL,
         las     = 2,
         cex     = 1.0,
         dev     = "screen",
         outfile = "BoxPlot",
         w       = 800,
         h       = 540,
         ...) 
{
   if (debug.xps()) print("------plotBoxplot------")

   ## check for correct class
   if (!extends(class(x), "ProcesSet")) {
      stop(paste(sQuote("x"), "is not derived from class", sQuote("ProcesSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   if (is.null(mar)) {
      bmar <- NULL;
   } else {
      bmar <- list(b=mar[1], cex=cex, w=w);
   }#if

   ## plot data
   boxplot(x,
           which   = which,
           size    = size,
           transfo = transfo,
           range   = range,
           names   = names,
           bmar    = bmar,
           las     = las,
           ...)

   if (dev != "screen") {
      dev.off();
   }#if
}#plotBoxplot

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotCall" <-
function(x,
         beside  = TRUE,
         names   = "namepart",
         col     = c("red","green","blue"),
         legend  = c("P","M","A"),
         ylim    = c(0,100),
         ylab    = "detection call [%]",
         las     = 2,
         dev     = "screen",
         outfile = "CallPlot",
         w       = 800,
         h       = 540,
         ...) 
{
   if (debug.xps()) print("------plotCall------")

   ## check for correct class
   if (!is(x, "CallTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("CallTreeSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot data
   callplot(x,
            beside = beside,
            names  = names,
            col    = col,
            legend = legend,
            ylim   = ylim,
            ylab   = ylab,
            las    = las,
            ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotCall

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotCOI" <-
function(x,
         type    = c("pos", "neg"),
         qualopt = "raw",
         radius  = 0.5,
         linecol = "gray70",
         visible = TRUE,
         dev     = "screen",
         outfile = "CenterOfIntensityPlot",
         w       = 540,
         h       = 540,
         ...)
{
   if (debug.xps()) print("------plotCOI------")

   ## check for correct class
   if (!is(x, "QualTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("QualTreeSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot data
   coiplot(x,
           type    = type,
           qualopt = qualopt,
           radius  = radius,
           linecol = linecol,
           visible = visible,
           ...) 

   if (dev != "screen") {
      dev.off();
   }#if
}#plotCOI

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotCorr" <-
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
         dev        = "screen",
         outfile    = "CorrelationPlot",
         w          = 540,
         h          = 540,
         ...) 
{
   if (debug.xps()) print("------plotCorr------")

   ## check for correct class
   if (!is(x, "ExprTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("ExprTreeSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot data
   corplot(x,
           which      = which,
           transfo    = transfo,
           method     = method,
           col        = col,
           names      = names,
           sort       = sort,
           reverse    = reverse,
           bmar       = bmar,
           add.legend = add.legend,
            ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotCorr

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotDensity" <-
function(x,
         which      = "",
         size       = 0,
         transfo    = log2,
         ylab       = "density",
         xlab       = "log intensity",
         names      = "namepart",
         type       = "l",
         col        = 1:6,
         lty        = 1:5,
         add.legend = FALSE,
         dev        = "screen",
         outfile    = "DensityPlot",
         w          = 540,
         h          = 540,
         verbose    = TRUE,
         ...) 
{
   if (debug.xps()) print("------plotDensity------")

   ## check for correct class
   if (!extends(class(x), "ProcesSet")) {
      stop(paste(sQuote("x"), "is not derived from class", sQuote("ProcesSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot data
   hist(x,
        which      = which,
        size       = size,
        transfo    = transfo,
        ylab       = ylab,
        xlab       = xlab,
        names      = names,
        type       = type,
        col        = col,
        lty        = lty,
        add.legend = add.legend,
        verbose    = verbose,
         ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotDensity

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotImage" <-
function(x,
         type    = character(),
         qualopt = c("raw", "adjusted", "normalized"),
         transfo = log2,
         col     = NULL,
         names   = character(),
         dev     = "screen",
         outfile = "Image",
         w       = 800,
         h       = 800,
         verbose = TRUE,
         ...) 
{
   if (debug.xps()) print("------plotImage------");

   qualopt <- match.arg(qualopt);
   if (names[1] == "namepart") {
      stop("<names='namepart'> is not possible for this function");
   }#if

   ## check for correct class
   if (is (x, "DataTreeSet")) {
      if (is.na(match(type, c("intensity", "background")))) {
         stop(paste(sQuote(type), "is not a valid residual option"));
      }#if

      if (names[1] == "*") names <- unlist(treeNames(x));

      bg <- FALSE;
   } else if (is (x, "ExprTreeSet")) {
      if (is.na(match(type, c("intensity", "background")))) {
         stop(paste(sQuote(type), "is not a valid data option"));
      }#if

      if (type == "intensity") {
         if (names[1] == "*") names <- getTreeNames(rootFile(x), treetype=ADJTYPE[1]);

         bg <- FALSE;
      } else {
         if (names[1] == "*") {
            exten <- unique(extenPart(getTreeNames(rootFile(x))));
            id    <- match(exten, BGDTYPE);
            id    <- id[!is.na(id)];
            names <- getTreeNames(rootFile(x), treetype=BGDTYPE[id]);
         }#if

         bg <- TRUE;
      }#if
   } else if (is (x, "QualTreeSet")) {
      if (is.na(match(type, RESIOPT))) {
         stop(paste(sQuote(type), "is not a valid residual option"));
      }#if

      if (names[1] == "*") names <- getTreeNames(rootFile(x), treetype=QUATYPE[4]);

      names <- names[grep(qualopt, names)];
   } else {
      stop(paste(sQuote("x"), "must be one of class <DataTreeSet, QualTreeSet>"));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if
   oldpar <- par(no.readonly=TRUE);

   ## layout for images
   size <- function(x){c(ceiling(sqrt(x)), round(sqrt(x)))}
   nr   <- size(length(names))[1];
   nc   <- size(length(names))[2];
   nfig <- layout(matrix(c(1:(nr*nc)), nr, nc, byrow=TRUE), respect=TRUE);

   ## images
   layout.show(nfig);
   for (i in 1:length(names)) {
      if (verbose) cat("drawing image", i, "of", length(names), "...\r");

      if (is (x, "DataTreeSet") | is(x, "ExprTreeSet")) {
         image(x,
               bg         = bg,
               transfo    = transfo,
               col        = col,
               names      = names[i],
               xlab       = "",
               ylab       = "",
               add.legend = FALSE,
               ...);
      } else if (is (x, "QualTreeSet")) {
         image(x,
               type       = type,
               qualopt    = qualopt,
               transfo    = transfo,
               col        = col,
               names      = names[i],
               xlab       = "",
               ylab       = "",
               add.legend = FALSE,
               ...);
      }#if
   }#for
   if (verbose) cat("finished drawing", length(names), "images.   \n");

   par(oldpar);
   if (dev != "screen") {
      dev.off();
   }#if
}#plotImage

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotIntensity2GC" <-
function(x,
         treename, 
         which   = "",
         transfo = log2,
         range   = 0,
         col     = c("lightblue", "darkblue"),
         dev     = "screen",
         outfile = "Intensity2GCPlot",
         w       = 540,
         h       = 540,
         ...) 
{
   if (debug.xps()) print("------plotIntensity2GC------")

   ## check for correct class
   if (!is(x, "DataTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("DataTreeSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot intensity vs GC content
   intensity2GCplot(x,
                    treename = treename, 
                    which    = which,
                    transfo  = transfo,
                    range    = range,
                    col      = col,
                    ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotIntensity2GC

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotMA" <-
function(x,
         transfo = log2,
         method  = "median",
         names   = "namepart",
         ylim    = c(-6,6),
         xlab    = "A",
         ylab    = "M",
         pch     = '.',
         mar     = c(3,3,2,1),
         dev     = "screen",
         outfile = "MvAPlot",
         w       = 540,
         h       = 540,
         ...) 
{
   if (debug.xps()) print("------plotMA------")

   ## check for correct class
   if (!is (x, "ExprTreeSet")) {
      stop(paste(sQuote("x"), "is not  class", sQuote("ExprTreeSet")));
   }#if

   ## get expression levels from data
   ds <- validData(x);

   if (is.function(transfo)) ds <- transfo(ds);
   if (method == "median")   mn <- apply(ds, 1, median)
   else                      mn <- rowMeans(ds);

   if (is.null(names))              names <- colnames(ds)
   else if (names[1] == "namepart") names <- namePart(colnames(ds))
   else                             ds    <- ds[, names, drop=F];

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if
   oldpar <- par(no.readonly=TRUE);

   size <- function(x){c(ceiling(sqrt(x)), round(sqrt(x)))}
   nc <-size(ncol(ds))[1];
   nr <-size(ncol(ds))[2];
#   nc <- nr <- ceiling(sqrt(ncol(ds)));
   nfig   <- layout(matrix(c(1:(nr*nc)), nr, nc, byrow=TRUE), respect=TRUE);

   ## images
   layout.show(nfig);
   for (i in 1:ncol(ds)) {
      plot(x    = mn,
           y    = (ds[,i] - mn),
           main = names[i],
           ylim = ylim,
           xlab = xlab,
           ylab = ylab,
           pch  = pch,
           ...);
      abline(h=0, col=gray(0.5));
   }#for

   par(oldpar);
   if (dev != "screen") {
      dev.off();
   }#if
}#plotMA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotMAD" <-
function(x,
         which      = "UnitName",
         transfo    = log2,
         col        = NULL,
         names      = "namepart",
         sort       = FALSE,
         bmar       = NULL,
         add.legend = FALSE,
         dev        = "screen",
         outfile    = "MADPlot",
         w          = 540,
         h          = 540,
         ...) 
{
   if (debug.xps()) print("------plotMAD------")

   ## check for correct class
   if (!is(x, "ExprTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("ExprTreeSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot data
   madplot(x,
           which      = which,
           transfo    = transfo,
           col        = col,
           names      = names,
           sort       = sort,
           bmar       = bmar,
           add.legend = add.legend,
           ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotMAD

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotNUSE" <-
function(x,
         which    = "UnitName",
         size     = 0,
         range    = 0,
         names    = "namepart",
         main     = "NUSE Plot",
         ylim     = c(0.8,1.2),
         las      = 2,
         add.line = TRUE,
         outline  = FALSE,
         dev      = "screen",
         outfile  = "NUSEPlot",
         w        = 800,
         h        = 540,
         ...) 
{
   if (debug.xps()) print("------plotNUSE------")

   ## check for correct class
   if (!extends(class(x), "ProcesSet")) {
      stop(paste(sQuote("x"), "is not derived from class", sQuote("ProcesSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot data
   if (is (x, "ExprTreeSet")) {
      nuseplot(x,
               which    = which,
               size     = size,
               range    = range,
               names    = names,
               main     = main,
               ylim     = ylim,
               las      = las,
               add.line = add.line,
               outline  = outline,
               ...);
   } else if (is (x, "QualTreeSet")) {
      nuseplot(x,
               range    = range,
               names    = names,
               main     = main,
               ylim     = ylim,
               las      = las,
               add.line = add.line,
               ...);
   }#if

   if (dev != "screen") {
      dev.off();
   }#if
}#plotNUSE

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotPCA" <-
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
         dev        = "screen",
         outfile    = "PCAPlot",
         w          = 540,
         h          = 540,
         ...) 
{
   if (debug.xps()) print("------plotPCA------")

   ## check for correct class
   if (!is(x, "ExprTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("ExprTreeSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot data
   pcaplot(x,
           which      = which,
           transfo    = transfo,
           method     = method,
           groups     = groups,
           screeplot  = screeplot,
           squarepca  = squarepca,
           pcs        = pcs,
           add.labels = add.labels,
           add.legend = add.legend,
           col        = col,
           names      = names,
           as.list    = as.list,
           ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotPCA

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotPM" <-
function(x,
         which   = "",
         size    = 0,
         transfo = NULL,
         method  = mean,
         names   = "namepart",
         beside  = TRUE,
         col     = c("red", "blue"),
         legend  = c("PM","MM"),
         las     = 2,
         ylab    = "mean intensities",
         dev     = "screen",
         outfile = "PMPlot",
         w       = 540,
         h       = 540,
         ...) 
{
   if (debug.xps()) print("------plotPM------")

   ## check for correct class
   if (!is(x, "DataTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("DataTreeSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot data
   pmplot(x,
          which   = which,
          size    = size,
          transfo = transfo,
          method  = method,
          names   = names,
          beside  = beside,
          col     = col,
          legend  = legend,
          las     = las,
          ylab    = ylab,
          ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotPM

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotProbeset" <-
function(x,
         unitID,
         unittype   = "transcript",
         which      = "pm",
         transfo    = log2,
         names      = "namepart",
         ylim       = NULL,
         col        = 1:6,
         lty        = 1:5,
         add.legend = FALSE,
         dev        = "screen",
         outfile    = "ProbesetPlot",
         w          = 540,
         h          = 540,
         ...) 
{
   if (debug.xps()) print("------plotProbeset------")

   ## check for correct class
   if (!is(x, "DataTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("DataTreeSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot probeset
   probesetplot(x,
        unitID     = unitID,
        unittype   = unittype,
        which      = which,
        transfo    = transfo,
        names      = names,
        ylim       = ylim,
        col        = col,
        lty        = lty,
        add.legend = add.legend,
         ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotProbeset

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotRLE" <-
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
         dev      = "screen",
         outfile  = "RLEPlot",
         w        = 800,
         h        = 540,
         verbose  = TRUE,
         ...) 
{
   if (debug.xps()) print("------plotRLE------")

   ## check for correct class
   if (!extends(class(x), "ProcesSet")) {
      stop(paste(sQuote("x"), "is not derived from class", sQuote("ProcesSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot data
   if (is (x, "ExprTreeSet")) {
      rleplot(x,
              which    = which,
              size     = size,
              range    = range,
              names    = names,
              main     = main,
              ylim     = ylim,
              las      = las,
              add.line = add.line,
              outline  = outline,
              ...);
   } else if (is (x, "QualTreeSet")) {
      rleplot(x,
              range    = range,
              names    = names,
              main     = main,
              ylim     = ylim,
              las      = las,
              add.line = add.line,
              ...);
   }#if

   if (dev != "screen") {
      dev.off();
   }#if
}#plotRLE

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

"plotVolcano" <-
function(x,
         labels      = "",
         p.value     = "pval",
         mask        = FALSE,
         show.cutoff = TRUE,
         cex.text    = 0.7,
         col.text    = "blue",
         col.cutoff  = "grey",
         xlim        = NULL,
         xlab        = "Log2(Fold-Change)",
         ylab        = "-Log10(P-Value)",
         pch         = '.',
         dev         = "screen",
         outfile     = "VolcanoPlot",
         w           = 540,
         h           = 540,
         ...) 
{
   if (debug.xps()) print("------plotVolcano------")

   ## check for correct class
   if (!is(x, "AnalysisTreeSet")) {
      stop(paste(sQuote("x"), "is not class", sQuote("AnalysisTreeSet")));
   }#if

   ## add extension to outfile
   outfile <- paste(outfile, dev, sep=".");

   ## output device
   if (dev == "screen") {
      x11();
   } else if (dev == "jpeg") {
      jpeg(filename=outfile, width=w, height=h);
   } else if (dev == "png") {
      png(filename=outfile, width=w, height=h);
   } else if (dev == "pdf") {
      pdf(file=outfile, width=w, height=h);
   } else if (dev == "ps") {
      postscript(file=outfile, width=w, height=h);
   } else {
      stop(paste("unknown device dev=", sQuote(dev)));
   }#if

   ## plot data
   volcanoplot(x,
               labels      = labels,
               p.value     = p.value,
               mask        = mask,
               show.cutoff = show.cutoff,
               cex.text    = cex.text,
               col.text    = col.text,
               col.cutoff  = col.cutoff,
               xlim        = xlim,
               xlab        = xlab,
               ylab        = ylab,
               pch         = pch,
               ...);

   if (dev != "screen") {
      dev.off();
   }#if
}#plotVolcano

#------------------------------------------------------------------------------#
