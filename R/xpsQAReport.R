#------------------------------------------------------------------------------#
# xpsQAReport: create quality control report
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
xpsQAReport <- 
function(xps.data,
         xps.expr    = NULL,
         xps.call    = NULL,
         xps.qual    = NULL,
         dataset     = character(0), 
         title       = "Quality Report",
         date        = "October, 2011",
         author      = "Christian Stratowa",
         outdir      = file.path(getwd(), "QAReport"),
         add.pseudo  = FALSE,
         overwrite   = FALSE,
         verbose     = TRUE,
         ...) 
{
   ## check for presence of data
   if (!is.null(xps.data)) {
      chipname <- chipName(schemeSet(xps.data));
      chiptype <- chipType(schemeSet(xps.data));
      numtrees <- xps.data@numtrees;
   } else if (!is.null(xps.expr)) { 
      chipname <- chipName(schemeSet(xps.expr));
      chiptype <- chipType(schemeSet(xps.expr));
      numtrees <- xps.expr@numtrees;
   } else if (!is.null(xps.call)) { 
      chipname <- chipName(schemeSet(xps.call));
      chiptype <- chipType(schemeSet(xps.call));
      numtrees <- xps.call@numtrees;
   } else if (!is.null(xps.qual)) { 
      chipname <- chipName(schemeSet(xps.qual));
      chiptype <- chipType(schemeSet(xps.qual));
      numtrees <- xps.qual@numtrees;
   } else { 
      stop("at least one of <xps.data, xps.expr, xps.call, xps.qual> must exist");
   }#if

   ## directory containing parts of QAReport.Rnw
   indir <- file.path(path.package("xps"), "QC");

   ## create directory containing final QAReport.Rnw
   if (file.exists(outdir)) {
      if (overwrite) unlink(outdir, recursive=TRUE)
      else           stop("QAReport does already exist.");
   }#if
   if (!dir.create(outdir)) 
      stop("could not create report directory");
   if (!dir.create(file.path(outdir, "png"))) 
      stop("could not create report subdirectory 'png'");
   docdir <- file.path(outdir, "png");

   ## optional pagebreak
   QCp <- readLines(file.path(indir, "QC.break.Rnw"));


   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   # QAReport.Rnw part containing title and introduction
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   QCb <- readLines(file.path(indir, "QC.begin.Rnw"));

   ## replace title, date, author
   QCb <- sub("@TITLE@",  title,  QCb);
   QCb <- sub("@DATE@",   date,   QCb);
   QCb <- sub("@AUTHOR@", author, QCb);

   ## dataset info
   QCb <- sub("@DATASET@",  dataset,  QCb);
   QCb <- sub("@NUMTREES@", numtrees, QCb);
   QCb <- sub("@CHIPNAME@", chipname, QCb);
   QCb <- sub("@CHIPTYPE@", chiptype, QCb);

   ## need to replace "_" with "\_"
   QCb <- gsub("_","\\\\_", QCb);

   write(QCb, file.path(docdir, "QAReport.Rnw"));


   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   # QAReport.Rnw part containing raw data
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if (!is.null(xps.data)) {
      if (verbose) print("quality assessment of DataTreeSet...");    
      if (!is(xps.data, "DataTreeSet")) {
         stop(paste(sQuote("xps.data"), "is not a class", sQuote("DataTreeSet")));
      }#if

      if (chiptype == "GeneChip") which <- "pm"
      else                        which <- "core";

      ## hist
      plotDensity(xps.data, 
                  which      = which,
                  add.legend = TRUE,
                  dev        = "png",
                  outfile    = file.path(docdir, "DensityPlotData"),
                  w          = 540,
                  h          = 540,
                  verbose    = verbose);

      ## boxplot
      plotBoxplot(xps.data,
                  which   = "userinfo:fIntenQuant",
                  dev     = "png",
                  outfile = file.path(docdir, "BoxplotData"),
                  w       = 540,
                  h       = 540);

      QCd <- readLines(file.path(indir, "QC.data.Rnw"));
      QCd <- sub("@WHICH@", which, QCd);

      write(QCd, file.path(docdir, "QAReport.Rnw"), append=TRUE);
      write(QCp, file.path(docdir, "QAReport.Rnw"), append=TRUE);
   }#if


   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   # QAReport.Rnw part containing expression data
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if (!is.null(xps.expr)) {
      if (verbose) print("quality assessment of ExprTreeSet...");    
      if (!is(xps.expr, "ExprTreeSet")) {
         stop(paste(sQuote("xps.expr"), "is not a class", sQuote("ExprTreeSet")));
      }#if

      exprtype <- toupper(xps.expr@exprtype);

      ## hist
      plotDensity(xps.expr, 
                  add.legend = TRUE,
                  dev        = "png",
                  outfile    = file.path(docdir, "DensityPlotExpr"),
                  w          = 540,
                  h          = 540,
                  verbose    = verbose);

      ## boxplot
      plotBoxplot(xps.expr,
                  which   = "userinfo:fLevelQuant",
                  dev     = "png",
                  outfile = file.path(docdir, "BoxplotExpr"),
                  w       = 540,
                  h       = 540);

      ## correlation plot
      plotCorr(xps.expr,
               sort       = TRUE,
               add.legend = TRUE,
               dev        = "png",
               outfile    = file.path(docdir, "CorrplotExpr"));

      ## MAD plot
      plotMAD(xps.expr,
             sort       = TRUE,
             add.legend = TRUE,
             dev        = "png",
             outfile    = file.path(docdir, "MADplotExpr"));

      ## PCA plot
      plotPCA(xps.expr,
             add.labels = TRUE,
             dev        = "png",
             outfile    = file.path(docdir, "PCAplotExpr"));

      QCx <- readLines(file.path(indir, "QC.expr.Rnw"));
      QCx <- sub("@EXPRTYPE@", exprtype, QCx);

      write(QCx, file.path(docdir, "QAReport.Rnw"), append=TRUE);
      write(QCp, file.path(docdir, "QAReport.Rnw"), append=TRUE);
   }#if


   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   # QAReport.Rnw part containing detection call data
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if (!is.null(xps.call)) {
      if (verbose) print("quality assessment of CallTreeSet...");    
      if (!is(xps.call, "CallTreeSet")) {
         stop(paste(sQuote("xps.call"), "is not a class", sQuote("CallTreeSet")));
      }#if

      calltype <- toupper(xps.call@calltype);
      treetype <- extenPart(treeNames(xps.call));

      ## percent present call
      call <- treeInfo(xps.call,
                       treetype = treetype,
                       varlist  = "userinfo:fPcAbsent:fPcMarginal:fPcPresent");
      write.table(t(call), file.path(docdir, "PrecentCall.txt"), sep="\t", col.names=NA);

      QCc <- readLines(file.path(indir, "QC.call.Rnw"));
      QCc <- sub("@CHIPNAME@", chipname, QCc);
      QCc <- sub("@DETCALL@",  calltype, QCc);
      QCc <- gsub("_","\\\\_", QCc);

      write(QCc, file.path(docdir, "QAReport.Rnw"), append=TRUE);
      write(QCp, file.path(docdir, "QAReport.Rnw"), append=TRUE);
   }#if


   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   # QAReport.Rnw part containing additional quality control data
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if (!is.null(xps.qual)) {
      if (verbose) print("quality assessment of QualTreeSet...");    
      if (!is(xps.qual, "QualTreeSet")) {
         stop(paste(sQuote("xps.qual"), "is not a class", sQuote("QualTreeSet")));
      }#if

      ## RNA degradation
      rnadeg <- AffyRNAdeg(xps.qual);

      ## plot RNA degradation
      png(filename = file.path(docdir, "RNADegradationPlot.png"), width=540, height=540);
      plotAffyRNAdeg(rnadeg, add.legend=TRUE);
      dev.off()

      ## plot slope of RNA degradation
      png(filename = file.path(docdir, "RNADegradationSlope.png"), width=540, height=440);
      plotAffyRNAdeg(rnadeg, summary=TRUE);
      dev.off();

      ## border plots
      plotBorder(xps.qual,
                 type    = "pos",
                 dev     = "png",
                 outfile = file.path(docdir, "BorderPlotPos"),
                 w       = 540,
                 h       = 500);

      plotBorder(xps.qual,
                 type    = "neg",
                 dev     = "png",
                 outfile = file.path(docdir, "BorderPlotNeg"),
                 w       = 540,
                 h       = 500);

      ## Center-of-Intensity plots
      plotCOI(xps.qual,
              dev     = "png",
              outfile = file.path(docdir, "COIPlot"),
              w       = 540,
              h       = 400);

      ## NUSE plot
      plotNUSE(xps.qual,
               names   = "namepart:normalized",
               dev     = "png",
               outfile = file.path(docdir, "NUSEPlot"),
               w       = 540,
               h       = 540);

      ## RLE plot
      plotRLE(xps.qual,
               names   = "namepart:normalized",
               dev     = "png",
               outfile = file.path(docdir, "RLEPlot"),
               w       = 540,
               h       = 540);

      QCq <- readLines(file.path(indir, "QC.qual.Rnw"));

      write(QCq, file.path(docdir, "QAReport.Rnw"), append=TRUE);
      write(QCp, file.path(docdir, "QAReport.Rnw"), append=TRUE);
   }#if


   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   # QAReport.Rnw part containing pseudo images
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if (!is.null(xps.qual) && add.pseudo) {
      if (verbose) print("pseudo images...");    

      treenames <- getTreeNames(rootFile(xps.qual), treetype=QUATYPE[4]);
      treenames <- treenames[grep("raw", treenames)];
      numtrees  <- length(treenames);

      QCpt <- readLines(file.path(indir, "QC.pseudo.txt.Rnw"));
      QCpt <- sub("@NUMTREES@", numtrees, QCpt);

      write(QCpt, file.path(docdir, "QAReport.Rnw"), append=TRUE);

      ## Pseudo-images in sets of four
      if (is.logical(add.pseudo)) n <- ceiling(numtrees/4)
      else                        n <- numtrees;
      for (i in 1:n) {
         if (is.logical(add.pseudo)) names <- treenames[(i*4-3):(i*4)]
         else                        names <- treenames[i];

         plotImage(xps.qual,
                   type    = "weights",
                   qualopt = "raw",
                   names   = names,
                   dev     = "png",
                   outfile = file.path(docdir, paste("PseudoImage", i, sep="")),
                   w       = 640,
                   h       = 680,
                   verbose = TRUE);

         QCpf <- readLines(file.path(indir, "QC.pseudo.fig.Rnw"));
         QCpf <- sub("@PSEUDOIMAGE@", paste("PseudoImage", i, sep=""), QCpf);

         write(QCpf, file.path(docdir, "QAReport.Rnw"), append=TRUE);
         write(QCp,  file.path(docdir, "QAReport.Rnw"), append=TRUE);
      }#for
   }#if


   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   # QAReport.Rnw part containing end of document
   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   QCe <- readLines(file.path(indir, "QC.end.Rnw"));
   QCe <- sub("@DATASET@",  dataset,  QCe);
   QCe <- gsub("_","\\\\_", QCe);

   write(QCe, file.path(docdir, "QAReport.Rnw"), append=TRUE);

   ## build vignette QC.pdf
   if (requireNamespace("tools")) {
      olddir <- getwd();
      setwd(docdir);

      tryCatch({Sweave("QAReport.Rnw");
                setwd(outdir);
                tools::texi2pdf(file.path(docdir, "QAReport.tex"), clean=TRUE)
               },
               finally = setwd(olddir)
              );
   }#if
}#xpsQAReport

#------------------------------------------------------------------------------#
