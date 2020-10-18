#------------------------------------------------------------------------------#
# function to correct certain Affymetrix annotation files, see:
# https://www.stat.math.ethz.ch/pipermail/bioconductor/2009-August/029049.html

# need to open binary file to prevent conversion to CRLF on WinXP, see:
# http://tolstoy.newcastle.edu.au/R/help/03b/5893.html
#------------------------------------------------------------------------------#
"updateAnnotation" <- function(infile, outfile, probeset, skip, eol="\n") {
   ## read header and probesets
   cat("reading", infile, "...\n");
   header <- readLines(infile, n=skip);
   annot  <- read.csv(infile, colClasses="character", comment.char="", skip=skip);

   ## delete probeset
   for (i in 1:length(probeset)) {
      line  <- which(annot[,"probeset_id"] == probeset[i]);
      if (length(line) > 0) {
         cat("deleting line", line, "for probeset", probeset[i], "...\n");
         annot <- annot[-line,];
      }#if
   }#for

   ## write header and append probesets
   cat("writing", outfile, "...\n");
   file <- file(outfile, "wb")
   writeLines(header, con=file, sep=eol);
   write.table(annot, file=file, append=TRUE, sep=",", eol=eol, row.names=FALSE);
   close(file) 
}#updateAnnotation

#------------------------------------------------------------------------------#
"deleteNegControlFromAffxControl" <- function(infile, outfile, eol="\n") {
   ifile <- file(infile, "r");
   ofile <- file(outfile, "wb");

   line <- "";
   while (length(line) > 0) {
      line <- readLines(ifile, n=1);

      skip <- grep("control->affx",line) & grep("neg_control",line);
      if (length(skip) && skip) next;

      writeLines(line, con=ofile, sep=eol);
   }#while

   close(ifile);
   close(ofile);
}#deleteNegControlFromAffxControl


#------------------------------------------------------------------------------#
