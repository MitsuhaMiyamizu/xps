\name{export.filter}
\alias{export.filter}
\title{Export filter data as text files}
\description{
  Export data from classes \code{\link{FilterTreeSet}} or \code{\link{AnalysisTreeSet}} to \code{outfile}.
}
\usage{
export.filter(xps.fltr, treename = "*", treetype = character(0), varlist = "*", outfile = character(0), sep = "\t", as.dataframe = FALSE, verbose = TRUE)
}
\arguments{
  \item{xps.fltr}{an object of type \code{\link{FilterTreeSet}} or \code{\link{AnalysisTreeSet}}.}
  \item{treename}{tree name to export.}
  \item{treetype}{type of tree(s) to export, 'pfr', 'ufr' or 'stt'.}
  \item{varlist}{names of tree leaves to export.}
  \item{outfile}{name of output file.}
  \item{sep}{column separator}
  \item{as.dataframe}{if \code{TRUE} a data.frame will be returned.}
  \item{verbose}{logical, if \code{TRUE} print status information.}
}
\details{
 Export data from classes \code{\link{FilterTreeSet}}, or \code{\link{AnalysisTreeSet}}
 to \code{outfile}.

 Parameter \code{varlist} lists the parameters to export: \cr
 - parameters are separated by ":", e.g. \code{varlist="fUnitName:fFlag"}. \cr
 - for \code{varlist="*"} all valid parameters will be exported.

 For class \code{FilterTreeSet} the following \code{varlist} parameters are valid:
 \tabular{lll}{
  \tab \code{fUnitName}:\tab unit name (probeset ID). \cr
  \tab \code{fFlag}:\tab mask.
 }

 For class \code{AnalysisTreeSet} \code{varlist} can contain annotation
 parameters and parameters of the resulting data. \cr
  Following \code{varlist} annotation parameters are valid:
 \tabular{lll}{
  \tab \code{fUnitName}:\tab unit name (probeset ID). \cr
  \tab \code{fTranscriptID}:\tab transcript\_id (probeset ID). \cr
  \tab \code{fName}:\tab gene name. \cr
  \tab \code{fSymbol}:\tab gene symbol. \cr
  \tab \code{fAccession}:\tab mRNA accession such as Refseq ID. \cr
  \tab \code{fEntrezID}:\tab entrez ID. \cr
  \tab \code{fChromosome}:\tab chromosome. \cr
  \tab \code{fStart}:\tab start position. \cr
  \tab \code{fStop}:\tab stop position. \cr
  \tab \code{fStrand}:\tab strand on chromosome. \cr
  \tab \code{fCytoBand}:\tab cytoband.
 }
 For class \code{AnalysisTreeSet} the following \code{varlist} parameters are valid:
 \tabular{lll}{
  \tab \code{mn1}:\tab mean of group 1. \cr
  \tab \code{mn2}:\tab mean of group 2. \cr
  \tab \code{fc}:\tab fold-change fc=mn2/mn1. \cr
  \tab \code{se}:\tab standard error. \cr
  \tab \code{df}:\tab degree of freedom. \cr
  \tab \code{stat}:\tab t-statistic. \cr
  \tab \code{pval}:\tab p-value. \cr
  \tab \code{nper}:\tab number of permutations. \cr
  \tab \code{pcha}:\tab p-chance. \cr
  \tab \code{padj}:\tab adjusted p-value. \cr
  \tab \code{flag}:\tab flag. \cr
  \tab \code{mask}:\tab only rows with \code{flag=1} will be exported.
 }
}
\value{
If \code{as.dataframe} is \code{TRUE}, the data will be imported into the
current R session as \code{data.frame}. Otherwise, \code{NULL} will be returned.
}
\author{Christian Stratowa}
\seealso{\code{\link{export-methods}}}
\keyword{manip}
