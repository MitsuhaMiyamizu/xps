\name{uniTest-methods}
\docType{methods}
\alias{uniTest-methods}
\alias{uniTest}
\alias{uniTest<-}
\title{A Two-Group Unitest}
\description{
 Unitest performs a a two group uni-test such as the \code{t.test} on each row of the
 expression dataframe. \cr
 The Unitest returns a dataframe containing the results of the test.

\emph{Usage}

\code{uniTest(object)} \cr
\code{uniTest(object, value)<-}
}
\arguments{
  \item{object}{object of class \code{UniFilter}.}
  \item{value}{character vector \code{c(type, alternative, correction, numperm, mu, paired, conflevel, varequ)}}
}
\details{
 The method \code{uniTest} initializes the following parameters:
 \tabular{lll}{
  \tab \code{type}:\tab a character string specifying the type of test: currently \code{"t.test"} (default) or \code{"normal.test"}. \cr
  \tab \code{alternative}:\tab a character string specifying the alternative hypothesis,
        must be one of \code{"two.sided"} (default), \code{"greater"} or \code{"less"}. \cr
  \tab \code{correction}:\tab a correction to adjust p-values for multiple comparisons: \cr
  \tab                   \tab \code{correction="none"}: no correction (default). \cr
  \tab                   \tab \code{correction="bonferroni"}: Bonferroni correction. \cr
  \tab                   \tab \code{correction="BH" or "fdr"}: correction for false discovery rate (Benjamini & Hochberg). \cr
  \tab                   \tab \code{correction="BY"}: correction for false discovery rate (Benjamini & Yekutieli). \cr
  \tab                   \tab \code{correction="hochberg"}: Hochberg correction. \cr
  \tab                   \tab \code{correction="holm"}: Holm correction. \cr
  \tab                   \tab \code{correction="wy"}: Westfall-Young step-down adjusted p-chance (E.Manduchi). \cr
  \tab \code{numperm}:\tab optional number of permutations  used to determine p-chance (default is 0). \cr
  \tab \code{mu}:\tab a number indicating the true value of the difference in means for a two sample test (default is 0). \cr
  \tab \code{paired}:\tab a logical indicating whether you want a paired uni-test (default is FALSE). \cr
  \tab \code{conflevel}:\tab confidence level of the interval (default is 0.95). \cr
  \tab \code{varequ}:\tab a logical variable indicating whether to treat the two variances as being
        equal. If \code{TRUE} then the pooled variance is used to estimate the variance otherwise
        the Welch (or Satterthwaite) approximation to the degrees of freedom is used (default is FALSE).
 }
}
\value{
  An initialized \code{\linkS4class{UniFilter}} object.
}
\author{Christian Stratowa}
\references{
  Benjamini, Y., and Hochberg, Y. (1995).
  Controlling the false discovery rate: a practical and powerful
  approach to multiple testing.
  \emph{Journal of the Royal Statistical Society Series} B, \bold{57}, 289--300.

  Benjamini, Y., and Yekutieli, D. (2001).
  The control of the false discovery rate in multiple testing under dependency.
  \emph{Annals of Statistics} \bold{29}, 1165--1188.


  Holm, S. (1979).
  A simple sequentially rejective multiple test procedure.
  \emph{Scandinavian Journal of Statistics}, \bold{6}, 65--70.

  Westfall P.H. and Young S.S. (1993)
  Resampling-based multiple testing:examples and methods for p-value adjustment.
  \emph{Wiley series in probability and mathematical statistics}; Wiley.

  Dudoit S., Yang Y.H., Callow M.J., Speed T.P.  (2000) 
  Statistical methods for identifying differentially expressed genes in replicated cDNA microarray experiments. 
  \emph{Technical report} \bold{578}; UC Berkeley.

  Manduchi E. (2000)
  Software: tpWY, see: \url{http://www.cbil.upenn.edu/tpWY/}
}
\examples{
unifltr <- UniFilter()
uniTest(unifltr) <- c("t.test","two.sided","none",0,0.0,FALSE,0.98,TRUE)
str(unifltr)
}
\keyword{methods}
