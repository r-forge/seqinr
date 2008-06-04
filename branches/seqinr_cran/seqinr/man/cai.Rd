\name{cai}
\alias{cai}
\title{Codon Adaptation Index (Sharp and Li 1987)}
\description{
Adapted from the documentation of the CAI function in the
program \code{codonW} writen by John Peden:
CAI is a measurement
of the relative adaptiveness of the codon usage of a gene towards the
codon usage of highly expressed genes. The relative adaptiveness (w) of
each codon is the ratio of the usage of each codon, to that of the most
abundant codon for the same amino acid. The CAI
index is defined as the geometric mean of these relative adaptiveness
values. Non-synonymous codons and termination codons (genetic code
dependent) are excluded. To aid computation, the CAI is calculated as
using a natural log summation, To prevent a codon having a relative
adaptiveness value of zero, which could result in a CAI of zero;
these codons have fitness of zero (<.0001) are adjusted to 0.01.
}
\usage{
  cai(seq, w, numcode = 1, zero.threshold = 0.0001, zero.to = 0.01)
}
\arguments{
  \item{seq}{a coding sequence as a vector of single characters}
  \item{w}{a vector for the relative adaptiveness of each codon}
  \item{numcode}{the genetic code number as in \code{\link{translate}}}
  \item{zero.threshold}{a value in \code{w} below this threshold is
                      considered as zero}
  \item{zero.to}{a value considered as zero in \code{w} is forced to
               this value}
}
\details{
  More doc needed here
}
\value{
  A single numerical value for the CAI.
}
\references{
  
  \code{citation("seqinr")}.
  
  The program codonW used here for comparison is available at 
  \url{http://codonw.sourceforge.net/}.
}
\seealso{
  \code{\link{caitab}} for some \code{w} values
}
\author{J.R. Lobry}
\examples{
#
# How to reproduce the results obtained with the C program codonW
# version 1.4.4 writen by John Peden. We use here the "input.dat"
# test file from codonW (Saccharomyces cerevisiae).
#
  inputdatfile <- system.file("sequences/input.dat", package = "seqinr")
  input <- read.fasta(file = inputdatfile) # read the FASTA file
#
# Import results obtained with codonW
#
  scucofile <- system.file("sequences/scuco.txt", package = "seqinr")
  scuco.res <- read.table(scucofile, header = TRUE) # read codonW result file
#
# Use w for Saccharomyces cerevisiae
#
  data(caitab)
  w <- caitab$sc
#
# Compute CAI and compare results:
#
  cai.res <- sapply(input, cai, w = w)
  plot(cai.res, scuco.res$CAI,
    main = "Comparison of seqinR and codonW results",
    xlab = "CAI from seqinR",
    ylab = "CAI from codonW",
    las = 1)
  abline(c(0,1))
}
