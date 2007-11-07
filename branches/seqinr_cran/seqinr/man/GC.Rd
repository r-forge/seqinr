\name{G+C Content}
\alias{GC}
\alias{GC1}
\alias{GC2}
\alias{GC3}
\alias{GCpos}
\title{Calculates the fractional G+C content of nucleic acid sequences.}
\description{
  Calculates the fraction of G+C bases of the input nucleic acid
  sequence(s). It reads in nucleic acid sequences, sums the number of
  'g' and 'c' bases and writes out the result as the fraction (in the
  interval 0.0 to 1.0) to the total number of 'a', 'c', 'g' and 't' bases. 
  Global G+C content \code{GC}, G+C in the first position of the codon bases
  \code{GC1}, G+C in the second position of the codon bases
  \code{GC2}, and G+C in the third position of the codon bases
  \code{GC3} can be computed. All functions can take ambiguous bases
  into account when requested.
}
\usage{
GC(seq, forceToLower = TRUE, exact = FALSE, NA.GC = NA, oldGC = FALSE)
GC1(seq, frame = 0, ...)
GC2(seq, frame = 0, ...)
GC3(seq, frame = 0, ...)
GCpos(seq, pos, frame = 0, ...)
}
\arguments{
  \item{seq}{a nucleic acid sequence as a vector of single characters}
  \item{frame}{for coding sequences, an integer (0, 1, 2) giving the frame}
  \item{forceToLower}{logical defaulting to \code{TRUE}: force sequence
  	characters in lower-case. Turn this to \code{FALSE} to save time
  	if your sequence is already in lower-case}
  \item{exact}{logical defaulting to \code{FALSE}: should ambiguous bases taken
  	into account when computing the G+C content (see details)}
  \item{NA.GC}{what should be returned when the GC is impossible to
    compute from data, for instance with NNNNNNN. This behaviour could
    be different when argument \code{exact} is \code{TRUE}, for instance
    the G+C content of WWSS is \code{NA} by default, but is 0.5 when
   \code{exact} is set to \code{TRUE}}
  \item{...}{arguments passed to the function \code{GC}}
  \item{pos}{for coding sequences, the codon position (1, 2, 3) that should be
    taken into account to compute the G+C content}
  \item{oldGC}{logical defaulting to \code{FALSE}: should the GC content computed
    as in seqinR <= 1.0-6, that is as the sum of 'g' and 'c' bases divided by
    the length of the sequence. As from seqinR >= 1.1-3, this argument is
    deprecated and a warning is issued.}
}
\details{
	When \code{exact} is set to \code{TRUE} the G+C content is estimated 
	with ambiguous bases taken into account. Note that this is time expensive.
	A first pass is made on non-ambiguous bases to estimate the probabilities
	of the four bases in the sequence. They are then used to weight the
	contributions of ambiguous bases to the G+C content. Let note nx
	the total number of base 'x' in the sequence. For instance
	suppose that there are nb bases 'b'. 'b' stands for "not a", that
	is for 'c', 'g' or 't'. The contribution of 'b' bases to the GC base
	count will be:
	
	nb*(nc + ng)/(nc + ng + nt)
	
	The contribution of 'b' bases to the AT base count will be:
	
	nb*nt/(nc + ng + nt)
	
	All ambiguous bases contributions to the AT and GC counts are weighted
	is similar way and then the G+C content is computed as ngc/(nat + ngc).
}
\value{
  \code{GC} returns the fraction of G+C (in [0,1]) as a numeric vector of length one.
  \code{GCpos} returns GC at position \code{pos}.
  \code{GC1}, \code{GC2}, \code{GC3} are wrappers for \code{GCpos} with the
    argument \code{pos} set to 1, 2, and 3, respectively.
  \code{NA} is returned when \code{seq} is \code{NA}.
  \code{NA.GC} defaulting to \code{NA} is returned when the G+C content 
    can not be computed from data.
}
\references{
  \code{citation("seqinr")}.
  
  The program codonW used here for comparison is available at 
  \url{http://codonw.sourceforge.net/}.
}
\seealso{You can use \code{\link{s2c}} to convert a string into a vetor of single
character and \code{\link{tolower}} to convert upper-case characters into
lower-case characters. Do not confuse with \code{\link{gc}} for garbage collection.
}
\author{D. Charif and L. Palmeira and J.R. Lobry}
\examples{
   mysequence <- s2c("agtctggggggccccttttaagtagatagatagctagtcgta")
   GC(mysequence)  # 0.4761905
   GC1(mysequence) # 0.6428571
   GC2(mysequence) # 0.3571429
   GC3(mysequence) # 0.4285714
#
# With upper-case characters:
#
  myUCsequence <- s2c("GGGGGGGGGA")
  GC(myUCsequence) # 0.9
#
# With ambiguous bases:
#
  GC(s2c("acgt")) # 0.5
  GC(s2c("acgtssss")) # 0.5
  GC(s2c("acgtssss"), exact = TRUE) # 0.75
#
# Missing data:
#
  stopifnot(is.na(GC(s2c("NNNN"))))
  stopifnot(is.na(GC(s2c("NNNN"), exact = TRUE)))
  stopifnot(is.na(GC(s2c("WWSS"))))
  stopifnot(GC(s2c("WWSS"), exact = TRUE) == 0.5)
#
# Coding sequences tests:
#
  cdstest <- s2c("ATGATG")
  stopifnot(GC3(cdstest) == 1)
  stopifnot(GC2(cdstest) == 0)
  stopifnot(GC1(cdstest) == 0)
#
# How to reproduce the results obtained with the C program codonW
# version 1.4.4 writen by John Peden. We use here the "input.dat"
# test file from codonW (there are no ambiguous base in these
# sequences).
#
  inputdatfile <- system.file("sequences/input.dat", package = "seqinr")
  input <- read.fasta(file = inputdatfile) # read the FASTA file
  inputoutfile <- system.file("sequences/input.out", package = "seqinr")
  input.res <- read.table(inputoutfile, header = TRUE) # read codonW result file

# remove stop codon before computing G+C content (as in codonW)

  GC.codonW <- function(dnaseq, ...){
  	 GC(dnaseq[seq_len(length(dnaseq) - 3)], ...)
  }
  input.gc <- sapply(input, GC.codonW, forceToLower = FALSE)
  max(abs(input.gc - input.res$GC)) # 0.0004946237

  plot(x = input.gc, y = input.res$GC, las = 1,
  xlab = "Results with GC()", ylab = "Results from codonW",
  main = "Comparison of G+C content results")
  abline(c(0, 1), col = "red")
  legend("topleft", inset = 0.01, legend = "y = x", lty = 1, col = "red")
}
\keyword{manip}
