\name{EXP}
\alias{EXP}
\docType{data}
\title{ Vectors of coefficients to compute linear forms }
\description{
This dataset is used to compute linear forms on codon frequencies:
if \code{codfreq} is a vector of codon frequencies then
\code{drop(freq \%*\% EXP$CG3)} will return for instance the G+C content
in third codon positions. Base order is the lexical order: a,
c, g, t (or u).
}

\usage{data(EXP)}

\format{ List of 42 vectors of coefficients
\describe{
  \item{A   }{num [1:4] 1 0 0 0 }
  \item{A3  }{num [1:64] 1 0 0 0 1 0 0 0 1 0 ...}
  \item{AGZ }{num [1:64] 0 0 0 0 0 0 0 0 1 0 ...}
  \item{ARG }{num [1:64] 0 0 0 0 0 0 0 0 1 0 ...}
  \item{AU3 }{num [1:64] 1 0 0 1 1 0 0 1 1 0 ...}
  \item{BC  }{num [1:64] 0 1 0 0 0 0 0 0 0 0 ...}
  \item{C   }{num [1:4] 0 1 0 0}
  \item{C3  }{num [1:64] 0 1 0 0 0 1 0 0 0 1 ...}
  \item{CAI }{num [1:64]  0.00  0.00 -1.37 -2.98 -2.58 ...}
  \item{CG  }{num [1:4] 0 1 1 0}
  \item{CG1 }{num [1:64] 0 0 0 0 0 0 0 0 0 0 ...}
  \item{CG12}{num [1:64] 0 0 0 0 0.5 0.5 0.5 0.5 0.5 0.5 ...}
  \item{CG2 }{num [1:64] 0 0 0 0 1 1 1 1 1 1 ...}
  \item{CG3 }{num [1:64] 0 1 1 0 0 1 1 0 0 1 ...}
  \item{CGN }{num [1:64] 0 0 0 0 0 0 0 0 0 0 ...}
  \item{CGQ }{num [1:64] 0 0 0 0 0 1 1 0 0 0 ...}
  \item{CU  }{num [1:64] 0 1 0 1 0 0 0 0 0 0 ...}
  \item{DEN }{num [1:64] 0 1 0 1 0 0 0 0 0 0 ...}
  \item{ECM }{num [1:64] 1 1 0 1 0 1 0 1 0 0 ...}
  \item{ECO }{num [1:64] 17 29 17 29 71 21 71 21 143 68 ...}
  \item{F1  }{num [1:64]  1.026  0.239  1.026  0.239 -0.097 ...}
  \item{FD  }{num [1:64] 1 1 1 1 1 1 1 1 1 0 ...}
  \item{FN  }{num [1:64] 1 1 0 0 0 1 0 1 0 0 ...}
  \item{G   }{num [1:4] 0 0 1 0}
  \item{G3  }{num [1:64] 0 0 1 0 0 0 1 0 0 0 ...}
  \item{GP  }{num [1:65] 6.77 4.56 6.77 4.56 3.63 ...}
  \item{HI  }{num [1:64]  0.00  0.00 -2.91 -0.99 -1.03 ...}
  \item{KD  }{num [1:64] -3.9 -3.5 -3.9 -3.5 -0.7 -0.7 -0.7 -0.7 -4.5 -0.8 ...}
  \item{LEC }{num [1:64] 1 0 1 0 1 1 0 1 1 0 ...}
  \item{LEV }{num [1:64] 2.86 0 1.27 0 1.1 1.1 0 1.1 1.11 0 ...}
  \item{MG  }{num [1:64]  0.000 -0.740 -1.236  0.000 -0.168 ...}
  \item{MG2 }{num [1:64]  0.000 -0.740 -1.236  0.000 -0.168 ...}
  \item{NUM }{num [1:64] 0 1 0 0 0 0 0 0 0 0 ...}
  \item{OP  }{num [1:64] 0.17 0.29 0.17 0.29 0.21 0.21 0.21 0.21 0.19 0.42 ...}
  \item{Q   }{num [1:64] 0 0 0 0 1 1 1 1 0 0 ...}
  \item{Q3  }{num [1:64] 0 0 0 0 0 1 1 0 0 0 ...}
  \item{QA3 }{num [1:64] 0 0 0 0 1 0 0 0 0 0 ...}
  \item{QC3 }{num [1:64] 0 0 0 0 0 1 0 0 0 0 ...}
  \item{QT  }{num [1:64] 0 0 0 0 1 1 1 1 0 0 ...}
  \item{SUB }{num [1:64]  0.00  0.00 -1.99 -1.26  0.00 ...}
  \item{U   }{num [1:4] 0 0 0 1}
  \item{U3  }{num [1:64] 0 0 0 1 0 0 0 1 0 0 ...}
  }
}
\details{
Note that it's better to work directly directly at the amino-acid level
when computing linear forms on amino-acid frequencies so as to have
a single coefficient vector. For instance EXP$KD to compute the Kyte
and Doolittle hydrophaty index from codon frequencies is valid only
for the standard genetic code.
}
\source{
ANALSEQ EXPFILEs for command EXP.\cr
\url{http://biomserv.univ-lyon1.fr/doclogi/docanals/manuel.html}
}
\references{
  To be completed...
\describe{
  \item{AGZ }{ ? }
  \item{ARG }{ ? }
  \item{AU3 }{ ? }
  \item{BC  }{ ? }
  \item{CAI }{ ? }
  \item{CU  }{ ? }
  \item{DEN }{ ? }
  \item{ECM }{ ? }
  \item{ECO }{ ? }
  \item{F1  }{From Table 2 in Lobry, J.R., Gautier, C. (1994) Hydrophobicity,
expressivity and aromaticity are the major trends of amino-acid usage in
999 \emph{Escherichia coli} chromosome-encode genes. \emph{Nucleic Acids
Research},\bold{22}:3174-3180.}
  \item{FD  }{ ? }
  \item{FN  }{ ? }
  \item{GP  }{ ? }
  \item{HI  }{ ? }
  \item{KD  }{Kyte, J., Doolittle, R.F. (1982) A simple method for displaying
the hydropathic character of a protein. \emph{J. Mol. Biol.},\bold{157}
:105-132.}
  \item{LEC }{ ? }
  \item{LEV }{ ? }
  \item{MG  }{ ? }
  \item{MG2 }{ ? }
  \item{NUM }{ ? }
  \item{OP  }{ ? }
  \item{SUB }{ ? }
  }
}
\examples{
data(EXP)
}
\keyword{datasets}
