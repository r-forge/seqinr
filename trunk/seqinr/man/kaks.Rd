\name{kaks}
\alias{kaks}
\title{ to Get an Estimation of Ka and Ks }
\description{
  Ks and Ka  are respectively the number of substitutions per synonymous site and per nonsynonymous site between two protein-coding genes. The ratio of nonsynonymous (Ka) to synonymous (Ks) nucleotide substitution rates is an indicator of selective pressures on genes. A ratio significantly greater than 1 indicates positive selective pressure. A ratio around 1 indicates either neutral evolution at the protein level or an averaging of sites under positive and negative selective pressures. A ratio less than 1 indicates pressures to conserve protein sequence. This function compute the Ka and Ks values for a set of aligned sequences and gives the associated variance matrix. 
}
\usage{
kaks(x)
}
\arguments{
  \item{x}{ An object of class \code{alignment} }
}
\details{
 
}
\value{
  \item{ ks }{ matrix of ks values }
  \item{ ka }{ matrix of ka values }
  \item{ vks }{ variance matrix of ks }
  \item{ vka }{ variance matrix of ka }	
}
\references{
Li WH. (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. 
\emph{J Mol Evol.}, \bold{Jan};36(1):96-9.\cr 
The C programm implementing this method was provided by M. Gouy.
}
\author{ D. Charif }
\seealso{\code{\link{read.alignment}} }
\examples{
 s = read.alignment(File=system.file("sequences/test.phylip", package = "seqinr"),format="phylip")
 kaks(s)	
}
\keyword{ manip }
