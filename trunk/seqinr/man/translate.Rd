\name{translate}
\alias{translate}
\title{ Translate nucleic acid sequences }
\description{
   This function translate nucleic acid sequences to the corresponding peptide sequence. It can translate in any of the 3 forward or three reverse sense frames. It can translate using the standard (Universal) genetic code and also with a selection of non-standard codes. 
}
\usage{
translate(seq, frame = 0, sens = "F", numcode = 1)
}
\arguments{
  \item{seq}{ an object of class \code{seq}. }
  \item{frame}{ Frame(s) (0,1,2) to translate. By default the frame \code{0} is used. }
  \item{sens}{ Sense to translate: \code{F} for forward sense and \code{R} for revers sense. }
  \item{numcode}{ The number of the code to use. By default the standard genetic code is used. }
}
}

\details{

 The following genetic codes are described here. The number preceding each code corresponds to \code{numcode}. 
\arguments{	
  \item{1}  standard
  \item{2}  vertebrate.mitochondrial
  \item{3}  yeast.mitochondrial
  \item{4}  protozoan.mitochondrial+mycoplasma
  \item{5}  invertebrate.mitochondrial
  \item{6}  ciliate+dasycladaceal
  \item{9}  echinoderm+flatworm.mitochondrial
  \item{10} euplotid
  \item{11} bacterial+plantplastid
  \item{12} alternativeyeast
  \item{13} ascidian.mitochondrial
  \item{14} alternativeflatworm.mitochondrial
  \item{15} blepharism
  \item{16} chlorophycean.mitochondrial
  \item{21} trematode.mitochondrial
  \item{22} scenedesmus.mitochondrial
  \item{23} hraustochytrium.mitochondria
} 
}
\value{
  \code{translate} returns a vector of character containig the peptide sequence in the standard one-letter IUPAC code. Termination (STOP) codons are translated as the character '*'.
}
\references{
    The genetic codes have been taken from the ncbi taxonomy database: \url{ http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c }. Last update October 05, 2000.}
 
\author{ D. CHARIF }

\examples{

a="atggatgccaaggtcgtcgccgtgctggccctcgtgctggctgcgctctgcctcagcgacgggaaaccggtcagcctgagctacagatgtccttgccgattcttcgagagccacgttgccagagccaacgtcaagcatctcaaaatcctcaacactcccaactgtgcccttcagatcgtggcgaggctgaagaacaacaatagacaagtgtgcattgacccgaaattgaagtggattcaggaatacctggagaaagctttaaacaagaggttcaagatgtga"
a=toupper(s2c(a))
translate(a)
translate(a,3,"R",6)

\keyword{ manip }

