\name{SEQINR.UTIL}
\alias{SEQINR.UTIL}
\docType{data}
\title{some utility data for seqinr}
\description{
  This data set gives the genetics code, the name of each codon, the IUPAC one-letter code for aminoacids and the physico-chemical class of amino acid.                                                                                                              
}
\usage{data(SEQINR.UTIL)}
\format{
 	\code{SEQINR.UTIL} is a list containing the 3 following objects:  
	\enumerate{	
 	\item{CODES.NCBI}{is a data frame containing the genetics code : The standard ('Universal') genetic code with a selection of non-standard codes. } 
 	\item{CODON.AA}{is a three columns data frame. The first column is a factor containing the codon. The second column is a factor giving the aminoacids names for each codon. The last column is a factor giving the IUPAC one-letter code for aminoacids}
 	\item{AA.PROPERTY}{is a list giving the physico-chemical class of amino acid. The differents classes are the following one : Tiny, Small, Aliphatic, Aromatic, Non.polar, Polar, Charged, Basic, Acidic }
 } 
}
\source{
  Data prepared by D.Charif \email{charif@biomserv.univ-lyon1.fr}.\cr The genetic codes have been taken from the ncbi taxonomy database: \url{http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c }. Last update October 05, 2000.\cr The IUPAC one-letter code for aminoacids is descibed at: \url{http://www.chem.qmul.ac.uk/iupac/AminoAcid/ }. 
}
\examples{
data(SEQINR.UTIL)
}
\keyword{datasets}
