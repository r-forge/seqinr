\name{mkbank}
\alias{mkbank}

\title{ To create your own ACNUC data base }
\description{
  \code{mkbank} create an ACNUC data base which can be request using the \code{choosebank} and \code{query} functions.
}
\usage{
mkbank(bankname, fileseq, bank="genbank", format="flatfiles")
}
\arguments{
  \item{bankname}{is the name of the bank to create}
  \item{fileseq}{A quoted string of chars giving the name of the file which contains the sequences. This file must be in the current working directory. \cr
  If the sequences come from:\cr
  -genbank, the .seq extension is needed \emph{example : toto.seq}\cr
  -embl or swissprot, the .dat extension is needed \emph{example : toto.dat}\cr
}
  \item{bank}{is a quoted string giving the bank from which come the sequences: \code{genbank}, \code{embl}, \code{swissprot} or \code{nbrf}}
  \item{format}{is the format of the sequences file(s). If ACNUC is interfaced with the flat-files as distributed by genbank, embl or swissprot, the parameters \code{flatfiles} must be choosen. Otherwise if ACNUC is interfaced with the GCG package, the parameters \code{gcg} must be choosen.}
}
\details{
  \code{mkbank} will first create two directories in the current working directory: \cr
  The directory \file{flatfiles} which contains the file(s) of sequences.\cr 
  The directory \file{index} which contains the index files of the ACNUC data base. \cr 

Secondly, the initf program will create an empty ACNUC data base. The following index files are created in the directory \file{index}:\cr\cr
BIBLIO\cr
ACCESS\cr 
LOCUS\cr
SUBSEQ\cr
AUTHOR\cr 
LONGL\cr 
SHORTL\cr 
SMJYT\cr
SPECIES\cr
TEXT\cr
MERES\cr
EXTRACT\cr

Thirdly, the gbemgener program will generate the ACNUC data base and fill the index files.
}

\value{
  An ACNUC data base is created in the current working directory. During the programs execution, some comments will appear.    
}

\references{\url{http://pbil.univ-lyon1.fr/databases/acnuc.html}\cr
Gouy, M., Milleret, F., Mugnier, C., Jacobzone, M., Gautier,C. (1984) ACNUC: a nucleic acid sequence data base and analysis system. 
\emph{Nucl. Acids Res.}, \bold{12}:121-127.\cr }
\author{ D. CHARIF }
\note{ Most of the documentation was imported from ACNUC help files written by Manolo Gouy }

\seealso{\code{\link{query}},\code{\link{choosebank}}}
\examples{

# Warning: you will create a database structured under ACNUC in your working directory. 	

  \dontrun{fileseq = system.file("sequences/felis.seq", package = "seqinr")}  
 \dontrun{mkbank("felis",fileseq)}
  \dontrun{choosebank("felis")}
  \dontrun{toto=2}

}
\keyword{ utilities }

