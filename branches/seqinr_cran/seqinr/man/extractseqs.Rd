\name{extractseqs}
\alias{extractseqs}
\alias{exseq}
\title{To extract the sequences information of a sequence or a list of sequence in different formats}
\description{
To extract the sequences information of a sequence or a list of sequence in different formats
}

\usage{
extractseqs(lrankseqnum,socket = "auto", format="fasta",operation="simple",zlib = FALSE, verbose = FALSE, npaquets = -1, nzlines=1000)
exseq(lrankseqnum,socket = "auto", format="fasta",operation="simple",zlib = FALSE, verbose = FALSE, npaquets = -1, nzlines=1000)
}

\arguments{
  \item{lrankseqnum}{the name of list on server, or the name of a sequence}
  \item{socket}{a socket of class connection and sockconn returned by \code{choosebank}. 
   Default value (auto) means
   that the socket will be set to to the socket component of the banknameSocket variable. }
  \item{format}{the format of output.Can be \code{acnuc}, \code{fasta}, or \code{flat} }
  \item{operation}{the type of extraction. Can be  \code{simple}, \code{translate}, \code{fragment}, \code{feature} or \code{region}  }
  \item{zlib}{if \code{TRUE}, compression mode is on}
  \item{verbose}{if \code{TRUE}, verbose mode is on}
  \item{npaquets}{number of paquets}
  \item{nzlines}{number of line in zlib mode}
}

\details{
 Under development.
}

}
\value{
  Sequence data.
}
\references{ 
\code{citation("seqinr")}
}
\author{S.Penel}
\seealso{
  \code{\link{choosebank}}, 
  \code{\link{query}}
   \code{\link{getlistrank}}
}
}
\examples{
 \dontrun{
 # Need internet connection
 choosebank("swissprot")
 query("MyListName", "k=globin")
 MyList.fasta <- exseq("MyListName", verbose=T)
 #
 
 }
}
\keyword{utilities}
