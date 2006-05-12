\name{write.fasta}
\alias{write.fasta}
\title{ Write file in fasta format }
\description{
  Writes sequences to a file in Fasta format.
}
\usage{write.fasta(sequences,names,nbchar=60,file.out)}
\arguments{
  \item{sequences}{ A DNA or protein sequence (i.e. a
    character vector) or a list of  sequences (i.e. a
    list of character vectors). }
  \item{names}{ The name(s) of the sequences. }
  \item{nbchar}{ The number of characters per line (default: 60) }
  \item{file.out}{ The name of the output file. }
}
\value{
 NULL
}
\references{
  \code{citation("seqinr")}
}
\author{ Anamaria Necsulea }
\seealso{ \code{\link{read.fasta}} }
\examples{
# read sequences from a Fasta file #
s = read.fasta(File=system.file("sequences/ortho.fasta",package =
"seqinr"))
# select only third codon positions #
t = lapply(s, function(x) x[seq(from=3,to=length(x), by=3)])
# write the modified sequences to a file #
write.fasta(sequences=t, names=names(s),nbchar=80,file.out="test.fasta")
}
\keyword{ utilities }
