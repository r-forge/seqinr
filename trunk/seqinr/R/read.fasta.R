read.fasta <- function(File = system.file("sequences/ct.fasta", package = "seqinr"))
{
  seq <- scan(File,what=character(),sep="\t") 
  ind <- c(grep(">",seq),length(seq)+1)
  nomseq <- as.list(c(rep(0,(length(ind)-1))))
  sequences <- as.list(rep(0,length(ind)-1))
  lclass <- as.list(rep(0,length(ind)-1))
  for(i in 1:(length(ind)-1))
  {
    nomseq[[i]]<-unlist(strsplit(seq[ind[i]]," "))[1]
    sequences[[i]]<-unlist(strsplit(seq[(ind[i]+1):(ind[i+1]-1)],split=NULL))
  }
  nomseq = lapply(nomseq,substr,2,nchar(nomseq))
  names(sequences)=nomseq
  for(i in 1:length(sequences)) lclass[[i]]=sequences[i]
  lapply(lclass,initSeqFasta)
}


