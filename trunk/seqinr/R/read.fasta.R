read.fasta <- function(Fichier = system.file("sequences/ct.fasta", package = "package"))
{
  seq <- scan(Fichier,what=character(),sep="\t") 
  ind <- c(grep(">",seq),length(seq)+1)
  nomseq <- as.list(c(rep(0,(length(ind)-1))))
  sequences <- as.list(rep(0,length(ind)-1))
  for(i in 1:(length(ind)-1))
  {
    nomseq[[i]]<-unlist(strsplit(seq[ind[i]]," "))[1]
    sequences[[i]]<-unlist(strsplit(seq[(ind[i]+1):(ind[i+1]-1)],split=NULL))
  }
  names(sequences)<-nomseq
  sequences
}


