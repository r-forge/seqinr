read.fasta <- function(File = system.file("sequences/ct.fasta", package = "seqinr"))
{
  seq <- scan(File,what=character(),sep="\t") 
  ind <- c(grep(">",seq),length(seq)+1)
  nomseq <- as.list(c(rep(0,(length(ind)-1))))
  sequences <- as.list(rep(0,length(ind)-1))
  for(i in 1:(length(ind)-1))
  {
    nomseq[[i]]<-unlist(strsplit(seq[ind[i]]," "))[1]
    nomseq[[i]]=substr(nomseq[[i]],2,nchar(nomseq))	
    sequences[[i]]<-unlist(strsplit(seq[(ind[i]+1):(ind[i+1]-1)],split=NULL))
    attributes(sequences[[i]])=list(name=nomseq[[i]])
  }
  names(sequences)=nomseq
  if(deftype(sequences[[1]])=="AA") lapply(sequences,initSeqFastaAA) 
  else lapply(sequences,initSeqFastadna)
}


