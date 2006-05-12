write.fasta<-function(sequences, names, nbchar=60, file.out){

  write.oneseq<-function(sequence, name, nbchar, file.out){
    outfile=file(description=file.out,open="a")
    writeLines(paste(">",name,sep=""),outfile)
    l=length(sequence)
    q=floor(l/nbchar)
    r=l-nbchar*q
    if(q>0){
      sapply(1:q, function(x) writeLines(paste(sequence[(nbchar*(x-1)+1):(nbchar*x)],collapse="",sep=""),outfile))
      
    }
    if(r>0){
      writeLines(paste(sequence[(nbchar*q+1):l],collapse="",sep=""),outfile)
    }
    close(outfile)
  }
  
  if(!is.list(sequences)){
    write.oneseq(sequence=sequences, name=names, nbchar=nbchar, file.out=file.out)
  }
  else{
    sapply(1:length(sequences), function(x) write.oneseq(sequence=as.character(sequences[[x]]), name=names[x], nbchar=nbchar, file.out=file.out))
   }
  return(NULL);
}
