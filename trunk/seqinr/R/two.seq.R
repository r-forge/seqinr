two.seq<-function(seq,phase=0){
	seq=toupper(unlist(strsplit(seq,"")))
	a<-seq(phase+1,length(seq),2)
	b<-seq(phase+2,length(seq),2)	
	paste(seq[a],seq[b],sep="")
}

