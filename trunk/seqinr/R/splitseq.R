splitseq<-function(seq,frame=0,word=3){
	if(word==1) seq[frame:length(seq)]
	l<-(floor((length(seq)-frame)/word)*word)+frame
	sq1=seq[seq(frame+1,l,word)]
	for(i in 2:word) sq1=paste(sq1,seq[seq(frame+i,l,word)],sep="")
	sq1
}

