uco<-function(seq,phase=0,freq=FALSE,as.data.frame=FALSE){
	sq<-splitseq(seq,phase,nombre=3)
	eff<-table(factor(sq,levels=CODON))
	if(freq == TRUE) eff<-round(eff/(floor(length(seq)/3)),4)
	if(as.data.frame == TRUE) eff=as.data.frame(cbind(as.character(AA),as.character(CODON),as.vector(eff)))
	eff
} 


