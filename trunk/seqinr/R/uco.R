uco<-function( seq, frame = 0, freq = FALSE, as.data.frame = FALSE){
	sq<-splitseq(seq,frame)
	eff<-table(factor(sq,levels=CODON))
	if(freq == TRUE) eff<-round(eff/(floor(length(seq)/3)),4)
	if(as.data.frame == TRUE){
	eff=as.data.frame(cbind(as.character(AA),as.character(CODON),as.vector(eff)))
	if(freq == TRUE) names(eff)=c("aa","codon","freq")
	else names(eff)=c("aa","codon","eff")
	}
	eff
} 


