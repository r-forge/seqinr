uco.freq<-function(seq,phase=0){
	seq=toupper(unlist(strsplit(seq,"")))
	sq<-codon.seq(seq,phase)
	l<-floor(length(seq)/3)
	freq<-round(table(factor(sq,levels=codon))/l,4)
	t<-as.data.frame(cbind(as.character(AA),as.character(codon),as.vector(freq)))
}

uco.table<-function(seq,phase=0){
	seq=toupper(unlist(strsplit(seq,"")))
	sq<-codon.seq(seq,phase)
	eff<-table(factor(sq,levels=codon))
	as.data.frame(cbind(as.character(AA),as.character(codon),as.vector(eff)))
}

uco<-function(seq,phase=0){
	seq=toupper(unlist(strsplit(seq,"")))
	sq<-codon.seq(seq,phase)
	eff<-table(factor(sq,levels=codon))
}
