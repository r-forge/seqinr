count<-function(seq,nombre){
	seq=toupper(unlist(strsplit(seq,"")))
	if (nombre==1){
	table(factor(seq,levels=levels(as.factor(alphabet(1)))))
}
	else if (nombre==2){
	s<-paste(seq,seq[2:length(seq)],sep="")
	table(factor(s,levels=levels(as.factor(alphabet(2)))))
}
	else if (nombre==3){
	s<-paste(seq,seq[2:length(seq)],seq[3:length(seq)],sep="")
	table(factor(s,levels=levels(as.factor(alphabet(3)))))
}
	else if (nombre==4){
	s<-paste(seq,seq[2:length(seq)],seq[3:length(seq)],seq[4:length(seq)],sep="")
	table(factor(s,levels=levels(as.factor(alphabet(4)))))
}
	else if (nombre==5){
	s<-paste(seq,seq[2:length(seq)],seq[3:length(seq)],seq[4:length(seq)],seq[5:length(seq)],sep="")
	table(factor(s,levels=levels(as.factor(alphabet(5)))))
}
}
