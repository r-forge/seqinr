count<-function(seq,nombre){
	s=seq
	if (nombre==1){
	table(factor(seq,levels=levels(as.factor(alphabet(1)))))
	}
	else{
	for(i in 2:nombre) s=paste(s,seq[i:length(seq)],sep="")
	table(factor(s,levels=levels(as.factor(alphabet(nombre)))))
	}
    }