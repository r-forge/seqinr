comp<-function(seq){
	seq<-replace(seq,seq=="A","tmp1")
	seq<-replace(seq,seq=="C","tmp2")
	seq<-replace(seq,seq=="G","C")
	seq<-replace(seq,seq=="T","A")
	seq<-replace(seq,seq=="tmp1","T")
	seq<-replace(seq,seq=="tmp2","G")
	seq
}

