alphabet<-function(al)
{
al1<-c("A","C","G","T")
if(al==1) al1
else
{
alx<-array(rep(0,4^(al)),rep(4,al))
alx<-array(c(paste(alphabet(al-1),array("A",rep(4,(al-1))),sep=""),paste(alphabet(al-1),array("C",rep(4,(al-1))),sep=""),paste(alphabet(al-1),array("G",rep(4,(al-1))),sep=""),paste(alphabet(al-1),array("T",rep(4,(al-1))),sep="")),rep(4,al),dimnames=as.list(rep(list(al1),al))
)
alx
}
}

