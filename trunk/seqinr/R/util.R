########################
# char to string
########################

c2s <- function( chars = c("m","e","r","g","e","d") )
{
  return( paste( chars, collapse = "" ) )
}

###########################
# string to char
############################

s2c <- function(string){
	return(.Call("s2c",string))
}


###########################
# Conversion of the numeric encoding of a DNA sequence into
# a vector of chars
############################

n2s <- function(nseq, levels = c("a", "c", "g", "t"), base4 = TRUE)
{
  if( base4 )
    levels[nseq + 1]
  else
    levels[nseq]
}

###############################
# simple numerical encoding of a DNA sequence that by default
# is independent of locale.
###############################

s2n <- function(seq, levels = c("a", "c", "g", "t"), base4 = TRUE, ... )
{
  if( base4 )
    codes(factor(seq, levels = levels , ...) ) - 1
  else
    codes(factor(seq, levels = levels , ...) )
}

################################
# GC.percent
#################################

GC = function(seq)
{
	c=count(seq,1)
	cc=(c[2]+c[3])/sum(c)
	return(as.vector(cc))
}

plotGC=function(sequence,fenetre){
	m=matrix(c(1:2),2,1)
	layout(m)
	l=as.list(splitseq(sequence,word=fenetre))
	ll=lapply(l,s2c)
	lll=lapply(ll,GC)
	ss=unlist(lll)
	hist(ss,xlab=("taux de GC"),main="histogramme du taux de GC")
	s=scale(ss)
	x=seq(0,length(ss)*fenetre-fenetre,fenetre)
	plot(x,s,type="n",xlab="bases",ylab="taux de GC",main="Evolution du taux de GC le long d'une sequence")
	abline(h=0)
	y2=s[which(s>0)]
	x2=x[which(s>0)]
	y3=s[which(s<0)]
	x3=x[which(s<0)]
	lines(x2,y2,col="green3",type="h")
	lines(x3,y3,col="violet",type="h")	
	box(lty='137', col = 'red')
}


##########################################
# Conversion one-letter code to 3-letters code for amino-acids
##########################################

aaa <- function( aa )
{
  aa1 <- s2c("*ACDEFGHIKLMNPQRSTVWY")
  aa3 <- c("Stp", "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile",
           "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr",
           "Val", "Trp", "Tyr")
  convert <- function( x )
  {
    if( all( x != aa1 ) )
    { 
      warning("Unknown one letter code for aminoacid")
      return( NA )
    }
    else
    {
      return( aa3[which( x == aa1 )] )
    }
  }
  return( as.vector(unlist(sapply( aa, convert ) ) ) )
}


#########################################
# revers a sequence
#######################################

invers<-function(seq)
{
	seq=toupper(s2c)
	rev(seq)
}

##########################################
#complement a sequences
###########################################

comp<-function(seq){
	seq<-replace(seq,seq=="a","tmp1")
	seq<-replace(seq,seq=="c","tmp2")
	seq<-replace(seq,seq=="g","c")
	seq<-replace(seq,seq=="t","a")
	seq<-replace(seq,seq=="tmp1","t")
	seq<-replace(seq,seq=="tmp2","g")
	seq
}


#########################################################
# donne tous les mots de longueur x sur l'alphabet {a,c,g,t} 
#########################################################

alphabet<-function(al)
{
al1<-s2c("acgt")
if(al==1) al1
else
{
alx<-array(rep(0,4^(al)),rep(4,al))
alx<-array(c(paste(alphabet(al-1),array("a",rep(4,(al-1))),sep=""),paste(alphabet(al-1),array("c",rep(4,(al-1))),sep=""),paste(alphabet(al-1),array("g",rep(4,(al-1))),sep=""),paste(alphabet(al-1),array("t",rep(4,(al-1))),sep="")),rep(4,al),dimnames=as.list(rep(list(al1),al))
)
alx
}
}

