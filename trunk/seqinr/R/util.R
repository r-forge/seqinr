#
# char to string
#
c2s <- function( chars = c("m","e","r","g","e","d") )
{
  return( paste( chars, collapse = "" ) )
}
#
# Conversion of the numeric encoding of a DNA sequence into
# a vector of chars
#
n2s <- function(nseq, levels = c("a", "c", "g", "t"), base4 = TRUE)
{
  if( base4 )
    levels[nseq + 1]
  else
    levels[nseq]
}
#
# simple numerical encoding of a DNA sequence that by default
# is independent of locale.
#
s2n <- function(seq, levels = c("a", "c", "g", "t"), base4 = TRUE, ... )
{
  if( base4 )
    codes(factor(seq, levels = levels , ...) ) - 1
  else
    codes(factor(seq, levels = levels , ...) )
}
GC = function(seq)
{
c=count(seq,1)
cc=(c[2]+c[3])/sum(c)
return(cc)
}#
# Conversion one-letter code to 3-letters code for amino-acids
#
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

invers<-function(seq)
{
	seq=toupper(s2c)
	rev(seq)
}

comp<-function(seq){
	seq<-replace(seq,seq=="A","tmp1")
	seq<-replace(seq,seq=="C","tmp2")
	seq<-replace(seq,seq=="G","C")
	seq<-replace(seq,seq=="T","A")
	seq<-replace(seq,seq=="tmp1","T")
	seq<-replace(seq,seq=="tmp2","G")
	seq
}

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


