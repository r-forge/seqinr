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

s2n <- function(seq, levels = c("a", "c", "g", "t"), base4 = TRUE)
{
  if( base4 )
    unclass(factor(seq, levels = levels ) ) - 1
  else
    unclass(factor(seq, levels = levels ) )
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

######################
# GC3		     #
######################

GC3 = function(seq){
	sequence <- splitseq( seq, 3)
	codons <- words(length = 3, alphabet = s2c("acgt"))	
	eff=table(factor( sequence , levels=codons))
	f =round(eff/(floor(length(seq)/3)),4)
	return(as.vector(f %*% EXP$CG3))
	}

######################
# GC2		     #
######################


GC2 = function(seq){
	sequence <- splitseq( seq, 3)
	codons <- words(length = 3, alphabet = s2c("acgt"))	
	eff=table(factor( sequence , levels=codons))
	f =round(eff/(floor(length(seq)/3)),4)
	return(as.vector(f %*% EXP$CG2))
	}
