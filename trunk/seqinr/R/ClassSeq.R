	#####################################################################################
	# 		classes de séquences						    #
	#  toutes les classes doivent avoir exactement la même interface à savoir:	    #
	#  une fonction initNomClasse qui retourne une instance de la classe		    #
	#  des spécialisations des fonctions                                                #
	#	     getSequence(seq) retourne un vecteur de char	     		    #
	#            getFrag(seq,begin,end) retourne un vecteur de caractères               #
	#            getLength(seq) retourne un "entier"                                    #
	#            getName(seq) retourne une chaîne                                       #
	#            getProp(seq) retourne une liste nommée                                 #
	#	     getAnnot(seq,nl) reourne un vecteur de string                          #
	#	     Translate(seq) retourne un vecteur de char                             #
	#####################################################################################


getSequence.default = function(object){
	if(length(object) == 1) object=s2c(object)	
 	xx = tolower(object)
 	if(length(grep("[acgtu]",xx)) != length(xx)) stop("Biological sequence is needed !")	
 	else return(xx)		
}

getFrag.default = function(object,begin,end){ 
	if(length(object) == 1) object=s2c(object)	
 	xx = tolower(object)
 	if(length(grep("[acgtu]",xx)) != length(xx)) stop("Biological sequence is needed !")
 	if(begin>length(xx) || end>length(xx) || begin>end) stop("borns are not correct")	
 	else return(xx[begin:end])		
}

getLength.default = function(object){
	if(length(object) == 1) object=s2c(object)	
 	xx = tolower(object)
 	if(length(grep("[acgtu]",xx)) != length(xx)) stop("Biological sequence is needed !")
 	return(length(xx))
}

getName.default = function(object){
 	stop("no name")
}

getProp.default = function(object){
 	stop("no property")
}

getAnnot.default = function(object,nbl){ 
 	stop("no annotation for this sequence")
}

Translate.default = function(seq,frame=0, sens= "F", numcode=1){
	translate(seq,frame,sens,numcode)
}

##################################################################

getFrag =  function(object,begin,end) {
	if(! inherits(object,c("SeqFastadna","SeqFastaAA","SeqAcnucLocal","SeqAcnucWeb","SeqFrag"))) { getFrag.default(object,begin,end) }
	else UseMethod("getFrag")
}

getSequence = function(object){
	if(! inherits(object,c("SeqFastadna","SeqFastaAA","SeqAcnucLocal","SeqAcnucWeb","SeqFrag"))) {getSequence.default(object)}
	else UseMethod("getSequence")
}


getLength =  function(object) {
	if(! inherits(object,c("SeqFastadna","SeqFastaAA","SeqAcnucLocal","SeqAcnucWeb","SeqFrag"))) {getLength.default(object)}
	else UseMethod("getLength")
}

getName =  function(object) {
	if(! inherits(object,c("SeqFastadna","SeqFastaAA","SeqAcnucLocal","SeqAcnucWeb","SeqFrag"))) {getName.default(object)}
	else UseMethod("getName")
}

getProp =  function(object) {
	if(! inherits(object,c("SeqFastadna","SeqFastaAA","SeqAcnucLocal","SeqAcnucWeb","SeqFrag"))) {getProp.default(object)}
	else UseMethod("getProp")
}


getAnnot = function(object,nbl) {
	if(! inherits(object,c("SeqFastadna","SeqFastaAA","SeqAcnucLocal","SeqAcnucWeb","SeqFrag"))) {getAnnot.default(object,nbl)}
	else UseMethod("getAnnot")
}

Translate = function(seq,frame=0, sens= "F", numcode=1){
	if(! inherits(seq,c("SeqFastadna","SeqFastaAA","SeqAcnucLocal","SeqAcnucWeb","SeqFrag"))) {Translate.default(seq,frame=0, sens= "F", numcode=1)}
	else UseMethod("Translate")
}



	########################################################################################################
	#		Classe de sequence SeqFastadna et ses méthodes:                                        #
	#	La classe de séquence SeqFasta pour les séquences résultants de la lecture d'un fichier au     # 
	#	format fasta.                                                                                  #
	########################################################################################################

	##################################################################################################
	# as.SeqFasta sera appelée au moment de la lecture d'un fichier au format fasta par read.fasta() #
	##################################################################################################

as.SeqFastadna = function(object){
	class(object)="SeqFastadna"	
        return(object)
        }

is.SeqFastadna = function(object){
	inherits(object,"SeqFastadna")
}

getSequence.SeqFastadna = function(object){
	return(object)
	}

getFrag.SeqFastadna = function(object, begin, end){
	if(end > getLength(object)) stop("invalid end")	
	newSeq = object[begin:end]
	newSeq = as.SeqFrag(newSeq,begin,end,compl=T,name=getName(object))
	return(newSeq)
	}

getLength.SeqFastadna = function(object){
	return(length(object))
	}

getName.SeqFastadna = function(object){
	return(attr(object,"name"))
}

getProp.SeqFastadna = function(object){
	return(list(seqtype="DNA"))
}

getAnnot.SeqFastadna = function(object){
	return(attr(object,"Annot"))
}

summary.SeqFastadna = function(object){
	compo=count(object,1)
	return(list(composition=compo,GC=GC(object)))
}

Translate.SeqFastadna =  function(object, frame=0, sens= "F", numcode=1){
	translate(object, frame=0, sens= "F", numcode=1)
}
	


	


	###############################################################################
	#		Classe de sequences SeqFastaAA et ses méthodes:               #
	###############################################################################

as.SeqFastaAA = function(object){
	class(object)="SeqFastaAA"	
        return(object)
        }

is.SeqFastaAA = function(object){
	inherits(object,"SeqFastaAA")
}

getSequence.SeqFastaAA = function(object){
	return(object)
	}


getFrag.SeqFastaAA = function(object, begin, end){
	if(end > getLength(object)) stop("invalid end")	
	newSeq = object[begin:end]
	newSeq = as.SeqFrag(newSeq,begin,end,compl=T,name=getName(object))
	return(newSeq)
	}


getLength.SeqFastaAA = function(object){
	return(length(object))
	}


getName.SeqFastaAA = function(object){
	return(attr(object,"name"))
}


getProp.SeqFastaAA = function(object){
	return(list(seqtype="AA"))
}

getAnnot.SeqFastaAA = function(object,nbl){
	return(attr(object,"Annot"))
}

summary.SeqFastaAA = function(object){
	compo=table(factor(object, levels = levels(SEQINR.UTIL$CODON.AA$L)))
	return(list(composition=compo/getLength(object),AA.Property=AAprop(object)))
}




	############################################################################################
	#	
	#		Classe de sequences SeqAcnucLocal et ses méthodes: 			   #
	# La classe de séquence SeqAcnucLocal pour les séquences provenant                         #
	# d' une requete dans les banques structurées sous ACNUC                                   #
	# Cette classe ne contiendra pas la sequence car le but est de ne pas la stocker.          #
	# On stockera le mnemo et le nom de la banque d'où elle provient.                          #
	#                                                                                          #
	############################################################################################


	###################################################################################################
	# l'initialisation se fera au moment de l'appel à la fonction getreq. On attribuera à chaque mnemo#
	# de la liste résultant de la requete la classe SeqAcnucLocal					  #	
	###################################################################################################



as.SeqAcnucLocal = function(object){
	class(object)="SeqAcnucLocal"
	return(object)
}


is.SeqAcnucLocal = function(object){

	inherits(object, "SeqAcnucLocal")
}


getFrag.SeqAcnucLocal = function(object,begin,end){
	b = getLength(object)
	if((end > b) || (begin > b)) stop("born out of limits")
	else{  
	s = .Call("getseq2",object,begin,(end-begin+1))
	seq = s2c(s)
	return(as.SeqFrag(seq,begin,end,compl=T,name=getName(object)))
	}
}


getSequence.SeqAcnucLocal = function(object){
	return(getseq(object,as.string=F))
	}

getName.SeqAcnucLocal = function(object){
	return(as.character(object))
}

getLength.SeqAcnucLocal = function(object){
	return(getAttribut(object)[[1]])
	}

getProp.SeqAcnucLocal = function(object){
	return(getAttribut(object)[2:3])
}

getAnnot.SeqAcnucLocal = function(object,nbl){
	return(getAnnots(object,nbl))
}


Translate.SeqAcnucLocal = function(object){
	seq = translateCDS(object)
	return(s2c(seq))
}


summary.SeqAcnucLocal = function(object){
	s=getSequence(object)
 	return(list(name=getName(object),GC.percent=GC(s),base.count=count(s,1)))
	}



	######################
	# fonctions annexes  #
	######################



AAprop = function(object){
	s=table(factor(object, levels = levels(SEQINR.UTIL$CODON.AA$L)))
	t=sum(s)
	list(Tiny=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Tiny)])/t,Small=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Small)])/t,Aliphatic=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Aliphatic)])/t,Aromatic=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Aromatic)])/t,Non.polar=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Non.polar)])/t,Polar=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Polar)])/t,Charged=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Charged)])/t,Basic=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Basic)])/t,Acidic=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Acidic)])/t)
}


####################################################################################################
#												   #
#	Classe de Sequences SeqAcnucWeb                                                            #
#												   #			
####################################################################################################





as.SeqAcnucWeb = function( object, socket=F ){

	class(object)="SeqAcnucWeb"
	attributes(object)=list(class="SeqAcnucWeb",socket=socket)
	object
}


is.SeqAcnucWeb = function( x ){	
	inherits(x,"SeqAcnucWeb")
}



getSequence.SeqAcnucWeb = function(){
	b=getLength( object )
	getsequence.socket(attr(object,"socket"),object,start=1,length=b)
}



getFrag.SeqAcnucWeb = function(object ,begin, end ){

	b = getLength(object)
	if((end > b) || (begin > b)) stop("born out of limits")  
	bb=end-begin+1
	newSeq = getsequence.socket(attr(object,"socket"),object,start=begin,length=bb)
	newSeq = as.SeqFrag(newSeq,begin=begin,end=end,compl=T,name=getName(object))
	return(newSeq)
}



getName.SeqAcnucWeb = function(object ){	

	return( as.character(object) )

}

getLength.SeqAcnucWeb = function( object ){

	return( getAttribut.socket(attr(object,"socket"),object)[[1]] )

}


getProp.SeqAcnucWeb = function(object){

	return( getAttribut.socket(attr(object,"socket"),object) ) 

}

getAnnot.SeqAcnucWeb = function(object, nbl ){
		
	return( readAnnots.socket( socket= attr(object,"socket"),name=object, nl=nbl) ) 

}


Translate.SeqAcnucWeb = function(object,frame=0, sens= "F", numcode=1){
	translate(object, frame=0, sens= "F", numcode=1)
}


	############################################################################
	#		Classe de sequences SeqFrag et ses méthodes:               #
	############################################################################




as.SeqFrag = function(object,begin,end,compl=F,name="frag"){
	if(compl){ attr(object,"seqMother") = name }
	else attr(object,"seqMother") = getName(seq)
        attr(object,"begin") = begin
	attr(object,"end") = end
	class(object) = "SeqFrag"
        return(object)
        }

is.SeqFrag = function(object){
	inherits(object,"SeqFrag")
}


getSequence.SeqFrag = function(object){
	return(object)
	}


getFrag.SeqFrag = function(object,begin,end){
        if((end<begin) || (end>getLength(object)))  stop("invalid end")
        newBegin = attr(object,"begin")+begin-1
        newEnd = attr(object,"begin")+end-1
	newSeq = object[begin:end]
        newSeq = as.SeqFrag(object=newSeq,begin=newBegin,end=newEnd,compl=T,name=getName(object))
	return(newSeq)
        }

getLength.SeqFrag = function(object){
	return(attr(object,"end")-(attr(object,"begin")+1))
}

getName.SeqFrag = function(object){
	return(attr(object,"seqMother"))
}

getprop.SeqFrag = function(object){
	return(list())
}


Translate.SeqFrag = function(seq, frame=0, sens= "F", numcode=1){
	translate(seq, frame=0, sens= "F", numcode=1)
}



	