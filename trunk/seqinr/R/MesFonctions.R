# BEWARE functions whith name ending in "C" must never
# be used. Their existence or prototype is in no way warranted in futur
# versions



# transforme une base en numero
# si le caractere est différent de acgt retourne -1

 base2num = function(x) {
r=grep(x,c("a","c","g","t"))
if(length(r)==0) r=0;
return(r-1)
}

#transforme un numero en base

num2base = function(n) {
if(any( 0:3==n)) return(c("a","c","g","t")[n+1])
return("n")
}

#transforme un vecteur de bases en vecteur numerique

seq2num = function(x) {
n=length(x)
r=1:n
for(i in 1:n) r[i]=base2num(x[i])
return(r)
}

#alphabet pour une longeur donné (auteur Delphine Charif)

alphabet<-function(al)
{
al1<-c("a","c","g","t")
if(al==1) al1
else
{
alx<-array(rep(0,4^(al)),rep(4,al))
alx<-array(c(paste(alphabet(al-1),array("a",rep(4,(al-1))),sep=""),paste(alphabet(al-1),array("c",rep(4,(al-1))),sep=""),paste(alphabet(al-1),array("g",rep(4,(al-1))),sep=""),paste(alphabet(al-1),array("t",rep(4,(al-1))),sep="")),rep(4,al),dimnames=as.list(rep(list(al1),al))
)
alx
}
}

# base4 donne la décomposition en base 4 d'un nombre retourné dans un vecteur de longueur l

base4 = function(x,l){
r=rep(0,l)
k=l
while(x>0) {
r[k]=x%%4
k=k-1
x=x%/%4
}
return(r)
}

# rbase4 fonction réciproque de base4

rbase4 = function(x) {
        r=0
        b=1
        for(i in length(x):1){
                r=r+b*x[i]
                b=b*4
                }
        return(r)
        }


#motif donne le motif de longueur l et de numero n

motif=function(n,l){
v=base4(n,l)
s=""
for(i in 1:l)s=paste(s,num2base(v[i]),sep="",collapse="")
return(s)
}


# frequences des motifs de longueur donné l, avec un pas de p, en debutant à d
# la sequence est sous forme d'un vecteur numérique
# permet de calculer:
#               la fréquence des bases
#               la fréquence des base en positin III des codons
#               l'usage du code ....

freqC = function(s,length,step,begin) { # MUST NOT BE CALLED
        colonne=list()
        x=rep(0,4^length)
        debuts=s(begin,length(s)-length,by=step)
        for(i in debuts){
                n=rbase4(seq[i:(i+l-1)])
                x[n+1]=x[n+1]+1 # R numerote betement a partir de 1 et non de 0
                }
        for(i in 1:(4^length)){
                colonne=c(colonne,motif(i-1,l))
                }
        names(x)=colonne
        return(as.data.frame(t(x)))
        }
        
        
#la meme mais utilisant une séquence

freqSeq = function(seq,length,step,begin){
        x=freqC(seq2num(getFrag(seq)),length,step,begin)
        rownames(x)=getName(seq)
        return(x)
        }
        
#la meme mais utilisant une liste de séquences


freqLSeq = function(listSeq,length,step,begin) {
        n=length(listSeq)
        x=freq(seq2num(getFrag(listSeq[[1]])),length,step,begin)
        rownames(x)=getName(listSeq[[1]])
        for(i in 2:n){
                 x1=freq(seq2num(getFrag(listSeq[[i]])),length,step,begin)
                 rownames(x1)=getName(listSeq[[i]])
                 x=rbind(x,x1)
                 }
        return(x)
        }

######################################################################
#classes de séquences
# toutes les classes doivent aoir exactement la même interface à savoir:
# une fonction initNomClasse qui retourne une instance de la classe
# des spécialisations des fonctions
#                    getFrag(seq,begin=0,end=getLength(seq)) retourne un tableau de caractères
#                    getLength(seq) retourne un "entier"
#                    getName(seq) retourne une chaîne
#                    getProp(seq) retourne une liste nommée
######################################################################
getFragC=function(x,y,z){return(NULL)}
getLengthC=function(x){return(0)}
getNameC=function(x){return(NULL)}
getPropC=function(x){return(list())}

getFrag =  function(x,y,z) {
if(is.null(attr(x,"class"))) {getFragC(x,y,z)}
else UseMethod("getFrag")
}

getLength =  function(x,y,z) {
if(is.null(attr(x,"class"))) {getLengthC(x,y,z)}
else UseMethod("getLength")
}

getName =  function(x,y,z) {
if(is.null(attr(x,"class"))) {getNameC(x,y,z)}
else UseMethod("getName")
}

getProp =  function(x,y,z) {
if(is.null(attr(x,"class"))) {getPropC(x,y,z)}
else UseMethod("getProp")
}



######################################
# Un exemple très simple (mais pas inutile!!!)
#   Classe de sequence SeqSimple
####################################

initSeqSimple = function(s,name="seq") {
        r=list(seq=s,name=name)
        class(r)="SeqSimple"
        return(r)
        }

getFrag.SeqSimple=function(seq,begin=1,end=getLength(seq)){
        if(end<begin)return(NULL)
        return(s2c(substr(seq$seq,begin,end)))
        }
        
getLength.SeqSimple = function(seq){return(nchar(seq$seq))}

getName.SeqSimple = function(seq){return(seq$name)}

getprop.SeqSimple = function(seq){return(list())}

#######################
#fin de la classe SeqSimple
#######################

########################################
# Ce qui touche aux sous-sequences
#  Classe SeqFrag : un fragment d'une séquence
#  Classe SeqDaughter : la concaténation de plusieurs séquences (par exemple des fragments)
########################################
#BEWARE seq must be a character string representing the name of the variable
initSeqFrag=function(seq,begin,end,name="Frag",compl=FALSE,num=0){
        if(compl){
                nom=paste(name,num,sep="")
                }
        else nom=name
        seq=as.symbol(seq)
        r=list(seqmother=seq,begin=begin,end=end,name=nom)
        class(r)="SeqFrag"
        return(r)
        }
getFrag.SeqFrag=function(seq,begin=1,end=getLength(seq)) {
        if(end<begin)return(NULL)
        s=eval(as.symbol(seq$seqmother),envir=globalenv())
        bm=seq$begin-1
        b=bm+begin
        em=bm+end
        e=min(em,s$end)
        return(getFrag(s,b,e))
        }

getLength.SeqFrag = function(seq){return(seq$end-seq$begin+1)}

getName.SeqFrag = function(seq){return(seq$name)}

getprop.SeqFrag = function(seq){return(list())}

############ Fin classe SeqFrag  ####################"

seqDC=function(seq,i){ # MUST NEVER BE CALLED
        if(class(seq)!="SeqDaughter") return(NULL)
        if(i<=seq$nbFrag)return(seq[1][[1]][i])
        else return(NULL)
        }

initSeqDaughter=function(name,...) {
        l=list(...)
        l=c(l,name=name)
        nbFrag=length(l[1][[1]])
        l=c(l,nbFrag=nbFrag)
        class(l)="SeqDaugther"
        length=0
        for(i in 1:l$nbFrag) length=length+getLength(SeqDC(l,i))
        l=c(l,length=length)
        return(l)
        }

getLength.SeqDaugther = function(seq){return(seq$length)}

getName.seqDaugther = function(seq){return(seq$name)}

getFrag.seqDaugther = function(seq,begin=1,end=getLength(seq)){
        if(end<begin)return(NULL)
        if(begin>getLength(seq))return(NULL)
        if(end>getLength(seq))end=getLength(seq)
        if(begin<1)begin=1
        k=0
# looking for first implied fragment
        i=1
        for(i in 1:seq$nbFrag) {
                k=k+getLength(seqDC(seq,i))
                if(begin<k)break
                }
        b=getLength(seqDC(seq,i))-k+begin #the begin in the current frag
        e=getLength(seqDC(seq,i))-k+end
        e1=min(e,getLength(seqDC(seq,i))) #the end in the current frag
        frag=getFrag(seqDC(seq,i),b,e1)
        e=e-e1
        while(e>0){
                i=i+1
                e1=min(e,getLength(seqDC(seq,i)))
                frag=c(frag,getFrag(seqDC(seq,i),1,e1))
                e=e-e1
                }
        return(frag)
        }
############################ end Class SeqDaughter ###########################



# La plupart des calculs sur l'usage des codons et des bases sont des formes lineaires
# de frequences de motifs. Les fonction suivantes ont pour objectif de gerer cette situation
# dans le contexte d'un resultat de freqLSeq. 

####################################################################
# une classe Forme lineaire (fl):

initFl = function(coef,name="fl",comment="") {
      l=list(coef=coef,name=name,comment=comment)
      class(l)="fl"
      return(l)
      }
# des exemples utiles
GCb=initFl(c(0,1,1,0),name="C+G",comment="compute G+C from bases frequencies")
GC3=rep(c(0,1,1,0),16)
GC3=initFl(GC3,name="CGIII",comment="compute G+C frequencies in position III of codons")
Len=initFl(rep(1,64),name="length")
# la fonction de base calcFl
# premier parametre une forme lineaire ou une liste de formes lineaires
# deuxieme parametre une data.frame
#retour une data.frame de meme lignes que la precedante et comme colonne la ou les
#formes lineaires.
 ErrorFl1="Error in calcFl: bad parameter type"
calcFl = function(fl,x){
      if(class(fl)=="fl") return(calcFlC(fl,x))
      if(!is.list(fl)) return(ErrorFl1)
      r=calcFlC(fl[[1]],x)
      for(i in 2:length(fl)){
            r=rbind(r,calcFlC(fl[[i]],x))
            }
      return(r)
      }
calcFlC=function(fl,x) { # MUST NEVER BE CALLED
      if(class(fl) != "fl") return(ERRORFl1)
      r= apply(fl$coef*t(x),2,sum)
      r=as.data.frame(r)
      names(r)=fl$name
      rownames(r)=rownames(x)
      return(r)
      }


##################################################################
# A simmple way to study a function along a slidding window
# is to define a list of sequence each of them corresponding to one window
# SeqFrag class allows to do that without copying the original sequence
# slidingWind(seq,begin=1,winlen=90,step=30) constructs this list
# BEWARE if original sequence is modified unpredicted results could result
##################################################################
#BEWARE seq must be a character string representing the name of the variable
slidingWin = function(seq,begin=1,winlen=90,step=30) {
        s=eval(as.symbol(seq),envir=globalenv())
        print(s)
        print(class(s))
        if(begin>getLength(s))return(NULL)
        if(winlen>getLength(s))return(NULL)
        n=(getLength(s)-winlen-begin+1)%/%step;
        if(n==0)return(NULL)
        l=list()
        b=begin
        for(i in 1:n){
        l=c(l,list(initSeqFrag(seq,begin=b,end=b+winlen-1,name=getName(s),compl=TRUE,num=i)))
        b=b+step
        }
        return(l)
        }
      
#############################################################
# Here is the begining of random generation of sequences
# All these methods generates instances of SeqSimple
# Available methods are:
# rMarkov0(length,proba=c(0.25,0.25,0.25,0.25),name="simMark0")
#############################################################

rMarkov0 = function(length,proba=c(0.25,0.25,0.25,0.25),name="simMark0") {
        x=sample(c("a","c","g","t"),length,replace=TRUE,prob=proba)
        s=initSeqSimple(c2s(x),name=name)
        return(s)
        }

# une fonction pour les utilisateurs naifs
aidefreqLseq="calcule les frequences des bases, des doublets, des codons, ... dans une liste de sequences"
aidecalcFl="peut calculer G+C, GCIII, Fop, ... a partir des resultats de freqLSeq"
aidewindow="génère une liste de séquences en faisant glisser une fenêtre sur une séquence donnée"
aide = function(q="") {
	if(q=="") {
	      print(paste("freqLSeq(listSeq,length,step,begin):",aidefreqLseq),quote=FALSE)
	      print(paste("calcFl(fl,x):",aidecalcFl),quote=FALSE)
              print(paste("slidingWin(seq,begin=1,winlen=90,step=30):",aidewindow),quote=FALSE)
	      }
	}


