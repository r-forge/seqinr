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

freq = function(s,l,p,d) {
        colonne=list()
        x=rep(0,4^l)
        debuts=seq(d,length(s)-l,by=p)
        for(i in debuts){
                n=rbase4(s[i:(i+l-1)])
                x[n+1]=x[n+1]+1 # R numerote betement a partir de 1 et non de 0
                }
        for(i in 1:(4^l)){
                colonne=c(colonne,motif(i-1,l))
                }
        names(x)=colonne
        return(as.data.frame(t(x)))
        }
        
        
#la meme mais utilisant une séquence

freqSeq = function(s,l,p,d){
        x=freq(seq2num(getFrag(s)),l,p,d)
        rownames(x)=getName(s)
        return(x)
        }
        
#la meme mais utilisant une liste de séquences


freqLSeq = function(liste,l,p,d) {
        n=length(liste)
        x=freq(seq2num(getFrag(liste[[1]])),l,p,d)
        rownames(x)=getName(liste[[1]])
        for(i in 2:n){
                 x1=freq(seq2num(getFrag(liste[[i]])),l,p,d)
                 rownames(x1)=getName(liste[[i]])
                 x=rbind(x,x1)
                 }
        return(x)
        }

######################################################################
#classes de séquences
# toutes les classes doivent aoir exactement la même interface à savoir:
# une fonction initNomClasse qui retourne une instance de la classe
# des spécialisations des fonctions
#                    getFrag(seq,debut=0,fin=getLength(seq)) retourne un tableau de caractères
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

getFrag.SeqSimple=function(seq,debut=1,fin=getLength(seq)){
        return(s2c(substr(seq$seq,debut,fin)))
        }
        
getLength.SeqSimple = function(seq){return(nchar(seq$seq))}

getName.SeqSimple = function(seq){return(seq$name)}

getprop.SeqSimple = function(seq){return(list())}

#######################
#fin de la classe SeqSimple
#######################






