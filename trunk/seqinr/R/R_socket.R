###################################################################################################
# Ces fonctions permettent d'interroger les banques structurées sous ACNUC au travers des sockets.#
###################################################################################################



choosebank.socket = function(bank = NA ,host = "pbil.univ-lyon1.fr", port = 5558){

	# ouverture d'un "client socket" sur le serveur pbil et sur le port 5558
	socket = socketConnection( host = host, port = port, server = F, blocking=T)
	rep1 = readLines(socket, n=1)

	# Si pas de banques spécifiées: liste des banques
	if(is.na(bank)){
	writeLines("knowndbs",socket, sep = "\n")
	rep = readLines(socket, n=1)
	res = readLines(socket, n=as.numeric(parser.socket(rep)))
	res = sapply(res,function(x){
		 pos=grep(" ",s2c(x))
		 substr(x,1,(pos[1]-1))
		})
	bank=as.vector(res[ menu(res)])
	}

	# Ouverture de la banque
	request = paste("acnucopen&db=",bank,sep="") 
	writeLines( request, socket, sep = "\n")
	rep2 = readLines(socket, n=1) 
        res = parser.socket(rep2)

	# vérification du bon fonctionnement de la connection et ouverture de la banque
	if(res[1] != "0"){
		print("bank name incorrect")
		rm(socket)
	}
	else{ 
		assign("banknameSocket",bank,.GlobalEnv)
		return(list(socket=socket, bankname = bank, totseqs = res[3], totspecs=res[4], totkeys=res[5]))
	}
}	
	


closebank.socket = function( socket ){
	
	# fermeture de la banque acnuc

	writeLines( "acnucclose", socket, sep="\n" )
	rep = readLines( socket ,n=1)
	if(parser.socket(rep)=="0") print("OK ACNUC SOCKET STOPED")
	close(socket)
	rm(socket)
}


parser.socket = function(p)
{
	p1=s2c(p)
	b=c(which(p1=="="))
  	a=c(which(p1=="&"),length(p1)+1)
	o=character(length(a))
	for(i in 1:length(a)) {o[i]=substr(p,b[i]+1,a[i]-1) }
	return(o)
}
	


getSequence.socket = function( socket, name, start, length){
	
	request2 = paste("gfrag&name=", name,"&start=", start, "&length=", length, sep= "")
	writeLines( request2, socket, sep="\n" )
	s = readLines(socket,n=1)

	if(length(s)==0){
		 print("invalid sequence name")
		}
	else{		
		s = s2c(s)
		sequence = s[(grep("&",s)+1):length(s)]
		return(sequence)
	}
}
	


getAttribut.socket = function( socket, name){
	
	# Récupération des attributs d'une séquence
	
	request = paste( "isenum&name=",name, sep="")
	writeLines( request, socket, sep="\n")
	res = readLines( socket, n=1 )
	p=parser.socket(res)
	l=list( length = as.numeric(p[2]),frame = as.numeric(p[3]), gencode =as.numeric(p[4]) )
	return(l)
}


readAnnots.socket = function( socket, name, nl){

	request= paste( "read_annots&name=", name, "&nl=", nl, sep= "")
	annots=character(nl)
	for(i in 1:nl){
	   writeLines( request , socket, sep="\n")
	   annots[i] =  readLines( socket , n=1 )
	}
	return(annots)
}




getNumber.socket = function( socket, name){
	
	request = paste( "isenum&name=",name, sep="")
	writeLines( request, socket, sep="\n")
	s = readLines(socket,n=1)
	return(parser.socket(s)[1])
}

query.socket = function (socket, listname, query) 
{
    writeLines("prep_requete", socket, sep = "\n")
    readLines(socket, n = 1)
    request = paste("proc_requete&query=\"", query, "\"&name=\"", 
        listname, "\"", sep = "")
    writeLines(request, socket, sep = "\n")
    res = readLines(socket, n = 1)
    p = parser.socket(res)
    if (as.numeric(p[1]) != 0) 
        stop(paste("invalid request:", p[2], sep = ""))
    lrank = p[2]
    first = 1
    liste = character(as.numeric(p[3]))
    for (i in 1:length(liste)) {
        writeLines(paste("nexteltinlist&lrank=", lrank, "&first=", 
            first, "&type=SQ", sep = ""), socket, sep = "\n")
        res = readLines(socket, n = 1)
        r = parser.socket(res)
        first = r[1]
        liste[i] = r[2]
    }
    liste = lapply(liste, function(x){substring(x,2,nchar(x)-1)})
    liste = lapply(liste, as.SeqAcnucWeb, socket)	
    result = list(call = match.call(), name = listname, req = as.list(liste), 
        socket = socket)
    class(result) = c("qaw")
    assign(listname, result, env = .GlobalEnv)
    print(result)
}


print.qaw <- function(x, ...)
{

	cat("\n")
	cat("\n$socket: ")
	print(x$socket)
	cat("\n$banque: ")
	cat(get("banknameSocket",env=.GlobalEnv))
	cat("\n$call: ")
	print(x$call)
	cat("$name: ")
	print(x$name)
	cat("\n")
	sumry <- array("", c(1, 4), list(1, c("list", "length", "mode", "content")))
	sumry[1, ] <- c("$req",length(x$req),"character","sequences")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
}



getKeyword.socket <- function( socket, name){

         writeLines(paste("isenum&name=",name,sep=""),socket,sep="\n")
         res = readLines( socket , n=1 )
         number = parser.socket(res)[1] 

         writeLines(paste("readsub&num=",number,sep=""),socket,sep="\n")
         res2 = readLines( socket , n=1 ) 
         rr = parser.socket(res2)
         
         writeLines(paste("readshrt&num=",rr[7],sep=""),socket,sep="\n")
         res3 = readLines( socket , n=1 ) 
         
	 p1=s2c(res3)
	 b=c(which(p1=="="))
  	 a=c(which(p1=="&"))
         d=c(which(p1==","),length(p1)+1)
	 o=character(length(a))
         o[1]=substr(res3,a[2]+1,d[1]-1)
	 s = seq(2,length(a)*2,by=2)
         for(i in 1:(length(s)-1)){o[i+1] = substr(res3,d[s[i]]+1,d[s[i]+1]-1)} 

	 lapply(o,function(x){
          	writeLines(paste("readkey&num=",x,sep=""),socket,sep="\n")	
	        res4 = readLines( socket , n=1 ) 
                parser.socket(res4)[2]
})

} 


getLocation.socket <- function( socket, name){

	 writeLines(paste("isenum&name=",name,sep=""),socket,sep="\n")
         res = readLines( socket , n=1 )
         number = parser.socket(res)[1] 
	 
         writeLines(paste("readsub&num=",number,sep=""),socket,sep="\n")
         res2 = readLines( socket , n=1 ) 
         rr = parser.socket(res2)
  
	 # Test si subsequence           
     
	 l=list()	
         if(as.numeric(rr[5]) != 0) stop(" pas un CDS : sequence mere\n")
         else {
		i=1
 	 	writeLines(paste("readext&num=",rr[6],sep=""),socket,sep="\n")		
		res3 = readLines( socket , n=1 )
		r = parser.socket(res3)
		l[[i]] = c(r[3],r[4])
		n=r[5] 
	}
        while(as.numeric(n) != 0){
		i=i+1
		writeLines(paste("readext&num=",n,sep=""),socket,sep="\n")	 
		res4 = readLines( socket , n=1 )
		rrr = parser.socket(res4)
		l[[i]] = c(rrr[3],rrr[4])
		n=rrr[5]
  		}
	return(l)
}	


getType.socket = function(socket){
	writeLines( "readfirstrec&type=SMJ",socket, sep="\n" ) 
	s = readLines(socket,n=1)
	rep = parser.socket(s)
	if(rep[1]!="0") stop("erreur")
	rep = as.numeric(rep)
	writeLines( paste("readsmj&num=",10,"&nl=",rep[2]-10,sep=""), socket, sep="\n" ) 
	ss = readLines(socket,n=rep[2]-9)
	occ = grep("name=\"04",ss)
	h = ss[occ]
	return(lapply(h,function(x){ c(noquote(parser.socket(x))[2],noquote(parser.socket(x))[4])  }))
}


plot.SeqAcnucWeb = function(name){

	if(! inherits(name,c("SeqAcnucWeb"))) stop("Sequence of class SeqAcnucWeb is needed")
	socket = attr(name,"socket")
	par(mfrow=c(2,1))
	q = paste("me n=",name,sep="")
	query.socket(socket,listname= "me",query = q)
	l=getLength(me$req[[1]])
	x=c(0,l+(1/10)*l)
	y = c(0,15)
	plot(x,y,ann=FALSE,type="n",axes=FALSE)
	axis(1, col.axis = "blue")
	title(main= paste("Physical position (base) of the subsequences","\n","on the parent sequence", me$req[[1]], sep=" "),font.main=3, col.main="blue")
	mtext(paste("(length of the parent sequence = ", l, ")", sep=""), cex = 0.7, font = 3)

	if(me$req[[1]] != name){
		cat("It is a subsequence\n")
		p = getLocation.socket(socket,name)
		kk=lapply(p,function(x){rect(x[1],0,x[2],1,col="red", border="red" )})
		plot(c(0,l),c(0,10),type="n",axes=FALSE,ann=FALSE)
		title("Legend",font.main=4)
		legend(9,legend=name,fill="red",bg="cornsilk",ncol = 1)
	}

	else{ 
		q = paste("fi n=",name,sep="")
		query.socket(socket = socket, listname = "filles", query = q )
		type = getType.socket(socket)
		t = lapply(type,function(x){ substring(x[1],4,nchar(x[1])-1) } )
		cat("It is a parent sequence\n")
		f = 1
		cou = 0
		rap=numeric()
		nb=numeric()

	while( f <= length(type) ){
		cat("Enter the number of the group of subsequences that you wish to plot on the parent sequence ?\n")
		cat("Enter 0 to stop\n\n")

	noret = lapply(type,function(x){ 
	cat(paste(x[2],"\n",sep=""),"\n")})
	f=scan()
	if(f == 0 || f>length(type) ) break
	if(sum(as.numeric(f == rap)) == 1){
	cat("This group of subsequences has already been ploted !")
	cou=cou-1
	} 
	cou = cou+1
	q=paste("filles et t=",t[[f]],sep="")
	query.socket(socket = socket, listname = "tmp", query = q ) 

	if(is.na(tmp$req[[1]]) || tmp$req[[1]] == name ){ 
	cat("There is no",type[[f]][2],"on this parent sequence\n")
	cou=cou-1
	}
	else{
		u = lapply(tmp$req,getLocation )
		lapply(u, function(x){ lapply(x,function(x){rect(x[1],0+cou,x[2],1+cou,col=cou, border=cou )})})
		rap[cou]=f
		nb[cou]=length(u)
	}
	}
	if(length(rap)==0) stop("There is no type for this sequence")
	plot(c(0,l),c(0,10),type="n",axes=FALSE,ann=FALSE)
	title("Legend",font.main=4)
	legend(9,legend=paste(t[rap],"(",nb,")",sep=""),fill=c(1:cou),bg="cornsilk",ncol = 4)
	par(mfrow=c(1,1))
	}
}
