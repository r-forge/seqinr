choosebank <- function( bank = "demo")
{

#
# Test if a bank is already in use
#

 if( exists("bankname",envir=globalenv()) && (get("bankname",envir=globalenv()) != bank ))
{ 
# stop("An ACNUC data base is already in use, you must change the data base using changebank()")
#
.C("Racnucclose")
}
else
{	
if( exists("bankname",envir=globalenv()) && get("bankname",envir=globalenv()) == bank )
{ 
stop("this data base is already in use")
}
}
 
#
# Define a function to print something when we are unable
# to locate the requested bank
#
  failure.warning <- function( bank )
  {
    warning(paste("Unable to get environment variable ", bank))
    warning("seting ACNUC bank to local demo in sequences/entero")
  }
#
# Set defaults to local bank demo:
#
  assign("bankname", bank, .GlobalEnv)	
  acnuc <- system.file("sequences/entero", package = "seqinr")
  gcgacnuc <- acnuc
#
# Check argument:
#
  if( missing( bank ) )
  {
    warning("bank argument is missing")
    warning("seting ACNUC bank to local demo in sequences/entero")
  }
#
# Try to find bank location:
#
  else
  {
    ad <- Sys.getenv( bank )
    if( ad == "" ) # We were unable to get bank location
    {
      if( Sys.info()["user"] == "ADE-User" & Sys.info()["nodename"] == "pbil" )
      {
        #
        # seqinr was run from Rweb at the URL:
        # http://pbil.univ-lyon1.fr/Rweb/Rweb.general.html
        #
        banklocs <- readLines("/misc/pub/banques_acnuc")
        targetbank <- grep( bank, banklocs )
        if( length( targetbank ) == 0 )
        {
          failure.warning( bank )
        }
        else
        {
           line <- banklocs[targetbank]
           words <- unlist( strsplit(line, split = " ") )
           words <- words[nchar(words)>0]
           acnuc <- substr(words[3], 2, nchar(words[3]))
           gcgacnuc <- substr(words[4], 1, nchar(words[4]) - 1)
        }
      }
      else
      {
        failure.warning( bank )
      }
    }
    else # We were able to get bank location
    {
      var <- unlist(strsplit(ad, " "))
      acnuc <- var[1]
      gcgacnuc <- var[2]
    }
  }
#
# Open ACNUC database:
#
  .C("Racnucopen", acnuc, gcgacnuc)
#
# print ACNUC database content
#
  helpfile <- paste(acnuc, "/HELP", sep = "")
  if( file.exists( helpfile ) )
  {
    help <- readLines(helpfile, n = 100)
    words <- unlist(strsplit(help[1], split = " "))
    words <- words[nchar(words)>0]
    nlines <- as.numeric(words[length(words)])
    for( i in 2:(nlines+1) )
      cat( help[i],"\n" )
  }
  else
  {
    warning(paste("HELP file not found in folder",acnuc))
  }
}


choixbanque <- choosebank # Just an alias for compatibility


query <- function(listname,request)
{	

#
# Test if a database is open
#
	if(exists("bankname",envir=globalenv()) == FALSE){
	stop("you must first choose a database using choosebank()") 
	}
#	
# call to the C funcion getreq() to obtain the list of sequences name 
#
	liste<-.Call("getreq",listname,request)
	liste<-as.character(liste)
	liste<-lapply(liste,s2c)
	p<-function(m) paste(m[m!=" "],collapse="")
	liste<-lapply(liste,p)
	liste = lapply(liste,initSeqReq)	
	toto=list(call = match.call(),name=listname,req=liste)
	class(toto)=c("requete")
	assign(listname,toto,env = .GlobalEnv)
	print(toto)
}

getseq <- function(name, as.string = TRUE)
{
  x <- .Call("getseq", name)
  ifelse( as.string, return(x), return(s2c(x)) )
}

getseq2 <- function(name,B1,B2,as.string = TRUE)
{
 	if(as.string == TRUE){
 	.Call("getseq2",name,B1,B2)
	}
	else{
 	tata<-.Call("getseq2",name,B1,B2)
 	seq <- unlist(strsplit(tata,""))
  	toupper(seq)
	}
}


changebank <- function(bankname)
{
	.C("Racnucclose")
	choosebank(bankname)
}


print.requete <- function(x, ...)
{
	cat("\n")
	cat("class: ")
	cat(class(x))
	cat("\n$banque: ")
	cat(get("bankname",env=.GlobalEnv))
	cat("\n$call: ")
	print(x$call)
	cat("$name: ")
	print(x$name)
	cat("\n")
	sumry <- array("", c(1, 4), list(1, c("list", "length", 
        "mode", "content")))
	sumry[1, ] <- c("$req",length(x$req),"character","mnemoniques")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
}

translateCDS <- function(name)
{
	.Call("translateCDS",name)
}

getKeyword <- function(name)
{
	k = .Call("getKey",name)
	k = as.character(k)
	k = as.list(k)
	k =strsplit(k," ")
	p<-function(m) paste(m[m!=""],collapse=" ")
	k = lapply(k,p)
	k = unlist(k)
	k
}


getExon <- function(name)
{
	a=.Call("getExon",name)
	deb=a[seq(1,length(a),2)]
	fin=a[seq(2,length(a),2)]	
	b=as.list(rep(0,(length(a)/2)))
	for(i in 1:(length(a)/2)) b[[i]]=c(deb[i],fin[i])
	b
}


getAttribut <- function(name)
{
	res = .Call("getAttribut",name)
	att = list(lseq=res[1],frame=res[2],gencode=res[3])		
	return(att)
}


getAnnots <-function(name,ligne)
{

	res = .Call("getAnnots",name,ligne)
	res = as.character(res)
	res
}


