# ==> acnucopen&db=xxxxx
# <== code=2	missing db= argument
# code=3	if no database with that name is known by the server
# code=4	if database is currently unavailable
# code=5	if a database is currently opened and has not been closed
# code=6&challenge=xx	if the database requires password authorization, server sends a challenge to client.
#	==> reply=xx	authorization data sent by client to server that must be the MD5 digest of
#		the string "challenge:dbname:md5-pw" where md5-pw is the MD5 digest of the password.
#	<== code=6	when authorization failed
# code=0&type=[GENBANK|EMBL|SWISSPROT|NBRF]&totseqs=xx&totspecs=xx&totkeys=xx
# &ACC_LENGTH=xx&L_MNEMO=xx&WIDTH_KW=xx&WIDTH_SP=xx&WIDTH_SMJ=xx&WIDTH_AUT=xx&WIDTH_BIB=xx&lrtxt=xx&SUBINLNG=xx
# Initiates remote access to an acnuc database.
# The db= argument identifies the target database by a logical name 
# that can be any dbname returned by the knowndbs command, or taken from the 1st column of this table,
# or the name of a database requiring password authorization.
# type : the type of database that was opened.
# totseqs, totspec, totkey : total number of seqs, species, keywords in opened database.
# ACC_LENGTH, L_MNEMO, WIDTH_KW, WIDTH_SP, WIDTH_SM, WIDTH_AUT, WIDTH_BIB, lrtxt, SUBINLNG: max lengths of 
# record keys in database

acnucopen <- function(db, socket, challenge = NA){
  #
  # Check arguments:
  #
  if(!is.na(challenge) ) stop("password protection not implemented yet")
  #
  # Build request:
  #
  request <- paste("acnucopen&db=", db, sep = "")
  writeLines(request, socket, sep = "\n")
  answerFromServer <- readLines(socket, n = 1)
  #
  # Check that there is an answer from server:
  #
  if(length(answerFromServer) == 0){
    stop("Empty answer from server")
  }
  res <- parser.socket(answerFromServer)
  #
  # Check that no error is returned:
  #
  if(res[1] != "0"){
    if( res[1] == "1" ){
      stop("unrecognized command") # should not happen
    }
    if( res[1] == "2" ){
      stop("missing db= argument") # should not happen
    }
    if( res[1] == "3" ){
      stop(paste("Database with name -->", db, "<-- is not known by server.\n", sep = ""))
    }
    if( res[1] == "4" ){
      stop(paste("Database with name -->", db, "<-- is currently off for maintenance, please try again later.\n", sep = ""))
    }
    if( res[1] == "5" ){
      stop("A database is currently opened and has not been closed.\n")
    }
    if( res[1] == "6" ){
      stop(paste("Database with name -->", db, "<-- is protected by a password (unimplemented).\n", sep = ""))
    }
    stop("I don't know what this error code means for acnucopen, please contact package maintener.\n")
  }
  return(list(type = res[2], totseqs = as.numeric(res[3]), totspecs = as.numeric(res[4]),
              totkeys = as.numeric(res[5]), ACC_LENGTH = as.numeric(res[6]),
              L_MNEMO = as.numeric(res[7]), WIDTH_KW = as.numeric(res[8]),
              WIDTH_SP = as.numeric(res[9]), WIDTH_SMJ = as.numeric(res[10]),
              WIDTH_AUT = as.numeric(res[11]), WIDTH_BIB = as.numeric(res[12]),
              lrtxt = as.numeric(res[13]), SUBINLNG = as.numeric(res[14])))
}

