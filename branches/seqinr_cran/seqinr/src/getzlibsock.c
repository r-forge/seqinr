#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <zlib.h>		/* needs to be before Rconnections.h */
#include "Rconnections.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <string.h>

void *prepare_sock_gz_r(int ncon );
char *z_read_sock(void *v);
int close_sock_gz_r(void *v);
/*static void *extract_opaque = NULL;*/

#define R_EOF	-1	
#define MAXESSAIS 1000

SEXP getzlibsock(SEXP sock, SEXP nmax, SEXP debug)
{

  void *extract_opaque;
  int debugon;
  int testc;
  SEXP essai,essai2; 
  int numsoc;
  SEXP numsocinR;
  
  /*------- */
  SEXP ans = R_NilValue, ans2;
  int i,j, n, nn, nnn, ok, warn, nread, c;
  int itest;
  Rconnection con = NULL;
  Rboolean wasopen;
  
  /*-------*/
  char *test;
  char *res,*megares;
  char *mot[30];
  
  
  /*a modifier?*/
  debugon = INTEGER_VALUE(debug);
  n=INTEGER_VALUE(nmax);
		  
  if(!inherits(sock, "connection")) {
  	Rprintf("Error!\n\n'con' is not a connection");}
  else{
  	if (debugon)
		Rprintf("'con' is a connection...\n");
/***	con = getConnection(asInteger(sock));
	wasopen = con->isopen;
	if(!wasopen) {
		if(!con->open(con)) 
			Rprintf("Error!\n\ncannot open the connection");
		}
	else { ***/
	/* for a non-blocking connection, more input may have become available, so re-position */
	/***
		if (debugon)
			Rprintf("I can open the connection...\n");
		if(con->canseek && !con->blocking)
	    		con->seek(con, con->seek(con, -1, 1, 1), 1, 1);
    		}
 	con->incomplete = FALSE;	
***/
	/*numsoc = asInteger(sock)+1;*/
	numsoc = asInteger(sock);
	if (debugon)
   		printf("Socket number is %d....\n",numsoc);
   	extract_opaque=prepare_sock_gz_r(numsoc);
   	res=(char *)malloc(500*sizeof(char)); 
   	res=z_read_sock(extract_opaque);
/*AJOUT PATCHE CRADO*/
	itest=0;
	while ( res == NULL){ 
		res=z_read_sock(extract_opaque);
		itest++;
		if (debugon)
			printf("#");
		if (itest> MAXESSAIS) {
			printf("Socket error!\n");
			printf("No answer from socket after %d trials!\n",itest);
			ans2 = allocVector(STRSXP, 1);
			SET_STRING_ELT(ans2, 0,mkChar("No answer from socket."));
			PROTECT(ans = ans2);
    			UNPROTECT(1);
			testc=close_sock_gz_r(extract_opaque);
			return ans ;
			}	
		}
	if( res != NULL){ 
		if (debugon)
   			printf("-->%s\n",res);
   		if (strcmp(res,"code=0") == 0) {
			if (debugon)
   			printf("Socket answer is ok %s(%d)\n",res, strlen(res));
			} 
    		nn = (n < 0) ? 1000 : n; /* initially allocate space for 1000 lines */
    		nnn = (n < 0) ? INT_MAX : n;
    		PROTECT(ans = allocVector(STRSXP, nn));
		nread=0;
		/*nnn=50;*/
		if (debugon)
    			printf("n=%d, nn=%d,nnn=%d\n",n,nn, nnn);
		while (res != NULL)  {
			if (nread >=nnn) {
				printf("Increasing memory...\n");
	    			ans2 = allocVector(STRSXP, 2*nn);
	    			for(i = 0; i < nn; i++)
					SET_STRING_ELT(ans2, i, STRING_ELT(ans, i));
	    			nn *= 2;
	    			nnn=nn;
	    			UNPROTECT(1); /* old ans */
	    			PROTECT(ans = ans2);
				};
			if  (strcmp(res,"\033count=1") != 0){
				SET_STRING_ELT(ans, nread, mkChar(res));
				nread++;
				}
			res=z_read_sock(extract_opaque);
    			}
		if (debugon)
			printf("Number of lines of data : %d\n",nread-2);
    		ans2 = allocVector(STRSXP, nread-2);
		j=0; /* Je remplis en evitanat la premiere et la derniere ligne*/
    		for(i = 1; i < nread-1; i++){
			SET_STRING_ELT(ans2, j, STRING_ELT(ans, i));
			j++;
			}
	    	UNPROTECT(1); /* old ans */
	    	PROTECT(ans = ans2);
    		UNPROTECT(1);
		if (debugon)
     			printf("Ok, everything is fine!\n");
		} 
	else {
		printf("Socket error!\n");
		ans2 = allocVector(STRSXP, 1);
		SET_STRING_ELT(ans2, 0,mkChar("Empty answer from socket."));
		PROTECT(ans = ans2);
    		UNPROTECT(1);
		
		}
    	testc=close_sock_gz_r(extract_opaque);
	}

 /*** if(!wasopen) con->close(con);***/
  return ans ;
}
