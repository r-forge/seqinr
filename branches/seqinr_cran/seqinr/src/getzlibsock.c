#include <R.h>
#include <Rdefines.h>


SEXP getzlibsock(SEXP sock , SEXP debug)
{


  int debugon;
  FILE *fsocket;



  SEXP rka;
  SEXP rks;
  SEXP rvka;
  SEXP rvks;
  SEXP res;
/*  SEXP lsequtil; The effective number of sites used, to be implemented */

  debugon = INTEGER_VALUE(debug);
 /* totseqs = INTEGER_VALUE(nbseq);*/
  prepare_sock_gz_r(fsocket);
   
  if(debugon) Rprintf("C> mode degug is on at C level\n");

}
