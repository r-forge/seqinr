#include <R.h>
#include <Rdefines.h>


SEXP getzlibsock(SEXP FILE , SEXP debug)
{


  int debugon;



  SEXP rka;
  SEXP rks;
  SEXP rvka;
  SEXP rvks;
  SEXP res;
/*  SEXP lsequtil; The effective number of sites used, to be implemented */

  debugon = INTEGER_VALUE(debug);
 /* totseqs = INTEGER_VALUE(nbseq);*/
   
  if(debugon) Rprintf("C> mode degug is on at C level\n");

}
