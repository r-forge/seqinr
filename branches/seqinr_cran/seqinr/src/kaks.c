#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <math.h>
#include <Rinternals.h>
#include <R.h>
#include <Rdefines.h>

#define maxnseqs 100
  

/********************************************************************/
/****************************** LWL93 *******************************/
/********************************************************************/

/* Le programme lwl93 calcule les taus de substitutions synonymes */
/* et non-synonymes au sens de Li (1993, J.M.E.36. 96-99) entre */
/* toutes les paires de sequences d'un fichier mase. La sortie est */
/* un ou deux fichier .num contenant les Ka et/ou Ks avec/sans leurs */
/* variances. */

int code_mt;

int readmaseseqs(char *, char **, char **, char **, int);
void reresh(char **, int, int);
void prefastlwl(float **, float **, float **, float **, float **, float **, float **, float **, float **, float **);
int fastlwl(char **, int, int, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **, float **);



SEXP kaks(SEXP sequences, SEXP nbseq)
{

  char **seq;
  float *tl0[64], *tl1[64], *tl2[64], *tti0[64], *tti1[64], *tti2[64], *ttv0[64], *ttv1[64], *ttv2[64];
  char  buff[40];
  int i, j, ii, jj, kk, k, l, lg, totseqs, lgseq, *aa, m, n,  *listbranche[maxnseqs], option2;
  float *ka[100], *ks[100], *bootka[100], *bootks[100], *rl[21], a, b, *vka[100], *vks[100];
  FILE *nuc;
  char *fname; 
  
 float mat[19][19] = {{.382, .382, .343, .382, .382, .382, .382, .128, .040, .128, .040, .128, .040, .128, .040, .128, .343, .128, .040 }, 
		     { .382, .382, .128, .343, .343, .343, .343, .128, .040, .128, .040, .128, .040, .128, .040, .128, .128, .040, .040 }, 
		     { .343, .128, .343, .382, .382, .382, .343, .128, .040, .128, .128, .343, .128, .343, .128, .343, .343, .128, .040 },
		     { .382, .343, .382, .343, .343, .343, .343, .343, .040, .343, .343, .382, .343, .382, .343, .382, .382, .382, .343 },
		     { .382, .343, .382, .343, .382, .382, .382, .343, .040, .343, .128, .343, .128, .128, .128, .343, .343, .128, .040 },
		     { .382, .343, .382, .343, .382, .382, .382, .343, .040, .343, .128, .343, .128, .128, .040, .128, .128, .128, .040 },
		     { .382, .343, .343, .343, .382, .382, .382, .343, .040, .343, .128, .343, .128, .128, .128, .128, .343, .128, .040 },
		     { .128, .128, .128, .343, .343, .343, .343, .343, .040, .343, .128, .343, .128, .343, .128, .343, .343, .128, .040 },
		     { .040, .040, .040, .040, .040, .040, .040, .040, .040, .382, .382, .382, .343, .343, .343, .128, .128, .343, .128 },
		     { .128, .128, .128, .343, .343, .343, .343, .343, .382, .040, .040, .128, .128, .040, .128, .040, .040, .040, .040 },
		     { .040, .040, .128, .343, .128, .128, .128, .128, .382, .040, .343, .343, .343, .343, .128, .128, .128, .128, .128 },
		     { .128, .128, .343, .382, .343, .343, .343, .343, .382, .128, .343, .343, .343, .343, .343, .128, .128, .343, .343 },
		     { .040, .040, .128, .343, .128, .128, .128, .128, .343, .128, .343, .343, .343, .382, .343, .343, .343, .343, .343 },
		     { .128, .128, .343, .382, .128, .128, .128, .343, .343, .040, .343, .343, .382, .343, .382, .128, .128, .343, .343 },
		     { .040, .040, .128, .343, .128, .040, .128, .128, .343, .128, .128, .343, .343, .382, .382, .343, .382, .382, .343 },
		     { .128, .128, .343, .382, .343, .128, .128, .343, .128, .040, .128, .128, .343, .128, .343, .343, .343, .382, .382 },
		     { .343, .128, .343, .382, .343, .128, .343, .343, .128, .040, .128, .128, .343, .128, .382, .343, .382, .343, .128 },
		     { .128, .040, .128, .382, .128, .128, .128, .128, .343, .040, .128, .343, .343, .343, .382, .382, .343, .343, .343 },
		     {.040, .040, .040, .343, .040, .040, .040, .040, .128, .040, .128, .343, .343, .343, .343, .382, .128, .343, .382 }};




  SEXP rka;
  SEXP rks;
  SEXP rvka;
  SEXP rvks;
  SEXP toto;
  SEXP res;

   totseqs = INTEGER_VALUE(nbseq);
  
   seq = (char **)malloc(totseqs*sizeof(char *));
   
   for(i=0;i<totseqs;i++){
      seq[i] = CHAR(STRING_ELT(sequences,i));
   }
   
   lgseq = strlen(seq[0]);

   PROTECT(res=allocVector(VECSXP,4));
   PROTECT(rka=NEW_NUMERIC(totseqs*totseqs));
   PROTECT(rks=NEW_NUMERIC(totseqs*totseqs));
   PROTECT(rvka=NEW_NUMERIC(totseqs*totseqs));
   PROTECT(rvks=NEW_NUMERIC(totseqs*totseqs));
   

	for (i = 0; i < 64; i++) {
		if ((tl0[i] = (float *) malloc(64 * sizeof(float))) == NULL) {
			error("Pas assez de memoire\n");
		}
		if ((tl1[i] = (float *) malloc(64 * sizeof(float))) == NULL) {
			error("Pas assez de memoire\n");
		}
		if ((tl2[i] = (float *) malloc(64 * sizeof(float))) == NULL) {
			error("Pas assez de memoire\n");
		}
		if ((tti0[i] = (float *) malloc(64 * sizeof(float))) == NULL) {
			error("Pas assez de memoire\n");
		}
		if ((tti1[i] = (float *) malloc(64 * sizeof(float))) == NULL) {
			error("Pas assez de memoire\n");
		}
		if ((tti2[i] = (float *) malloc(64 * sizeof(float))) == NULL) {
			error("Pas assez de memoire\n");
		}
		if ((ttv0[i] = (float *) malloc(64 * sizeof(float))) == NULL) {
			error("Pas assez de memoire\n");
		}
		if ((ttv1[i] = (float *) malloc(64 * sizeof(float))) == NULL) {
			error("Pas assez de memoire\n");
		}
		if ((ttv2[i] = (float *) malloc(64 * sizeof(float))) == NULL) {
			error("Pas assez de memoire\n");
		}
	}


	for (i = 0; i < 21; i++)
		rl[i] = (float *) malloc(21 * sizeof(float));

	for (i = 2; i < 21; i++) {
		for (j = 1; j < i; j++) {
			*(rl[i] + j) = mat[j-1][i-2] ;
		}
	}

	for (i = 1; i <= 20; i++) {
		*(rl[i] + i) = 1.0;
		for (j = i + 1; j <= 20; j++)
			*(rl[i] + j) = *(rl[j] + i);
	}


	lgseq = strlen(seq[0]);

	for (i=0;i<totseqs;i++){
		for(j=0;j<lgseq;j++){
			if ((*(seq[i]+j)!='A') && (*(seq[i]+j)!='G') && (*(seq[i]+j)!='C') && (*(seq[i]+j)!='T') && (*(seq[i]+j)!='-') ) {
				if (j%3==0) {
					*(seq[i]+j)='-';
					*(seq[i]+j+1)='-';
					*(seq[i]+j+2)='-';
				}

				if (j%3==1) {
					*(seq[i]+j)='-';
					*(seq[i]+j+1)='-';
					*(seq[i]+j-1)='-';
				}

				if (j%3==2) {
					*(seq[i]+j)='-';
					*(seq[i]+j-1)='-';
					*(seq[i]+j-2)='-';
				}
			}
		}
	}


	reresh(seq,totseqs,0);
	
	
	for (i = 0; i <= totseqs; i++) {
		ka[i] = (float *) malloc((totseqs + 1) * sizeof(float));
		vka[i] = (float *) malloc((totseqs + 1) * sizeof(float));
		ks[i] = (float *) malloc((totseqs + 1) * sizeof(float));
		vks[i] = (float *) malloc((totseqs + 1) * sizeof(float));
	}

	lgseq = strlen(seq[0]);

	prefastlwl(rl, tl0, tl1, tl2, tti0, tti1, tti2, ttv0, ttv1, ttv2);
	fastlwl(seq, totseqs, lgseq, ka, ks, tti0, tti1, tti2, ttv0, ttv1, ttv2, tl0, tl1, tl2, vka, vks);

	
    /********************************************************************************/
    /* Remplissage de l'objet R (matrice de taille nb_seq * nb_seq  avec ks       */
    /********************************************************************************/

	n=0;
			
       	for(i=0;i<totseqs-1;i++){
	  for(j=i+1;j<totseqs;j++){
		    REAL(rks)[n+j-i]=ks[i][j];
		  }
		  n=n+totseqs+1;
	}
	

  /********************************************************************************/
    /* Remplissage de l'objet R (matrice de taille nb_seq * nb_seq  avec ka       */
    /********************************************************************************/

	n=0;
			
       	for(i=0;i<totseqs-1;i++){
	  for(j=i+1;j<totseqs;j++){
		    REAL(rka)[n+j-i]=ka[i][j];
		  }
		  n=n+totseqs+1;
	}
	
  /********************************************************************************/
    /* Remplissage de l'objet R (matrice de taille nb_seq * nb_seq  avec vks       */
    /********************************************************************************/

	n=0;
			
       	for(i=0;i<totseqs-1;i++){
	  for(j=i+1;j<totseqs;j++){
		    REAL(rvks)[n+j-i]=vks[i][j];
		  }
		  n=n+totseqs+1;
	}
	

      /********************************************************************************/
	/* Remplissage de l'objet R (matrice de taille nb_seq * nb_seq  avec vka       */
    /********************************************************************************/

	n=0;
			
       	for(i=0;i<totseqs-1;i++){
	  for(j=i+1;j<totseqs;j++){
		    REAL(rvka)[n+j-i]=vka[i][j];
		  }
		  n=n+totseqs+1;
	}
	

	 SET_ELEMENT(res,0,rka);
	 SET_ELEMENT(res,1,rks);
	 SET_ELEMENT(res,2,rvks);
	 SET_ELEMENT(res,3,rvka);

	UNPROTECT(5);
	return(res);
}


int fastlwl(char **seq, int nbseq, int lgseq, float **ka, float **ks, float **tti0, float **tti1, float **tti2, float **ttv0, float **ttv1, float **ttv2, float **tl0, float **tl1, float **tl2, float **vka, float **vks)
{

	const float     trois = 3.0;
	float           l[3], a[3], b[3], p[3], q[3], ti[3], tv[3], cc[3],
	                aaa[3], bb[3], flgseq, va[3], vb[3],es1,es2,es3,es4,es5,es6,es7,es8;
	char            ci1, ci2, ci3, cj1, cj2, cj3, cod1[3], cod2[3];
	int             i, j, ii, jj, nbdiff, cat, pos[3], num1, num2,
	                sat, sat1, sat2;

	sat = sat1 = sat2 = 2;
	flgseq = (float) lgseq;
	if (flgseq / trois != lgseq / 3) {
		printf("Nombre de nt non multiple de trois.\n");
		return;
	}
	for (i = 0; i < nbseq - 1; i++) {
		for (j = i + 1; j < nbseq; j++) {
			l[0] = l[1] = l[2] = 0;
			ti[0] = ti[1] = ti[2] = tv[0] = tv[1] = tv[2] = 0;

	     		for (ii = 0; ii < lgseq / 3; ii++) {
				cod1[0] = *(seq[i] + 3 * ii);
				cod1[1] = *(seq[i] + 3 * ii + 1);
				cod1[2] = *(seq[i] + 3 * ii + 2);
				cod2[0] = *(seq[j] + 3 * ii);
				cod2[1] = *(seq[j] + 3 * ii + 1);
				cod2[2] = *(seq[j] + 3 * ii + 2);
				num1 = num(cod1);
				num2 = num(cod2);
				l[0] += *(tl0[num1] + num2);
				l[1] += *(tl1[num1] + num2);
				l[2] += *(tl2[num1] + num2);
				ti[0] += *(tti0[num1] + num2);
				ti[1] += *(tti1[num1] + num2);
				ti[2] += *(tti2[num1] + num2);
				tv[0] += *(ttv0[num1] + num2);
				tv[1] += *(ttv1[num1] + num2);
				tv[2] += *(ttv2[num1] + num2);
			}



			for (ii = 0; ii < 3; ii++) {
				p[ii] = ti[ii] / l[ii];
				q[ii] = tv[ii] / l[ii];
				aaa[ii] = 1 / (1 - 2 * p[ii] - q[ii]);
				bb[ii] = 1 / (1 - 2 * q[ii]);
				cc[ii] = (aaa[ii] + bb[ii]) / 2;

				if (bb[ii] <= 0) {
					b[ii] = 10;
				} else
					b[ii] = 0.5 * (float) log(bb[ii]);

				if ((aaa[ii] <= 0) || (bb[ii] <= 0)) {
					a[ii] = 10;
				} else
					a[ii] = 0.5 * (float) log(aaa[ii]) - 0.25 * log(bb[ii]);




es1=aaa[ii] * aaa[ii] * p[ii] + cc[ii] * cc[ii] * q[ii];
es2=(aaa[ii] * p[ii] + cc[ii] * q[ii]) * ( aaa[ii] * p[ii] + cc[ii] * q[ii]);

				va[ii] = (aaa[ii] * aaa[ii] * p[ii] + cc[ii] * cc[ii] * q[ii] - (aaa[ii] * p[ii] + cc[ii] * q[ii]) * ( aaa[ii] * p[ii] + cc[ii] * q[ii])) / l[ii];


				vb[ii] = bb[ii] * bb[ii] * q[ii] * (1 - q[ii]) / l[ii];


			}



			if ((a[1] != 10) && (a[2] != 10) && (b[2] != 10)){
				ks[i][j] = (l[1] * a[1] + l[2] * a[2]) / (l[2] + l[1]) + b[2];
				vks[i][j] = (l[1] * l[1]  * va[1] + l[2] * l[2] * va[2]) /  ((l[1] + l[2]) * (l[1]+l[2])) + vb[2] - bb[2] * q[2] * (2 * aaa[2] * p[2] - cc[2] * (1 - q[2]))/(l[1]+l[2]);
			}
			else {
				sat1 = 1;
				vks[i][j]=ks[i][j] = 9.999999;

			}

			if ((a[0] != 10) && (b[0] != 10) && (b[1] != 10)){
				ka[i][j] = a[0] + (l[0] * b[0] + l[1] * b[1]) / (l[0] + l[1]);
				vka[i][j] = (l[0] * l[0]  * vb[0] + l[1] * l[1] * vb[1]) /  ((l[1] + l[0]) * (l[1]+l[0])) + va[0] - bb[0] * q[0] * (2 * aaa[0] * p[0] - cc[0] * (1 - q[0]))/(l[1]+l[0]);
	}

		else {
				vka[i][j]=ka[i][j] = 9.999999;
				sat2 = 1;
			}


		}
	}
	if (sat1 == 1)
		sat = 1;
	if (sat2 == 1)
		sat = 0;
	return sat;
}



int num(char *cod)
{
	int             n1, n2, n3;

	n1 = n2 = n3 = 0;
	if (cod[0] == 'C')
		n1 = 1;
	if (cod[1] == 'C')
		n2 = 1;
	if (cod[2] == 'C')
		n3 = 1;
	if (cod[0] == 'G')
		n1 = 2;
	if (cod[1] == 'G')
		n2 = 2;
	if (cod[2] == 'G')
		n3 = 2;
	if (cod[0] == 'T')
		n1 = 3;
	if (cod[1] == 'T')
		n2 = 3;
	if (cod[2] == 'T')
		n3 = 3;

	return 16 * n1 + 4 * n2 + n3;
}





int catsite(char c1, char c2, char c3, int i) {

	/* renvoie 0 si le site i du codon c1c2c3 est non degenere */
	/* 1                                  2-fold degenerate */
	/* 2                                  4-fold degenerate */

	if (i == 3) {
	   if( !code_mt ) {
		if ( (c1 == 'A') && (c2 == 'T') && (c3 == 'G'))
			return 0;
		if ( (c1 == 'T') && (c2 == 'G') && (c3 == 'A'))
			return 0;
		if ( (c1 == 'T') && (c2 == 'G') && (c3 == 'G'))
			return 0;
		}
		if (c2 == 'C')
			return 2;
		if ((c1 == 'C') && (c2 == 'T'))
			return 2;
		if ((c1 == 'G') && (c2 == 'T'))
			return 2;
		if ((c1 == 'G') && (c2 == 'G'))
			return 2;
		if ((c1 == 'C') && (c2 == 'G'))
			return 2;
		return 1;
	}
	else if (i == 1) {
		if ((c1 == 'C') && (c2 == 'T') && (c3 == 'A'))
			return 1;
		if ((c1 == 'C') && (c2 == 'T') && (c3 == 'G'))
			return 1;
		if ((c1 == 'T') && (c2 == 'T') && (c3 == 'A'))
			return 1;
		if ((c1 == 'T') && (c2 == 'T') && (c3 == 'G'))
			return 1;
	   if( !code_mt ) {
		if ((c1 == 'A') && (c2 == 'G') && (c3 == 'A'))
			return 1;
		if ((c1 == 'A') && (c2 == 'G') && (c3 == 'G'))
			return 1;
		if ((c1 == 'C') && (c2 == 'G') && (c3 == 'A'))
			return 1;
		if ((c1 == 'C') && (c2 == 'G') && (c3 == 'G'))
			return 1;
		}
		return 0;
	}
	return 0;
}


char transf(char nt1, char nt2)
{
	if (nt1 == nt2) {
		printf("Same nt, patate.\n");
		return 'S';
	}
	if ((nt1 == 'A') && (nt2 == 'C'))
		return 'v';
	if ((nt1 == 'A') && (nt2 == 'G'))
		return 'i';
	if ((nt1 == 'A') && (nt2 == 'T'))
		return 'v';
	if ((nt1 == 'G') && (nt2 == 'C'))
		return 'v';
	if ((nt1 == 'G') && (nt2 == 'T'))
		return 'v';
	if ((nt1 == 'C') && (nt2 == 'T'))
		return 'i';
	if ((nt1 == 'C') && (nt2 == 'A'))
		return 'v';
	if ((nt1 == 'G') && (nt2 == 'A'))
		return 'i';
	if ((nt1 == 'T') && (nt2 == 'A'))
		return 'v';
	if ((nt1 == 'C') && (nt2 == 'G'))
		return 'v';
	if ((nt1 == 'T') && (nt2 == 'G'))
		return 'v';
	if ((nt1 == 'T') && (nt2 == 'C'))
		return 'i';

	printf("Error\n%c, %c\n", nt1, nt2);
	return 'E';
}


void titv1(char *cod1, char *cod2, float poids, float *ti, float *tv, float* l)
{
	int             i, j, jj;
	char            a, b, ci1, ci2, ci3, cj1, cj2, cj3;
	char            transf(char, char);
	float		poids2 = poids/2;


	ci1 = cod1[0];
	ci2 = cod1[1];
	ci3 = cod1[2];
	cj1 = cod2[0];
	cj2 = cod2[1];
	cj3 = cod2[2];




	
	for (i = 0; i <= 2; i++)
		if (cod1[i] != cod2[i]) {

			l[catsite(ci1, ci2, ci3, i + 1)]+=0.5 * poids;
			l[catsite(cj1, cj2, cj3, i + 1)]+=0.5 * poids;

			a = cod1[i];
			b = cod2[i];
			if (transf(a, b) == 'i') {
				ti[catsite(ci1, ci2, ci3, i + 1)] += 0.5 * poids;
				ti[catsite(cj1, cj2, cj3, i + 1)] += 0.5 * poids;
			} else {
				tv[catsite(ci1, ci2, ci3, i + 1)] += 0.5 * poids;
				tv[catsite(cj1, cj2, cj3, i + 1)] += 0.5 * poids;
			}

			if( code_mt ) continue;  /* il n'y a plus les pb de TI non-syno et de TV syno avec code_mt ! */
		
	if (((ci2 == 'T') && (cj2 == 'T')) || ((ci2 == 'G') && (cj2 == 'G'))) { /* T ou G ensemble en pos 2 des 2 codons */


if (i==0){ /* pos 1 */	
		/* tous ces cas sont des transitions en un site 2-fold non-syno pour le code universel:
il faut les enlever du comptage des TI 2-fold (ti[1]) et les ajouter au comptage des TV 2-fold (tv[1])
pour le code_mt ce sont des sites non dege qui ont ete traites simplement comme il faut */
		if ((ci1 == 'C') && (ci2 == 'G') && (ci3 == 'A') && (cj1 == 'T') && (cj2 == 'G') && (cj3 == 'A')) {
			ti[1] -= 0.5 * poids; /* CGA / TGA */
			tv[1] += 0.5 * poids;
		}
		if ((ci1 == 'C') && (ci2 == 'G') && (ci3 == 'G') && (cj1 == 'T') && (cj2 == 'G') && (cj3 == 'G')) {
			ti[1] -= 0.5 * poids; /* CGG / TGG */
			tv[1] += 0.5 * poids;
		}
		if ((ci1 == 'A') && (ci2 == 'G') && (ci3 == 'G') && (cj1 == 'G') && (cj2 == 'G') && (cj3 == 'G')) {
			ti[1] -= 0.5 * poids; /* AGG / GGG */
			tv[1] += 0.5 * poids;
		}
		if ((ci1 == 'A') && (ci2 == 'G') && (ci3 == 'A') && (cj1 == 'G') && (cj2 == 'G') && (cj3 == 'A')) {
			ti[1] -= 0.5 * poids; /* AGA / GGA */
			tv[1] += 0.5 * poids;
		}
		if ((ci1 == 'T') && (ci2 == 'G') && (ci3 == 'A') && (cj1 == 'C') && (cj2 == 'G') && (cj3 == 'A')) {
			ti[1] -= 0.5 * poids; /* TGA / CGA */
			tv[1] += 0.5 * poids;
		}
		if ((ci1 == 'T') && (ci2 == 'G') && (ci3 == 'G') && (cj1 == 'C') && (cj2 == 'G') && (cj3 == 'G')) {
			ti[1] -= 0.5 * poids; /* TGG / CGG */
			tv[1] += 0.5 * poids;
		}
		if ((ci1 == 'G') && (ci2 == 'G') && (ci3 == 'G') && (cj1 == 'A') && (cj2 == 'G') && (cj3 == 'G')) {
			ti[1] -= 0.5 * poids; /* GGG / AGG */
			tv[1] += 0.5 * poids;
		}
		if ((ci1 == 'G') && (ci2 == 'G') && (ci3 == 'A') && (cj1 == 'A') && (cj2 == 'G') && (cj3 == 'A')) {
			ti[1] -= 0.5 * poids; /* GGA / AGA */
			tv[1] += 0.5 * poids;
		}


/* tous ces cas sont 
code universel: TV syno en sites 2-fold il faut les enlever du comptage des TV 2-fold (tv[1]) et les ajouter au comptage des TI 2-fold (ti[1])
code_mt: TV non syno en site non dege qui ont ete correctement comptes
*/
		if ((ci1 == 'C') && (ci2 == 'G') && (ci3 == 'A') && (cj1 == 'A') && (cj2 == 'G') && (cj3 == 'A')) {
			tv[1] -= poids; /* CGA / AGA : TV syno code univ, non code mt */
			ti[1] += poids;
		}
		if ((ci1 == 'A') && (ci2 == 'G') && (ci3 == 'A') && (cj1 == 'C') && (cj2 == 'G') && (cj3 == 'A')) {
			tv[1] -= poids; /* AGA / CGA : TV syno code univ, non code mt */
			ti[1] += poids;
		}
		if ((ci1 == 'C') && (ci2 == 'G') && (ci3 == 'G') && (cj1 == 'A') && (cj2 == 'G') && (cj3 == 'G')) {
			tv[1] -= poids; /* CGG / AGG : TV syno code univ, non code mt */
			ti[1] += poids;
		}
		if ((ci1 == 'A') && (ci2 == 'G') && (ci3 == 'G') && (cj1 == 'C') && (cj2 == 'G') && (cj3 == 'G')) {
			tv[1] -= poids; /* AGG / CGG : TV syno code univ, non code mt */
			ti[1] += poids;
		}
}

if (i==2){ /* pos 3 */	
/* tous ces cas sont
code universel: des TV syno en site 2-fold il faut les enlever des TV 2-fold (iv[1]) et ajouter aux TI 2-fold (ti[1])
code_mt: ce sont des TV non syno en site 2-fold qui int ete comptees normalement
*/
		if ((ci1 == 'A') && (ci2 == 'T') && (ci3 == 'A') && (cj1 == 'A') && (cj2 == 'T') && (cj3 == 'T')) {
			tv[1] -= poids; /* TV ATA / ATT : syno code univ, non code mt */
			ti[1] += poids;
		}
		if ((ci1 == 'A') && (ci2 == 'T') && (ci3 == 'T') && (cj1 == 'A') && (cj2 == 'T') && (cj3 == 'A')) {
			tv[1] -= poids; /* TV ATT / ATA : syno code univ, non code mt */
			ti[1] += poids;
		}
		if ((ci1 == 'A') && (ci2 == 'T') && (ci3 == 'A') && (cj1 == 'A') && (cj2 == 'T') && (cj3 == 'C')) {
			tv[1] -= poids; /* TV ATA / ATC : syno code univ, non code mt */
			ti[1] += poids;
		}
		if ((ci1 == 'A') && (ci2 == 'T') && (ci3 == 'C') && (cj1 == 'A') && (cj2 == 'T') && (cj3 == 'A')) {
			tv[1] -= poids; /* TV ATC / ATA : syno code univ, non code mt */
			ti[1] += poids;
		}



/* ces 2 cas sont
code universel: des TI non syno en site 2-fold il faut les enlever des TI 2-fold (ti[1]) et les ajouter aux TV 2-fold (tv[1])
code_mt: des TI syno en site 2-fold qui ont ete comptees normalement
*/
		if ((ci1 == 'A') && (ci2 == 'T') && (ci3 == 'A') && (cj1 == 'A') && (cj2 == 'T') && (cj3 == 'G')) {
			ti[1] -= 0.5 * poids; /* TI ATA / ATG : non syno code univ, syno code mt */
			tv[1] += 0.5 * poids;
		}
		if ((ci1 == 'A') && (ci2 == 'T') && (ci3 == 'G') && (cj1 == 'A') && (cj2 == 'T') && (cj3 == 'A')) {
			ti[1] -= 0.5 * poids; /* TI ATG / ATA : non syno code univ, syno code mt */
			tv[1] += 0.5 * poids;
		}
}
		}

	}
}


void titv2(char *cod1, char *cod2, float *ti, float *tv, float* l, int *aa, float **rl, int* pos)
{

	char            codint1[3], codint2[3];
	int             i, j, n, a1, a2, a3, a4, b1, b2, b3,b4;
	float           l1, l2, p1, p2,f1,f2,f3,f4;
	void            titv1(char *, char *, float, float *, float *,float*);


	memcpy(codint1, cod1, 3);
	memcpy(codint2, cod1, 3);
	for (i = 0; i < 2; i++) {
		if (cod1[i] != cod2[i])
			codint1[i] = cod2[i];
		if (cod1[i] != cod2[i])
			break;
	}
	for (j = i + 1; j <= 2; j++) {
		if (cod1[j] != cod2[j])
			codint2[j] = cod2[j];
		if (cod1[j] != cod2[j])
			break;
	}

	
	l1 = *(rl[aa[num(cod1)]] + aa[num(codint1)]) * *(rl[aa[num(codint1)]] + aa[num(cod2)]);
	l2 = *(rl[aa[num(cod1)]] + aa[num(codint2)]) * *(rl[aa[num(codint2)]] + aa[num(cod2)]);

	p1 = l1 / (l1 + l2);
	p2 = 1 - p1;
	for (i=0;i<3;i++) if (pos[i]==0) n=i+1;
	l[catsite(cod1[0], cod1[1] ,cod1[2], n)]+=0.333333;
	l[catsite(cod2[0], cod2[1] ,cod2[2], n)]+=0.333333;
	l[catsite(codint1[0], codint1[1] ,codint1[2], n)]+=0.333333*p1;
	l[catsite(codint2[0], codint2[1] ,codint2[2], n)]+=0.333333*p2;
	titv1(cod1, codint1, p1, ti, tv,l);
	titv1(cod2, codint1, p1, ti, tv,l);
	titv1(cod1, codint2, p2, ti, tv,l);
	titv1(cod2, codint2, p2, ti, tv,l);


}

void titv3(char *cod1, char *cod2, float *ti, float *tv, float* l, int *aa, float **rl)
{

	char           *codint1[6], *codint2[6];
	int             i, j, ii,a,b,c,d,aaa,aab,aac,aad;
	float           like[6], p[6], somli, rlab, rlbc, rlcd;
	void            titv1(char *, char *, float, float *, float *, float*);
	int             num(char *);

	for (i = 0; i < 6; i++) {
		if ((codint1[i] = (char *) malloc(3 * sizeof(char))) == NULL) {
			printf("Erreur d'allocation\n");
			exit(1);
		}
		if ((codint2[i] = (char *) malloc(3 * sizeof(char))) == NULL) {
			printf("Erreur d'allocation\n");
			exit(1);
		}
	}
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3 ; j++)
			if (j != i) {
				if ((i == 0) || ((i == 1) && (j == 0))) {
					ii = 3 * i + j - 1;
				} else {
					ii = 3 * i + j - 2;
				}
				memcpy(codint1[ii], cod1, 3);
				*(codint1[ii] + i) = cod2[i];
				memcpy(codint2[ii], codint1[ii], 3);
				*(codint2[ii] + j) = cod2[j];
				a=num(cod1);
				b=num(codint1[ii]);
				c=num(codint2[ii]);
				d=num(cod2);
				aaa=aa[a];
				aab=aa[b];
				aac=aa[c];
				aad=aa[d];
				rlab=*(rl[aaa]+aab);
				rlbc=*(rl[aab]+aac);
				rlcd=*(rl[aac]+aad);
				like[ii] = rlab*rlbc*rlcd;
			}
	}

	somli = 0;
	for (i = 0; i < 6; i++)
		somli += like[i];
	for (i = 0; i < 6; i++) {
		p[i] = like[i] / somli;
		titv1(cod1, codint1[i], p[i], ti, tv,l);
		titv1(codint1[i], codint2[i], p[i], ti, tv,l);
		titv1(codint2[i], cod2, p[i], ti, tv,l);
	}


}



void prefastlwl(float **rl, float **tl0, float **tl1, float **tl2, float **tti0, float **tti1, float **tti2, float **ttv0, float **ttv1, float **ttv2)
{

	float           l[3], k[3], ti[3], tv[3],cc[3], aaa[3], bb[3], flgseq;
	char             cod1[3], cod2[3];
	int             i, j, ii, jj, nbdiff, cat, pos[3], aa[64], n1, n2, n3;
	void            titv2(char *, char *, float *, float *, float *, int *, float **, int *pos);
	void            titv3(char *, char *, float *, float *, float *, int *, float **);
	void            titv1(char *, char *, float, float *, float *, float *);
	float		 minrl;

/* code des acides amines:
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20    0
F W Y H M L I V P  C  A  G  T  S  Q  N  K  R  E  Q  stop
*/

	aa[0] = 17;/* aaa K */
	aa[1] = 16;/* aac N */
	aa[2] = 17;/* aag K */
	aa[3] = 16;/* aat N */
	aa[4] = 13;/* aca T */
	aa[5] = 13;/* acc T */
	aa[6] = 13;/* acg T */
	aa[7] = 13;/* act T */
if(code_mt)
	aa[8] = 0;/* aga * */
else
	aa[8] = 18;/* aga R */
	aa[9] = 14;/* agc S */
if(code_mt)
	aa[10] = 0;/* agg * */
else
	aa[10] = 18;/* agg R */
	aa[11] = 14;/* agt S */
if(code_mt)
	aa[12] = 5;/* ata M */
else
	aa[12] = 7;/* ata I */
	aa[13] = 7;/* atc I */
	aa[14] = 5;/* atg M */
	aa[15] = 7;/* att I */
	aa[16] = 15;
	aa[17] = 4;
	aa[18] = 15;
	aa[19] = 4;
	aa[20] = 9;
	aa[21] = 9;
	aa[22] = 9;
	aa[23] = 9;
	aa[24] = 18;
	aa[25] = 18;
	aa[26] = 18;
	aa[27] = 18;
	aa[28] = 6;
	aa[29] = 6;
	aa[30] = 6;
	aa[31] = 6;
	aa[32] = 19;
	aa[33] = 20;
	aa[34] = 19;
	aa[35] = 20;
	aa[36] = 11;
	aa[37] = 11;
	aa[38] = 11;
	aa[39] = 11;
	aa[40] = 12;
	aa[41] = 12;
	aa[42] = 12;
	aa[43] = 12;
	aa[44] = 8;
	aa[45] = 8;
	aa[46] = 8;
	aa[47] = 8;
	aa[48] = 0;/* taa * */
	aa[49] = 3;/* tac Y */
	aa[50] = 0;/* tag * */
	aa[51] = 3;/* tat Y */
	aa[52] = 14;/* tca S */
	aa[53] = 14;/* tcc S */
	aa[54] = 14;/* tcg S */
	aa[55] = 14;/* tct S */
if(code_mt)
	aa[56] = 2;/* tga W */
else
	aa[56] = 0;/* tga * */
	aa[57] = 10;/* tgc */
	aa[58] = 2;/* tgg W */
	aa[59] = 10;/* tgt */
	aa[60] = 6;/* tta */
	aa[61] = 1;/* ttc */
	aa[62] = 6;/* ttg */
	aa[63] = 1;/* ttt */

/* ajoute par M. Gouy */
/* calcul minrl = val minimale du tableau rl */
minrl=rl[1][1];
for(i=1; i<=20; i++)
	for(j=i+1; j<=20; j++)
		if(rl[i][j] < minrl ) minrl=rl[i][j];
/* chargement rl[0][i] et rl[i][0] avec minrl correspond a aa = stop */
	for(i= 0; i<=20; i++) rl[0][i] = rl[i][0] = minrl;


	for (i = 0; i < 63; i++) {

		for (j = i; j < 64; j++) {
	
	
			for(ii=0;ii<3;ii++){
				l[ii]=ti[ii]=tv[ii]=0;
			}


			n1 = i / 16;
			n2 = (i - 16 * n1) / 4;
			n3 = i - 16 * n1 - 4 * n2;
			cod1[0] = 'A';
			if (n1 == 1)
				cod1[0] = 'C';
			if (n1 == 2)
				cod1[0] = 'G';
			if (n1 == 3)
				cod1[0] = 'T';
			cod1[1] = 'A';
			if (n2 == 1)
				cod1[1] = 'C';
			if (n2 == 2)
				cod1[1] = 'G';
			if (n2 == 3)
				cod1[1] = 'T';
			cod1[2] = 'A';
			if (n3 == 1)
				cod1[2] = 'C';
			if (n3 == 2)
				cod1[2] = 'G';
			if (n3 == 3)
				cod1[2] = 'T';

			n1 = j / 16;
			n2 = (j - 16 * n1) / 4;
			n3 = j - 16 * n1 - 4 * n2;
			cod2[0] = 'A';
			if (n1 == 1)
				cod2[0] = 'C';
			if (n1 == 2)
				cod2[0] = 'G';
			if (n1 == 3)
				cod2[0] = 'T';
			cod2[1] = 'A';
			if (n2 == 1)
				cod2[1] = 'C';
			if (n2 == 2)
				cod2[1] = 'G';
			if (n2 == 3)
				cod2[1] = 'T';
			cod2[2] = 'A';
			if (n3 == 1)
				cod2[2] = 'C';
			if (n3 == 2)
				cod2[2] = 'G';
			if (n3 == 3)
				cod2[2] = 'T';




			nbdiff = 0;
			pos[0] = pos[1] = pos[2] = 0;
			if (cod1[0] != cod2[0]) {
				nbdiff++;
				pos[0] = 1;
			}
			if (cod1[1] != cod2[1]) {
				nbdiff++;
				pos[1] = 1;
			}
			if (cod1[2] != cod2[2]) {
				nbdiff++;
				pos[2] = 1;
			}
			if (nbdiff != 2)
				for (jj = 0; jj < 3; jj++)
					if (pos[jj] == 0) {
						l[catsite(cod1[0], cod1[1], cod1[2], jj + 1)] += 0.5;
						l[catsite(cod2[0], cod2[1], cod2[2], jj + 1)] += 0.5;
					}
			if (nbdiff == 1)
				titv1(cod1, cod2, 1.0, ti, tv, l);
			if (nbdiff == 2)
				titv2(cod1, cod2, ti, tv, l, aa, rl, pos);
			if (nbdiff == 3)
				titv3(cod1, cod2, ti, tv, l, aa, rl);
			
			*(tl0[i]+j)=*(tl0[j]+i)=l[0];
			*(tl1[i]+j)=*(tl1[j]+i)=l[1];
			*(tl2[i]+j)=*(tl2[j]+i)=l[2];
			*(tti0[i]+j)=*(tti0[j]+i)=ti[0];
			*(tti1[i]+j)=*(tti1[j]+i)=ti[1];
			*(tti2[i]+j)=*(tti2[j]+i)=ti[2];
			*(ttv0[i]+j)=*(ttv0[j]+i)=tv[0];
			*(ttv1[i]+j)=*(ttv1[j]+i)=tv[1];
			*(ttv2[i]+j)=*(ttv2[j]+i)=tv[2];

		}
	}
	
return;	

}


void reresh(char** seq, int nbseq, int option){

/* Si option = 0, toutes les positions avec au moins un gap sont eliminees */
	

  int lgseq, l, drapeau, i, j, k;
  char **seqref; 

   seqref=(char **)malloc(nbseq*sizeof(char *));
  
   lgseq=strlen(seq[1]);

   for(i=0;i<nbseq;i++)
     if ((seqref[i]=(char*)malloc(lgseq*sizeof(char)))==NULL){
       error("Erreur d'allocation");
     }

	l=-1;
	if (option==0){
		for(i=0;i<lgseq;i++){
			drapeau=0;
			for(j=0;j<nbseq;j++){
				if (*(seq[j]+i)=='-') drapeau=1;
			}
			if (drapeau==0){
				l++;
				for(k=0;k<nbseq;k++) *(seqref[k]+l)=*(seq[k]+i);
			}	
		}
	}
	else{
		for(i=0;i<lgseq;i++){
			drapeau=0;
			for(j=0;j<nbseq;j++){
				if (*(seq[j]+i)!='-') {
					drapeau=1;
					break;
				}
			}
			if (drapeau==1){
				l++;
				for(k=0;k<nbseq;k++) *(seqref[k]+l)=*(seq[k]+i);
			}		
		}
	}
	for(i=0;i<nbseq;i++){
		for (j=l+1;j<lgseq;j++) {
			*(seqref[i]+j)='\0';
		}
	}
	for (i=0;i<nbseq;i++) {
		for (j=0;j<lgseq;j++){
			*(seq[i]+j)=*(seqref[i]+j);
		}
	}

	
}
