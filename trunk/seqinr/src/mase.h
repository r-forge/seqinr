#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#define MAXNSEQS  256
#define MAXMNMASE 30
#define MAXSTRING 5000


struct SEQMASE
{
    char mn[MAXMNMASE];
    char *com;
    char *seq;
    int lg;
};

char *check_alloc(int nbrelt, int sizelt);
void rem_blank(char *string);
void free_mase(struct SEQMASE * aln, int nbsq);
