#include "dir_acnuc.h"
#include <limits.h>
#include <string.h>

#ifdef unix
#define MODULO_MINI 1000000
#else
#define MODULO_MINI 100000
#endif
#define PLUSGRAND(a,b) (a > b ? a : b)
#define CROIX "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"


/* interface publique */
/* initialisation du système : appeler avant tout ajout de numero d'accession
retour: NULL ssi erreur
*/
void *init_acchash(void);

/* recherche ou création d'accession number
name: acc sans espace terminal et avec 0 final
create: TRUE pour création si inexistant
        FALSE pour recherche sans création
retour: != 0 si trouvé ou créé
           0 si pas trouvé et non créé (création non demandée ou impossible cause mémoire)
*/
int find_cre_acc(char *name, int create);

/* fonctions privees */
static void accaddhash(int recnum);
static int acchashfind(char *nom);
static int acchash(char *buffer, int modulo);


/* external procedures */
extern void *acc_binary_tree;
int find_key(char *name, void *tree, int create_if_not_found, int **p_aux_data);


/* global data */
static int *acchashcode, *acctabnext, modulo, maxacctab;

/* pour le tester  /

int main(void)
{
void *p;
int total, i, j;
char nom[20];

dir_acnucopen("RO");
/* dir_set_mmap(kacc);
/
p = init_acchash();
if(p == NULL) exit(ERREUR);

total = read_first_rec(kacc, 0); 
for(i = 2; i <= total; i++) {
if(i % 50000 == 0) printf("i=%d\n", i);
	readacc(i);
	memcpy(nom, pacc->name, ACC_LENGTH);
	nom[ACC_LENGTH] = 0; trim_key(nom);
	j = find_cre_acc(nom, FALSE);
	if(j != i) {
		printf("pb");
		}
	}
return 0;
}

/  fin du test  */



int find_cre_acc(char *name, int create)
{
int num, total, acclastsorted;

if(acc_binary_tree != NULL) { /* old procedure */
	return find_key(name, acc_binary_tree, create, NULL);
	}

num = acchashfind(name);
if(num != 0 || ! create) return num;

total = read_first_rec(kacc, &acclastsorted); 
if(total >= maxacctab) {
	int plus, *p;
	plus = (maxacctab + 100) * 1.3;
	p = (int *)realloc(acctabnext, (plus + 1)*sizeof(int) );
	if(p == NULL) return 0;
	acctabnext = p;
	maxacctab = plus;
	}
total++;
padtosize(pacc->name, name, ACC_LENGTH);
pacc->plsub = 0;
writeacc(total);
write_first_rec(kacc, total, acclastsorted);
accaddhash(total);
return total;
}


void *init_acchash(void)
{
int total, i;

total = read_first_rec(kacc, NULL); 
modulo = PLUSGRAND(total/4, MODULO_MINI);
acchashcode = (int *) calloc(modulo , sizeof(int));
if(acchashcode ==  NULL) return NULL;
maxacctab = 1.5 * total + 100;
acctabnext = (int *) malloc( (maxacctab + 1) * sizeof(int));
if(acctabnext ==  NULL) return NULL;
for(i = 2; i <= total; i++) {
	readacc(i);
	if(memcmp(pacc->name, CROIX, ACC_LENGTH) == 0) continue;
	accaddhash(i);
	}
return (void *)acchashcode;
}


static void accaddhash(int recnum)
{
int num;

readacc(recnum);
pacc->name[ACC_LENGTH] = 0;
trim_key(pacc->name);
num = acchash(pacc->name, modulo);
if(acchashcode[num] == 0) {
	acchashcode[num] = recnum;
	acctabnext[recnum] = 0;
	}
else	{
	acctabnext[ recnum ] = acchashcode[num];
	acchashcode[num] = recnum;
	}
}
	

static int acchashfind(char *nom)
{
static char buffer[ACC_LENGTH+1];
int num;

num = acchash(nom, modulo);
num = acchashcode[num];
if(num == 0) return 0;
padtosize(buffer, nom, ACC_LENGTH);
while(num != 0) {
	readacc(num);
	if(strncmp(pacc->name, buffer, ACC_LENGTH) == 0) return num;
	num = acctabnext[num];
	} 
return 0;
}


static int acchash(char *buffer, int modulo)
{
/* inspire de celui de java
sans le % pour maintenir valeur positive le resultat est un peu meilleur
avec mnemos (et meilleur que CRC32) et un peu moins bon avec keyw
mais risque de depassement entier qui peut-etre peut arreter le programme 
*/
const int maxi = (INT_MAX - 256) / 37;
int c, h = 0;
    while( ( c = *(buffer++) ) != 0 ) {
	h = h * 37 + c;
	if(h >= maxi) h %= modulo;
	}
    return abs(h) % modulo;
}

