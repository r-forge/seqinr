#include "dir_acnuc.h"
#include <ctype.h>

#define TVALS 10
struct lvals {
	int val[TVALS];
	struct lvals *next;
	} ;
typedef struct lvals LVALS;

/* globals */
int *blist;
LVALS **smjlistes, **speclistes, **hostlistes, **keylistes, **sublistes,
	**acclistes, **biblistes, *meres, *disparues;


/* included routines */
void addval(LVALS **point, int val);
void processlng(DIR_FILE *kan, int numrec, int offset, LVALS *listvals);
void process_shrt(DIR_FILE *kan,int numrec, LVALS *listvals);
void updateklng(int totsmj, int totspec, int totkey, int totsub);
void updatekshrt(int totacc, int totbib);
void lvals_to_blist(LVALS *plist, int *blist, int *p_mini, int *pmaxi);
int proc_lng_list(int debut, LVALS *plist);
int proc_shrt_list(int debut, LVALS *plist);
/* repetition sous forme variante */
int mdlng(DIR_FILE *kan,int numrec, int offset, int val, int *new);
int addlng(int point, int val);
int suplng(int point, int val);
int mdshrt(DIR_FILE *kan, int numrec, int offset, int val, int *newplist);
int proc_mmap_options(int argc, char **argv);

extern void *mycalloc(int n, size_t taille);

#ifdef unix
int check_term(void);
void write_quick_meres(void);
#endif


int main(int argc, char *argv[])
{
FILE *list_file;
char fname[100], sname[100], *pcar, 
	*lng_buff, *shrt_buff, *sub_buff, *spec_buff, *key_buff, *loc_buff;
int l, isub, big_buffer_size;
int totsmj, totkey, totspec, totacc, totbib, totsub, use_mmap;

fname[0] = 0;
#ifdef unix
if(argc  >= 2) strcpy(fname, argv[1]);
#endif
if(fname[0] == 0) {
	printf("File containing list of sequences to suppress? ");
	gets(fname);
	}
list_file = fopen(fname ,"r");
if(list_file == NULL) {
	printf("File not found: %s\n", fname);
	exit(ERREUR);
	}

dir_acnucopen("WP");
use_mmap = proc_mmap_options(argc, argv);
if(use_mmap) {
	dir_set_mmap(kloc);
	dir_set_mmap(kshrt);
	dir_set_mmap(ksub);
	dir_set_mmap(kbib);
	dir_set_mmap(kacc);
	if( !(nbrf || swissprot) )dir_set_mmap(kext);
	dir_set_mmap(kkey);
	dir_set_mmap(kspec);
	dir_set_mmap(ksmj);
	}

totsmj = read_first_rec(ksmj,NULL);
totspec = read_first_rec(kspec,NULL);
totkey = read_first_rec(kkey,NULL);
totbib = read_first_rec(kbib, NULL);
totacc = read_first_rec(kacc, NULL);
totsub = read_first_rec(ksub, NULL);

sublistes = (LVALS **)mycalloc(totsub+1, sizeof(LVALS *));
smjlistes = (LVALS **)mycalloc(totsmj+1, sizeof(LVALS *));
speclistes = (LVALS **)mycalloc(totspec+1, sizeof(LVALS *));
hostlistes = (LVALS **)mycalloc(totspec+1, sizeof(LVALS *));
keylistes = (LVALS **)mycalloc(totkey+1, sizeof(LVALS *));
acclistes = (LVALS **)mycalloc(totacc+1, sizeof(LVALS *));
biblistes = (LVALS **)mycalloc(totbib+1, sizeof(LVALS *));
meres = (LVALS *)mycalloc(1, sizeof(LVALS));
disparues = (LVALS *)mycalloc(1, sizeof(LVALS));
blist = (int *)mycalloc(lenw, sizeof(int));

#ifdef unix
check_term(); /* possibilite d'interruption */
#endif

printf("Suppressing sequences listed in file %s\n", fname);

while(fgets(sname,sizeof(sname),list_file) != NULL) {
	l = strlen(sname);
	if(sname[l-1]=='\n') sname[l-1] = '\0';
	pcar = sname; while(*pcar != 0) { *pcar = toupper(*pcar); ++pcar; }
	isub = isenum(sname);
	if(isub==0)
		printf("Seq does not exist:%s\n",sname);
	else	{
		delseq(isub);
		printf("%s suppressed\n",sname);
		}
	fflush(stdout);
#ifdef unix
	if( check_term() ) {
		printf("Clean interruption\n");
		break;
		}
#endif
	}

printf("Updating short lists\n"); fflush(stdout);
updatekshrt(totacc, totbib);
if(use_mmap) {
	dir_set_normal(kshrt);
	dir_set_normal(kloc);
	dir_set_mmap(klng);
	}
printf("Updating long lists\n"); fflush(stdout);
updateklng(totsmj, totspec, totkey, totsub);

dir_acnucclose();
#ifdef unix
write_quick_meres();
#endif
printf("Normal end\n");
exit(0);
} /* end of main */


void addval(LVALS **debut, int val)
{
int i= -1;

while(++i < TVALS)
	if( (*debut)->val[i] == 0) break;
if(i >= TVALS) {
	LVALS *point;
	point = mycalloc(1, sizeof(LVALS));
	point->next = *debut;
	*debut = point;
	i = 0;
	}
(*debut)->val[i] = val;
}


void processlng(DIR_FILE *kan, int numrec, int offset, LVALS *listvals)
{
char *pbuffer;
int *point, deblng, vide, val;

if(kan == ksub) {
	readsub(numrec);
	point= &(psub->pext);
	deblng = abs(*point);
	pbuffer = (char *)psub;
	}
else if(kan == kspec || kan == kkey) {
	dir_read(kan,numrec,1,pspec);
	if(offset==2)
		point= &(pspec->plsub);
	else
		point= &(pspec->plhost);
	deblng= *point;
	pbuffer = (char *)pspec;
	}
else if(kan == ksmj) {
	readsmj(numrec);
	point= &(psmj->plong);
	deblng= *point;
	pbuffer = (char *)psmj;
	}
else if(kan == klng) {
	deblng = offset;
	}
if(deblng == 0) return;
vide = proc_lng_list(deblng, listvals);
if(vide) {
	if(kan != klng)	{
		*point = 0;
		dir_write(kan, numrec, 1, pbuffer);
		}
	else	{
		memset(plng, 0, lrlng);
		writelng(deblng);
		}
	}
}


void process_shrt(DIR_FILE *kan, int numrec, LVALS *listvals)
{
int deblist, vide;

if(kan == kbib) {
	readbib(numrec);
	deblist = pbib->plsub;
	}
else	{
	readacc(numrec);
	deblist = pacc->plsub;
	}
if(deblist == 0)return;
vide = proc_shrt_list(deblist, listvals);
if(vide) {
	if(kan == kbib) {
		pbib->plsub = 0;
		writebib(numrec);
		}
	else	{
		pacc->plsub = 0;
		writeacc(numrec);
		}
	}
}


void updateklng(int totsmj, int totspec, int totkey, int totsub)
{
int numrec, i, val, loop, deblng, totlng;
for(numrec = 2; numrec <= totsmj; numrec++) {
	if(smjlistes[numrec] != NULL)
		processlng(ksmj,numrec,1,smjlistes[numrec]);
	}
for(numrec = 2; numrec <= totspec; numrec++) {
	if(speclistes[numrec] != NULL)
		processlng(kspec,numrec,2,speclistes[numrec]);
	if(hostlistes[numrec] != NULL)
		processlng(kspec,numrec,6,hostlistes[numrec]);
	}
for(numrec = 2; numrec <= totkey; numrec++) {
	if(keylistes[numrec] != NULL)
		processlng(kkey,numrec,2,keylistes[numrec]);
	}
for(numrec = 2; numrec <= totsub; numrec++) {
	if(sublistes[numrec] != NULL)
		processlng(ksub,numrec,3,sublistes[numrec]);
	}
/* liste des meres */
processlng(klng, 0, 2, meres);
/* liste des disparues */
lngbit(3, blist);
while(disparues != NULL) {
	for(i = 0; i<TVALS; i++) 
		if(disparues->val[i]) bit1(blist,disparues->val[i]);
	disparues = disparues->next;
	}
val = 0; loop = -1; deblng = 3;
totlng = read_first_rec(klng,NULL);
while((val = irbit(blist,val,nseq)) != 0) {
	if(++loop >= SUBINLNG) {
		plng->next = ++totlng;
		dir_write(klng, deblng, 1, plng);
		deblng = plng->next;
		loop = 0;
		}
	plng->sub[loop] = val;
	}
while(++loop < SUBINLNG) plng->sub[loop] = 0;
plng->next = 0;
dir_write(klng, deblng, 1, plng);
if(deblng == totlng) write_first_rec(klng, totlng, 0);
}


void updatekshrt(int totacc, int totbib)
{
int numrec;
for(numrec = 2; numrec <= totbib; numrec++) {
	if(biblistes[numrec] != NULL)
		process_shrt(kbib, numrec, biblistes[numrec]);
	}
for(numrec = 2; numrec <= totacc; numrec++) {
	if(acclistes[numrec] != NULL)
		process_shrt(kacc, numrec, acclistes[numrec]);
	}
}


int mdlng(DIR_FILE *kan,int numrec, int offset, int val, int *new)
{
if(offset > 0) {
	printf("Logic error in mdlng\n");
	exit(ERREUR);
	}
offset = abs(offset);
if(kan == ksmj) {
	if(smjlistes[numrec]==NULL) {
		smjlistes[numrec] = (LVALS *)mycalloc(1, sizeof(LVALS));
		}
	addval( &(smjlistes[numrec]), val);
	}
else if(kan == kkey) {
	if(keylistes[numrec]==NULL) {
		keylistes[numrec] = (LVALS *)mycalloc(1, sizeof(LVALS));
		}
	addval( &(keylistes[numrec]), val);
	}
else if(kan == kspec) {
	if(offset == 2) {
		if(speclistes[numrec] == NULL) {
			speclistes[numrec] = (LVALS *)mycalloc(1,sizeof(LVALS));
			}
		addval( &(speclistes[numrec]), val);
		}
	else {
		if(hostlistes[numrec] == NULL) {
			hostlistes[numrec] = (LVALS *)mycalloc(1,sizeof(LVALS));
			}
		addval( &(hostlistes[numrec]), val);
		}
	}
else if(kan == ksub) {
	if(sublistes[numrec] == NULL) {
		sublistes[numrec] = (LVALS *)mycalloc(1,sizeof(LVALS));
		}
	addval( &(sublistes[numrec]), val);
	}

}


int addlng(int point, int val)
{
if(point != 3) {
	printf("Logic error in addlng\n");
	exit(ERREUR);
	}
addval( &disparues, val);
}


int suplng(int point, int val)
{
if(point != 2) {
	printf("Logic error in suplng\n");
	exit(ERREUR);
	}
addval( &meres, val);
}


int mdshrt(DIR_FILE *kan, int numrec, int offset, int val, int *newplist)
{
if(kan == kbib) {
	if(biblistes[numrec] == NULL) {
		biblistes[numrec] = (LVALS *)mycalloc(1, sizeof(LVALS));
		}
	addval( &(biblistes[numrec]), val);
	}
else if(kan == kacc) {
	if(acclistes[numrec] == NULL) {
		acclistes[numrec] = (LVALS *)mycalloc(1,sizeof(LVALS));
		}
	addval( &(acclistes[numrec]), val);
	}
else 	{
	printf("Logic error in mdshrt\n");
	exit(ERREUR);
	}
}


void lvals_to_blist(LVALS *plist, int *blist, int *p_mini, int *p_maxi)
{
int num, mini, maxi, val, mot;

mini = -1;
while(plist != NULL) {
	for(num = 0; num < TVALS; num++) {
		if( (val = plist->val[num]) == 0) break;
		mot = (val - 1) / lmot;
		if(mini == -1) {
			mini = maxi = mot;
			blist[mot] = 0;
			}
		else if( mot > maxi) {
			while(maxi < mot) blist[++maxi] = 0;
			}
		else if(mot < mini) {
			while(mini > mot) blist[--mini] = 0;
			}
		bit1(blist, val);
		}
	plist = plist->next;
	}
*p_mini = mini * lmot + 1; 
*p_maxi = (maxi + 1) * lmot;
}


int proc_lng_list(int debut, LVALS *plist)
{
static struct rlng aux;
int mini, maxi, val, newdebut, lastaux, vide, num;

lvals_to_blist(plist, blist, &mini, &maxi);
lastaux = -1; newdebut = debut; vide = TRUE;
memset(&aux, 0, SUBINLNG * sizeof(int) );
while(debut != 0) {
	readlng(debut);
	if(newdebut == debut) aux.next = plng->next;
	for(num = 0; num < SUBINLNG; num++) {
		if( (val = plng->sub[num]) == 0) continue;
		if( val < mini || val > maxi || !testbit(blist, val) ) {
			if(lastaux >= SUBINLNG - 1) {
				if( dir_write(klng, newdebut, 1, &aux) ) 
					dir_writeerr(klng, newdebut);
				newdebut = aux.next;
				if(newdebut == debut)
					aux.next = plng->next;
				else
					if(dir_read(klng, newdebut, 1, &aux) != 1) 
						dir_readerr(klng, newdebut);
				memset(&aux, 0, SUBINLNG * sizeof(int) );
				lastaux = -1;
				}
			aux.sub[++lastaux] = val;
			vide = FALSE;
			}
		}
	debut = plng->next;
	}
if( !vide ) {
	aux.next = 0;
	if( dir_write(klng, newdebut, 1, &aux) ) dir_writeerr(klng, newdebut);
	}
return vide;
}


int proc_shrt_list(int debut, LVALS *plist)
{
int mini, maxi, val, next, newdebut, preval;

lvals_to_blist(plist, blist, &mini, &maxi);
readshrt(debut);
val = pshrt->val;
if(val >= mini && val <= maxi && testbit(blist, val) ) {
	do	{
		newdebut = pshrt->next;
		if(newdebut == 0) return TRUE;
		readshrt(newdebut);
		val = pshrt->val;
		}
	while(val >= mini && val <= maxi && testbit(blist, val) );
	writeshrt(debut);
	}
while(TRUE) {
	readshrt(debut);
	preval = pshrt->val;
	next = pshrt->next;
	if(next == 0) break;
	readshrt(next);
	val = pshrt->val;
	if(val >= mini && val <= maxi && testbit(blist, val) ) {
		pshrt->val = preval;
		writeshrt(debut);
		}
	else	{
		debut = next;
		}
	}
return FALSE;
}


int proc_mmap_options(int argc, char **argv)
{
int num = 0;
int use_mmap = FALSE;

#ifdef unix
while (++num < argc) {
	if(strcmp(argv[num], "-mmap") == 0) use_mmap = TRUE;
	}
#endif
return use_mmap;
}
