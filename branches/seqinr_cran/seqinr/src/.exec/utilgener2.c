#include "dir_acnuc.h"

typedef struct {
	off_t fpos;
	off_t gcg_fpos;
	int modified;
	char date[9];
	} s_to_do_extra;
	
typedef struct _fille_request {
	char *fille;
	long pannots;
	int div;
	int type;
	struct _fille_request *next;
	} fille_request;

#define TVALS 5
typedef struct lvals {
	int val[TVALS];
	struct lvals *next;
	} LVALS;


/* local functions */
void *calc_meres_to_do(FILE *in_flat, int current_div, FILE *seq_file,
	FILE *log_file);
static off_t find_in_seq_file(FILE *seq_file, char *name, int *erreur);
off_t fpos_from_to_do(void *loop_to_do);
off_t gcg_fpos_from_to_do(void *loop_to_do);
int compare_date_and_length(char *line, size_t lline, int *p_oldseq, 
	char *seqname, char *newdate, FILE *in_flat, int current_div, 
	off_t addr_info, off_t addr_seq);
int compare_dates(char *nouvelle, char *ancienne);
void calc_gbk_date(char *gbk, char *newdate);
void remove_modified_seqs(void *dynlist_meres_to_do, 
	void **p_dynlist_mere_with_request, FILE *log_file);
static int fast_remove(int nsub, void *dynlist_mere_with_request,
	LVALS **smjlistes, LVALS **speclistes, LVALS **hostlistes, 
	LVALS **keylistes, LVALS **sublistes, 
	LVALS **acclistes, LVALS **biblistes,
	LVALS **p_meres, LVALS **p_disparues);
static int delsubseq(int isub, int mere, 
	LVALS **smjlistes, LVALS **keylistes, LVALS **sublistes, 
	LVALS **p_disparues);
static void free_after_fast_remove(LVALS **smjlistes, LVALS **speclistes, 
	LVALS **hostlistes, 
	LVALS **keylistes, LVALS **sublistes, 
	LVALS **acclistes, LVALS **biblistes, LVALS *meres, LVALS *disparues);
static void free_lval_tab(LVALS **tab, int total);
static void free_lvals(LVALS *list);
static LVALS *addval(LVALS *debut, int val);
static void processlng(DIR_FILE *kan, int numrec, int offset, LVALS *listvals, 
	int *auxliste);
static void process_shrt(DIR_FILE *kan, int numrec, LVALS *listvals, 
	int *auxlist);
static void updateklng(int totsmj, int totspec, int totkey, int totsub, 
	LVALS **smjlistes, LVALS **speclistes, LVALS **hostlistes,
	LVALS **keylistes, LVALS **sublistes,
	LVALS *meres, LVALS *disparues, int *auxliste);
static void updatekshrt(int totacc, int totbib, LVALS **biblistes, 
	LVALS **acclistes, int *auxliste);
static void lvals_to_blist(LVALS *plist, int *blist, int *p_mini, int *p_maxi);
static int proc_lng_list(int debut, LVALS *plist, int *auxlist);
static int proc_shrt_list(int debut, LVALS *plist, int *auxlist);
void sup_fille_request_from_mere(char *nom_mere, 
	void *dynlist_mere_with_request);
int add_fille_request(int fille, void *dynlist_mere_with_request);
void process_fille_request(fille_request *p_request, char *location, 
	char *qualifiers, int max_loc_qual, int *p_fille_num, 
	int *p_tot_filles_cre, FILE *log_file);
void process_genetic_code(int mere);


/* external prototypes */
void *init_dynlist(void);
char *find_dynlist(char *name, void *tree, int create_if_not_found);
int sizeof_dynlist(void *arbre);
int remove_dynlist(char *name, void *tree);
char *next_dynlist(void *arbre, void **p_noeud);
int optimize_dynlist(void *tree); 
void free_dynlist(void *tree);
char *find_g_dynlist(char *name, void *tree, int create_if_not_found, void **node);
void *dynlist_get_extra(void *node);
void dynlist_set_extra(void *node, void *extra);
int trim_key(char *);
int read_addr_loc_qualif(long faddr, int div, 
	char *location, int maxlocat, char *type, char *qualif, int maxqualif);
int mkextr(char *location, int seq_num, int *p_fille_num, char *name,
	int smj_type, char *feat_name, long info_addr, int div, char *qualif);
void goffset(int point,int *div,int *offset);
#ifdef unix
int check_term(void); /* gestion de la possibilite d'interruption propre */
#endif


void *calc_meres_to_do(FILE *in_flat, int current_div, FILE *seq_file,
	FILE *log_file)
/* return NULL if not enough memeory */
{
void *meres_to_do, *node;
char ligne[200], seqname[L_MNEMO + 1], newdate[9];
int old_seq, newer;
s_to_do_extra *extra, *old_extra;
off_t addr_info, gcg_addr;

if( (meres_to_do = init_dynlist() ) == NULL) return NULL;
while( TRUE) {
	addr_info = ftello(in_flat);
	if( fgets(ligne, sizeof(ligne), in_flat) == NULL ) break;
	if( ( genbank && strncmp(ligne, "LOCUS ", 6) != 0 ) ||
		( embl && strncmp(ligne, "ID", 2) != 0) ) continue;
	if( !flat_format ) { /* cas GCG */
		char *p, *q; int erreur;
		if(genbank) p = ligne + 5;
		else p =  ligne + 2;
		while( *p == ' ') p++;
		q = strchr(p, ' ');
		memcpy(seqname, p, q - p); seqname[q - p] = 0;
		gcg_addr = find_in_seq_file(seq_file, seqname, &erreur);
		if(erreur) {
			fprintf(log_file, 
			"Error in synchro betwen .ref and .seq files for %s\n"
			"Bank is left unchanged\n",
			seqname);
			exit(ERREUR);
			}
		}
	newer = compare_date_and_length(ligne, sizeof(ligne), &old_seq, 
		seqname, newdate, in_flat, current_div, addr_info, gcg_addr);
	if( flat_format ) { /* sauter tout le reste de la seq, sauf si GCG */
		do	fgets(ligne, sizeof(ligne), in_flat);
		while(strncmp(ligne, "//", 2) != 0);
		}
	/* seq deja connue et non plus recente */
	if(old_seq != 0 && !newer) continue; 
	if( find_g_dynlist(seqname, meres_to_do, TRUE, &node) == NULL) 
		return NULL;
	extra = (s_to_do_extra *)malloc(sizeof(s_to_do_extra));
	if(extra == NULL) return NULL;
	extra->fpos = addr_info;
	if( !flat_format ) extra->gcg_fpos = gcg_addr;
	memcpy(extra->date, newdate, 9);
	if(old_seq == 0) { /* nouvelle sequence */
		extra->modified = FALSE;
		}
	else	{ /* seq d'acnuc plus ancienne */
		extra->modified = TRUE;
		}
	old_extra = dynlist_get_extra(node);
	if(old_extra != NULL) { /* deja seq de meme nom connue ds dynlist */
		if( compare_dates(newdate, old_extra->date) >= 0) {
			free(old_extra); /* nouvelle est plus recente */
			dynlist_set_extra(node, extra);
			}
		else free(extra); /* precedente plus recente */
		}
	else /* pas de seq deja connue de ce nom */
		dynlist_set_extra(node, extra);
	}
return meres_to_do;
}


static off_t find_in_seq_file(FILE *seq_file, char *name, int *erreur)
/* chercher la sequence name ou name_0 dans fichier .seq au format GCG */
{
int l;
off_t debnuc;
static char ligne[100];
char target1[20], target2[20];

/* chercher "nom " ou "nom_0" car on peut avoir nom, nom_0, nom_00 */
strcpy(target1, name); strcat(target1, " "); 
strcpy(target2, name); strcat(target2, "_0");
l = strlen(target1);
do	{
	debnuc = ftello(seq_file);
	if( fgets(ligne, sizeof(ligne), seq_file) == NULL ) {
		*erreur = TRUE;
		return 0;
		}
	}
while(strncmp(ligne, ">>>>", 4) != 0 || 
	( strncmp(ligne + 4, target1, l) != 0 && strncmp(ligne + 4, target2, l + 1) != 0 ) );
*erreur = FALSE;
return debnuc;
}


off_t fpos_from_to_do(void *loop_to_do)
{
s_to_do_extra *extra;
extra = dynlist_get_extra(loop_to_do);
return extra->fpos;
}


off_t gcg_fpos_from_to_do(void *loop_to_do)
{
s_to_do_extra *extra;
extra = dynlist_get_extra(loop_to_do);
return extra->gcg_fpos;
}


int compare_date_and_length(char *line, size_t lline, int *p_oldseq, 
	char *seqname, char *newdate, FILE *in_flat, int current_div, 
	off_t addr_info, off_t addr_seq)
{
char *p, *q, *gbk_date, *gbk_name;
int oldseq, length, compar;

if(genbank) {
	length = decode_locus_rec(line, &gbk_name, NULL, NULL, NULL, 
				&gbk_date);
	strcpy(seqname, gbk_name);
	}
else {
	p = line + 2;
	while(*p == ' ') p++;
	q = strchr(p, ' '); *q = 0; strcpy(seqname, p); *q = ' ';
	}
majuscules(seqname);
oldseq = isenum(seqname);
*p_oldseq = oldseq;
if(oldseq == 0) { /* nouvelle seq */
	return FALSE;
	}
else	{
	/* recherche des dates et longueurs */
	if(genbank) {
		calc_gbk_date(gbk_date, newdate);
		}
	else if(embl || swissprot) {
		if(embl) { 
			p=strchr(line, ';');
			p = strchr(p + 1, ';');
			p = strchr(p + 1, ';');
			sscanf(p + 1, "%d", &length);
			}
		else	{
			p=strchr(line, ';');
			p = strchr(p + 1, ';');
			sscanf(p + 1, "%d", &length);
			}
		/* recherche de la 2eme ligne DT */
		do	fgets(line, lline, in_flat);
		while( strncmp(line,"DT",2) != 0 && strncmp(line,"SQ",2) != 0 );
		strcpy(newdate,"xxxx");
		if(strncmp(line,"DT",2) == 0) {
			fgets(line, lline, in_flat);
			if(swissprot) fgets(line, lline, in_flat);
			if(strncmp(line, "DT", 2) == 0) 
				calc_gbk_date(line+5, newdate);
			}
		}
	readsub(oldseq);
	readloc(psub->plinf);
	if(newdate[0] == 'x') { /* date absente des annots */
		memcpy(newdate, ploc->date, sizeof(ploc->date));
		compar = 0;
		}
	else	compar = compare_dates(ploc->date, newdate);
	if(compar < 0) return TRUE; /* nouvelle seq plus recente que dans indexs */
	else if(compar == 0) { /* meme date */
		int old_div, old_gcg_seq; long old_info;
		seq_to_annots(oldseq, &old_info, &old_div);
		if(old_div != current_div) /* autre division que dans indexs */
			return TRUE;
		else	{ /* meme division dans indexs */
			if(flat_format) {
				if(old_info == addr_info && psub->length == length)
					return FALSE; /* meme addr et meme long */
				else if(old_info > addr_info) {
					/* seq indexee plus loin: verifier si ok */
					p = read_annots(old_info, old_div);
					return (p == NULL || strstr(p, seqname) == NULL);
					}
				else /* seq indexee avant dans fichier */
					return TRUE;
				}
			else	{ /* cas GCG: get pointer to seq in indexs */
				if(big_annots) old_gcg_seq = ploc->pnuc;
				else	goffset( ploc->pnuc, &old_div, &old_gcg_seq);
				/* TRUE des que pointeurs differents */
				return (old_info != addr_info || old_gcg_seq != addr_seq );
				}
			}
		}
	return FALSE;
	}
}


int compare_dates(char *nouvelle, char *ancienne)
{
int test, old_year, new_year;

test = strncmp(nouvelle + 6, ancienne + 6, 2);
if(test != 0) {
	sscanf(nouvelle + 6, "%2d", &new_year);
	sscanf(ancienne + 6, "%2d", &old_year);
	if(new_year >= 50) new_year += 1900;
	else new_year += 2000;
	if(old_year >= 50) old_year += 1900;
	else old_year += 2000;
	test = new_year - old_year;
	return test;
	}
test = strncmp(nouvelle, ancienne, 2);
if(test != 0) return test;
return strncmp(nouvelle + 3, ancienne + 3, 2);
}


void calc_gbk_date(char *gbk, char *newdate)
/* input gbk pointe sur 21-JUL-1994
retour: newdate rempli avec 07/21/94
*/
{
int jour, num_mois, an;
char *p;
static char mois[]="JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC";

sscanf(gbk,"%d",&jour);
sscanf(gbk+9,"%d",&an);
gbk[6] = 0;
p = strstr(mois,gbk+3);
gbk[6] = '-';
if(p!=NULL)  num_mois = (p-mois)/3 + 1;
else 	num_mois = 0;
sprintf(newdate,"%2.2d/%2.2d/%2.2d",num_mois,jour,an);
return;
}


void remove_modified_seqs(void *dynlist_meres_to_do, 
	void **p_dynlist_mere_with_request, FILE *log_file)
{
void *loop;
s_to_do_extra *extra;
int status;
char *p;
int totsmj, totkey, totspec, totacc, totbib, totsub, *auxliste, total,
	interrupted;
LVALS **smjlistes, **speclistes, **hostlistes, **keylistes, **sublistes,
	**acclistes, **biblistes, *meres, *disparues;

/* initialisations */
totsmj = read_first_rec(ksmj,NULL);
totspec = read_first_rec(kspec,NULL);
totkey = read_first_rec(kkey,NULL);
totbib = read_first_rec(kbib, NULL);
totacc = read_first_rec(kacc, NULL);
totsub = read_first_rec(ksub, NULL);

sublistes = (LVALS **)calloc(totsub+1, sizeof(LVALS *));
smjlistes = (LVALS **)calloc(totsmj+1, sizeof(LVALS *));
speclistes = (LVALS **)calloc(totspec+1, sizeof(LVALS *));
hostlistes = (LVALS **)calloc(totspec+1, sizeof(LVALS *));
keylistes = (LVALS **)calloc(totkey+1, sizeof(LVALS *));
acclistes = (LVALS **)calloc(totacc+1, sizeof(LVALS *));
biblistes = (LVALS **)calloc(totbib+1, sizeof(LVALS *));
meres = (LVALS *)calloc(1, sizeof(LVALS));
disparues = (LVALS *)calloc(1, sizeof(LVALS));
auxliste = (int *)malloc(lenw * sizeof(int));
if(sublistes == NULL || smjlistes == NULL || speclistes == NULL || 
	hostlistes == NULL || keylistes == NULL || acclistes == NULL ||
	biblistes == NULL || meres == NULL || disparues == NULL ||
	auxliste == NULL) {
		fprintf(log_file, "Not enough memory to delete modified seqs\n"
			"Base is left clean\n");
		exit(ERREUR);
		}

status = total = 0;
interrupted = FALSE;
*p_dynlist_mere_with_request = init_dynlist();
fprintf(log_file, "Removing updated sequences\n");
loop = NULL;
while( (p = next_dynlist(dynlist_meres_to_do, &loop) ) != NULL ) {
	extra = dynlist_get_extra(loop);
	if(!extra->modified) continue;
	fprintf(log_file, "%s  %s\n", p, extra->date);
	status = fast_remove(isenum(p), *p_dynlist_mere_with_request,
		smjlistes, speclistes, hostlistes, keylistes, sublistes, 
		acclistes, biblistes, &meres, &disparues);
	total++;
	if(status) break;
#ifdef unix
	if( check_term() ) { /* tester si interruption demandee */
		fprintf(log_file, "Clean program interruption\n");
		fflush(log_file);
		interrupted = TRUE;
		break; 
		}
#endif
	}
	
/* operations de fin d'effacement */
if(total > 0) {
	fprintf(log_file, "Updating short lists\n"); fflush(log_file);
	updatekshrt(totacc, totbib, biblistes, acclistes, auxliste);
	fprintf(log_file, "Updating long lists\n"); fflush(log_file);
	updateklng(totsmj, totspec, totkey, totsub, smjlistes, speclistes, 
		hostlistes, keylistes, sublistes, meres, disparues, auxliste);
	}
if(status != 0) {
	fprintf(log_file, "Not enough memory to delete modified seqs\n"
			"Base is left in bad state!\n");
	exit(ERREUR);
	}
/* liberation de memoire */
free(auxliste);
free_after_fast_remove(smjlistes, speclistes, hostlistes, keylistes, sublistes, 
	acclistes, biblistes, meres, disparues);
fprintf(log_file, "%d modified sequences removed\n", total); fflush(log_file);
if(interrupted) exit(0);
}


static void free_after_fast_remove(LVALS **smjlistes, LVALS **speclistes, 
	LVALS **hostlistes, 
	LVALS **keylistes, LVALS **sublistes, 
	LVALS **acclistes, LVALS **biblistes, LVALS *meres, LVALS *disparues)
{
free_lvals(meres);
free_lvals(disparues);
free_lval_tab( smjlistes, read_first_rec(ksmj, NULL) );
free_lval_tab( speclistes, read_first_rec(kspec, NULL) );
free_lval_tab( hostlistes, read_first_rec(kspec, NULL) );
free_lval_tab( keylistes, read_first_rec(kkey, NULL) );
free_lval_tab( sublistes, read_first_rec(ksub, NULL) );
free_lval_tab( acclistes, read_first_rec(kacc, NULL) );
free_lval_tab( biblistes, read_first_rec(kbib, NULL) );
}


static void free_lval_tab(LVALS **tab, int total)
{
int i;
LVALS *list, *next;

for(i = 0; i <= total; i++) free_lvals(tab[i]);
free(tab);
}


static void free_lvals(LVALS *list)
{
LVALS *next;

while(list != NULL) {
	next = list->next;
	free(list);
	list = next;
	}
}


static LVALS *addval(LVALS *debut, int val)
{
int i= -1;

if(debut == NULL) debut = calloc(1, sizeof(LVALS));
if(debut == NULL) return NULL;
while(++i < TVALS)
	if( debut->val[i] == 0) break;
if(i >= TVALS) {
	LVALS *point;
	point = calloc(1, sizeof(LVALS));
	if(point == NULL) return NULL;
	point->next = debut;
	debut = point;
	i = 0;
	}
debut->val[i] = val;
return debut;
}


static void processlng(DIR_FILE *kan, int numrec, int offset, LVALS *listvals, 
	int *auxliste)
{
void *pbuffer;
int *point, deblng, vide, val;

if(kan == ksub) {
	readsub(numrec);
	point= &(psub->pext);
	deblng = abs(*point);
	pbuffer = psub;
	}
else if(kan == kspec || kan == kkey) {
	dir_read(kan,numrec,1,pspec);
	if(offset==2)
		point= &(pspec->plsub);
	else
		point= &(pspec->plhost);
	deblng= *point;
	pbuffer = pspec;
	}
else if(kan == ksmj) {
	readsmj(numrec);
	point= &(psmj->plong);
	deblng= *point;
	pbuffer = psmj;
	}
else if(kan == klng) {
	deblng = offset;
	}
if(deblng == 0) return;
vide = proc_lng_list(deblng, listvals, auxliste);
if(vide && kan != klng){
	*point = 0;
	dir_write(kan, numrec, 1, pbuffer);
	}
}


static void process_shrt(DIR_FILE *kan, int numrec, LVALS *listvals, 
	int *auxlist)
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
vide = proc_shrt_list(deblist, listvals, auxlist);
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


static void updateklng(int totsmj, int totspec, int totkey, int totsub, 
	LVALS ** smjlistes, LVALS **speclistes, LVALS **hostlistes,
	LVALS **keylistes, LVALS **sublistes,
	LVALS *meres, LVALS *disparues, int *auxliste)
{
int numrec, i, val, loop, deblng, totlng;
for(numrec = 2; numrec <= totsmj; numrec++) {
	if(smjlistes[numrec] != NULL)
		processlng(ksmj,numrec,1,smjlistes[numrec], auxliste);
	}
for(numrec = 2; numrec <= totspec; numrec++) {
	if(speclistes[numrec] != NULL)
		processlng(kspec,numrec,2,speclistes[numrec], auxliste);
	if(hostlistes[numrec] != NULL)
		processlng(kspec,numrec,6,hostlistes[numrec], auxliste);
	}
for(numrec = 2; numrec <= totkey; numrec++) {
	if(keylistes[numrec] != NULL)
		processlng(kkey,numrec,2,keylistes[numrec], auxliste);
	}
for(numrec = 2; numrec <= totsub; numrec++) {
	if(sublistes[numrec] != NULL)
		processlng(ksub,numrec,3,sublistes[numrec], auxliste);
	}
/* liste des meres */
processlng(klng, 0, 2, meres, auxliste);
/* liste des disparues */
lngbit(3, auxliste);
while(disparues != NULL) {
	for(i = 0; i<TVALS; i++) 
		if(disparues->val[i]) bit1(auxliste, disparues->val[i]);
	disparues = disparues->next;
	}
val = 0; loop = -1; deblng = 3;
totlng = read_first_rec(klng,NULL);
while((val = irbit(auxliste,val,nseq)) != 0) {
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


static void updatekshrt(int totacc, int totbib, LVALS **biblistes, 
	LVALS **acclistes, int *auxliste)
{
int numrec;
for(numrec = 2; numrec <= totbib; numrec++) {
	if(biblistes[numrec] != NULL)
		process_shrt(kbib, numrec, biblistes[numrec], auxliste);
	}
for(numrec = 2; numrec <= totacc; numrec++) {
	if(acclistes[numrec] != NULL)
		process_shrt(kacc, numrec, acclistes[numrec], auxliste);
	}
}


static void lvals_to_blist(LVALS *plist, int *blist, int *p_mini, int *p_maxi)
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


static int proc_lng_list(int debut, LVALS *plist, int *auxlist)
{
static struct rlng aux;
int mini, maxi, val, newdebut, lastaux, vide, num;

lvals_to_blist(plist, auxlist, &mini, &maxi);
lastaux = -1; newdebut = debut; vide = TRUE;
memset(&aux, 0, SUBINLNG * sizeof(int) );
while(debut != 0) {
	readlng(debut);
	if(newdebut == debut) aux.next = plng->next;
	for(num = 0; num < SUBINLNG; num++) {
		if( (val = plng->sub[num]) == 0) continue;
		if( val < mini || val > maxi || !testbit(auxlist, val) ) {
			if(lastaux >= SUBINLNG - 1) {
				if( dir_write(klng, newdebut, 1, &aux) ) 
					dir_writeerr(klng, newdebut);
				newdebut = aux.next;
				if(newdebut == debut)
					aux.next = plng->next;
				else if(dir_read(klng,newdebut,1,&aux) != 1)
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
aux.next = 0;
if( dir_write(klng, newdebut, 1, &aux) ) dir_writeerr(klng, newdebut);
return vide;
}


static int proc_shrt_list(int debut, LVALS *plist, int *auxlist)
{
int mini, maxi, val, next, newdebut, preval;

lvals_to_blist(plist, auxlist, &mini, &maxi);
readshrt(debut);
val = pshrt->val;
if(val >= mini && val <= maxi && testbit(auxlist, val) ) {
	do	{
		newdebut = pshrt->next;
		if(newdebut == 0) return TRUE;
		readshrt(newdebut);
		val = pshrt->val;
		}
	while(val >= mini && val <= maxi && testbit(auxlist, val) );
	writeshrt(debut);
	}
while(TRUE) {
	readshrt(debut);
	preval = pshrt->val;
	next = pshrt->next;
	if(next == 0) break;
	readshrt(next);
	val = pshrt->val;
	if(val >= mini && val <= maxi && testbit(auxlist, val) ) {
		pshrt->val = preval;
		writeshrt(debut);
		}
	else	{
		debut = next;
		}
	}
return FALSE;
}


#define addval_and_test(a, b) if( (a = addval(a, b)) == NULL ) return 1

static int delsubseq(int isub, int mere, 
	LVALS **smjlistes, LVALS **keylistes, LVALS **sublistes, 
	LVALS **p_disparues)
/* pour traiter la partie subseq de l'effacement d'une sequence
isub=num de la seq
mere=num de la mere qui est aussi effacee si on efface une fille a cause de 
     sa mere, =0 si on efface une fille ou une mere
*/
{
struct rsub rsub2, *psub2= &rsub2;
int valeur, point;

dir_read(ksub,isub,1,psub2);
addval_and_test(*p_disparues, isub);
suphsh(isub, ksub);
memset(psub, 0, lrsub);
memset(psub->name, 'x', L_MNEMO);
writesub(isub);
if(psub2->type) {
	addval_and_test(smjlistes[psub2->type], isub);
	}
point=psub2->plkey;
while(point) {
	readshrt(point);
	valeur = pshrt->val;
	point = pshrt->next;
	addval_and_test(keylistes[valeur], isub);
	}
point=psub2->pext;
if(point>0) {
	int autremere=1, premere;
	char name2[L_MNEMO + 1], *pos, *fin = name2 + L_MNEMO - 1;
	memcpy(name2, psub2->name, L_MNEMO);
	name2[L_MNEMO] = 0;
	pos = strchr(name2,'.');
	if(pos != NULL) {
		while(pos<=fin) *(pos++)=' ';
		}
	premere=0;
	while(point) {
		readext(point);
		point=pext->next;
		if(pext->mere==premere)continue;
		premere= pext->mere;
		if(premere != mere) {
			addval_and_test(sublistes[premere], isub);
			}
		if(autremere) {
			readsub(premere);
			autremere = strncmp(psub->name,name2,16);
			}
		}
	if(autremere) {
		autremere = isenum(name2);
		if(autremere != mere) {
			addval_and_test(sublistes[autremere], isub);
			}
		}
	}	
return 0;
} /* end of delsubseq */


static int fast_remove(int nsub, void *dynlist_mere_with_request,
	LVALS **smjlistes, LVALS **speclistes, LVALS **hostlistes, 
	LVALS **keylistes, LVALS **sublistes, 
	LVALS **acclistes, LVALS **biblistes,
	LVALS **p_meres, LVALS **p_disparues)
{
struct rloc rloc2;
struct rlng rlng2;
int locus, list_fi, i, point, fille, l_mere;
char mere_dot[L_MNEMO + 2];

if(nsub == 0) return 0;
readsub(nsub);
memcpy(mere_dot, psub->name, L_MNEMO); mere_dot[L_MNEMO] = 0;
trim_key(mere_dot); l_mere = strlen(mere_dot);
strcpy(mere_dot + l_mere, "."); l_mere++;
locus = psub->plinf;
list_fi = - psub->pext;
readloc(locus);
if(ploc->molec) {
	addval_and_test(smjlistes[ploc->molec], nsub);
	}
if(ploc->stat) {
	addval_and_test(smjlistes[ploc->stat], nsub);
	}
if(ploc->org) {
	addval_and_test(smjlistes[ploc->org], nsub);
	}
addval_and_test(*p_meres, nsub);
point = ploc->placc;
while(point) {
	readshrt(point);
	point = pshrt->next;
	addval_and_test(acclistes[pshrt->val], nsub);
	}
if( ( nbrf || swissprot ) && ploc->spec != 0) {
	point = ploc->spec;
	while(point) {
		readshrt(point);
		point=pshrt->next;
		addval_and_test(speclistes[pshrt->val], nsub);
		}
	}
else if (ploc->spec != 0) {
	addval_and_test(speclistes[ploc->spec], nsub);
	}
if(ploc->host) {
	addval_and_test(hostlistes[ploc->host], nsub);
	}
point = ploc->plref;
while(point) {
	int nbib;
	readshrt(point);
	point = pshrt->next;
	nbib = pshrt->val;
	addval_and_test(biblistes[nbib], nsub);
	readbib(nbib);
	if(pbib->j) {
		addval_and_test(smjlistes[pbib->j], nsub);
		}
	if(pbib->y) {
		addval_and_test(smjlistes[pbib->y], nsub);
		}
	}
if(ploc->bef) {
	dir_read(kloc,ploc->bef,1,&rloc2);
	rloc2.next = 0;
	dir_write(kloc,ploc->bef,1,&rloc2);
	}
if(ploc->next) {
	dir_read(kloc,ploc->next,1,&rloc2);
	rloc2.bef = 0;
	dir_write(kloc,ploc->next,1,&rloc2);
	}
memset(&rloc2, 0, lrloc);
memset(rloc2.date, 'x', sizeof(rloc2.date));
dir_write(kloc, locus, 1, &rloc2);

while(list_fi) {
	dir_read(klng, list_fi, 1, &rlng2);
	list_fi = rlng2.next;
	for(i = 0; i < SUBINLNG; i++ ) {
		fille = rlng2.sub[i];
		if(fille == 0) continue;
		readsub(fille);
		if(psub->length == 0) continue; /* seq deja effacee avant */
		if(strncmp(mere_dot, psub->name, l_mere) != 0) {
			add_fille_request(fille, dynlist_mere_with_request);
			}
		if( delsubseq(fille, nsub, 
			smjlistes, keylistes, sublistes, p_disparues)) return 1;
		}
	}
if( delsubseq(nsub, 0, smjlistes, keylistes, sublistes, p_disparues) ) return 1;
point = read_first_rec(ksub, NULL);
write_first_rec(ksub, point, 0);
return 0; 
} /* end of fast_remove */


void sup_fille_request_from_mere(char *nom_mere, 
	void *dynlist_mere_with_request)
{
fille_request *request, *next;
void *node;

if( find_g_dynlist(nom_mere, dynlist_mere_with_request, FALSE, &node)
	== NULL ) return;
request = dynlist_get_extra(node);
while(request != NULL) {
printf("Enleve %s pannots=%ld div=%d\n", request->fille, request->pannots, request->div);
	next = request->next;
	free(request->fille);
	free(request);
	request = next;
	}
remove_dynlist(nom_mere, dynlist_mere_with_request);
}


int add_fille_request(int fille, void *dynlist_mere_with_request)
/* returns TRUE if memory error, FALSE otherwise */
{
fille_request *p_request, *elt;
char *p, *q;
void *node;

p_request = malloc(sizeof(fille_request));
if(p_request != NULL) p_request->fille = malloc(L_MNEMO + 1);
if(p_request == NULL || p_request->fille == NULL) return TRUE;
readsub(fille);
memcpy(p_request->fille, psub->name, L_MNEMO); p_request->fille[L_MNEMO] = 0;
trim_key(p_request->fille);
seq_to_annots(fille, &(p_request->pannots), &(p_request->div));
p_request->type = psub->type;
p_request->next = NULL;
/* ajouter mere principale de cette fille a dynlist_mere_with_request */
p = strchr(p_request->fille, '.'); *p = 0;
q = find_g_dynlist(p_request->fille, dynlist_mere_with_request, TRUE, &node);
*p = '.';
if( q == NULL) return TRUE;
elt = dynlist_get_extra(node);
if(elt == NULL)
	dynlist_set_extra(node, p_request);
else	{
	while(elt->next != NULL) elt = elt->next;
	elt->next = p_request;
	}
printf("Ajoute fille=%s annots=%ld div=%d type=%d\n",
	p_request->fille, p_request->pannots,
	p_request->div, p_request->type);
return FALSE;
}


void process_fille_request(fille_request *p_request, char *location, 
	char *qualifiers, int max_loc_qual, int *p_fille_num, 
	int *p_tot_filles_cre, FILE *log_file)
{
char *p, feat_name[20], nom_fille[L_MNEMO + 1];
int num_mere, erreur, created, point, i, trouve;
fille_request debut;

debut.next = p_request;
p_request = &debut;
while(p_request->next != NULL) {
	p_request = p_request->next;
	strcpy(nom_fille, p_request->fille);
	p = strchr(nom_fille, '.'); *p = 0;
	num_mere = isenum(nom_fille);
	if(num_mere == 0) continue;
	/* cette adresse est-elle deja traitee par une fille ? */
	readsub(num_mere);
	point = abs(psub->pext); trouve = FALSE;
	while(!trouve && point != 0) {
		readlng(point); point = plng->next;
		for(i = 0; i < SUBINLNG; i++) {
			if(plng->sub[i] == 0) break;
			readsub(plng->sub[i]);
			readshrt(psub->plinf);
			if(pshrt->val == p_request->pannots) {
				trouve = TRUE;
				break;
				}
			}
		}
	if(trouve) continue;
	/* attempt creating new sequence */
	erreur = read_addr_loc_qualif(p_request->pannots, 
		p_request->div, location, 
		max_loc_qual, feat_name, qualifiers, max_loc_qual);
	if(erreur) continue;
	strcpy(nom_fille, p_request->fille);
	created = mkextr(location, num_mere, p_fille_num, nom_fille, 
		p_request->type, 
		feat_name, p_request->pannots, p_request->div, qualifiers);
	if(created) {
		write_first_rec(ksub, *p_fille_num, 0);
		++(*p_tot_filles_cre);
		readsub(*p_fille_num); 
		memcpy(nom_fille, psub->name, L_MNEMO); nom_fille[L_MNEMO] = 0;
		fprintf(log_file, "Created %s\n", nom_fille);
		}
	}
}


void process_genetic_code(int mere)
/* met le code genetique de toutes les filles CDS de mere
n'est utile que pour les CDS qui n'ont pas /transl_table
et que si le libelle de l'espece indique le code de l'espece
*/
{
char libel[lrtxt + 1], *p, *q, *debut;
int mtgc, gc, is_mt, point, i, fille, current_gc, phase;
static int cds_type = 0;

if( !(genbank || embl) ) return;
if(cds_type == 0) cds_type = fcode(ksmj, "04CDS", 20);
readsub(mere);
readloc(psub->plinf);
readspec(ploc->spec);
if(pspec->libel == 0) return;
readtxt(pspec->libel);
memcpy(libel, ptxt, lrtxt); libel[lrtxt] = 0;
debut = libel;
gc = mtgc = 0; /* standard code */
while( (p = strchr(debut, '|')) != NULL) {
	*p = 0;
	if( ( q = strstr(debut, "mtgc:")) != NULL) 
		sscanf(q + 5, "%d", &mtgc);
	else if( ( q = strstr(debut, "gc:")) != NULL) 
		sscanf(q + 3, "%d", &gc);
	debut = p + 1;
	}
if(gc == 0 && mtgc == 0) return;
is_mt = FALSE;
if(ploc->org != 0) {
	readsmj(ploc->org);
	if(strncmp(psmj->name, "05MITOCHONDRION ", 16) == 0 ||
	   strncmp(psmj->name, "05KINETOPLAST ", 14) == 0) is_mt = TRUE;
	}
current_gc = (is_mt ? mtgc : gc);
if(current_gc == 0) return;
if(psub->type == cds_type) {
	readsub(mere);
	phase = (psub->phase % 100);
	psub->phase = 100 * current_gc + phase;
	writesub(mere);
	}
point = - psub->pext;
while(point > 0) {
	readlng(point); point = plng->next;
	for(i = 0; i < SUBINLNG; i++) {
		if( (fille = plng->sub[i]) == 0) break;
		readsub(fille);
		if(psub->type != cds_type) continue;
		phase = (psub->phase % 100);
		psub->phase = 100 * current_gc + phase;
		writesub(fille);
		}
	}
}
