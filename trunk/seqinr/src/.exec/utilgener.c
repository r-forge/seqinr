#include "dir_acnuc.h"
#include <ctype.h>

#define MEM_ONE_REC_TYPE 10000
#define MAX_LOC_QUAL 100000

/* structures - typedefs */
	/* gestion de memorisation des requetes d'accession */
struct _acc_request {
	long inf;
	int div;
	int mere;
	int type;
	char acc[ACC_LENGTH + 1];
	char nom[L_MNEMO + 1];
	};

	/* memorisation des fragments d'une fille */
#define MAX_FRAGS 500
struct _fille_frags { 
	int nfmax;  /* maximum # of fragments in a subseq */
	int nfrags; /* # of fragments in subseq (from 1) */
	struct _one_frag {
		int debut; /* beginning of fragment */
		int fin; /* end of fragment */
		int pnuc; /* pointer to nucleot for that fragment */
		int mere; /* # of parent seq of that fragment	*/
		} one_frag[MAX_FRAGS];
	};
typedef void (*type_extra_function)(char *, int);


/* included functions */
int mkextr(char *location, int seq_num, int *p_fille_num, char *name,
	int smj_type, char *feat_name, long info_addr, int div, char *qualif);
char *make_unique_seq_name(char *oldname);
char *add_gene_extension(int isub, char *extension);
static int proc_locat_fragment(char *debut, char *fin,
	int mere, int *part5, int *part3, int *fromlocus, char *access);
int find_label(char *label, char *access, long *pannot, int *div, int *isub);
char *next_word(char *ligne, char *word, char *separs, int max_word, 
	int *erreur);
void add_liste_meres(int seqnum);
void close_if_open(void);
#ifdef unix
void proc_mmap_options(int argc, char **argv);
#endif
int *get_smj_extra(int smjnum);
int *get_keyw_extra(int kwnum);
int *get_spec_extra(int specnum);
int *get_bib_extra(int bibnum);
void fast_add_seq_to_spec(int specnum, int seqnum);
void fast_add_seq_to_keyw(int kwnum, int seqnum);
static void add_to_spec_kw_lng( DIR_FILE *fichier, int num, int seqnum, 
	int *extra);
void add_to_smjlng( int smjnum, int seqnum, int *extra);
char *nextpar(char *debut);
void connect_locus_bib( int locnum, int seqnum, int bibnum, int *extra);
int process_year_bib(char *annee, int **extra);
int next_start_seq(FILE *addr_file, char *gcgacnuc_value, char *ligne, 
	size_t lligne, off_t *p_addr, int *p_div);
char *proc_classif(char *ligne, size_t lligne);
int extend_acc_requests(void);


/* external prototypes */
char *get_qualif_value(char *p);
int get_acnuc_gc_number(int ncbi_gc);
int find_key(char *name, void *tree, int create_if_not_found,
	int **p_aux_data);
int poffset(int div, int offset);
int read_addr_loc_qualif(long faddr, int div, 
	char *location, int maxlocat, char *type, char *qualif, int maxqualif);
int find_cre_acc(char *name, int create);


/* global variables */
struct _acc_request *acc_request = NULL;
int tot_request =  0;
struct _fille_frags fille_frags;

/* extern variables */
extern FILE *in_flat, *log_file;
extern void *smj_binary_tree, *bib_binary_tree;
extern int non_chargee;
extern char *gcgname[];
extern type_extra_function extra_function;


int mkextr(char *location, int seq_num, int *p_fille_num, char *name,
	int smj_type, char *feat_name, long info_addr, int div, char *qualif)
/* fabrication d'une fille a partir d'une location
seq_num: numero de la mere
*p_fille_num: rendu avec num de la fille fabriquee (inchange si pas de fille faite)
name: nom de la fille a creer
smj_type: numero du type de la fille
feat_name: feature key
info_addr + div: file_addr + division de la feature entry traitee
qualif: les qualifiers deja lus
rendu TRUE si fille creee, FALSE sinon en particulier si fille == mere
*/
{
static char err_text[9][40] = {"OK", "unused",
 	"Syntax error", "Reference to other database",
	"Label not found", "Primary accession # not found",
	"Increase parameter MAX_FRAGS",
	"Endpoints go beyond parent sequence",
	"Increase parameter MAX_LOC_QUAL" 
	};
char *p, *debut, *fin, *q, *virgule, access[20],
	*qualif_value;
int lpart5, lpart3, fromlocus, compl, erreur, in_join,
	num, totext, fille_num, prev_mere, longueur, *extra,
	done, numkey, created, info_fille;
static int first = TRUE;
static int npart5, npart3;

if(first) {
	first = FALSE;
	npart5 = crekeyword("PARTIAL", "5'-PARTIAL");
	npart3 = crekeyword("PARTIAL", "3'-PARTIAL");
	}
created = FALSE;
/* if location contains order( or group( or one-of( or replace( or ^
 or "sequence" do not process further
*/
if(strstr(location,"ORDER(") != NULL ||
	strstr(location,"GROUP(") != NULL ||
	strstr(location,"ONE-OF(") != NULL ||
	strstr(location,"ONE_OF(") != NULL ||
	strstr(location,"REPLACE(") != NULL ||
	strchr(location,'^') != NULL ||
	strchr(location,'"') != NULL ) return FALSE;
/* if location contains approximate endpoint, i.e., a single dot,
 do not process further
*/
p = location;
while ( ( p = strchr(p, '.') ) != NULL) {
	p++;
	if (*p == 0) break;
	if( *p != '.' ) {/* mais accepter .1234567890: */
		while(isdigit(*p)) p++;
		if(*p != ':') return FALSE;
		}
	p++;
	}
lpart5 = lpart3 = fromlocus = compl = in_join = FALSE;
debut = location; fin = location + strlen(location) - 1;
fille_frags.nfrags = 0;
fille_frags.nfmax = MAX_FRAGS;
if(strncmp(location, "COMPLEMENT(JOIN(", 16) == 0) {
	compl = TRUE;
	debut += 11;
	}
if(strncmp(debut, "JOIN(", 5) == 0) {
	in_join = TRUE;
	fin = nextpar(debut + 4);
	debut += 5;
/* remove redundant JOIN() within a previous JOIN() */
	while( (p = strstr(debut, "JOIN(")) != NULL &&
		p < fin) {
		q = nextpar(p + 4);
		memset(p, ' ', 5);
		*q = ' ';
		}
	}
else	{
	virgule = fin + 1;
	}
do	{
	if(in_join) {
		virgule = debut;
		do	{
			virgule++;
			if(*virgule == '(') virgule = nextpar(virgule) + 1;
			}
		while(*virgule != ',' && *virgule != ')' && *virgule != 0);
		}
	erreur = proc_locat_fragment(debut, virgule - 1, seq_num, 
		&lpart5, &lpart3, &fromlocus, access);
	if(erreur) break;
	if(in_join) {
		debut = virgule + 1;
		}
	}
while (*virgule != 0 && *virgule != ')');
if(erreur != 0) {
    if(erreur == 5) { /* memorize an accession # request */
		for(num = 0; num < tot_request; num++) {
			if( strcmp(acc_request[num].acc, access) == 0 &&
				acc_request[num].inf == info_addr && 
				acc_request[num].div == div) break;
			}
		if(num >= tot_request) {
			if(extend_acc_requests()) {
				fprintf(log_file, 
					"Warning: not enough memory for new accession request\n");
				}
			else	{
				strcpy(acc_request[tot_request].acc, access);
				acc_request[tot_request].inf = info_addr;
				acc_request[tot_request].div = div;
				acc_request[tot_request].mere = seq_num;
				strcpy(acc_request[tot_request].nom, name);
				acc_request[tot_request].type = smj_type;
				tot_request++;
				}
			}
    	}
    else {
    	fprintf(log_file, "Warning: cannot process feature: %s %s\n%s\n",
    		feat_name, err_text[erreur], location);
    	}
    return FALSE;
    }

/* check if subseq clears an accession # request */
do	{
	done = FALSE;
	for(num = 0; num < tot_request; num++) {
		if(acc_request[num].inf == info_addr && 
					acc_request[num].div == div) {
			int i;
			for(i = num + 1; i < tot_request; i++) 
				acc_request[i - 1] = acc_request[i];
			tot_request--;
			done = TRUE;
			}
		}
	}
while(done);
/* check if subseq is identical to mother seq */
if(fille_frags.nfrags == 1 && fille_frags.one_frag[0].debut == 1) {
	readsub(seq_num);
	if(fille_frags.one_frag[0].fin == psub->length) {
/* subseq same as locus seq */
		fille_num = seq_num;
		goto suite;
		}
	}
/* do not create subseqs shorter than 3 nucleotides */
if(fille_frags.nfrags == 1 && 
	abs(fille_frags.one_frag[0].fin - fille_frags.one_frag[0].debut) <= 3)
		return FALSE;
/* ici une vraie sous-sequence */
(*p_fille_num)++;
fille_num = *p_fille_num;
created = TRUE;
/* calcul de la longueur */
prev_mere = 0;
longueur = 0;
for(num=0; num < fille_frags.nfrags; num++) {
	if(fille_frags.one_frag[num].mere != prev_mere) {
/* ajouter a liste des filles de ses meres */
		prev_mere = fille_frags.one_frag[num].mere;
		mdlng(ksub, prev_mere, 3, fille_num, NULL);
		}
	longueur += abs(fille_frags.one_frag[num].fin - 
			fille_frags.one_frag[num].debut) + 1;
	}
/* si ne contient pas de morceau de seq # seq_num, mettre en plus en fille de seq_num */
if(!fromlocus) mdlng(ksub, seq_num, 3, fille_num, NULL);
if(compl) {
	int num2, tmp;
/* traitement du complement( externe */
	if(lpart5 != lpart3) {
		lpart5 = !lpart5;
		lpart3 = !lpart3;
		}
	/* inverser l'ordre des exons et de deb et fin */
	for(num = 0; num < fille_frags.nfrags / 2; num++) {
		num2 = fille_frags.nfrags - num - 1;
		tmp = fille_frags.one_frag[num].debut;
		fille_frags.one_frag[num].debut = 
			fille_frags.one_frag[num2].fin;
		fille_frags.one_frag[num2].fin = tmp;
		tmp = fille_frags.one_frag[num].fin;
		fille_frags.one_frag[num].fin = 
			fille_frags.one_frag[num2].debut;
		fille_frags.one_frag[num2].debut = tmp;
		tmp = fille_frags.one_frag[num].pnuc;
		fille_frags.one_frag[num].pnuc = 
			fille_frags.one_frag[num2].pnuc;
		fille_frags.one_frag[num2].pnuc = tmp;
		tmp = fille_frags.one_frag[num].mere;
		fille_frags.one_frag[num].mere = 
			fille_frags.one_frag[num2].mere;
		fille_frags.one_frag[num2].mere = tmp;
		}
	if( fille_frags.nfrags % 2 == 1 ) { 
		/* inverser deb et fin pour exon central */
		num = fille_frags.nfrags / 2;
		tmp = fille_frags.one_frag[num].debut;
		fille_frags.one_frag[num].debut = fille_frags.one_frag[num].fin;
		fille_frags.one_frag[num].fin = tmp;		
		}
	}
/* verifier que name n'est pas deja utilise (peut arriver dans cas particuliers)
sinon en creer un nouveau */
if(isenum(name) != 0) p = make_unique_seq_name(name);
else	p = name;
memset(psub, 0, lrsub);
padtosize(psub->name, p, L_MNEMO);
psub->length = longueur;
totext = read_first_rec(kext, NULL);
psub->pext = totext + 1;
writesub(fille_num);
if( !compl ) 
	compl = (fille_frags.one_frag[0].debut > fille_frags.one_frag[0].fin);
/* traitement liste infos */
if( ! big_annots ) info_fille = poffset(div, (int) info_addr);
else	info_fille = (int) info_addr;
mdshrt(ksub, fille_num, 5, info_fille, NULL);
/* ecriture fichier EXTRACT */
for(num = 0; num < fille_frags.nfrags; num++) {
	totext++;
	pext->deb = fille_frags.one_frag[num].debut;
	pext->fin = fille_frags.one_frag[num].fin;
	pext->pnuc = fille_frags.one_frag[num].pnuc;
	pext->mere = fille_frags.one_frag[num].mere;
	pext->next = ( num == fille_frags.nfrags - 1 ? 0 : totext + 1 );
	writeext(totext);
	}
write_first_rec(kext, totext, 0);
/* connecter au hashing des mnemos */
addhsh(fille_num, ksub);

suite:
/* connect seq <--> type */
readsub(fille_num);
if(psub->type == 0) { 
	/* cas rare ou plusieurs fois fille = mere avec types differents! */
	psub->type = smj_type;
	writesub(fille_num);
	extra = get_smj_extra(smj_type);
	add_to_smjlng(smj_type, fille_num, extra);
	}
/* associer avec mots-cles partiels */
if(lpart5) {
	erreur = mdshrt(ksub, fille_num, 4, npart5, NULL);
	if(erreur != 2) fast_add_seq_to_keyw(npart5, fille_num);
	}
if(lpart3) {
	erreur = mdshrt(ksub, fille_num, 4, npart3, NULL);
	if(erreur != 2) fast_add_seq_to_keyw(npart3, fille_num);
	}
/* traitement des qualifiers */
/* qualifier /EC_NUMBER */
p = strstr(qualif, "/EC_NUMBER=");
if( p != NULL && (qualif_value = get_qualif_value(p) ) != NULL ) {
	numkey = crekeyword("EC_NUMBERS", qualif_value);
	erreur = mdshrt(ksub, fille_num, 4, numkey, NULL);
	if(erreur != 2) fast_add_seq_to_keyw(numkey, fille_num);
	}
/* qualifier /EVIDENCE */
p = strstr(qualif, "/EVIDENCE=");
if( p != NULL && (qualif_value = get_qualif_value(p) ) != NULL ) {
	numkey = crekeyword(NULL, qualif_value);
	erreur = mdshrt(ksub, fille_num, 4, numkey, NULL);
	if(erreur != 2) fast_add_seq_to_keyw(numkey, fille_num);
	}
/* qualifier /PARTIAL */
p = strstr(qualif, "/PARTIAL");
if( p != NULL && !( lpart5 || lpart3 ) ) {
	numkey = crekeyword(NULL, "PARTIAL");
	erreur = mdshrt(ksub, fille_num, 4, numkey, NULL);
	if(erreur != 2) fast_add_seq_to_keyw(numkey, fille_num);
	}
/* qualifier /PSEUDO */
p = strstr(qualif, "/PSEUDO");
if( p != NULL ) {
	numkey = crekeyword(NULL, "PSEUDO");
	erreur = mdshrt(ksub, fille_num, 4, numkey, NULL);
	if(erreur != 2) fast_add_seq_to_keyw(numkey, fille_num);
	}
/* qualifier /PRODUCT */
p = strstr(qualif, "/PRODUCT=");
if( p != NULL && (qualif_value = get_qualif_value(p) ) != NULL ) {
	numkey = crekeyword("PRODUCTS", qualif_value);
	erreur = mdshrt(ksub, fille_num, 4, numkey, NULL);
	if(erreur != 2) fast_add_seq_to_keyw(numkey, fille_num);
	}
/* qualifier  /codon_start= */
	/* traitement de AA x at y pas fait. Est-il encore utile ?????
	*/
p = strstr(qualif, "/CODON_START=");
if( p != NULL) {
	int phase;
	qualif_value = get_qualif_value(p);
	phase = 0;
	if(qualif_value != NULL) sscanf(qualif_value, "%d", &phase);
/* should say only 1,2,or 3 as a frame identification */
	if(phase >= 1 && phase <= 3) {
		phase--;
		readsub(fille_num);
		psub->phase = phase;
		writesub(fille_num);
		}
	}
/* qualifier  /gene= or /standard_name= or (/nomgen= for EMBL) */
p = strstr(qualif, "/GENE=");
if( p == NULL) p = strstr(qualif, "/STANDARD_NAME=");
if( p == NULL) p = strstr(qualif, "/NOMGEN=");
if( p != NULL && (qualif_value = get_qualif_value(p) ) != NULL ) {
	numkey = crekeyword("GENETIC NAMES", qualif_value);
	erreur = mdshrt(ksub, fille_num, 4, numkey, NULL);
	if(erreur != 2) fast_add_seq_to_keyw(numkey, fille_num);
	add_gene_extension(fille_num, qualif_value);
	}
/* qualifier /TRANSL_TABLE */
p = strstr(qualif, "/TRANSL_TABLE=");
if( p != NULL ) {
	int gc = 1;
	sscanf(p + 14, "%d", &gc);
	gc = get_acnuc_gc_number(gc);
	readsub(fille_num);
	psub->phase = 100 * gc + (psub->phase % 100);
	writesub(fille_num);
	}

/* permettre traitement optionnel d'autres qualifiers */
if(extra_function != NULL) extra_function(qualif, fille_num);

return created;
} /* end of mkextr */


char *make_unique_seq_name(char *oldname)
{
static char newname[L_MNEMO + 10];
char *p;
int rank;

strcpy(newname, oldname);
p = strchr(newname, '.');
if(p == NULL) return oldname;
do p++;
while( *p != 0 && !isdigit(*p) );
if( *p == 0 ) return oldname;
sscanf(p, "%d", &rank);
do	{
	rank++;
	sprintf(p, "%d", rank);
	}
while(isenum(newname) != 0);
return newname;
}


char *add_gene_extension(int isub, char *extension)
/* changer l'extension du nom d'une fille avec extension
isub: numero de la fille
rend nouveau nom si nom change', NULL sinon
*/
{
char *p;
int l;
static char newname[ 2 * L_MNEMO];

readsub(isub);
if(psub->length == 0 || psub->pext <= 0) return NULL;
p = strchr(psub->name, '.');
if(p == NULL || p >= psub->name + L_MNEMO ) return NULL;
l = p - psub->name + 1;
if( l + strlen(extension) >= sizeof(newname) ) return NULL;
memcpy(newname, psub->name, l);
strcpy(newname + l, extension);
/* replace internal spaces by '-' */
while ( ( p = strchr(newname, ' ') ) != NULL ) *p = '-';
/* remove @ characters */
while ( ( p = strchr(newname, '@') ) != NULL && p >= newname + l ) *p = ' ';
/* remove [ characters */
while ( ( p = strchr(newname, '[') ) != NULL && p >= newname + l ) *p = ' ';
/* remove ] characters */
while ( ( p = strchr(newname, ']') ) != NULL && p >= newname + l ) *p = ' ';
compact(newname); 
if( (int) strlen(newname) > L_MNEMO) return NULL;
if(isenum(newname) != 0) return NULL;
padtosize(newname, newname, L_MNEMO);
suphsh(isub, ksub);
readsub(isub);
memcpy(psub->name, newname, L_MNEMO);
writesub(isub);
addhsh(isub, ksub);
return newname;
}


static int proc_locat_fragment(char *debut, char *fin,
	int mere, int *part5, int *part3, int *fromlocus, char *access)
/* traitement d'un exon d'une location, avec cas envoi avec access # ou avec
label.
debut, fin: location traitee, extremites de sa partie traitee,
mere: # mere traitee,
part5, part3: mis a TRUE si partiel
*fromlocus: mis a TRUE si contient un morceau de la mere traitee
access: rendu avec accession # cause de la non creation
valeur retournee: 0 si OK, sinon code d'erreur
*/
{
int compl, point, next, nacc, a, b, err, expr_a_b, sv_mere, div;
long pinf;
char *p, label[16];
static char location[MAX_LOC_QUAL];

sv_mere=mere;
boulot:
while( *debut==' ' && debut<=fin) debut++;
if(debut>fin) return 2;
if(strncmp(debut,"COMPLEMENT(",11)==0) {
	debut+=11;
	p=nextpar(debut-1)-1; if(p>fin) return 2;
	fin = p;
	compl=TRUE;
	}
else compl=FALSE;
/* recherche renvoi a autre banque */
if((p=strstr(debut,"::"))!=NULL && p<fin) return 3;
p=strchr(debut,':');
if(p!=NULL && p<fin) {
/* renvoi a autre seq par access # */
	memcpy(access,debut,p-debut);
	access[p-debut]=0;
	if( (debut = strchr(access,'.')) != NULL) *debut = 0;
	debut = p + 1;
	nacc = find_cre_acc(access, FALSE);
/* chercher une seq ayant cet acc # en primary ou seq unique */
	if(nacc == 0) return 5;
	readacc(nacc);
	point = pacc->plsub;
	if(point == 0) return 5;
	readshrt(point);  
	if(pshrt->next == 0) { /* une seule seq associee a access */
		mere = pshrt->val;
		}
	else	{ /* plusieurs seqs associees, chercher celle dont primary */
		while(TRUE) {
			if(point == 0) return 5;
			readshrt(point); next = pshrt->next; 
			mere = pshrt->val;
			readsub(mere); readloc(psub->plinf);
			readshrt(ploc->placc);
			if(pshrt->val == nacc) break;
			point = next;
			}
		}
	}
else	{
/* on reste dans la seq de la table de features */
	*access = 0;
	}
/* traitement partie location */
p = strstr(debut, ".."); if(p != NULL && p >= fin) p = NULL;
if( p!= NULL) 
	expr_a_b = TRUE;
else	{
	memcpy(label,debut,fin-debut+1); label[fin-debut+1]=0;
	if((p=strchr(label,'<'))!=NULL) *p=' ';
	if((p=strchr(label,'>'))!=NULL) *p=' ';
	err=sscanf(label,"%d",&a);
	if(err==1) { /* residu unique */
		*part5= *part3 = FALSE; /* ??? */
		b=a; expr_a_b=FALSE;
		}
	else	{ /* c'est un label */
		err=find_label(label,(*access==0 ? NULL : access), &pinf, &div,
				 &mere);
		if(err==1) return 5;
		else if(err==2) return 4;
		err = read_addr_loc_qualif(pinf, div, location, MAX_LOC_QUAL, 
				NULL, NULL, 0);
		if(err) return 8;
		debut = location;
		fin = debut + strlen(debut) - 1;
		goto boulot;
		}
	}
if(expr_a_b) {
	char *q;
	if(*debut=='<' || *debut=='>') {
		if(compl) *part3=TRUE;
		else *part5=TRUE;
		debut++;
		}
	/* traiter cas (pos1)..pos2 */
	if( (*debut == '(') &&
			(q = strchr(debut, ')') ) != NULL && q < p)
		q = debut + 1;
	else
		q = debut;
	err = sscanf(q,"%d",&a);
	if(err == 0) return 2;
	p += 2;
	if(*p=='<' || *p=='>') {
		if(compl) *part5=TRUE;
		else *part3=TRUE;
		p++;
		}
	/* traiter cas pos1..(pos2) */
	if( (*p == '(') &&
			(q = strchr(p, ')') ) != NULL && q <= fin)
		q = p + 1;
	else
		q = p;
	err = sscanf(q,"%d",&b);
	if(err == 0) return 2;
	}
if(fille_frags.nfrags >= fille_frags.nfmax) return 2;
++fille_frags.nfrags;
readsub(mere);
readloc(psub->plinf);
/* check fragment endpoints relatively to size of mother seq */
if(a<1 || b<1 || a>psub->length || b>psub->length) return 7;
if(a > b) { /* forcer que toujours a <= b */
	next=a; a=b; b=next;
	}
if(compl) {
	next=a; a=b; b=next;
	}
fille_frags.one_frag[fille_frags.nfrags - 1].debut = a;
fille_frags.one_frag[fille_frags.nfrags - 1].fin = b;
fille_frags.one_frag[fille_frags.nfrags - 1].pnuc = ploc->pnuc;
fille_frags.one_frag[fille_frags.nfrags - 1].mere = mere;
if(mere == sv_mere) *fromlocus = TRUE;
return 0;
}
/* end of proc_locat_fragment */


int find_label(char *label, char *access, long *pannot, int *div, int *isub)
/* recherche de la feature  access-label
si access==NULL, recherche locale ds la table de seq # *isub,
si trouve, *isub = # mere ou trouve et *pannot + *div = pointeur info debut features
        et rendu 0
si pas trouve, rendu 1 si accession manque  et 2 si label manque 
*/
{
char *p;
int nacc, trouve;
long point;

if(access==NULL)
	nacc=0;
else	{
	nacc = find_cre_acc(access, FALSE);
	if(nacc!=0) {
		readacc(nacc);
		nacc=pacc->plsub;
		}
	if(nacc==0) return 1;
	}
do	{
	if(nacc!=0) 	{
		readshrt(nacc);
		*isub=pshrt->val; nacc=pshrt->next;
		}
	seq_to_annots(*isub, pannot, div);
	read_annots(*pannot, *div);
	trouve=FALSE;
	do	{
		next_annots(pannot);
		if( (genbank && strncmp(pinfo->line,"BASE COUNT",10) == 0 ) || 
			(embl && strncmp(pinfo->line,"SQ ",3) == 0 ) ) break;
		trouve = (genbank && strncmp(pinfo->line,"FEATURES",8) == 0) || 
				(embl && strncmp(pinfo->line,"FH",2) == 0 );
		}
	while(!trouve);
	if(!trouve) break;
	point= *pannot;
	do	next_annots(&point);
	while( embl && strncmp(pinfo->line,"FH",2) == 0 );
	while( (genbank && strcmptrail(pinfo->line, 2, NULL, 0) == 0 ) || 
				(embl && strncmp(pinfo->line, "FT", 2) == 0) ) {
		*pannot = point;
		do	{
			majuscules(pinfo->line);
			p=strstr(pinfo->line,"/LABEL=");
			if(p != NULL) {
				if(*(p+7)==0) {
					next_annots(&point);
					majuscules(pinfo->line);
					p=pinfo->line+14;
					}
				if(strncmp(p+7,label,strlen(label)) == 0) 
					return 0;
				}
			next_annots(&point);
			}
		while ( (genbank && strcmptrail(pinfo->line,10,NULL,0) == 0 ) ||
			(embl && strcmptrail(pinfo->line,10,"FT",2) == 0) );
		}
	}
while(nacc!=0);
return 2;
}


char *next_word(char *ligne, char *word, char *separs, int max_word, 
	int *erreur)
/* trouve dans ligne le prochain mot selon separateurs de separ
et le met dans word
maxword = taille autorisee dans word
erreur rendu TRUE si word n'est pas assez grand, et rend word tronque
si pas trouve, fait word[0] = 0
rend pointeur vers debut prochain mot
rend NULL si pas de prochain mot
*/
{
static char ligne2[MEM_ONE_REC_TYPE];
char *p;
int l;

strcpy(ligne2, ligne);
p = strtok(ligne2, separs);
if( p != NULL) {
	l = strlen(p);
	if(l >= max_word) {
		if( erreur != NULL) *erreur = TRUE;
		l = max_word - 1;
		memcpy(word, p, l); word[l] = 0;
		} 
	else 	{
		if( erreur != NULL) *erreur = FALSE;
		memcpy(word, p, l + 1);
		}
	p = strtok(NULL, separs);
	}
else	word[0] = 0;
return (p == NULL ? NULL : p - ligne2 + ligne);
}


void add_liste_meres(int seqnum)
{
static int debut = 2;

readlng(debut);
while(plng->next != 0) {
	debut = plng->next;
	readlng(debut);
	}
addlng(debut, seqnum);
}


void close_if_open(void)
{
if(ksub != NULL) dir_acnucclose();
}


#ifdef unix
void proc_mmap_options(int argc, char **argv)
{
int num = 0;
while (++num < argc) {
	if(strcmp(argv[num], "-mmap") != 0) continue;
	num++;
	if(strcmp(argv[num], "kshrt") == 0) dir_set_mmap(kshrt);
	else if(strcmp(argv[num], "klng") == 0) dir_set_mmap(klng);
	else if(strcmp(argv[num], "ksub") == 0) dir_set_mmap(ksub);
	else if(strcmp(argv[num], "kloc") == 0) dir_set_mmap(kloc);
	else if(strcmp(argv[num], "kspec") == 0) dir_set_mmap(kspec);
	else if(strcmp(argv[num], "kkey") == 0) dir_set_mmap(kkey);
	else if(strcmp(argv[num], "kacc") == 0) dir_set_mmap(kacc);
	else if(strcmp(argv[num], "kaut") == 0) dir_set_mmap(kaut);
	else if(strcmp(argv[num], "kbib") == 0) dir_set_mmap(kbib);
	else if(strcmp(argv[num], "ksmj") == 0) dir_set_mmap(ksmj);
	}
}
#endif


int *get_smj_extra(int smjnum)
{
static char code[sizeof(psmj->name) + 1];
int *extra;
const int width = sizeof(psmj->name);

readsmj(smjnum); 
memcpy(code, psmj->name, width); code[width] = 0; trim_key(code);
smjnum = find_key(code, smj_binary_tree, FALSE, &extra);
if( smjnum == 0) extra = NULL;
return extra;
}


int *get_keyw_extra(int kwnum)
{
static int *extra = NULL;
static int maxkey = 0;
int totkey, *p_extra;

totkey = read_first_rec(kkey, NULL);
if( totkey > maxkey ) {
	int *tmp, new_maxkey;
	new_maxkey = totkey + 1000;
	tmp = calloc( new_maxkey , sizeof(int) );
	if(tmp != NULL) {
		memcpy(tmp, extra, maxkey * sizeof(int) );
		if(extra != NULL) free(extra);
		extra = tmp;
		maxkey = new_maxkey;
		}
	}
if(kwnum > maxkey)
	p_extra = NULL;
else
	p_extra = extra + kwnum - 1;
return p_extra;
}


int *get_spec_extra(int specnum)
{
static int *extra = NULL;
static int maxspec = 0;
int totspec, *p_extra;

totspec = read_first_rec(kspec, NULL);
if( totspec > maxspec ) {
	int *tmp, new_maxspec;
	new_maxspec = totspec + 200;
	tmp = calloc( new_maxspec , sizeof(int) );
	if(tmp != NULL) {
		memcpy(tmp, extra, maxspec * sizeof(int) );
		if(extra != NULL) free(extra);
		extra = tmp;
		maxspec = new_maxspec;
		}
	}
if(specnum > maxspec)
	p_extra = NULL;
else
	p_extra = extra + specnum - 1;
return p_extra;
}


int *get_bib_extra(int bibnum)
{
static char code[sizeof(pbib->name) + 1];
int *extra;
const int width = sizeof(pbib->name);

readbib(bibnum); 
memcpy(code, pbib->name, width); code[width] = 0; trim_key(code);
bibnum = find_key(code, bib_binary_tree, FALSE, &extra);
if( bibnum == 0) extra = NULL;
return extra;
}


void fast_add_seq_to_spec(int specnum, int seqnum)
{
int *p_extra;

p_extra = get_spec_extra(specnum);
add_to_spec_kw_lng(kspec, specnum, seqnum, p_extra);
}


void fast_add_seq_to_keyw(int kwnum, int seqnum)
{
int *p_extra;

p_extra = get_keyw_extra(kwnum);
add_to_spec_kw_lng(kkey, kwnum, seqnum, p_extra);
}


static void add_to_spec_kw_lng( DIR_FILE *fichier, int num, int seqnum, 
	int *extra)
/* ajouter seqnum a liste longue de sp/kw #num avec memorisation dans *extra 
de la fin courante de cette liste
*/
{
int debut;
static struct rspec buffer;

if(extra == NULL) { /* faire un ajout normal dans la liste longue */
	mdlng(fichier, num, 2, seqnum, NULL);
	return;
	}
if(*extra == 0) {
	if(dir_read(fichier, num, 1, &buffer) != 1) dir_readerr(fichier,num);
	if(buffer.plsub == 0) {
		debut = read_first_rec(klng, NULL) + 1;
		buffer.plsub = debut;
		addlng(debut, seqnum);
		if( dir_write(fichier, num, 1, &buffer) ) 
			dir_writeerr(fichier,num);
		*extra = debut;
		readshrt(buffer.desc); 
		/* mettre debut liste desc avec signe < 0 */
		if(pshrt->val > 0) {
			pshrt->val = - pshrt->val;
			writeshrt(buffer.desc);
			}
		return;
		}
	else debut = buffer.plsub;
	}
else debut = *extra;
readlng(debut);
while(plng->next != 0) {
	debut = plng->next;
	readlng(debut);
	}
*extra = debut;
addlng(debut, seqnum);
}


void add_to_smjlng( int smjnum, int seqnum, int *extra)
/* ajouter seqnum a liste longue de smjnum avec memorisation dans *extra 
de la fin courante de cette liste
*/
{
int debut;

if(extra == NULL) { /* faire un ajout normal dans la liste longue */
	mdlng(ksmj, smjnum, 1, seqnum, NULL);
	return;
	}
if(*extra == 0) {
	readsmj(smjnum);
	if(psmj->plong == 0) {
		debut = read_first_rec(klng, NULL) + 1;
		psmj->plong = debut;
		addlng(debut, seqnum);
		writesmj(smjnum);
		*extra = debut;
		return;
		}
	else debut = psmj->plong;
	}
else debut = *extra;
readlng(debut);
while(plng->next != 0) {
	debut = plng->next;
	readlng(debut);
	}
*extra = debut;
addlng(debut, seqnum);
}



char *nextpar(char *debut)
/* rend pos de la ) correspondant a la ( qui est en *debut
rend NULL si pas trouvee
*/
{
while( *(++debut) != 0) {
	if(*debut == '(') {
		debut = nextpar(debut);
		if(debut == NULL) return NULL;
		}
	else if(*debut == ')') {
		return debut;
		}
	}
return NULL;
}


void connect_locus_bib( int locnum, int seqnum, int bibnum, int *extra)
/* attacher locus-locnum-sequence-seqnuma reference bibnum
en utilisant et mettant a jour son extra data pointee par extra
*/
{
int debut;

mdshrt(kloc, locnum, 8, bibnum, NULL);

if( *extra == 0) {
	readbib(bibnum);
	if(pbib->plsub == 0) {
		readshrt(1);
		debut = pshrt->val + 1;
		pbib->plsub = debut;
		writebib(bibnum);
		addshrt(debut, seqnum);
		*extra = debut;
		return;
		}
	else debut = pbib->plsub;
	}
else	debut = *extra;
readshrt(debut);
while(pshrt->next != 0) {
	debut = pshrt->next;
	readshrt(debut);
	}
addshrt(debut, seqnum);
readshrt(debut);
if(pshrt->next != 0) debut = pshrt->next;
*extra = debut;
}


int process_year_bib(char *annee, int **extra)
{
static char aux[21];
int y_number, tot_txt, l;

sprintf(aux, "03%s", annee); majuscules(aux);
y_number = find_key(aux, smj_binary_tree, TRUE, extra);
if(y_number == 0) {
	fprintf(log_file, "cannot create %s\n", aux);
	non_chargee = TRUE;
	return 0;
	}
readsmj(y_number);
if(psmj->libel == 0) {
	tot_txt = read_first_rec(ktxt, NULL);
	tot_txt++;
	sprintf(ptxt, "published in %s", annee);
	l = strlen(ptxt); memset(ptxt + l, ' ', lrtxt - l);
	writetxt(tot_txt);
	write_first_rec(ktxt, tot_txt, 0);
	psmj->libel = tot_txt;
	writesmj(y_number);
	}
return y_number;
}


int next_start_seq(FILE *addr_file, char *gcgacnuc_value, char *ligne, 
	size_t lligne, off_t *p_addr, int *p_div)
{
char name[100];
int l;

if( fgets(name, sizeof(name), addr_file) == NULL) return FALSE;
/* enlever \n et espaces terminaux */
l = strlen(name) - 1;
name[l] = 0; while(name[--l] == ' ') name[l] = 0;
if(strncmp(name, "//", 2) == 0) {
	if(in_flat != NULL) fclose(in_flat);
	in_flat = NULL;
	return next_start_seq(addr_file, gcgacnuc_value, ligne, lligne, p_addr, p_div);
	}
else if(in_flat == NULL) {
	char flatname[200]; int div;
	strcpy(flatname, gcgacnuc_value);
	strcat(flatname, name);
	if(genbank) strcat(flatname, ".seq");
	else	strcat(flatname, ".dat");
	in_flat = fopen(flatname, "r");
	if(in_flat == NULL) {
		fprintf(log_file, "ERROR: cannot open sequence file %s\n", flatname);
		return FALSE;
		}
	/* recherche de la division a charger */
	for(div = 0; div <= divisions; div++) {
		if(strcmp(gcgname[div], name) == 0) break;
		}
	if(div > divisions) {
		fprintf(log_file, "ERROR: unknown division %s\n", name);
		return FALSE;
		}
	*p_div = div;
	return next_start_seq(addr_file, gcgacnuc_value, ligne, lligne, p_addr, p_div);
	}
{ unsigned long myul; sscanf(name, "%lu", &myul); *p_addr = myul; }
fseeko(in_flat, *p_addr, SEEK_SET);
fgets(ligne, lligne, in_flat);
if( ( genbank && strncmp(ligne, "LOCUS", 5) != 0 ) ||
	( (embl || swissprot) && strncmp(ligne, "ID", 2) != 0 ) ) {
	fprintf(log_file, "Warning: trouble with line %s of address file\n", 
		name);
	return next_start_seq(addr_file, gcgacnuc_value, ligne, lligne, p_addr, 
		p_div);
	}
return TRUE;
}


char *proc_classif(char *ligne, size_t lligne)
/* traitement de la partie classification des especes
Le dernier taxon connu est repere, les taxons suivants sont crees
au depart, ligne: 1ere des lignes de classification
au retour, ligne: ligne de type suivant
valeur rendue: nom du taxon ou classer cette espece
*/
{
char *p, *q, *debut, *fin, *last_pos;
const int width = sizeof(pspec->name);
static char buffer[1500], classif[1500], unclassif[] = "UNCLASSIFIED";
int num, total, nspec, last_known, l, doit;

/* lire toutes les lignes de classification car peuvent etre coupees
entre lignes dans GenBank */
debut = buffer; fin = buffer + sizeof(buffer) - 1; doit = TRUE;
while( *ligne == ' ' || 
	strncmp(ligne, "OC ", 3) == 0 || strncmp(ligne, "XX", 2) == 0  ) {
	if( *ligne == 'X' ) {
		fgets(ligne, lligne, in_flat);
		continue;
		}
	p = ligne + 2;
	while( *p == ' ') p++;
	p--;
	q = p + strlen(p) - 2;
	while( q > p && *q == ' ' ) q--;
	*(q + 1) = 0;
	l = q - p + 1;
	if(debut + l >= fin) doit = FALSE;
	if(doit) memcpy(debut, p, l); debut += l;
	fgets(ligne, lligne, in_flat);
	}
/* pas de classification ou trop longue */
if(debut == buffer || !doit) return unclassif; 
if(*(debut - 1) == '.') debut--;
*debut = 0;
/* memoriser tous les niveaux dans des chaines successives */
q = strtok(buffer, ";"); debut = classif; total = 0; 
fin = debut + sizeof(classif);
while(q != NULL) {
	while( *q == ' ') q++;
	l = strlen(q);
	if( l > width ) {
		l = width; q[width] = 0;
		}
	majuscules(q);
	if(debut + l + 1 >= fin ) return unclassif;
	strcpy(debut, q);
	total++; debut += l + 1;
	q = strtok(NULL, ";");
	}
if(total == 0) return unclassif;
/* recherche du dernier point connu de la classification */
debut = classif; last_pos = NULL;
for( num = 0; num < total; num++) {
	nspec = iknum(debut, kspec);
	if(nspec != 0) { last_known = num; last_pos = debut; }
	debut += strlen(debut) + 1;
	}
if(last_pos == NULL) { /* rien de connu dans la classification */
	crespecies(NULL, classif);
	last_pos = classif; last_known = 0;
	}
/* classer les derniers niveaux nouveaux */
debut = last_pos; num = last_known;
while( num < total - 1) {
	p = debut + strlen(debut) + 1;
	crespecies(debut, p);
	debut = p;
	num++;
	}
/* rendre la fin de la classification */
return debut;
}


int extend_acc_requests(void)
{
#define SLICE_REQUEST 100
static int max_request =  0;
struct _acc_request *new_tab;
int new_max_request;

if(tot_request < max_request) return FALSE;
new_max_request = max_request + SLICE_REQUEST;
new_tab = (struct _acc_request *)malloc(new_max_request * 
	sizeof(struct _acc_request) );
if(new_tab == NULL) return TRUE;
if(tot_request > 0) memcpy(new_tab, acc_request, 
	tot_request * sizeof(struct _acc_request) );
if(acc_request != NULL) free(acc_request);
acc_request = new_tab;
max_request = new_max_request;
return FALSE;
#undef SLICE_REQUEST
}

