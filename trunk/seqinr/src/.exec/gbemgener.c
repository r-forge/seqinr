/* 
gbemgener.c   generation de donnees GenBank ou EMBL sous acnuc
*/

#include "dir_acnuc.h"
#include <ctype.h>
#include <time.h>

#define MAX_LOC_QUAL 100000 /* attention!! aussi dans utilgener.c */
#define MEM_ONE_REC_TYPE 10000 /* attention!! aussi dans utilgener.c */

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

typedef struct {
	char *code;
	int number;
	} s_journal_code;
	
typedef void (*type_extra_function)(char *, int);

/* global variables */
FILE *in_flat, *log_file;
int seq_num, loc_num, non_chargee;
void *acc_binary_tree, *smj_binary_tree, *bib_binary_tree, *aut_binary_tree;
void *dynlist_journals;

/* local functions */
int process_1_sequence(int *p_status, int last_sub, int *p_last_loc, 
	char *ligne, size_t lligne, int use_address, int *p_tot_read, 
	off_t addr_info, off_t gcg_addr, int current_div,
	int tot_types, int *type_ranks, char **type_codes, char **type_abbrevs, 
	int *type_num,
	char *location, char *qualifiers, int type_locus, int *extra_locus,
	int current_rel, int *p_tot_charge, int *p_tot_filles_cre,
	FILE *seq_file, char *err_file_name);
int process_locus(char *ligne, int *p_last_sub, int last_loc, char *seq_name);
int process_id(char *ligne, int *p_last_sub, int last_loc, char *seq_name);
void process_dt(char *ligne, size_t lligne);
off_t calc_addr_seq(char *ligne, size_t lligne);
int process_accession(char *ligne, size_t taille);
void complete_acc_request(int loc_num, int *p_fille_num,
	char *location, char *qualifiers, int max_loc_qual);
void process_os(char *ligne, size_t lligne);
void extract_spec_name(char *ligne, char *outname);
char *organelle_from_line(char *ligne);
void process_source(char *ligne, size_t lligne);
static char *read_one_rec_type(char *ligne, size_t lligne);
void process_keywords(char *ligne, size_t lligne);
int process_features(FILE *in_flat, char *ligne, size_t lligne, 
	int tot_types, int *type_ranks, char **type_codes, char **type_abbrevs, 
	int *type_num,
	int current_div, char *location, char *qualifiers, int max_loc_qual);
void process_ft_key(int *type_ranks, char **type_codes,
			char **type_abbrevs, int *type_num, int type,
			char *location, int maxlocat, char *nmere, 
			int *p_fille_num,
			char *feat_name, long info_addr, int div, 
			char *qualifiers);
void process_exon(char *exon_qualifs, char *location, char *qualifiers, 
	int max_loc_qual,
	int tot_types, char **type_codes, char **type_abbrevs, int *type_ranks,
	int *p_fille_num);
int read_flat_loc_qualif(FILE *in_flat, char *ligne, size_t lligne, 
	off_t *next_feat_addr,
	char *location, int maxlocat, char *type, char *qualif, int maxqualif);
void process_reference(char *ligne, size_t lligne, int numref);
void process_journal(char *ligne, size_t lligne, char *author_data, int numref);
static char *read_author_lines(char *ligne, size_t lligne);
void process_authors(char *author_data, int num_bib, int first_aut_only);
char *next_auth_name(char *auth_data, char *auth_name, int max_name,
	int first_aut_only);
void format_auteur(char *auteur);
void process_book(char *journal_data, char *author_data);
void process_thesis(char *journal_data, char *author_data);
void process_unpublished(char *journal_data, char *author_data, int numref);
void process_patent(char *journal_data, char *author_data);
void set_reference(char *reference, 
	char *code_j, int num_j, int *extra_journal, 
	char *annee, char *author_data, int first_aut_only);
void *load_j_code_libel(void);
void remove_current_seq(int mere, int locus);
void remove_subseq_record(int num);
void process_origin(char *ligne);
void prep_chromosome(char *p, int seq);


/* external prototypes */
int mkextr(char *location, int seq_num, int *p_fille_num, char *name,
	int smj_type, char *feat_name, long info_addr, int div, char *qualif);
int find_label(char *label, char *access, long *pannot, int *div, int *isub);
char *next_word(char *ligne, char *word, char *separs, int max_word, 
	int *erreur);
void close_if_open(void);
int *get_smj_extra(int smjnum);
int *get_keyw_extra(int kwnum);
int *get_spec_extra(int specnum);
int *get_bib_extra(int bibnum);
void fast_add_seq_to_spec(int specnum, int seqnum);
void fast_add_seq_to_keyw(int kwnum, int seqnum);
void add_to_smjlng( int smjnum, int seqnum, int *extra);
int trim_key(char *name);
int process_year_bib(char *annee, int **extra);
void connect_locus_bib( int locnum, int seqnum, int bibnum, int *extra);
int next_start_seq(FILE *addr_file, char *gcgacnuc_value, char *ligne, 
	size_t lligne, off_t *p_addr, int *p_div);
int find_key(char *name, void *tree, int create_if_not_found,
	int **p_aux_data);
void *load_index_bt(DIR_FILE *k, int char_width, int use_aux_data);
int poffset(int div, int offset);
void goffset(int point,int *div,int *offset);
void cre_new_division(char *name);
void update_div_size(int numdiv);
int loadtypes( int **type_ranks, char ***type_codes, char ***type_abbrevs);
int read_addr_loc_qualif(long faddr, int div, 
	char *location, int maxlocat, char *type, char *qualif, int maxqualif);
char *get_qualif_value(char *p);
char *proc_classif(char *ligne, size_t lligne);
#ifdef unix
void proc_mmap_options(int argc, char **argv);
int check_term(void); /* gestion de la possibilite d'interruption propre */
void write_quick_meres(void);
#endif
void *init_dynlist(void);
char *find_g_dynlist(char *name, void *tree, int create_if_not_found, 
	void **node);
char *next_dynlist(void *arbre, void **p_noeud);
void *dynlist_get_extra(void *node);
void dynlist_set_extra(void *node, void *extra);
void *calc_meres_to_do(FILE *in_flat, int current_div, FILE *seq_file,
	FILE *log_file);
off_t fpos_from_to_do(void *loop_to_do);
void remove_modified_seqs(void *dynlist_meres_to_do, 
	void **p_dynlist_mere_with_request, FILE *log_file);
void sup_fille_request_from_mere(char *nom_mere, 
	void *dynlist_mere_with_request);
void process_fille_request(void *p_request, char *location, char *qualifiers,
	int max_loc_qual, int *p_fille_num, int *p_tot_filles_cre, 
	FILE *log_file);
void add_liste_meres(int seqnum);
void process_genetic_code(int mere);
off_t gcg_fpos_from_to_do(void *loop_to_do);
void *init_acchash(void);
int find_cre_acc(char *name, int create);

/* external globals */
extern struct _acc_request  *acc_request;
extern int tot_request;
extern int acc_length_on_disk;
extern char *gcgname[];
extern int divisions;


/* permettre traitement optionnel d'autres qualifiers 
en compilant avec -DEXTRA_QUALIF 
*/
#ifdef EXTRA_QUALIF 

extern void process_extra_qualif(char *, int);
type_extra_function extra_function = process_extra_qualif;

/* modele de fonction a ecrire dans un fichier separe:
#include "dir_acnuc.h"
char *get_qualif_value(char *);
int crekeyword(char *, char *);
void fast_add_seq_to_keyw(int, int);

void process_extra_qualif(char *qualif, int fille_num)
{
char *p, *qualif_value;
int numkey, erreur;

p = strstr(qualif, "/MYQUALIF=");
if( p != NULL && (qualif_value = get_qualif_value(p) ) != NULL ) {
	numkey = crekeyword("MYKEYWORD", qualif_value);
	erreur = mdshrt(ksub, fille_num, 4, numkey, NULL);
	if(erreur != 2) fast_add_seq_to_keyw(numkey, fille_num);
	}
}
*/

#else

type_extra_function extra_function = NULL;

#endif



int main(int argc, char *argv[])
{
static char flat_name[100], ligne[100], division_name[21], 
	location[MAX_LOC_QUAL], qualifiers[MAX_LOC_QUAL];
int current_div, tot_read, tot_charge, 
	last_sub, last_loc, tot_filles_cre, status,
	type_locus, *extra_locus, current_rel, use_address;
int tot_types, *type_ranks, *type_num;
char **type_codes, **type_abbrevs;
time_t heure_debut, heure_fin;
FILE *addr_file, *seq_file;
char argument[200], gcgacnuc_value[200], err_file_name[200], 
	seqname[L_MNEMO + 1];
void *dynlist_meres_to_do = NULL, *loop_to_do, *request;
void *dynlist_mere_with_request = NULL;
off_t addr_info, gcg_addr;

#ifdef unix
if(argc < 3 || ( argv[1][0] != 'a' && argv[1][0] != 'd' ) ) {
	fprintf(stderr, "Usage: gbemgener a name-of-address-file [-mmap name]\n"
		"or\n"
		        "       gbemgener d name-of-division [-mmap name]\n");
	exit(ERREUR);
	}
use_address = (argv[1][0] == 'a');
strcpy(argument, argv[2]);

#else

printf("address (a) or division (d)? ");
gets(argument);
use_address = (*argument == 'a' || *argument == 'A');
if(use_address) printf("name of address file? ");
else	printf("division name? ");
gets(argument);

#endif

log_file = stdout;
time(&heure_debut);
fprintf(log_file, "Program started at %s\n", asctime(localtime(&heure_debut)));
dir_acnucopen("WP");
#ifdef unix
proc_mmap_options(argc, argv);
/* preparer possibilite d'interruption propre */
check_term();
#endif
atexit(close_if_open);

if( ! ( genbank || embl ) ) {
	fprintf(log_file, "ERROR: Should open a GenBank or EMBL data base\n");	
	exit(ERREUR);
	}
if( use_address && !flat_format ) {
	fprintf(log_file, "ERROR: cannot use address file for GCG format\n");	
	exit(ERREUR);
	}
	{
	char *p;
	p = prepare_env_var("gcgacnuc");
	if(p == NULL) {
		fprintf(log_file, "ERROR: gcgacnuc is not defined\n");	
		exit(ERREUR);
		}
	strcpy(gcgacnuc_value, p);
	}
if(use_address) {
	char *p;
	addr_file = fopen(argument, "r");
	if(addr_file == NULL) {
		fprintf(log_file, "ERROR: Cannot open address file %s\n", 
			argument);	
		exit(ERREUR);
		}
	in_flat = NULL;
	strcpy(err_file_name, argument);
	p = strchr(err_file_name, '.');
	if(p == NULL) p = err_file_name + strlen(err_file_name);
	strcpy(p, ".err");
	}
else	{
	strcpy(division_name, argument);

	/* ouverture du fichier de la division a charger */
	strcpy(flat_name, gcgacnuc_value);
	strcat(flat_name, division_name);
	if(!flat_format) strcat(flat_name, ".ref");
	else if(genbank) strcat(flat_name, ".seq");
	else strcat(flat_name, ".dat");
	in_flat = fopen(flat_name, "r");
	if(in_flat == NULL) {
		fprintf(log_file, "ERROR: cannot open flat file %s\n", 
			flat_name);	
		exit(ERREUR);
		}
	if(!flat_format) { /* ouverture fichier .seq de gcg */
		strcpy(flat_name + strlen(flat_name) - 4, ".seq");
		seq_file = fopen(flat_name, "r");
		if(seq_file == NULL) {
			fprintf(log_file, "ERROR: cannot open .seq file %s\n", 
				flat_name);	
			exit(ERREUR);
			}
		}
		
	/* recherche ou creation de la division a charger */
	for(current_div = 0; current_div <= divisions; current_div++) {
		if(strcmp(gcgname[current_div], division_name) == 0) break;
		}
	if(current_div > divisions) {
		if(divisions >= 0) update_div_size(divisions);
		cre_new_division(division_name);
		fprintf(log_file, "New division created: %s\n", division_name);
		current_div = divisions;
		}
	else if(current_div == divisions) update_div_size(divisions);
	dir_acnucflush();
	/* calculer listes des seqs de la division nouvelles ou modifiees */
	if( (dynlist_meres_to_do = 
			calc_meres_to_do(in_flat, current_div, seq_file, 
				log_file) ) == NULL) {
		fprintf(log_file, "Not enough memory in calc_meres_to_do\n");
		exit(ERREUR);
		}
	/* supprimer les seqs moins recentes dans acnuc que dans division */
	remove_modified_seqs(dynlist_meres_to_do, 
		&dynlist_mere_with_request, log_file);
	}

/* chargement des arbres binaires */
acc_binary_tree = NULL; /* this characterizes the new procedure: hash for acc */
time(&heure_debut);
fprintf(log_file, "Loading acc nos started at %s", 
	asctime(localtime(&heure_debut)) ); fflush(log_file);
if(init_acchash() == NULL) {
	fprintf(log_file, 
		"ERROR: not enough memory to hash accession numbers\n");	
	exit(ERREUR);
	}
time(&heure_debut);
fprintf(log_file, "finished at %s", asctime(localtime(&heure_debut)) ); 
fflush(log_file);

smj_binary_tree = load_index_bt(ksmj, sizeof(psmj->name), TRUE);
if(smj_binary_tree == NULL) {
	fprintf(log_file, "ERROR: not enough memory to load SMJYT\n");	
	exit(ERREUR);
	}
dynlist_journals = load_j_code_libel();

bib_binary_tree = load_index_bt(kbib, sizeof(pbib->name), TRUE);
if(bib_binary_tree == NULL) {
	fprintf(log_file, "ERROR: not enough memory to load BIBLIO\n");	
	exit(ERREUR);
	}

aut_binary_tree = load_index_bt(kaut, sizeof(paut->name), FALSE);
if(aut_binary_tree == NULL) {
	fprintf(log_file, "ERROR: not enough memory to load AUTHOR\n");	
	exit(ERREUR);
	}

/* chargement des types connus */
tot_types = loadtypes( &type_ranks, &type_codes, &type_abbrevs);
type_num = malloc(tot_types * sizeof(int));
if(type_num == NULL) {
	fprintf(log_file, "ERROR: not enough memory to load type_num\n");	
	exit(ERREUR);
	}
/* memoriser type LOCUS */
if(genbank) type_locus = find_key("04LOCUS", smj_binary_tree, TRUE, 
					&extra_locus);
else type_locus = find_key("04ID", smj_binary_tree, TRUE, &extra_locus);
if(type_locus == 0) {
	fprintf(log_file, "cannot create 04LOCUS / 04ID\n");
	exit(ERREUR);
	}
if(genbank) {
/* prepare for release number: look for first desc of keyw RELEASE NUMBERS
and put its number in current_rel, or 0 if not found. 
*/
	current_rel = iknum("RELEASE NUMBERS", kkey);
	if(current_rel != 0) {
		readkey(current_rel);
		readshrt(pkey->desc);
		if(pshrt->next != 0) {
			readshrt(pshrt->next);
			readshrt(pshrt->val);
			current_rel = abs(pshrt->val);
			}
		else	current_rel = 0;
		}
	}
last_sub = read_first_rec(ksub, NULL);
last_loc = read_first_rec(kloc, NULL);
tot_read = tot_charge = tot_filles_cre = 0;
time(&heure_debut);
fprintf(log_file, "Sequence loading started at %s\n", asctime(localtime(&heure_debut)) );

loop_to_do = NULL;
while( TRUE ) {
	char *p;
	
	if(use_address) {
		if( ! next_start_seq(addr_file, gcgacnuc_value, ligne, 
			sizeof(ligne), &addr_info, &current_div ) ) break;
		}
	else	{
		if( (p = next_dynlist(dynlist_meres_to_do, &loop_to_do) ) 
			== NULL) break;
		addr_info = fpos_from_to_do(loop_to_do);
		if(!flat_format) gcg_addr = gcg_fpos_from_to_do(loop_to_do);
		fseeko(in_flat, addr_info, SEEK_SET);
		fgets(ligne, sizeof(ligne), in_flat);
		strcpy(seqname, p);
		sup_fille_request_from_mere(seqname, dynlist_mere_with_request);
		}
	last_sub = process_1_sequence(&status, last_sub, &last_loc, 
		ligne, sizeof(ligne),
		use_address, &tot_read, addr_info, gcg_addr, current_div,
		tot_types, type_ranks, type_codes, type_abbrevs, type_num,
		location, qualifiers, type_locus, extra_locus,
		current_rel, &tot_charge, &tot_filles_cre,
		seq_file, err_file_name);
	fflush(log_file);
	if(status != 0) break;
#ifdef unix
	if( check_term() ) { /* tester si interruption demandee */
		fprintf(log_file, "Clean program interruption\n");
		break; 
		}
#endif
	} /* end of loop over all seqs */
/* traiter les filles supprimees et peut-etre non re-chargees */
loop_to_do = NULL;
while( next_dynlist(dynlist_mere_with_request, &loop_to_do) != NULL) {
	request = dynlist_get_extra(loop_to_do);
	process_fille_request(request, location,
		qualifiers, MAX_LOC_QUAL, &last_sub, &tot_filles_cre, log_file);
	}
	
if(tot_request > 0) {
/* unsatisfied accession # requests */
	int num;
	fprintf(log_file, 
	"\nWarning: %d accession # requests remain unsatisfied; list follows\n",
		tot_request);
	for(num = 0; num < tot_request; num++) {
		readsub(acc_request[num].mere);
		read_addr_loc_qualif(acc_request[num].inf, acc_request[num].div,
			location, MAX_LOC_QUAL, NULL, NULL, 0);
		fprintf(log_file, "%.16s calls %s in %s\n",
			psub->name, acc_request[num].acc, location);
		}
	}
dir_acnucclose();
if(use_address) fclose(addr_file);
if( !flat_format ) fclose(seq_file);
time(&heure_fin);
#ifdef unix
if(tot_charge > 0) write_quick_meres();
#endif	
fprintf(log_file, "Program finished at %s\n", asctime(localtime(&heure_fin)) );
fprintf(log_file, "\nlues=%d  chargees=%d  difference=%d\n"
	"created subsequences=%d  seqs/second=%.2f\n", 
	tot_read, tot_charge, tot_read - tot_charge, tot_filles_cre,
	tot_read / (float)(heure_fin - heure_debut) );

return 0;
}  /* end of main */


int process_1_sequence(int *p_status, int last_sub, int *p_last_loc, 
	char *ligne, size_t lligne,
	int use_address, int *p_tot_read, off_t addr_info, off_t gcg_addr, 
	int current_div, int tot_types, int *type_ranks, char **type_codes, 
	char **type_abbrevs, int *type_num,
	char *location, char *qualifiers, int type_locus, int *extra_locus,
	int current_rel, int *p_tot_charge, int *p_tot_filles_cre,
	FILE *seq_file, char *err_file_name)
{
int old_last_sub, numref, old_tot_request, match_pending;
char seq_name[L_MNEMO + 1];
off_t addr_seq;

old_last_sub = last_sub;
if(genbank) seq_num = process_locus(ligne, &last_sub, *p_last_loc, 
						seq_name);
else seq_num = process_id(ligne, &last_sub, *p_last_loc, seq_name);
if(use_address || seq_num != 0) 
	fprintf(log_file, "---->%s\n", seq_name); fflush(log_file);
if(seq_num == 0) {
	if(use_address) fprintf(log_file, "already in database\n");
	*p_status = 0;
	return old_last_sub;
	}
(*p_tot_read)++;
old_tot_request = tot_request; match_pending = FALSE;
if( flat_format ) addr_seq = calc_addr_seq(ligne, lligne);
else	addr_seq = gcg_addr; /* cas GCG */
	
if( !big_annots ){ 
	/* check addresses not too large for encoding */
	int coded, div, addr2, addr;
	addr = (unsigned)( addr_info > addr_seq ? addr_info : addr_seq );
	coded = poffset(current_div, addr);
	goffset(coded, &div, &addr2);
	if(div != current_div || addr2 != addr) {
		fprintf(log_file, 
			"ERROR: file too big for encoding address\n");
		*p_status = 1;
		return old_last_sub;
		}
	addr_info = poffset(current_div, (unsigned)addr_info);
	addr_seq = poffset(current_div, (unsigned)addr_seq);
	}
readloc(loc_num);
ploc->pinf = (unsigned)addr_info;
ploc->pnuc = (unsigned)addr_seq;
if(big_annots) ploc->div = current_div;
writeloc(loc_num);
addhsh(seq_num, ksub); /* connecter au hashing des mnemos */
numref = 0;
fgets(ligne, lligne, in_flat);
if(genbank) {
    do	{
	if(strncmp(ligne, "ACCESSION", 9) == 0) {
		match_pending = process_accession(ligne, lligne);
		}
	else if(strncmp(ligne, "KEYWORDS", 8) == 0) {
		process_keywords(ligne, lligne);
		}
	else if(strncmp(ligne, "SOURCE", 6) == 0) {
		process_source(ligne, lligne);
	}
	else if(strncmp(ligne, "REFERENCE ", 10) == 0) {
		numref++;
		process_reference(ligne, lligne, numref);
		}
	else if(strncmp(ligne, "FEATURES", 8) == 0) {
		last_sub = process_features(in_flat, ligne, 
			lligne, 
			tot_types, type_ranks, type_codes, type_abbrevs,
			type_num,
			current_div, location, qualifiers, 
			MAX_LOC_QUAL);
		}
	else 	{
		fgets(ligne, lligne, in_flat);
		}
	}
    while( (!non_chargee) && strncmp(ligne, "ORIGIN", 6) != 0);
    if(strncmp(ligne, "ORIGIN", 6) == 0) process_origin(ligne);
    }
else { /* embl */
    do	{
	if(strncmp(ligne, "AC", 2) == 0) {
		match_pending = process_accession(ligne, lligne);
		}
	else if(strncmp(ligne, "DT", 2) == 0) {
		process_dt(ligne, lligne);
		}
	else if(strncmp(ligne, "KW", 2) == 0) {
		process_keywords(ligne, lligne);
		}
	else if(strncmp(ligne, "OS", 2) == 0) {
		process_os(ligne, lligne);
		}
	else if(strncmp(ligne, "RN", 2) == 0) {
		numref++;
		process_reference(ligne, lligne, numref);
		}
	else if(strncmp(ligne, "FH", 2) == 0) {
		last_sub = process_features(in_flat, ligne, 
			lligne, 
			tot_types, type_ranks, type_codes, type_abbrevs,
			type_num,
			current_div, location, qualifiers, 
			MAX_LOC_QUAL);
		}
	else 	{
		fgets(ligne, lligne, in_flat);
		}
	}
    while( (!non_chargee) && strncmp(ligne, "SQ", 2) != 0);
    }
if( non_chargee ) {
/* ici defaire liens faits pour seq non chargee */
	fprintf(log_file, "this sequence has not been inserted\n");
	if(use_address) {
		static int pre_div = -1;
		static FILE *err_file = NULL;
		if(err_file == NULL) err_file = fopen(err_file_name, "w");
		if(err_file != NULL) {
			if(current_div != pre_div) {
				if(pre_div != -1) fprintf(err_file, "//\n");
				pre_div = current_div;
				fprintf(err_file, "%s\n", gcgname[current_div]);
				}
			fprintf(err_file, "%s\n", seq_name);
			fflush(err_file);
			}
		}
	remove_current_seq(seq_num, loc_num);
	tot_request = old_tot_request;
	last_sub = old_last_sub;
	}
else	{
	if(match_pending) complete_acc_request(loc_num, &last_sub,
			location, qualifiers, MAX_LOC_QUAL);
	add_liste_meres(seq_num);
	/* mettre type de la mere a LOCUS sauf si deja fait */
	readsub(seq_num);
	if(psub->type == 0) {
		add_to_smjlng(type_locus, seq_num, extra_locus);
		psub->type = type_locus;
		writesub(seq_num);
		}
	/* connecter avec numero de release */
	if(genbank && current_rel != 0) {
		int err;
		err = mdshrt(ksub, seq_num, 4, current_rel, NULL);
		if( err != 2 ) 
			fast_add_seq_to_keyw(current_rel, seq_num);
		}
	/* si pas d'espece, mettre UNKNOWN */
	readloc(loc_num);
	if(ploc->spec == 0) {
		int unknown_spec;
		unknown_spec = crespecies("UNCLASSIFIED", "UNKNOWN");
		fast_add_seq_to_spec(unknown_spec, seq_num);
		ploc->spec = unknown_spec;
		writeloc(loc_num);
		}
	/* mettre le code genetique des filles de type CDS */
	process_genetic_code(seq_num);
	*p_last_loc = loc_num;
	write_first_rec(ksub, last_sub, 0);
	write_first_rec(kloc, *p_last_loc, 0);
	(*p_tot_charge)++;
	*p_tot_filles_cre += (last_sub - old_last_sub - 1);
	}
*p_status = 0;
return last_sub;
}


int process_locus(char *ligne, int *p_last_sub, int last_loc, char *seq_name)
{
int num, longueur, *extra, num_mol, mois, numkey, err;
char *p;
static char molecule[20] = "01", acnuc_date[9], 
	flat_divname[20] = "DIVISION ";
static char month[]="JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC";
char *name, *molec, *div3, *date;
int is_circ;

longueur = decode_locus_rec(ligne, &name, &molec, &is_circ, &div3, &date);
name[L_MNEMO] = 0;
strcpy(seq_name, name);
num = isenum(seq_name);
if(num != 0) return 0;
(*p_last_sub)++;
seq_num = num = *p_last_sub;
loc_num = last_loc + 1;
non_chargee = FALSE;
memset(psub, 0, lrsub);
padtosize(psub->name, seq_name, L_MNEMO);
psub->length = longueur;
psub->plinf = loc_num;
writesub(num);
memset(ploc, 0, lrloc);
ploc->sub = num;
/* traitement de la date */
date[2] = date[6] = 0;
mois = (strstr(month, date + 3) - month)/3 + 1;
sprintf(acnuc_date, "%2.2d/%s/%s", mois, date, date + 9);
memcpy(ploc->date, acnuc_date, 8);
memcpy(ploc->date + 8, acnuc_date, 8);
/* traitement molecule */
if( molec[0] == ' ' || molec[0] == '-' ) 
	strcpy(molecule + 2, "DNA");
else 	strcpy(molecule + 2, molec);
majuscules(molecule + 2); trim_key(molecule);
num_mol = find_key(molecule, smj_binary_tree, TRUE, &extra);
if(num_mol == 0) {
	fprintf(log_file, "cannot create %s\n", molecule);
	non_chargee = TRUE;
	}
else	{
	add_to_smjlng(num_mol, seq_num, extra);
	ploc->molec = num_mol;
	}
/* traitement champ circular */
if( is_circ) {
	numkey = crekeyword(NULL, "CIRCULAR");
	err = mdshrt(ksub, seq_num, 4, numkey, NULL);
	if(err != 2) fast_add_seq_to_keyw(numkey, seq_num);
	}
/* traitement division */
if( div3[0] != ' ') {
	strcpy(flat_divname + 9, div3); flat_divname[12] = 0;
	numkey = crekeyword("DIVISION NAMES", flat_divname);
	err = mdshrt(ksub, seq_num, 4, numkey, NULL);
	if(err != 2) fast_add_seq_to_keyw(numkey, seq_num);
	}
writeloc(loc_num);
return num;
}


int process_id(char *ligne, int *p_last_sub, int last_loc, char *seq_name)
{
int num, longueur, numstat, *extra, num_mol, circular = FALSE;
char *p;
static char status[20] = "00", molecule[20] = "01", 
	flat_divname[20] = "DIVISION ";

majuscules(ligne);
p = next_word(ligne + 5, seq_name, "; ", L_MNEMO, NULL);
p = next_word(p, status + 2, "; ", sizeof(status) - 2, NULL);
if(strncmp(p, "CIRCULAR", 8) == 0) {
	circular = TRUE;
	p += 8; 
	while( *p == ' ') p++;
	}
p = next_word(p, molecule + 2, "; ", sizeof(molecule) - 2, NULL);
p = next_word(p, flat_divname + 9, "; ", sizeof(flat_divname) - 9, NULL);
sscanf(p, "%d", &longueur);
num = isenum(seq_name);
if(num != 0) return 0;
++(*p_last_sub);
seq_num = num = *p_last_sub;
loc_num = last_loc + 1;
non_chargee = FALSE;
memset(psub, 0, lrsub);
padtosize(psub->name, seq_name, L_MNEMO);
psub->length = longueur;
psub->plinf = loc_num;
writesub(num);
/* traitement status */
numstat = find_key(status, smj_binary_tree, TRUE, &extra);
if(numstat == 0) {
	fprintf(log_file, "cannot create %s\n", status);
	non_chargee = TRUE;
	}
else	add_to_smjlng( numstat, seq_num, extra);
/* traitement molecule */
if(strcmp(molecule + 2, "XXX") == 0) strcpy(molecule + 2, "DNA");
num_mol = find_key(molecule, smj_binary_tree, TRUE, &extra);
if(num_mol == 0) {
	fprintf(log_file, "cannot create %s\n", molecule);
	non_chargee = TRUE;
	}
else	{
	add_to_smjlng(num_mol, seq_num, extra);
	}
memset(ploc, 0, lrloc);
ploc->sub = num;
ploc->stat = numstat;
ploc->molec = num_mol;
writeloc(loc_num);
/* traitement circular */
if( circular ) {
	int err, numkey;
	numkey = crekeyword(NULL, "CIRCULAR");
	err = mdshrt(ksub, seq_num, 4, numkey, NULL);
	if(err != 2) fast_add_seq_to_keyw(numkey, seq_num);
	}
/* traitement division */
if( flat_divname[9] != ' ') {
	int err, numkey;
	numkey = crekeyword("DIVISION NAMES", flat_divname);
	err = mdshrt(ksub, seq_num, 4, numkey, NULL);
	if(err != 2) fast_add_seq_to_keyw(numkey, seq_num);
	}
return num;
}


void process_dt(char *ligne, size_t lligne)
{
char *p;
static char date[20], acnuc_date[9];
static char month[]="JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC";
static char rel_name[30] = "RELEASE ";
int mois, relnum, numkey;

fgets(ligne, lligne, in_flat);
/* traitement de la date */
next_word(ligne + 5, date, " ", sizeof(date), NULL);
date[2] = date[6] = date[11] = 0;
mois = (strstr(month, date + 3) - month)/3 + 1;
sprintf(acnuc_date, "%2.2d/%s/%s", mois, date, date + 9);
readloc(loc_num);
memcpy(ploc->date, acnuc_date, 8);
memcpy(ploc->date + 8, acnuc_date, 8);
writeloc(loc_num);

/* traitement des release numbers */
relnum  = -1;
p = strstr(ligne, "(Rel.");
if(p != NULL) sscanf(p + 5, "%d", &relnum);
if(relnum != -1) {
	int err;
	sprintf(rel_name + 8, "%.2d", relnum);
	numkey = crekeyword("RELEASE NUMBERS", rel_name);
	err = mdshrt(ksub, seq_num, 4, numkey, NULL);
	if(err != 2) fast_add_seq_to_keyw(numkey, seq_num);
	}
fgets(ligne, lligne, in_flat);
}


off_t calc_addr_seq(char *ligne, size_t lligne)
/* ok pour genbank et embl */
{
off_t addr_seq, current;

current = ftello(in_flat);
if(genbank) {
	do	fgets(ligne, lligne, in_flat);
	while(strncmp(ligne, "ORIGIN", 6) != 0);
	}
else	{
	do	fgets(ligne, lligne, in_flat);
	while(strncmp(ligne, "SQ", 2) != 0);
	}
addr_seq = ftello(in_flat);
fseeko(in_flat, current, SEEK_SET);
return addr_seq;
}


int process_accession(char *ligne, size_t lligne)
/* ok pour genbank et embl */
/* returns TRUE if an accession number is in list of pending acc 
*/
{
char accession[50], *p;
int numacc, marge, num_back = -1;

if(genbank) marge = 12;
else marge = 5;
do	{
	p = ligne + marge;
	do	{
		p = next_word(p, accession, ";\n ", sizeof(accession), NULL);
		if(strlen(accession) <= ACC_LENGTH)
			numacc = find_cre_acc(accession, TRUE);
		else	numacc = 0;
		if(numacc == 0) {
			fprintf(log_file, "cannot create %s\n", accession);
			non_chargee = TRUE;
			}
		else	{
			mdshrt(kacc, numacc, 1, seq_num, NULL);
			mdshrt(kloc, loc_num, 10, numacc, NULL);
			if(num_back == -1) {
				for(num_back = 0; num_back < tot_request;
							num_back++)
					if(strcmp(acc_request[num_back].acc, 
						accession) == 0) break;
				if(num_back >= tot_request) num_back = -1;
				}
			}
		}
	while(p != NULL);
	fgets(ligne, lligne, in_flat);
	}
while( ( genbank && *ligne == ' ' ) || 
	( embl && strncmp(ligne, "AC", 2) == 0 ) );
return (num_back != -1);
}


void complete_acc_request(int loc_num, int *p_fille_num,
	char *location, char *qualifiers, int max_loc_qual)
{
int next, pre, erreur, created, num_back;
static char access[ACC_LENGTH + 1], feat_name[20], nom_fille[L_MNEMO + 1];

readloc(loc_num);
next = ploc->placc;
while(next != 0) {
	readshrt(next); next = pshrt->next;
	readacc(pshrt->val);
	memcpy(access, pacc->name, ACC_LENGTH);
	access[ACC_LENGTH] = 0;
	trim_key(access);
	pre = 0;
	suite:
	for(num_back = pre; num_back < tot_request; num_back++) {
		if(strcmp(acc_request[num_back].acc, access) == 0) break;
		}
	if(num_back >= tot_request) continue;
	/* attempt creating new sequence */
	erreur = read_addr_loc_qualif(acc_request[num_back].inf, 
			acc_request[num_back].div, location, 
			max_loc_qual, feat_name, qualifiers, max_loc_qual);
	if(erreur) goto suite;
	strcpy(nom_fille, acc_request[num_back].nom);
	created = mkextr(location, acc_request[num_back].mere, p_fille_num, 
		nom_fille, acc_request[num_back].type, feat_name, 
		acc_request[num_back].inf, acc_request[num_back].div, 
		qualifiers);
/* if request not successful, goto next request for same accession # */
	if(created) pre = 0;
	else	pre = num_back + 1;
	goto suite;
	}
}


void process_os(char *ligne, size_t lligne)
{
static char spec_1[300], spec_2[100], organelle[100], org_code[50] = "05";
int numspec, num_org, l, *extra;
char *ascendant;

organelle[0] = 0; spec_2[0] = 0; ascendant = NULL;
extract_spec_name(ligne + 5, spec_1);
do	fgets(ligne, lligne, in_flat);
while (strncmp(ligne, "OS", 2) == 0 || strncmp(ligne, "XX", 2) == 0 );
/* traiter lignes de classification */
if(strncmp(ligne, "OC", 2) == 0) {
	if(iknum(spec_1, kspec) == 0) { /* new species, classify it */
		ascendant = proc_classif(ligne, lligne);
		}
	else	{ /* known species, skip classification data */
		ascendant = NULL;
		do	fgets(ligne, lligne, in_flat);
		while (strncmp(ligne, "OC", 2) == 0 || 
					strncmp(ligne, "XX", 2) == 0 );
		}
	}
if(strncmp(ligne, "OG", 2) == 0) {
	l = strlen(ligne) - 1;
	memcpy(organelle, ligne + 5, l - 5);
	organelle[l - 5] = 0;
	do	fgets(ligne, lligne, in_flat);
	while (strncmp(ligne, "OG", 2) == 0 || strncmp(ligne, "XX", 2) == 0 );
	}
if(strncmp(ligne, "OS", 2) == 0 ) { /* cas de 2 lignes OS */
	extract_spec_name(ligne + 5, spec_2);
	do	fgets(ligne, lligne, in_flat);
	while (*ligne == 'O' || strncmp(ligne, "XX", 2) == 0 );
	}
if(organelle[0] != 0) { 
	majuscules(organelle);
	num_org = TRUE;
	if(strncmp(organelle, "MITOCHONDRION", 13) == 0 ) 
		strcpy(org_code + 2, "MITOCHONDRION");
	else if(strncmp(organelle, "PLASTID", 7) == 0 ) 
		strcpy(org_code + 2, "CHLOROPLAST");
	else	num_org = FALSE;
	if(num_org) {
		num_org = find_key(org_code, smj_binary_tree, TRUE, &extra);
		if(num_org == 0) {
			fprintf(log_file, "cannot create %s\n", org_code);
			non_chargee = TRUE;
			}
		else	{
			add_to_smjlng(num_org, seq_num, extra);
			readloc(loc_num);
			ploc->org = num_org;
			writeloc(loc_num);
			}
		}
	}
if(spec_2[0] != 0) { /* cas avec deux especes */
	strcat(spec_1, " X ");
	strcat(spec_1, spec_2);
	crespecies("HYBRIDS", spec_1);
	ascendant = NULL;
	}
if(ascendant != NULL)
	numspec = crespecies(ascendant, spec_1);
else
	numspec = crespecies("UNCLASSIFIED", spec_1);
fast_add_seq_to_spec(numspec, seq_num);
readloc(loc_num);
ploc->spec = numspec;
writeloc(loc_num);
}


void extract_spec_name(char *ligne, char *outname)
{
char *p;
int l;

majuscules(ligne);
p = ligne + strlen(ligne) - 2;
while(*p == ' ') p--;
/* enlever le contenu de () finale sauf si INFLUENZA VIRUS apparait */
if(strstr(ligne, "VIRUS") == NULL || strstr(ligne, "INFLUENZA") == NULL) {
	if(*p == ')') {
		do p--;
		while( *p != '(' && p > ligne + 1);
		p--;
		}
	while(*p == ' ') p--;
	}
l = p - ligne + 1;
memcpy(outname, ligne, l); outname[l] = 0;
}


char *organelle_from_line(char *ligne)
{
int i;
char *p;
static char *list[] = {
	"MITOCHONDRION"/*0*/,"KINETOPLAST"/*1*/,
	"CHLOROPLAST"/*2*/,
	"PLASTID"/*3*/,"APICOPLAST","CHROMOPLAST","CYANELLE",
	"LEUCOPLAST","PROTOPLAST","NUCLEOMORPH"
	};
int num = sizeof(list)/sizeof(char *);
static int known[] = {0,1,2,3,3,3,3,3,3,3};
static char organelle[100];

strcpy(organelle, ligne + 12);
p = strtok(organelle, " ");
if(p == NULL) return NULL;
majuscules(organelle);
for(i = 0; i < num; i++) {
	if(strcmp(list[i], organelle) == 0) 
		return list[known[i]];
	}
return NULL;
}


void process_source(char *ligne, size_t lligne)
{
static char name[150], organ[20];
char *p, *q, *d_name, *ascendant;
int numspec, num_org, *extra, do_it = FALSE;

/* new format: organelle at start of SOURCE record */
p = organelle_from_line(ligne); 
if(p != NULL) { 
	strcpy(organ + 2, p); do_it = TRUE;
	}
do 	fgets(ligne, lligne, in_flat);
while(strncmp(ligne, "  ORGANISM", 10) != 0);
p = ligne + strlen(ligne) - 1; *(p--) = 0; /* enlever \n terminal */
/* enlever les espaces terminaux */
trim_key(ligne);
d_name = ligne + 12;
if(! do_it) {
	/* old format: organelle at start of ORGANISM record */
	q = organelle_from_line(ligne);
	if(q != NULL) { 
		strcpy(organ + 2, q); do_it = TRUE;
		d_name = strchr(ligne + 12, ' ');
		while(*d_name == ' ') d_name++;
		}
	}
memcpy(name, d_name, p - d_name + 2);
majuscules(name);
/* traiter lignes de classification */
if( ( numspec = iknum(name, kspec) ) == 0) { /* new species, classify it */
	fgets(ligne, lligne, in_flat);
	ascendant = proc_classif(ligne, lligne);
	numspec = crespecies(ascendant, name);
	}
else	{ /* known species, skip classification data */
	do fgets(ligne, lligne, in_flat);
	while(*ligne == ' ');
	readspec(numspec);
	while(pspec->syno > 0) {
		numspec = pspec->syno;
		readspec(numspec);
		}
	}
fast_add_seq_to_spec(numspec, seq_num);
readloc(loc_num);
ploc->spec = numspec;
writeloc(loc_num);
if(do_it) {
	memcpy(organ, "05", 2);
	num_org = find_key(organ, smj_binary_tree, TRUE, &extra);
	if(num_org == 0) {
		fprintf(log_file, "cannot create %s\n", name);
		non_chargee = TRUE;
		}
	else	{
		add_to_smjlng(num_org, seq_num, extra);
		readloc(loc_num);
		ploc->org = num_org;
		writeloc(loc_num);
		}
	}
}


static char *read_one_rec_type(char *ligne, size_t lligne)
/* lit la suite du meme sous-type 
ok pour GenBank et embl*/
{
static char name[MEM_ONE_REC_TYPE], rec_type[3];
char *d_name, *p;
int skip;

if(genbank) skip = 12;
else 	{
	skip = 5;
	memcpy(rec_type, ligne, 2);
	}
d_name = name;
do	{
	p = ligne + strlen(ligne) - 1; *p = 0; /* enlever \n terminal */
	/* enlever les espaces terminaux */
	while ( *(--p) == ' ') *p = 0; 
	if(p >= ligne + skip) {
		if( d_name - name + (p - ligne - skip + 1) >= sizeof(name) ) {
			fprintf(log_file, 
			"Warning: increase parameter MEM_ONE_REC_TYPE\n");
			non_chargee = TRUE;
			return name;
			}
		memcpy(d_name, ligne + skip, p - ligne - skip + 1);
		d_name += (p - ligne - skip + 1); 
	}
	*d_name = 0;
	fgets(ligne, lligne, in_flat);
	if( ( genbank && ( ligne[0] != ' ' || ligne[2] != ' ' ) ) ||
		( embl && strncmp(ligne, rec_type, 2) != 0 ) ) break;
	*(d_name++) = ' '; /* encore a lire, ajouter espace en fin de deja lu */
	}
while( TRUE );
return name;
}


void process_keywords(char *ligne, size_t lligne)
/* ok pour genbank et embl */
{
static char keyword[100];
char *p, *d_name, *q;
int numkey, err;

p = read_one_rec_type(ligne, lligne);
if( *p == 0 ) return;
do	{
	p = next_word(p, keyword, ";", sizeof(keyword), NULL);
	if(keyword[0] == 0) break;
	d_name = keyword;
	while( *d_name == ' ') d_name++;
	q = d_name + strlen(d_name);
	while ( *(q - 1) == ' ' ) *(--q) = 0;/* enlever les espaces terminaux */
	if( *(q - 1) == '.' ) *(--q) = 0; /* enlever le . terminal */
	if( q <= d_name ) continue;
	numkey = crekeyword(NULL, d_name);
	err = mdshrt(ksub, seq_num, 4, numkey, NULL);
	if(err != 2) fast_add_seq_to_keyw(numkey, seq_num);
	}
while ( p != NULL );
}


int process_features(FILE *in_flat, char *ligne, size_t lligne, 
	int tot_types, int *type_ranks, char **type_codes, char **type_abbrevs, int *type_num,
	int current_div, char *location, char *qualifiers, int max_loc_qual)
{
static char feat_name[20], nmere[L_MNEMO + 1];
int i, type, erreur, fille_num;
char *p;
off_t debut_ft_key, debut_next_entry;

fille_num = seq_num;
readsub(seq_num);
memcpy(nmere, psub->name, L_MNEMO); nmere[L_MNEMO] = 0;
memset(type_num, 0, tot_types * sizeof(int));
while(embl && strncmp(ligne, "FH", 2) == 0) {
	debut_ft_key = ftello(in_flat);
	fgets(ligne, lligne, in_flat);
	}
if(genbank) {
	debut_ft_key = ftello(in_flat);
	fgets(ligne, lligne, in_flat);
	}

while( ligne[0] == ' ' || strncmp(ligne, "FT", 2) == 0 ) {
	erreur = read_flat_loc_qualif(in_flat, ligne, lligne, &debut_next_entry,
		location, max_loc_qual, feat_name, qualifiers, max_loc_qual);
	if(erreur) {
		non_chargee = TRUE;
		fprintf(log_file, "Warning: increase parameter MAX_LOC_QUAL\n");
		fprintf(log_file, "%s %s\n", feat_name, location);
		continue;
		}
/* recherche d'une sous-seq a creer */
	for(type = 0; type < tot_types; type++) 
		if(strcmp(feat_name, type_codes[type] ) == 0) break;
	if(type < tot_types) {
		process_ft_key(type_ranks, type_codes,
				type_abbrevs, type_num, type,
				location, max_loc_qual, 
				nmere, &fille_num, feat_name, debut_ft_key, 
				current_div, qualifiers);
		}
	debut_ft_key = debut_next_entry;
/* for SOURCE: look for CHROMOSOME data */
	if(strcmp(feat_name, "SOURCE") == 0) {
		p = strstr(qualifiers, "/CHROMOSOME=");
		if(p != NULL) {
			p = get_qualif_value(p);
			if(p != NULL) prep_chromosome(p, seq_num);
			}
		}
	else if(strcmp(feat_name, "EXON") == 0 && 
			strstr(qualifiers, "/USEDIN=") != NULL) {
/* si exon tester et si presence de /usedin */
		char *exon_qualifs;
		int l;
	/* dupliquer la chaine de qualifiers */
		l = strlen(qualifiers);
		exon_qualifs = malloc(l + 1);
		if(exon_qualifs != NULL) {
			memcpy(exon_qualifs, qualifiers, l + 1);
			process_exon(exon_qualifs, location, qualifiers, 
				max_loc_qual,
				tot_types, type_codes, type_abbrevs, type_ranks,
				&fille_num);
			free(exon_qualifs);
			}
		else	{
			fprintf(log_file, 
				"Warning: not enough memory to process exon\n");
			}
		}
/* for all features but "-" accrocher le nom de feature en mot-cle */
	if(strcmp(feat_name, "-") != 0) {
		int numkey, err;
		if( strcmp("MISC_FEATURE", feat_name) != 0 )
			numkey = crekeyword("MISC_FEATURE", feat_name);
		else
			numkey = crekeyword(NULL, "MISC_FEATURE");
		err = mdshrt(ksub, seq_num, 4, numkey, NULL);
		if(err != 2) fast_add_seq_to_keyw(numkey, seq_num);
		}
	}
return fille_num;
}


void process_ft_key(int *type_ranks, char **type_codes,
			char **type_abbrevs, int *type_num, int type, 
			char *location, int maxlocat, char *nmere, 
			int *p_fille_num,
			char *feat_name, long info_addr, int div, char *qualif)
{
char name[L_MNEMO + 5], aux[10];

if(type_ranks[type] == 0) {
	non_chargee = TRUE;
	fprintf(log_file, "Error, type missing in file SMJYT: %s\n", type_codes[type]);
	return;
	}
/* fabrication du nom de seq */
sprintf(name, "%s.%s", nmere, type_abbrevs[type]);
compact(name);
++type_num[type];
sprintf(aux, "%d", type_num[type]);
if( (int)strlen(name) + (int)strlen(aux) > L_MNEMO) {
	fprintf(log_file, "Warning: subsequence name too long: %s%s\n", 
		name, aux);
	return;
	}
strcat(name, aux);
mkextr(location, seq_num, p_fille_num, name, type_ranks[type],
		feat_name, info_addr, div, qualif);
}


void process_exon(char *exon_qualifs, char *location, char *qualifiers, 
	int max_loc_qual,
	int tot_types, char **type_codes, char **type_abbrevs, int *type_ranks,
	int *p_fille_num)
{
char *p, *q, *value, access[20], label[20];
long pannot;
int div, mere_label, num, next, type_index, doit, cur_rank, new_rank, erreur;
			
/* tester si presence de /usedin et repeter car possible plusieurs fois */
while( (p = strstr(exon_qualifs, "/USEDIN=")) != NULL) {
	exon_qualifs = p + 1;
	/* essayer de fabriquer la seq composite indiquee par /usedin= */
	value = get_qualif_value(p);
	if(value == NULL) continue;
	if( (q = strchr(value, ':')) == NULL) continue;
	memcpy(access, value, q - value); access[q - value] = 0;
	strcpy(label, q + 1);
	/* recherche de la sequence appellee */
	num = find_label(label, access, &pannot, &div, &mere_label);
	if(num) continue;
	/* si label trouve dans meme table des features pas la peine de traiter 
	 car cela a ete deja fait ou sera fait plus loin */
	 if(mere_label == seq_num) continue;
/* si label trouve est dans la liste des acc # requests, utiliser ce mecanisme */
	for(num = 0; num < tot_request; num++)
		if(acc_request[num].inf == pannot && 
				acc_request[num].div == div) break;
	if(num < tot_request) continue;
	/* calcul du type de cette sequence */
	read_annots(pannot, div);
	num = 0;
	p = pinfo->line + 4; while(*(++p) != ' ') label[num++] = toupper(*p);
	label[num] = 0;
	for(type_index = 0; type_index <tot_types; type_index++)
		if(strcmp(label, type_codes[type_index]) == 0) break;
	if(type_index >= tot_types) continue;
	if(type_ranks[type_index] == 0)continue;
	/* calcul du nom de la sequence a creer */
	readsub(mere_label);
	next = abs(psub->pext);
	doit = TRUE;
	new_rank = 0;
	while(doit && next != 0) {
		readlng(next); next = plng->next;
		for(num = 0; num < SUBINLNG; num++) {
			long faddr; int div2;
			if(plng->sub[num] == 0) break;
	/* check whether the feature entry is already processed */
			seq_to_annots(plng->sub[num], &faddr, &div2);
			if(faddr == pannot && div2 == div) {
				doit = FALSE;
				break;
				}
			readsub(plng->sub[num]);
			p = strchr(psub->name, '.');
			if(p == NULL || p >= psub->name + L_MNEMO) continue;
			if(strncmp(p + 1, type_abbrevs[type_index], 2) != 0) 
				continue;
			cur_rank = -1;
			sscanf(p + 3, "%d", &cur_rank);
			if(cur_rank == -1) continue;
			if(cur_rank > new_rank) new_rank = cur_rank;
			}
		}
	if( ! doit ) continue;
	new_rank++;
	readsub(mere_label);
	memcpy(label, psub->name, L_MNEMO);
	num = L_MNEMO;
	while(label[num - 1] == ' ') num--;
	label[num] = 0;
	p = strchr(label, '.');
	if(p == NULL) {
		label[num++] = '.';
		p = label + num - 1;
		}
	sprintf(p + 1, "%s%d", type_abbrevs[type_index], new_rank);
	/* essayer de fabriquer la sous-seq */
	erreur = read_addr_loc_qualif(pannot, div, location, max_loc_qual, NULL,
				qualifiers, max_loc_qual);
	if(erreur) {
		fprintf(log_file, "Warning: increase parameter MAX_LOC_QUAL\n");
		fprintf(log_file, "%s %s\n", type_codes[type_index], location);
		}
	else	{
		mkextr( location, mere_label, p_fille_num, label,
				type_ranks[type_index], type_codes[type_index], 
				pannot, div, qualifiers);
		}
	}
}


int read_flat_loc_qualif(FILE *in_flat, char *ligne, size_t lligne, 
	off_t *next_feat_addr,
	char *location, int maxlocat, char *type, char *qualif, int maxqualif)
/* lecture de la location courante et de ses qualifiers dans fichier in_flat
ligne, lligne: pour lire dedans les lignes
!!!! A l'APPEL ligne contient la 1ere ligne d'une feature entry !!!!
!!!! AU RETOUR ligne contient la 1ere ligne de l'entry suivante !!!!
*next_feat_addr: rendu avec addr debut entry suivante
location: memoire pour mettre la location dedans
maxlocat: place memoire max disponible
type: (si != NULL) memoire pour mettre le type dedans
qualif: memoire pour metre les qualifiers
maxqualif: place memoire max disponible
valeur rendue: FALSE si ok; TRUE si erreur pas assez de memoire
*/
{
int i, lloc=0, l, erreur, next_feat, lu;
char *ps, *p;

erreur = FALSE; next_feat = FALSE;
maxlocat--;
location[0] = 0;
if(type != NULL) { /* memorisation du type */
	l = 0; ps = ligne + 4; 
	while(*(++ps) != ' ') type[l++] = toupper(*ps); type[l] = 0;
	}
while( (ps = strchr(ligne,'/')) == NULL ) {
	i = strlen(ligne);
	while(ligne[i - 1] =='\n' || ligne[i - 1] == ' ') i--;
	i -= 21;
	if(lloc+i >= maxlocat) {
		erreur = TRUE;
		i = maxlocat - lloc;
		}
	if(i>0) {
		memcpy(location+lloc, ligne + 21, i);
		lloc += i; location[lloc] = 0;
		}
	*next_feat_addr = ftello(in_flat);
	fgets(ligne, lligne, in_flat);
	ps = NULL;
	if(strcmptrail(ligne,20,NULL,0) != 0 && 
		strcmptrail(ligne,20,"FT",2) !=0 ) {
		next_feat = TRUE;
		break;
		}
	}
if(ps!=NULL) {
	i = ps - ligne + 1;
	if( strcmptrail(ligne+2, i-3, NULL, 0) != 0 && i > 22 ) {
		if(lloc+i-22 >= maxlocat) {
			erreur = TRUE;
			i = maxlocat - lloc + 22 - 1;
			}
		if(i > 22) {
			memcpy(location+lloc,ligne+21,i-22);
			lloc += (i-22); location[lloc] = 0;
			}
		}
	}
compact(location);
majuscules(location);

/* lecture des qualifiers */
qualif[0] = 0;
if (next_feat) return erreur;
maxqualif--;
l = 0;
do	{
	lu = strlen(ps);
	while( lu > 0 && ( *(ps + lu - 1) == '\n' || *(ps + lu - 1) == ' ' ) ) 
		lu--;
	if( l + lu + 1>= maxqualif ) {
		erreur = TRUE;
		lu = maxqualif - l - 1;
		}
	if(lu > 1) {
		memcpy(qualif + l, ps, lu);
		l += lu + 1;
		qualif[l - 1] = ' '; qualif[l] = 0;
		}
	do	{
		*next_feat_addr = ftello(in_flat);
		fgets(ligne, lligne, in_flat);
		}
	while(ligne[0] == '\n');
	ps = ligne + 21;
	}
while (strcmptrail(ligne, 21, NULL, 0) == 0 || 
			strcmptrail(ligne, 21, "FT", 2) == 0);
majuscules(qualif);
return erreur;
}


void process_reference(char *ligne, size_t lligne, int numref)
/* for both GenBank and EMBL formats */
{
char *p, *q, *author_lines;

author_lines = NULL;
do	{
	if( ( genbank && strncmp(ligne, "  AUTHORS", 9) == 0 ) ||
		( embl && strncmp(ligne, "RA", 2) == 0 ) ) {
		/* memoriser rubrique AUTHORS convertie en majuscules */
		author_lines = read_author_lines(ligne, lligne);
		}
	else if( ( genbank && strncmp(ligne, "  JOURNAL", 9) == 0 ) ||
		( embl && strncmp(ligne, "RL", 2) == 0 ) ) {
		process_journal(ligne, lligne, author_lines, numref);
		}
	else fgets(ligne, lligne, in_flat);
	}
while( ( genbank && *ligne == ' ') || 
	( embl &&  *ligne == 'R' && strncmp(ligne, "RN", 2) != 0 ) );
}


void process_journal(char *ligne, size_t lligne, char *author_data, int numref)
{
int is_book, is_thesis, number,
	tot_txt, num_bib, saved_addr, l, special_format,
	*extra_journal, is_patent, is_unpublished;
static char aux[50], volume[20], first_page[20], annee[20],
	journal[500], first_author[100], lc_ligne[sizeof(journal)];
char *p, *q, *colon, *name;
void *node;
s_journal_code *data;

is_book = is_thesis = is_patent = is_unpublished = FALSE;
name = read_one_rec_type(ligne, lligne);
/* memoriser etat avec minuscules */
memcpy(lc_ligne, name, sizeof(lc_ligne)); lc_ligne[sizeof(lc_ligne) - 1] = 0;
majuscules(name);
if(strncmp(name, "(IN)", 4) == 0) is_book = TRUE;
else if(strncmp(name, "THESIS", 6) == 0) is_thesis = TRUE;
else if( (genbank && strstr(name, "PATENT:") != NULL ) ||
	(embl && strncmp(name, "PATENT NUMBER", 13) == 0) ) is_patent = TRUE;
else if(strncmp(name, "SUBMITTED", 9) == 0 ||
	strncmp(name, "SUBMISSION", 10) == 0 ||
	strstr(name, "IN PRESS") != NULL ||
	strstr(name, "PERSONAL COMMUNICATION") != NULL ||
	strstr(name, "UNPUBLISHED") != NULL ) 
		is_unpublished = TRUE;

if( is_book ) {
	process_book(name, author_data);
	return;
	}
else if( is_thesis ) {
	process_thesis(name, author_data);
	return;
	}
else if( is_patent ) {
	process_patent(name, author_data);
	return;
	}
else if( is_unpublished ) {
	process_unpublished(name, author_data, numref);
	return;
	}

/* recherche derniere ',' separant journal + volume du reste */
if(genbank) colon = strrchr(name, ','); 
else colon = strchr(name, ':');
if(colon == NULL) return;

/* recherche volume = 1er champ numerique apres au moins un espace depuis name */
p = strchr(name, ' ');
if(p == NULL) return;
do	p++;
while( *p != 0 && strchr("0123456789", *p) == NULL );
if( *p != 0 ) {
	*(p - 1) = 0; /* fin de la partie nom quand volume commence */
	q = p;
	do	q++;
	while(strchr("0123456789", *q) != NULL);
	*q = 0;
	if(q > colon) return; /* JOURNAL mal forme */
	strcpy(volume, p);
	}
else	strcpy(volume, "0");

/* calcul first_page */
p = colon + 1;
while( *p != 0 && *p != '-' && *p != '(' ) p++;
if(*p == 0 || p - 1 - colon >= sizeof(first_page) ) {
	strcpy(first_page, "0");
	}
else	{
	memcpy(first_page, colon + 1, p - 1 - colon); 
	first_page[p - 1 - colon] = 0;
	compact(first_page);
	}

/* calcul annee */
p = strchr(p, '(');
if(p != NULL) q = strchr(p, ')');
if( p == NULL || q == NULL || q - 1 - p >= sizeof(annee) ) {
	strcpy(annee, "0"); /* annee mal formee */
	}
else	{
	memcpy(annee, p + 1, q - 1 - p); annee[q - 1 - p] = 0;
	if(q - 1 - p > 4) compact(annee);
	}

/* calcul du nom du journal */
l = strlen(name);
if(l >= sizeof(journal)) {
	l = sizeof(journal) - 1;
	name[l] = 0;
	}
strcpy(journal, name);
lc_ligne[l] = 0;
q = journal - 1;
while( *(++q) != 0) {
	if ( *q == '.' || *q == ',' ) *q = ' ';
	}
compact(journal); journal[lrtxt] = 0;
/* recherche du libel du journal */
if( find_g_dynlist(journal, dynlist_journals, TRUE, &node) == NULL) {
	fprintf(log_file, "cannot create new journal %s\n", name);
	non_chargee = TRUE;
	return;
	}
data = dynlist_get_extra(node);
if( data != NULL ) {
	extra_journal = get_smj_extra(data->number);
	}
else	{ /* creation nouveau journal */
	data = (s_journal_code *)malloc(sizeof(s_journal_code));
	if(data == NULL) {
		fprintf(log_file, "cannot create new journal %s\n", name);
		non_chargee = TRUE;
		return;
		}
	l = strlen(journal); if(l > 18) l = 18;
	data->code = malloc(l + 1);
	if( data->code == NULL) {
		fprintf(log_file, "not enough memory\n");
		non_chargee = TRUE;
		return;
		}
	memcpy(data->code, journal, l); 
	data->code[l] = 0;
	sprintf(aux, "02%s", data->code);
	number = find_key(aux, smj_binary_tree, TRUE, &extra_journal);
	if(number == 0) {
		fprintf(log_file, "cannot create %s\n", aux);
		non_chargee = TRUE;
		return;
		}
	data->number = number;
	dynlist_set_extra(node, data);
	tot_txt = read_first_rec(ktxt, NULL) + 1;
	l = strlen(lc_ligne);
	memset(ptxt, ' ', lrtxt);
	memcpy(ptxt, lc_ligne, (l > lrtxt ? lrtxt : l) );
	writetxt(tot_txt);
	write_first_rec(ktxt, tot_txt, 0);
	readsmj(number);
	psmj->libel = tot_txt;
	writesmj(number);
	}
/* creation de la reference */
special_format = ( strcmp(volume, "0") == 0 || strcmp(first_page, "0") == 0 );
/* cas JOURNAL 0:0-0(annee). */
if( special_format ) {
/* ici fabriquer JOURNAL/ANNEE/AUTEUR1 */
	if( author_data == NULL ) return; /* refer sans RA */
	next_auth_name(author_data, first_author, sizeof(first_author), TRUE);
	if( *first_author == 0 ) return; /* refer sans auteur */
	format_auteur(first_author);
	sprintf(name, "%s/%s/%s", data->code, annee, first_author);
	set_reference(name, NULL, data->number, extra_journal, annee, 
		first_author, TRUE);
	}
else	{
	sprintf(name, "%s/%s/%s", data->code, volume, first_page);
	set_reference(name, NULL, data->number, extra_journal, annee, 
		author_data, FALSE);
	}
}


static char *read_author_lines(char *ligne, size_t lligne)
/* lit et memorise toutes les lignes auteur et ajoute "," en fin de chaque ligne si absente
ok pour GenBank et embl*/
{
static char name[MEM_ONE_REC_TYPE], rec_type[3];
char *d_name, *p;
int skip;

if(genbank) skip = 11; /* et pas 12 pour que "and" en debut de ligne marche bien ! */
else 	{
	skip = 4;
	memcpy(rec_type, ligne, 2);
	}
d_name = name;
do	{
	p = ligne + strlen(ligne) - 1; *p = 0; /* enlever \n terminal */
	while( p > ligne && *(--p) == ' ' ) *p = 0;/*enlever espaces terminaux*/
	if( d_name - name + (p - ligne - skip + 1) >= sizeof(name) - 2 ) {
		fprintf(log_file, 
			"Warning: increase parameter MEM_ONE_REC_TYPE\n");
		non_chargee = TRUE;
		return name;
		}
	memcpy(d_name, ligne + skip, p - ligne - skip + 1);
	d_name += (p - ligne - skip + 1); *d_name = 0;
	fgets(ligne, lligne, in_flat);
	if( ( genbank && ( ligne[0] != ' ' || ligne[2] != ' ' ) ) ||
		( embl && strncmp(ligne, rec_type, 2) != 0 ) ) break;
	/* encore a lire, ajouter ' ,' en fin de deja lu */
	if( *(d_name - 1) != ',') {
		memcpy(d_name, " ,", 2);
		d_name += 2;
		} 
	}
while( TRUE );
majuscules(name);
return name;
}


void process_authors(char *author_data, int num_bib, int first_aut_only)
/* for both GenBank and EMBL formats 
   if first_aut_only is TRUE, author_data = name of 1st author
   else, author_data = text of author line(s)
*/
{
int num_aut, tot_old, num;
static char auteur[100];
static int *old_auteurs = NULL, max_old = 0;

/* traitement auteurs */
if(author_data == NULL) return; /* reference sans ligne AUTHORS */
/* lecture des auteurs existants de cette reference */
readbib(num_bib); num = pbib->plaut; tot_old = 0;
while(num != 0) {
	readshrt(num); num = pshrt->next;
	if(tot_old >= max_old) {
		int *tmp;
		max_old += 50;
		tmp = (int *)malloc(max_old * sizeof(int));
		if(tmp == NULL) {
			if(old_auteurs != NULL) free(old_auteurs); 
			max_old = tot_old = 0;
			break;
			}
		if(tot_old > 0) memcpy(tmp, old_auteurs, tot_old * sizeof(int));
		if(old_auteurs != NULL) free(old_auteurs);
		old_auteurs = tmp;
		}
	old_auteurs[tot_old++] = pshrt->val;
	}
do	{
	if(first_aut_only)
		strcpy(auteur, author_data);
	else	{
		author_data = next_auth_name(author_data, auteur, 
			sizeof(auteur), FALSE);
		if(auteur[0] == 0) continue;
		format_auteur(auteur);
		}
	num_aut = find_key(auteur, aut_binary_tree, TRUE, NULL);
	if(num_aut == 0) {
		fprintf(log_file, "cannot create author %s\n", auteur);
		non_chargee = TRUE;
		}
	else	{
/* auteur num_aut doit etre associe */
		for(num = 0; num < tot_old; num++)
			if(old_auteurs[num] == num_aut) break;
		if(num < tot_old) 
/* cet auteur est deja la */
			old_auteurs[num] = 0;
		else	{
/* cet auteur doit etre ajoute */
			mdshrt(kaut, num_aut, 1, num_bib, NULL);
			mdshrt(kbib, num_bib, 2, num_aut, NULL);
			}
		}
	}
while( (!first_aut_only) && author_data != NULL);
/* enlever les anciens auteurs qui ne sont plus cites */
for(num = 0; num < tot_old; num++) {
	if( (num_aut = old_auteurs[num]) == 0) continue;
	mdshrt(kaut, num_aut, -1, num_bib, NULL);
	mdshrt(kbib, num_bib, -2, num_aut, NULL);
	}
}


char *next_auth_name(char *auth_data, char *auth_name, int max_name,
	int first_aut_only)
/* for both 
GenBank   "Nom, I., Nom, I. and Nom. I."  
and EMBL "Nom I., Nom J.;" formats
search next author name from auth_data, and puts it in auth_name 
	(with max_name bound)
returns NULL if auth_data has been completely scanned,
	or pointer for next scan
after return, auth_name can be filled, or auth_name[0]=0 can be ignored
if first_aut_only == TRUE, AND is processed differently
*/
{
char *virg, *fin, *p, *old_fin;
int l;

auth_name[0] = 0;

encore:
while(*auth_data == ' ' || *auth_data == ',') auth_data++;
if(*auth_data == 0) return NULL;
virg = strchr(auth_data, ',');
if(virg == NULL) virg = strchr(auth_data, ';');
if(virg == NULL) virg = auth_data + strlen(auth_data);
fin = virg - 1;
if(first_aut_only) {
	p = strstr(auth_data, " AND ");
	if(p != NULL && p + 3 <= fin) fin = p - 1;
	}
/* on traite [auth_data - fin] */
p = strchr(auth_data, '.');
if(p != NULL && p < fin) {
/* ignorer les III 3RD II 2ND ET finaux apres un . */
	if(strncmp(fin - 3, " III", 4) == 0) fin -= 3;
	else if(strncmp(fin - 3, " 3RD", 4) == 0) fin -= 3;
	else if(strncmp(fin - 2, " II", 3) == 0) fin -= 2;
	else if(strncmp(fin - 3, " 2ND", 4) == 0) fin -= 3;
	else if(strncmp(fin - 2, " ET", 3) == 0) fin -= 2;
	}
/* ignorer tous les mots termines par . */
while(fin >= auth_data) {
	while(fin >= auth_data && *fin == ' ') fin--;
	if( fin >= auth_data && ( *fin == '.' || *(fin - 1) == '.') ) {
		fin--;
		if( first_aut_only ) old_fin = fin;
		while(fin >= auth_data && *fin != ' ') fin--;
		if( first_aut_only && fin < auth_data ) fin = old_fin;
		}
	else break;
	}
if( fin < auth_data ) return (*virg == 0 ? NULL : virg + 1);
/* name found
ignore part before "AND" for GenBank format
*/
p = strstr(auth_data - 1, " AND "); /* le " - 1" est important */
if( p != NULL && p + 3 <= fin) {
	return p + 5; /* retraiter apres le AND */
	}
/* recommencer apres le . s'il est en position 1 ou 2 */
if(*auth_data == '.') {
	auth_data++;
	goto encore;
	}
if(*(auth_data + 1) == '.') {
	auth_data += 2;
	goto encore;
	}
l = fin - auth_data + 1;
if( l >= max_name) l = max_name - 1;
memcpy(auth_name, auth_data, l); auth_name[l] = 0;
return (*virg == 0 ? NULL : virg + 1);
}


void format_auteur(char *auteur)
{
char *q;
int l;
const int width = sizeof(paut->name);

l = strlen(auteur);
if( l > width ) { /* tronquer a width */
	l = width; auteur[l] = 0;
	}
/* suppression des espaces terminaux */
q = auteur + l;	
while( q > auteur && *(--q) == ' ' )*q  = 0;
}


void process_book(char *journal_data, char *author_data)
{
char *p;
static char annee[5], auteur[50], lu[100];

if( author_data == NULL ) return; /* refer sans AUTHORS */
p = strstr(journal_data, "(19");
if(p == NULL) p = strstr(journal_data, "(20");
if(p == NULL) strcpy(annee, "0");
else	{
	memcpy(annee, p + 1, 4); annee[4] = 0;
	}
next_auth_name(author_data, auteur, sizeof(auteur), TRUE);
if( *auteur == 0 ) return; /* refer sans auteur */
format_auteur(auteur);
sprintf(lu, "BOOK/%s/%s", annee, auteur);
set_reference(lu, "02BOOK", 0, NULL, annee, auteur, TRUE);
}


void process_thesis(char *journal_data, char *author_data)
{
char *p;
static char annee[5], auteur[50], lu[100];

if( author_data == NULL ) return; /* refer sans RA */
p = strstr(journal_data, "(19");
if(p == NULL) p = strstr(journal_data, "(20");
if( p == NULL) strcpy(annee, "0");	
else 	{
	memcpy(annee, p + 1, 4); annee[4] = 0;
	}
next_auth_name(author_data, auteur, sizeof(auteur), TRUE);
if( *auteur == 0 ) return; /* refer sans auteur */
format_auteur(auteur);
sprintf(lu, "THESIS/%s/%s", annee, auteur);
set_reference(lu, "02THESIS", 0, NULL, annee, auteur, TRUE);
}


void process_unpublished(char *journal_data, char *author_data, int numref)
{
char *p;
static char annee[5], auteur[50], lu[100];

if(author_data == NULL) return; /* reference sans AUTHORS */
if(strstr(journal_data, "SUBMITTED") != NULL ||
	strstr(journal_data, "SUBMISSION") != NULL ) {
	p = strstr(journal_data, "-19");
	if(p == NULL) p = strstr(journal_data, "-20");
	if(p == NULL || *(p + 5) != ')') strcpy(annee, "0");
	else	memcpy(annee, p + 1, 4); annee[4] = 0;
	}
else	{
	p = strstr(journal_data, "(19");
	if(p == NULL) p = strstr(journal_data, "(20");
	if(p == NULL || *(p + 5) != ')') strcpy(annee, "0");
	else	memcpy(annee, p + 1, 4); annee[4] = 0;
	}
if(strcmp(annee, "0") == 0 && numref >= 2) return;

next_auth_name(author_data, auteur, sizeof(auteur), TRUE);
if( *auteur == 0 ) return; /* refer sans auteur */
format_auteur(auteur);
sprintf(lu, "UNPUBL/%s/%s", annee, auteur);
set_reference(lu, "02UNPUBL", 0, NULL, annee, auteur, TRUE);
}


void process_patent(char *journal_data, char *author_data)
{
char *p, *q;
static char annee[5], reference[200], number[100];

/* recherche du numero */
if(genbank) {
	p = strstr(journal_data, "PATENT:") + 7;
	q = p;
	do	q++;
	while(*q == ' ');
	q = strchr(q, ' ');
	do	q++;
	while(*q == ' ');
	q = strchr(q, ' ');
	}
else	{ /* embl */
	p = journal_data + 13;
	while(*p == ' ') p++;
	q = p;
	while( *q != '/' && *q != ','  && *q != 0 ) q++;
	}
if( q - p >= sizeof(number) ) q = p + sizeof(number) - 1;
memcpy(number, p, q - p ); number[ q - p ] = 0;
/* recherche de l'annee */
p = strstr(q, "-19");
if(p == NULL) p = strstr(q, "-20");
if( p != NULL) {
	memcpy(annee, p + 1, 4); annee[4] = 0;
	}
else	strcpy(annee, "0");
sprintf(reference, "PATENT/%s", number);
set_reference(reference, "02PATENT", 0, NULL, annee, author_data, FALSE);
}


void set_reference(char *reference, 
	char *code_j, int num_j, int *extra_journal, 
	char *annee, char *author_data, int first_aut_only)
/* if first_aut_only is TRUE, author_data = name of 1st author
   else, author_data = text of author line(s)
*/
{
int num_bib, num_year, *extra_year, deja_connu, *extra_bib;
const int maxref = sizeof(pbib->name);

compact(reference);
if( (int) strlen(reference) > maxref) reference[maxref] = 0;
num_bib = find_key(reference, bib_binary_tree, TRUE, &extra_bib); 
if(num_bib == 0) {
	fprintf(log_file, "cannot create %s\n", reference);
	non_chargee = TRUE;
	return;
	}
if(num_j == 0) {
	num_j = find_key(code_j, smj_binary_tree, TRUE, &extra_journal); 
	if(num_j == 0) {
		fprintf(log_file, "cannot create %s\n", code_j);
		non_chargee = TRUE;
		return;
		}
	}
num_year = process_year_bib(annee, &extra_year);
readbib(num_bib);
/* si refer a au moins une seq associee, ne pas refaire les auteurs */
deja_connu = ( pbib->plsub != 0 );
connect_locus_bib(loc_num, seq_num, num_bib, extra_bib);
readbib(num_bib);
pbib->j = num_j;
/* ne pas changer l'annee d'une reference deja avec des sequences */
if( ( pbib->y == 0 ) || !deja_connu )
	pbib->y = num_year;
else if( num_year != pbib->y ) {
	num_year = pbib->y;
	extra_year = get_smj_extra(num_year);
	}
writebib(num_bib);
add_to_smjlng(num_j, seq_num, extra_journal);
if(num_year != 0) add_to_smjlng(num_year, seq_num, extra_year);
/* traitement auteurs */
if(! deja_connu) process_authors(author_data, num_bib, first_aut_only);
}


void *load_j_code_libel(void)
{
int tot_smj, num, l;
char *p, libel[lrtxt + 1];
s_journal_code *data;
void *dynlist, *node;

dynlist = init_dynlist();
tot_smj = read_first_rec(ksmj, NULL);
for(num = 2; num <= tot_smj; num++) {
	readsmj(num);
	if(strncmp(psmj->name, "02", 2) != 0) continue;
	if(     strncmp(psmj->name, "02BOOK ", 7) == 0 ||
		strncmp(psmj->name, "02THESIS ", 9) == 0 ||
		strncmp(psmj->name, "02UNPUBL ", 9) == 0 ||
		strncmp(psmj->name, "02PATENT ", 9) == 0 ) continue;
	data = (s_journal_code *)malloc(sizeof(s_journal_code));
	if(data == NULL) goto memoire;
	data->code = (char *)malloc(19);
	if( data->code == NULL ) goto memoire;
	memcpy( data->code, psmj->name + 2, 18);
	data->code[18] = 0;
	trim_key(data->code);
	data->number = num;
	readtxt(psmj->libel);
	memcpy(libel, ptxt, lrtxt); libel[lrtxt] = 0;
	p = libel - 1;
	while( *(++p) != 0 ) {
		if(*p == '.' || *p == ',') *p = ' ';
		else	*p = toupper(*p);
		}
	compact(libel);	
	if( find_g_dynlist(libel, dynlist, TRUE, &node) == NULL) goto memoire;
	dynlist_set_extra(node, data);
	}
return dynlist;

memoire:
fprintf(log_file, "ERROR: not enough memory for load_j_code_libel\n");
exit(ERREUR);
}


void remove_current_seq(int mere, int locus)
/* 
efface la seq de rang mere ds SUBSEQ et locus ds LOCUS apres tentative de 
chargement
*/
{
int point, valeur, num_smj, *extra, journal, annee,
	num, fille, lngrec;
struct rlng rlng2;

/* traitement LOCUS */
readloc(locus);
/* traitement espece */
if( ploc->spec != 0 ) {
	if(swissprot || nbrf) {
		point = ploc->spec;
		while(point != 0) {
			readshrt(point); point = pshrt->next;
			mdlng(kspec, pshrt->val, -2, mere, NULL);
			extra = get_spec_extra(pshrt->val);
			if(extra != NULL) *extra = 0;
			}
		}
	else	{
		mdlng(kspec, ploc->spec, -2, mere, NULL);
		extra = get_spec_extra(ploc->spec);
		if(extra != NULL) *extra = 0;
		}
	}
/* traitement hote */
if( ploc->host != 0 ) {
	mdlng(kspec, ploc->host, -6, mere, NULL);
	}
point = ploc->plref; /* traitement references */
while( point != 0 ) {
	readshrt(point);
	valeur = pshrt->val;
	point = pshrt->next;
	mdshrt(kbib, valeur, -1, mere, NULL);
	extra = get_bib_extra(valeur);
	if( extra != NULL ) *extra = 0;
	readbib(valeur);
	journal = pbib->j;
	annee = pbib->y;
	if(journal != 0) {
		mdlng(ksmj, journal, -1, mere, NULL);
		extra = get_smj_extra(journal);
		if( extra != NULL ) *extra = 0;
		}
	if(annee != 0) {
		mdlng(ksmj, annee, -1, mere, NULL);
		extra = get_smj_extra(annee);
		if( extra != NULL ) *extra = 0;
		}
	}
/* traitement molecule */
if(ploc->molec != 0) {
	mdlng(ksmj, ploc->molec, -1, mere, NULL);
	extra = get_smj_extra(ploc->molec);
	if( extra != NULL ) *extra = 0;
	}
point = ploc->placc; /* traitement accession */
while( point != 0 ) {
	readshrt(point);
	point = pshrt->next;
	lngrec = pshrt->val;
	lngrec = mdshrt(kacc, lngrec, -1, mere, NULL);
	}
/* traitement status */
if(ploc->stat != 0) {
	mdlng(ksmj, ploc->stat, -1, mere, NULL);
	extra = get_smj_extra(ploc->stat);
	if( extra != NULL ) *extra = 0;
	}
/* traitement organelle */
if(ploc->org != 0) {
	mdlng(ksmj, ploc->org, -1, mere, NULL);
	extra = get_smj_extra(ploc->org);
	if( extra != NULL ) *extra = 0;
	}

/* traitement des filles */
readsub(mere);
lngrec = abs(psub->pext);
while(lngrec != 0) {
	dir_read(klng, lngrec, 1, &rlng2); lngrec = rlng2.next;
	for(num = 0; num < SUBINLNG; num++) {
		if( (fille = rlng2.sub[num]) == 0) break;
		remove_subseq_record(fille);
		}
	}

/* traitement SUBSEQ */
remove_subseq_record(mere);

/* effacement du record de LOCUS */
memset(ploc, 0, lrloc); writeloc(locus);
}


void remove_subseq_record(int num)
{
int point, *extra, pre, extrec, val;
char *p, autre_mere[L_MNEMO + 1];

readsub(num);
if(psub->length == 0) return;
point = psub->plkey; /* traitement mots-cles */
while( point != 0 ) {
	readshrt(point);
	val = pshrt->val;
	point = pshrt->next;
	mdlng(kkey, val, -2, num, NULL);
	extra = get_keyw_extra(val);
	if(extra != NULL) *extra = 0;
	}
/*traitement type */
if(psub->type != 0) {
	mdlng(ksmj, psub->type, -1, num, NULL);
	extra = get_smj_extra(psub->type);
	if( extra != NULL ) *extra = 0;
	}
/* destroy mere --> subseq for mothers of this subseq */
extrec = psub->pext;
if(extrec > 0) {
	p = strchr(psub->name, '.');
	if(p != NULL) {
		*p = 0;
		pre = isenum(psub->name);
		mdlng(ksub, pre, -3, num, NULL);
		}
	else	pre = 0;
	}
while(extrec > 0) {
	readext(extrec); extrec = pext->next;
	if(pext->mere != pre) {
		mdlng(ksub, pext->mere, -3, num, NULL);
		pre = pext->mere;
		}
	}
suphsh(num, ksub);
/* effacement du record de SUBSEQ */
memset(psub, 0, lrsub); memset(psub->name, 'x', L_MNEMO); writesub(num);
}


void process_origin(char *ligne)
{
char *p;

majuscules(ligne);
if( (p = strstr(ligne, "CHROMOSOME ")) == NULL) return;
prep_chromosome(p + 11, seq_num);
}


void prep_chromosome(char *p, int seq)
{
#define KEY_WIDTH sizeof(pkey->name)
char *q;
int num, erreur;
static char name[KEY_WIDTH + 1] = "CHROMOSOME ";

while(*p == ' ') p++;
if(strstr(p, "DEGREE") != NULL) return;
if(strncmp(p, "CHROMOSOME", 10) == 0) p += 10;
while(*p == ' ') p++;
/* ignorer apres tous ces separateurs: ' ,-;/().:' */
q = p;
while( *q != 0 ) {
	if (strchr(" ,-;/().:\n", *q) != NULL) {
		*q = 0;
		break;
		}
	q++;
	}
if( *p == 0 ) return;
/* prendre tous les mots "chiffres romains" en chiffres et XYWZ */
if(strchr("IVXYWZ123456789", *p) == NULL) return;
/* se limiter jusqu'a P et Q a cause de 8p3 ou Xq12 */
q = strchr(p, 'P');
if(q == NULL) q = strchr(p, 'Q');
if(q != NULL) *q = 0;
/* ignorer si KB ou BP dans nom du chromosome */
if(strstr(p, "KB") != NULL || strstr(p, "BP") != NULL)
		return;
if(strlen(p) > KEY_WIDTH - 11) *(p + KEY_WIDTH - 11) = 0;
strcpy(name + 11, p);
num = crekeyword("CHROMOSOMES", name);
erreur = mdshrt(ksub, seq, 4, num, NULL);
if(erreur != 2) fast_add_seq_to_keyw(num, seq);
#undef KEY_WIDTH
}




