#include "dir_acnuc.h"
#include <ctype.h>
#include <time.h>

#define MAX_JOURNALS 5000
#define MEM_ONE_REC_TYPE 10000

/* global variables */
FILE *in_flat, *log_file;
int seq_num, loc_num, last_loc, last_sub, non_chargee;
void *acc_binary_tree, *smj_binary_tree, *bib_binary_tree,
	*aut_binary_tree;
int tot_journals;
char *j_code[MAX_JOURNALS], *j_libel[MAX_JOURNALS];
int j_number[MAX_JOURNALS];

/* local functions */
int process_id(char *ligne, char *name);
void process_sq(char *ligne);
void process_ac(char *ligne, size_t taille);
void process_dr(char *ligne, size_t lligne);
void process_dt(char *ligne, size_t lligne);
void process_gn(char *ligne, size_t lligne);
void process_os(char *ligne, size_t lligne);
void process_oc(char *ligne, size_t lligne);
static int is_desc(char *pere, int numspec);
char *check_spec_with_sp(char *name);
void process_og(char *ligne, size_t lligne);
void process_og_plasmid(char *ligne, size_t lligne);
void process_kw(char *ligne, size_t lligne);
void process_de(char *ligne, size_t lligne);
static void process_de_part(char *name, char *p);
static char *process_terminal_parenth(char *name, char *p);
char *remove_terminal_parenth(char *name, char *p);
void process_cc(char *ligne, size_t lligne);
void process_ft(char *ligne, size_t lligne);
void process_prodom(char *ligne);
void process_subcell_loc(char *ligne, size_t lligne);
void process_gene_family(char *ligne, size_t lligne);
void process_r_lines(char *ligne, size_t lligne);
void process_rl(char *ligne, size_t lligne, off_t addr_ra);
void process_ra(off_t addr_ra, int num_bib, int first_aut_only);
void remove_initials(char *auteur);
void process_book(char *ligne, size_t lligne, off_t addr_ra);
void process_thesis(char *ligne, size_t lligne, off_t addr_ra);
void process_submitted(char *ligne, size_t lligne, off_t addr_ra);
void process_patent(char *ligne, size_t lligne, off_t addr_ra);
int load_j_code_libel(void);
int remove_current_seq(int use_address, char *err_file_name, int current_div);


/* external prototypes */
char *proc_classif(char *ligne, size_t lligne);
int find_key(char *name, void *tree, int create_if_not_found,
	int **p_aux_data);
void *load_index_bt(DIR_FILE *k, int char_width, int use_aux_data);
int poffset(int div, int offset);
char *next_word(char *ligne, char *word, char *separs, int max_word, 
	int *erreur);
void add_liste_meres(int seqnum);
void close_if_open(void);
void proc_mmap_options(int argc, char **argv);
int *get_smj_extra(int smjnum);
int *get_keyw_extra(int kwnum);
int *get_spec_extra(int specnum);
int *get_bib_extra(int bibnum);
void fast_add_seq_to_spec(int specnum, int seqnum);
void fast_add_seq_to_keyw(int kwnum, int seqnum);
void add_to_smjlng( int smjnum, int seqnum, int *extra);
char *nextpar(char *debut);
int trim_key(char *name);
int process_year_bib(char *annee, int **extra);
void connect_locus_bib( int locnum, int seqnum, int bibnum, int *extra);
int next_start_seq(FILE *addr_file, char *gcgacnuc_value, char *ligne, 
	size_t lligne, off_t *p_addr, int *p_div);
void cre_new_division(char *name);
void update_div_size(int numdiv);
void *init_acchash(void);
int find_cre_acc(char *name, int create);
#ifdef unix
int check_term(void);
void write_quick_meres(void);
#endif



int main(int argc, char *argv[])
{
char flat_name[100], ligne[100], division_name[21], argument[200],
	gcgacnuc_value[200], err_file_name[200], seq_name[L_MNEMO + 1];
int current_div, tot_read, tot_charge, num_div_smj, 
	div_size, tot_smj, use_address;
time_t heure_debut, heure_fin;
FILE *addr_file;
off_t addr_seq, addr_info;

#ifdef unix
if(argc < 3 || ( argv[1][0] != 'a' && argv[1][0] != 'd' ) ) {
	fprintf(stderr, "Usage: swgener  a  name-of-address-file [-mmap name]\n"
		"or\n"
		        "       swgener  d  name-of-division  [-mmap name]\n");
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
fprintf(log_file, "Program started at %s\n",asctime(localtime(&heure_debut)) );
dir_acnucopen("WP");
#ifdef unix
proc_mmap_options(argc, argv);
#endif
if(!swissprot) {
	fprintf(log_file, "ERROR: Should open a swissprot data base\n");	
	exit(ERREUR);
	}
if(!flat_format) { 
	fprintf(log_file, "ERROR: Cannot process a gcg data base\n");	
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
	strcpy(division_name, "06FLT");
	strcat(division_name, argument);

	/* ouverture du fichier de la division a charger */
	strcpy(flat_name, gcgacnuc_value);
	strcat(flat_name, division_name + 5);
	strcat(flat_name, ".dat");
	in_flat = fopen(flat_name, "r");
	if(in_flat == NULL) {
		fprintf(log_file, "ERROR: cannot open flat file %s\n", 
			flat_name);	
		exit(ERREUR);
		}		
	/* recherche ou creation de la division a charger */
	for(current_div = 0; current_div <= divisions; current_div++) {
		if(strcmp(gcgname[current_div], division_name + 5) == 0) break;
		}
	if(current_div > divisions) {
		if(divisions >= 0) update_div_size(divisions);
		cre_new_division(division_name + 5);
		fprintf(log_file, "New division created: %s\n", 
			division_name + 5);
		current_div = divisions;
		}
	else if(current_div == divisions) update_div_size(divisions);
	dir_acnucflush();
	}

/* chargement des arbres binaires */
acc_binary_tree = NULL; /* this characterizes the new procedure: hash for acc */
if(init_acchash() == NULL) {
	fprintf(log_file, 
		"ERROR: not enough memory to hash accession numbers\n");	
	exit(ERREUR);
	}

smj_binary_tree = load_index_bt(ksmj, sizeof(psmj->name), TRUE);
if(smj_binary_tree == NULL) {
	fprintf(log_file, "ERROR: not enough memory to load SMJYT\n");	
	exit(ERREUR);
	}
tot_journals = load_j_code_libel();
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

last_sub = read_first_rec(ksub, NULL);
last_loc = read_first_rec(kloc, NULL);
tot_read = tot_charge = 0;
time(&heure_debut);
fprintf(log_file, "Sequence loading started at %s\n", asctime(localtime(&heure_debut)) );

#ifdef unix
/* preparer possibilite d'interruption propre */
check_term();
#endif
atexit(close_if_open);

while( TRUE ) {
	if(use_address) {
		if( ! next_start_seq(addr_file, gcgacnuc_value, ligne, 
			sizeof(ligne), &addr_info, &current_div ) ) break;
		}
	else	{
		addr_info = ftello(in_flat);
		if( fgets(ligne, sizeof(ligne), in_flat) == NULL ) break;
		if( strncmp(ligne, "ID", 2) != 0 ) continue;
		}

	non_chargee = FALSE;
	seq_num = process_id(ligne, seq_name);
	if(use_address || seq_num != 0) 
		fprintf(log_file, "---->%s\n", seq_name); fflush(log_file);
	if(seq_num == 0) {
		if(use_address) fprintf(log_file, "already in database\n");
		continue;
		}
	tot_read++;
	fgets(ligne, sizeof(ligne), in_flat);
	do	{
		if(strncmp(ligne, "AC ", 3) == 0) {
			process_ac(ligne, sizeof(ligne));
			}
		else if(strncmp(ligne, "DT ", 3) == 0) {
			process_dt(ligne, sizeof(ligne));
			}
		else if(strncmp(ligne, "DE ", 3) == 0) {
			process_de(ligne, sizeof(ligne));
			}
		else if(strncmp(ligne, "GN ", 3) == 0) {
			process_gn(ligne, sizeof(ligne));
			}
		else if(strncmp(ligne, "OS ", 3) == 0) {
			process_os(ligne, sizeof(ligne));
			}
		else if(strncmp(ligne, "OG ", 3) == 0) {
			process_og(ligne, sizeof(ligne));
			}
		else if(strncmp(ligne, "OC ", 3) == 0) {
			process_oc(ligne, sizeof(ligne));
			}
		else if(strncmp(ligne, "RN ", 3) == 0) {
			process_r_lines(ligne, sizeof(ligne));
			}
		else if(strncmp(ligne, "CC ", 3) == 0) {
			process_cc(ligne, sizeof(ligne));
			}
		else if(strncmp(ligne, "DR ", 3) == 0) {
			process_dr(ligne, sizeof(ligne));
			}
		else if(strncmp(ligne, "KW ", 3) == 0) {
			process_kw(ligne, sizeof(ligne));
			}
		else if(strncmp(ligne, "FT ", 3) == 0) {
			process_ft(ligne, sizeof(ligne));
			}
		else 	{
			fgets(ligne, sizeof(ligne), in_flat);
			}
		}
	while( (!non_chargee) && strncmp(ligne, "SQ ", 3) != 0);
	if( non_chargee ) {
/* ici defaire liens faits pour seq non chargee */
		fprintf(log_file, "this sequence has not been inserted\n");
		if( remove_current_seq(use_address, err_file_name, current_div) ) break;
		}
	else	{
		process_sq(ligne);
		addr_seq = ftello(in_flat);
		readloc(loc_num);
		if(!big_annots) {
			addr_info = poffset(current_div, (unsigned)addr_info);
			addr_seq = poffset(current_div, (unsigned)addr_seq);
			}
		else
			ploc->div = current_div;
		ploc->pinf = (unsigned)addr_info;
		ploc->pnuc = (unsigned)addr_seq;
		writeloc(loc_num);
		last_sub = seq_num;
		last_loc = loc_num;
		addhsh( last_sub, ksub);
		add_liste_meres(seq_num);
		write_first_rec(ksub, last_sub, 0);
		write_first_rec(kloc, last_loc, 0);
		tot_charge++;
		}
	fflush(log_file);
#ifdef unix
	if( check_term() ) { /* tester si interruption demandee */
		fprintf(log_file, "Clean program interruption\n");
		break; 
		}
#endif
	} /* end of loop over all seqs */
dir_acnucclose();
if(use_address) fclose(addr_file);
time(&heure_fin);
#ifdef unix
if(tot_charge > 0) write_quick_meres();
#endif	
fprintf(log_file, "Program finished at %s\n", asctime(localtime(&heure_fin)) );
fprintf(log_file, "lues=%d  chargees=%d  difference=%d seqs/second=%.2f\n", 
	tot_read, tot_charge, tot_read - tot_charge, 
	tot_read / (float)(heure_fin - heure_debut) );
return 0;
}  /* end of main */


int process_id(char *ligne, char *name)
{
int num, longueur, numstat, *extra;
char *p;
static char status[20], molecule[20];

p = next_word(ligne + 5, name, "; ", L_MNEMO + 1, NULL);
p = next_word(p, status + 2, "; ", sizeof(status) - 2, NULL);
p = next_word(p, molecule, "; ", sizeof(molecule) - 2, NULL);
sscanf(p, "%d", &longueur);
num = isenum(name);
if(num != 0) {
	return 0;
	}
seq_num = num = last_sub + 1;
loc_num = last_loc + 1;
memset(psub, 0, lrsub);
padtosize(psub->name, name, L_MNEMO);
psub->length = longueur;
psub->plinf = loc_num;
writesub(num);
memcpy(status, "00", 2);
numstat = find_key(status, smj_binary_tree, TRUE, &extra);
if(numstat == 0) {
	fprintf(log_file, "cannot create %s\n", status);
	return 0;
	}
add_to_smjlng( numstat, seq_num, extra);
memset(ploc, 0, lrloc);
ploc->sub = num;
ploc->stat = numstat;
writeloc(loc_num);
return num;
}


void process_sq(char *ligne)
{
char *p;
int vlength;

if(strncmp(ligne, "SQ ", 3) != 0) return;
p = strchr(ligne + 5, ' ');
if(p == NULL) return;
sscanf(p, "%d", &vlength);
readsub(seq_num);
if(psub->length != vlength) {
	fprintf(log_file, "Warning: lengths in ID and SQ lines differ\n");
	psub->length = vlength;
	writesub(seq_num);
	}
}


void process_ac(char *ligne, size_t lligne)
{
char accession[50], *p;
int numacc;

do	{
	p = ligne + 5;
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
			}
		}
	while(p != NULL);
	fgets(ligne, lligne, in_flat);
	}
while(strncmp(ligne, "AC ", 3) == 0);
}


void process_dr(char *ligne, size_t lligne)
/* traiter les lignes DR   EMBL; xxxxx;
*/
{
char accession[50], *p;
int numacc = 0;

if(strncmp(ligne + 5, "EMBL;", 5) == 0) {
	p = ligne + 10;
	p = next_word(p, accession, ";\n ", sizeof(accession), NULL);
	if(accession[0] != 0 && strlen(accession) <= ACC_LENGTH ) 
		numacc = find_cre_acc(accession, TRUE);
	if(numacc != 0) {
		mdshrt(kacc, numacc, 1, seq_num, NULL);
		mdshrt(kloc, loc_num, 10, numacc, NULL);
		}

	p = next_word(p, accession, ";\n ", sizeof(accession), NULL);
	numacc = 0;
	if( ( p = strchr(accession, '.') ) != NULL) *p = 0;
	if(accession[0] != '-' && strlen(accession) <= ACC_LENGTH ) 
		numacc = find_cre_acc(accession, TRUE);
	if(numacc != 0) {
		mdshrt(kacc, numacc, 1, seq_num, NULL);
		mdshrt(kloc, loc_num, 10, numacc, NULL);
		}

	}
fgets(ligne, lligne, in_flat);
}


void process_dt(char *ligne, size_t lligne)
{
char *p;
static char date[20], acnuc_date[9];
static char month[]="JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC";
int mois, relnum, numkey;
enum {swiss, trembl, embl} db_type;
static char key_root[3][20] = {"SWISS-PROT RELEASES", "TREMBL RELEASES",
				"EMBL RELEASES"};
static char rel_name[3][40] = { "SP-RELEASE ", "TREMBL-RELEASE ",
				"EMBL-RELEASE " };
static const int l_name[3] = { 11, 15, 13 };

fgets(ligne, lligne, in_flat);
fgets(ligne, lligne, in_flat);
majuscules(ligne);
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
p = strchr(ligne, '(');
if(p != NULL) {
	if(strncmp(p, "(REL.", 5) == 0) {
		db_type = swiss;
		}
	else if(strncmp(p, "(TREMBLREL.", 11) == 0) {
		db_type = trembl;
		}
	else	{
		db_type = embl;
		}
	p = strchr(p, ' ');
	relnum  = -1;
	if(p != NULL) sscanf(p, "%d", &relnum);
	if(relnum != -1) {
		sprintf(rel_name[db_type] + l_name[db_type], "%.2d", relnum);
		numkey = crekeyword(key_root[db_type], rel_name[db_type]);
		fast_add_seq_to_keyw(numkey, seq_num);
		mdshrt(ksub, seq_num, 4, numkey, NULL);
		}
	}

fgets(ligne, lligne, in_flat);
}


void process_gn(char *ligne, size_t lligne)
{
char name[100], *p, *q;
int numkey;

/* enlever les espaces points et \n terminaux */
p = ligne + strlen(ligne);
while( *(p-1) == '\n' || *(p-1) == '.' || *(p - 1) == ' ' ) p--;
*p = 0;
/* enlever les () */
while ( (p = strchr(ligne, '(') ) != NULL ) *p = ' ';
while ( (p = strchr(ligne, ')') ) != NULL ) *p = ' ';
p = ligne;
while( (p = strstr(ligne, " OR ") ) != NULL ) memset(p + 1, '~', 2);
while( (p = strstr(ligne, " AND ") ) != NULL ) memset(p + 1, '~', 3);
p = ligne + 5;
do	{
	while ( *p == ' ' || *p == '~' ) p++;
	if( *p == 0) break;
	p = next_word(p, name, "~", sizeof(name), NULL);
	q = name + strlen(name);
	while( *(--q) == ' ' ) *q = 0;
	numkey = crekeyword("GENETIC NAMES", name);
	fast_add_seq_to_keyw(numkey, seq_num);
	mdshrt(ksub, seq_num, 4, numkey, NULL);
	}
while ( p != NULL);
fgets(ligne, lligne, in_flat);
}


void process_os(char *ligne, size_t lligne)
{
static char name[1000];
char *p, *d_name, *q;
int numspec, trouve;

d_name = name;
do	{
	p = ligne + strlen(ligne) - 1; *p = 0; /* enlever \n terminal */
	/* enlever les espaces et . terminaux */
	while ( *(--p) == ' ' ||  *p == '.') *p = 0; 
	memcpy(d_name, ligne + 5, p - ligne - 4);
	d_name += (p - ligne - 4); *d_name = 0;
	fgets(ligne, lligne, in_flat);
	if( strncmp(ligne, "OS ", 3) != 0 ) break;
	*(d_name++) = ' '; /* si encore OS, ajouter espace en fin de deja lu */
	}
while( TRUE );
majuscules(name);
p = d_name - 1; /* p vers dernier caractere */
d_name = name;
do	{ /* decouper selon ", " hors les () */
	q = d_name; trouve = FALSE;
	do	{
		if( *q == '(' ) {
			q = nextpar(q);
			if ( q == NULL ) break;
			}
		else if( *q == ',' && *(q+1) == ' ') {
			trouve = TRUE;
			break;
			}
		q++;
		}
	while( *q != 0 );
	if( trouve ) {
		*q = 0;
		p = q - 1;
		}
	else	q = NULL;
	/* enlever recursivement les parentheses terminales (nom anglais) */
	p = remove_terminal_parenth(d_name, p);
	while( *d_name == ' ') d_name++; /* enlever espaces initiaux */
	if(strncmp(d_name, "AND ", 4) == 0) { /* enlever AND initial */
		d_name += 4;
		while( *d_name == ' ') d_name++; /* enlever espaces initiaux */
		}
	/* tronquer a WIDTH_KS */
	if( p + 1 - d_name > WIDTH_KS ) d_name[WIDTH_KS] = 0; 
	d_name = check_spec_with_sp(d_name);
	numspec = crespecies("UNCLASSIFIED", d_name);
	fast_add_seq_to_spec(numspec, seq_num);
	mdshrt(kloc, loc_num, 6, numspec, NULL);
	if(q != NULL) {
		d_name = q + 2;
		while( *d_name == ' ') d_name++;
		p = d_name + strlen(d_name) - 1;
		}
	}
while( q != NULL);
}


void process_oc(char *ligne, size_t lligne)
{
int numspec, doit, npere, desc, next;
char *ascend, specnam[WIDTH_KS + 1];

doit = TRUE;
readloc(loc_num);
if(ploc->spec == 0) doit = FALSE;
if(doit) {
	readshrt(ploc->spec);
	numspec = pshrt->val;
	readspec(numspec);
	memcpy(specnam, pspec->name, WIDTH_KS); specnam[WIDTH_KS] = 0;
	trim_key(specnam);
	}
if(doit && !is_desc("UNCLASSIFIED", numspec) ) doit = FALSE;
if(doit) {
	ascend = proc_classif(ligne, lligne);
	if( strcmp(ascend, "UNCLASSIFIED") == 0) doit = FALSE;
	if( strcmp(ascend, specnam) == 0) doit = FALSE;
	}
if(doit) {
/* mettre numspec descendant de ascend */
	npere = iknum(ascend, kspec);
	readspec(npere);
	while(pspec->syno > 0) {
		npere = pspec->syno;
		readspec(npere);
		}
	readspec(numspec);
	desc = pspec->desc;
	readspec(npere);
	readshrt(pspec->desc);
	pshrt->val = - npere;
	writeshrt(pspec->desc);
	addshrt(pspec->desc, desc);
	if(pspec->plsub == 0 && pspec->plhost == 0) {
		readshrt(pspec->desc);
		pshrt->val = npere;
		writeshrt(pspec->desc);
		}
/* enlever numspec descendant de unclassified */
	npere = iknum("UNCLASSIFIED", kspec);
	readspec(npere);
	readshrt(pspec->desc);
	next = pshrt->next;
	doit = supshrt(next, desc);
	if(doit == 2) {
		readshrt(pspec->desc);
		pshrt->next = 0;
		writeshrt(pspec->desc);
		}
	}
while(strncmp(ligne, "OC", 2) == 0)
	fgets(ligne, lligne, in_flat);
}


static int is_desc(char *pere, int numspec)
{
int npere, next, val;

npere = iknum(pere, kspec);
if(npere == 0) return  FALSE;
readspec(numspec);
val = pspec->desc;
readspec(npere);
readshrt(pspec->desc);
next = pshrt->next;
while(next != 0) {
	readshrt(next);
	if(pshrt->val == val) return TRUE;
	next = pshrt->next;
	}
return FALSE;
}


char *check_spec_with_sp(char *name)
/* ajouter . si espece se termine par " SP" */
{
static char new_name[WIDTH_KS + 1];
int l;

l = strlen(name);
if(l >= WIDTH_KS || strcmp(name + l - 3, " SP") != 0) return name;
memcpy(new_name, name, l);
strcpy(new_name + l, ".");
return new_name;
}


void process_og(char *ligne, size_t lligne)
{
int num_org, do_it, *extra;
static char name[20] = "05";

majuscules(ligne);
if( strncmp(ligne + 5, "PLASMID", 7) == 0) {
	process_og_plasmid(ligne, lligne);
	return;
	}
do_it = FALSE;
if( strncmp(ligne + 5, "CHLOROPLAST", 11) == 0 ) {
	memcpy(name + 2, ligne + 5, 11);
	name[13] = 0;
	do_it = TRUE;
	}
else if( strncmp(ligne + 5, "CYANELLE", 8) == 0 ) {
	memcpy(name + 2, ligne + 5, 8);
	name[10] = 0;
	do_it = TRUE;
	}
else if( strncmp(ligne + 5, "MITOCHONDRION", 13) == 0 ) {
	memcpy(name + 2, ligne + 5, 13);
	name[15] = 0;
	do_it = TRUE;
	}
if(do_it) {
	num_org = find_key(name, smj_binary_tree, TRUE, &extra);
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
fgets(ligne, lligne, in_flat);
}


void process_og_plasmid(char *ligne, size_t lligne)
{
static char name[1000];
char *p, *d_name, *q;
int numspec, trouve;

d_name = name;
do	{
	p = ligne + strlen(ligne) - 1; *p = 0; /* enlever \n terminal */
	/* enlever les espaces et . terminaux */
	while ( *(--p) == ' ' ||  *p == '.') *p = 0; 
	memcpy(d_name, ligne + 5, p - ligne - 4);
	d_name += (p - ligne - 4); *d_name = 0;
	fgets(ligne, lligne, in_flat);
	if( strncmp(ligne, "OG ", 3) != 0 ) break;
	*(d_name++) = ' '; /* si encore OG, ajouter espace en fin de deja lu */
	}
while( TRUE );
p = d_name - 1; /* p vers dernier caractere */
d_name = name;
do	{ /* decouper selon ", " hors les () */
	q = strstr(d_name, ", ");
	if( q != NULL ) {
		*q = 0;
		p = q - 1;
		}
	else	q = NULL;
	/* enlever recursivement les parentheses terminales */
	p = remove_terminal_parenth(d_name, p);
	while( *d_name == ' ') d_name++; /* enlever espaces initiaux */
	if(strncmp(d_name, "AND ", 4) == 0) { /* enlever AND initial */
		d_name += 4;
		while( *d_name == ' ') d_name++; /* enlever espaces initiaux */
		}
	/* tronquer a WIDTH_KS */
	if( p + 1 - d_name > WIDTH_KS ) d_name[WIDTH_KS] = 0; 
	d_name = check_spec_with_sp(d_name);
	numspec = crespecies("PLASMIDS", d_name);
	fast_add_seq_to_spec(numspec, seq_num);
	mdshrt(kloc, loc_num, 6, numspec, NULL);
	if(q != NULL) {
		d_name = q + 2;
		while( *d_name == ' ') d_name++;
		p = d_name + strlen(d_name) - 1;
		}
	}
while( q != NULL);
}


void process_kw(char *ligne, size_t lligne)
{
static char name[100];
char *p, *d_name, *q;
int numkey;

p = ligne + 5;
do	{
	p = next_word(p, name, ";\n", sizeof(name), NULL);
	if(name[0] == 0) break;
	d_name = name;
	while( *d_name == ' ') d_name++;
	q = d_name + strlen(d_name);
	while ( *(q - 1) == ' ' ) *(--q) = 0;/* enlever les espaces terminaux */
	if( *(q - 1) == '.' ) *(--q) = 0; /* enlever le . terminal */
	if( q <= d_name ) continue;
	numkey = crekeyword(NULL, d_name);
	fast_add_seq_to_keyw(numkey, seq_num);
	mdshrt(ksub, seq_num, 4, numkey, NULL);
	}
while ( p != NULL );
fgets(ligne, lligne, in_flat);
}


void process_de(char *ligne, size_t lligne)
{
static char name[1000], ec_code[20];
char *p, *d_name, *q;
int numkey, l;

d_name = name; *name = 0;
do	{
	p = ligne + strlen(ligne) - 1; *p = 0; /* enlever \n terminal */
	while ( *(--p) == ' ' ) *p = 0; /* enlever les espaces terminaux */
	if(p - ligne - 4 > 0) {
		memcpy(d_name, ligne + 5, p - ligne - 4);
		d_name += (p - ligne - 4); *d_name = 0;
		}
	fgets(ligne, lligne, in_flat);
	if( strncmp(ligne, "DE ", 3) != 0 ) break;
	*(d_name++) = ' '; *d_name = 0; /* espace de liaison entre lignes */
	}
while( TRUE );
trim_key(name);
l =  strlen(name);
if(l == 0) return;
p = name + l - 1; /* p vers dernier caractere */
/* process fragment(s) */
majuscules(name);
if( ( q = strstr(name, "(FRAGMENT)") ) != NULL ||
		 ( q = strstr(name, "(FRAGMENTS)") ) != NULL ){
	numkey = crekeyword(NULL, "PARTIAL");
	fast_add_seq_to_keyw(numkey, seq_num);
	mdshrt(ksub, seq_num, 4, numkey, NULL);
	/* enlever de name tout le (FRAGMENTs)  */
	while ( *q != ')' ) *(q++) = ' ';
	*q = ' ';
	while (p > name && *p == ' ') p--;
	}
/* process EC NUMBERS (plusieurs) */
while( ( d_name = strstr(name, " (EC ") ) != NULL) { 
	next_word(d_name + 5, ec_code,  ") ", sizeof(ec_code), NULL );
	numkey = crekeyword("EC_NUMBERS", ec_code);
	fast_add_seq_to_keyw(numkey, seq_num);
	mdshrt(ksub, seq_num, 4, numkey, NULL);
	/* enlever de name tout le (EC xxx)  */
	q = d_name + 4;
	while( *q != ')' && *(q + 1) != 0 ) q++;
	q++; 
	memmove( d_name, q, p - q + 2 );
	p -= (q - d_name);
	}
if(p <= name) return;
process_de_part(name, p);
}


static void process_de_part(char *name, char *p)
/* name: debut partie a traiter
p: fin partie a traiter (  *(p+1) == 0 )
*/
{
int numkey;
char *d_name, *q, *fin;

/* traiter les () terminales */
p = process_terminal_parenth(name, p);
/* decouper selon   " / "   */
d_name = name;
do	{
	q = strstr(d_name, " / ");
	if( q == NULL ) 
		q = p + 1;
	*q = 0;	
	fin = process_terminal_parenth(d_name, q - 1);
	if( fin > d_name ) {
		*(fin + 1) = 0;
		if(strncmp(d_name, "CONTAINS: ", 10) == 0) {
			d_name += 10;
			while( *d_name == ' ') d_name++;
			}
		if( fin - 9 > d_name && 
				strncmp(fin - 9, " PRECURSOR", 10) == 0) {
			fin -= 10;
			*(fin + 1) = 0;
			}
		numkey = crekeyword(NULL, d_name);
		fast_add_seq_to_keyw(numkey, seq_num);
		mdshrt(ksub, seq_num, 4, numkey, NULL);
		}
	d_name = q + 3;
	if(d_name <= p) 
		while( *d_name == ' ') d_name++;
	}
while( d_name <= p );
}


static char *process_terminal_parenth(char *name, char *p)
/* name: debut partie a traiter
p: fin partie a traiter (  *(p+1) == 0 )
*/
{
int tot_par;
char *debut, *fin;

/* enlever les espaces ou . terminaux */
while ( *p == ' ' || *p =='.' ) *(p--) = 0; 
/* enlever les parentheses terminales (plusieurs et recursivement) */
while( *p == ')' ) {
	tot_par = 0;
	fin = p - 1;
	do	{
		if( *p == ')' ) tot_par++;
		if( *p == '(' ) tot_par--;
		p--;
		}
	while( tot_par > 0  && p >= name );
	*(p+1) = 0;
	debut = p + 2; *(fin + 1) = 0;
	process_de_part(debut, fin);
	/* enlever les espaces ou . terminaux */
	while ( *p == ' ' || *p == '.' ) *(p--) = 0; 
	}
return p;
}


char *remove_terminal_parenth(char *name, char *p)
/* name: debut partie a traiter
p: fin partie a traiter (  *(p+1) == 0 )
*/
{
int tot_par;

/* enlever les espaces ou . terminaux */
while ( *p == ' ' || *p =='.' ) *(p--) = 0; 
/* enlever les parentheses terminales (plusieurs et recursivement) */
while( *p == ')' ) {
	tot_par = 0;
	do	{
		if( *p == ')' ) tot_par++;
		if( *p == '(' ) tot_par--;
		p--;
		}
	while( tot_par > 0  && p >= name );
	*(p+1) = 0;
	/* enlever les espaces ou . terminaux */
	while ( *p == ' ' || *p == '.' ) *(p--) = 0; 
	}
return p;
}


void process_cc(char *ligne, size_t lligne)
{
do	{
	if(strncmp(ligne + 5, "-!- SUBCELLULAR LOCATION:", 25) == 0 )
			process_subcell_loc(ligne, lligne);
	else if(strncmp(ligne + 5, "-!- GENE_FAMILY:", 16) == 0 )
			process_gene_family(ligne, lligne);
	else fgets(ligne, lligne, in_flat);
	}
while( strncmp(ligne, "CC ", 3) == 0);
}


void process_subcell_loc(char *ligne, size_t lligne)
{
char *p, *q, *d_name;
static char keyword[300];
int numkey;

d_name = keyword;
do	{
	p = ligne + strlen(ligne) - 1; *p = 0; /* enlever \n terminal */
	/* enlever les espaces et . terminaux */
	while ( *(p-1) == '.' ||  *(p-1) == ' ') *(--p) = 0; 
	q = ligne + 5; while( *q == ' ') q++;
	memcpy(d_name, q, p - q);
	d_name += (p - q); *d_name = 0;
	fgets(ligne, lligne, in_flat);
	if( strcmptrail(ligne, 7, "CC", 2) != 0 ) break;
/* si encore CC subcell-loc, ajouter espace en fin de deja lu */
	*(d_name++) = ' '; 
	}
while( TRUE );
p = keyword + 25;
while(*p == ' ') p++;
*(p + 40) = 0; /* tronquer a 40 */
trim_key(p);
numkey = crekeyword("SUBCELLULAR LOCATION", p);
fast_add_seq_to_keyw(numkey, seq_num);
mdshrt(ksub, seq_num, 4, numkey, NULL);
}


void process_gene_family(char *ligne, size_t lligne)
{
char *p, *q, *d_name;
static char keyword[300];
int numkey;

p = ligne + 21;
while(*p == ' ') p++;
q = ligne + strlen(ligne) - 1; *(q--) = 0;
while(*q == ' ' || *q == '.') *(q--) = 0;
numkey = crekeyword("GENE_FAMILY", p);
fast_add_seq_to_keyw(numkey, seq_num);
mdshrt(ksub, seq_num, 4, numkey, NULL);
fgets(ligne, lligne, in_flat);
}


void process_ft(char *ligne, size_t lligne)
{
static char feat_name[20];
int numkey;

do	{
	if( ligne[5] != ' ') {
		next_word(ligne + 5, feat_name, " ", sizeof(feat_name), NULL);
		numkey = crekeyword("MISC_FEATURE", feat_name);
		fast_add_seq_to_keyw(numkey, seq_num);
		mdshrt(ksub, seq_num, 4, numkey, NULL);
		if(strcmp(feat_name, "PRODOM") == 0) process_prodom(ligne);
		}
	fgets(ligne, lligne, in_flat);
	}
while( strncmp(ligne, "FT ", 3) == 0);
}


void process_prodom(char *ligne)
{
static char feat_name[WIDTH_KS + 1];
int numkey;

next_word(ligne + 34, feat_name, " ", sizeof(feat_name), NULL);
if(feat_name[0] == 0) return;
numkey = crekeyword("PRODOM DOMAINS", feat_name);
fast_add_seq_to_keyw(numkey, seq_num);
mdshrt(ksub, seq_num, 4, numkey, NULL);
}


void process_r_lines(char *ligne, size_t lligne)
{
off_t addr_ra, current;

addr_ra = 0;
do	{
	current = ftello(in_flat);
	fgets(ligne, lligne, in_flat);
	if(strncmp(ligne, "RA ", 3) == 0 && addr_ra == 0) 
			addr_ra = current;
	else if(strncmp(ligne, "RL ", 3) == 0) {
		process_rl(ligne, lligne, addr_ra);
		/* passer le reste de RL s'il y en a */
		while(strncmp(ligne, "RL ", 3) == 0) {
			fgets(ligne, lligne, in_flat);
			}
		}
	}
while( *ligne == 'R' && strncmp(ligne, "RN ", 3) != 0);
}


void process_rl(char *ligne, size_t lligne, off_t addr_ra)
{
int is_journal, is_book, is_thesis, rank, trouve, number, y_number,
	tot_txt, num_bib, l, special_format, compute_authors,
	*extra_year, *extra_journal, is_patent, is_submitted, *extra_bib;
static char name[300], j_name[150], aux[150],
	volume[20], first_page[20], annee[5],
	journal[200], first_author[100];
char *p, *d_name, *q, *colon;
off_t saved_addr;

is_journal = TRUE; is_book = FALSE; is_thesis = FALSE; is_patent = FALSE;
is_submitted = FALSE;
memcpy(aux, ligne + 5, 20); aux[20] = 0; majuscules(aux);
if(strncmp(aux, "(IN)", 4) == 0 ) is_book = TRUE;
else if(strncmp(aux, "THESIS ", 7) == 0) is_thesis = TRUE;
else if(strncmp(aux, "PATENT NUMBER", 13) == 0) is_patent = TRUE;
else if(strncmp(aux, "SUBMITTED ", 10) == 0 ) is_submitted = TRUE;
else if(strncmp(aux, "UNPUBLISHED ", 12) == 0 ) is_journal = FALSE;

if( is_book ) {
	process_book(ligne, lligne, addr_ra);
	return;
	}
else if( is_thesis ) {
	process_thesis(ligne, lligne, addr_ra);
	return;
	}
else if( is_patent ) {
	process_patent(ligne, lligne, addr_ra);
	return;
	}
else if( is_submitted ) {
	process_submitted(ligne, lligne, addr_ra);
	return;
	}
else if( ! is_journal ) {
	do	{
		fgets(ligne, lligne, in_flat);
		}
	while(strncmp(ligne, "RL ", 3) == 0);
	return;
	}
/* chargement partie journal sur plusieurs lignes */
d_name = name;
do	{
	p = ligne + strlen(ligne) - 1; *p = 0; /* enlever \n terminal */
	while ( *(--p) == ' ' ) *p = 0; /* enlever les espaces terminaux */
	memcpy(d_name, ligne + 5, p - ligne - 4);
	d_name += (p - ligne - 4); *d_name = 0;
	fgets(ligne, lligne, in_flat);
	}
while( strncmp(ligne, "RL ", 3) == 0 && *p != '.' );
/* recherche ':' separant journal + volume du reste */
colon = strchr(name, ':'); 
do	{
	if( colon == NULL ) return;
	q = colon - 1;
	while( q > name && *q != ' ' ) q--;
	if(q <= name) colon = strchr(colon + 1, ':');
	}
while(q <= name);
*colon = 0; *q = 0;
strcpy(volume, q + 1);
/* calcul first_page */
p = colon + 1;
while( *p != '-' && *p != '(' ) p++;
memcpy(first_page, colon + 1, p - 1 - colon); first_page[p - 1 - colon] = 0;
compact(first_page);
/* calcul annee */
while( *p != '(' ) p++;
q = strchr(p, ')');
if( q == NULL) return; /* annee mal formee */
memcpy(annee, p + 1, q - 1 - p); annee[q - 1 - p] = 0;
if(strncmp(annee, "19", 2) != 0 && strncmp(annee, "20", 2) != 0) 
				strcpy(annee, "0");

/* calcul du nom du journal */
strcpy(journal, name);
q = journal - 1;
while( *(++q) != 0) {
	*q = toupper(*q);
	if ( *q == '.') *q = ' ';
	}
compact(journal);
/* recherche du libel du journal */
trouve = FALSE;
for(rank = 0; rank < tot_journals; rank++) {
	if(strcmp(journal, j_libel[rank]) == 0) {
		trouve = TRUE;
		break;
		}
	}
if( trouve ) {
	sprintf(aux, "02%s", j_code[rank]);
	number = find_key(aux, smj_binary_tree, FALSE, &extra_journal);
	if(number != j_number[rank]) {
		fprintf(log_file, "error in j_number for %s\n", aux);
		non_chargee = TRUE;
		return;
		}
	}
else	{ /* creation nouveau journal */
	if(tot_journals >= MAX_JOURNALS) {
		fprintf(log_file, "cannot create new journal %s\n", name);
		non_chargee = TRUE;
		return;
		}
	l = strlen(journal); if(l > 18) l = 18;
	j_code[tot_journals] = malloc(l + 1);
	if( j_code[tot_journals] == NULL) {
		fprintf(log_file, "not enough memory\n");
		non_chargee = TRUE;
		return;
		}
	memcpy(j_code[tot_journals], journal, l); 
	j_code[tot_journals][l] = 0;
	l = strlen(journal); if(l > lrtxt) l = lrtxt;
	j_libel[tot_journals] = malloc(l + 1);
	if( j_libel[tot_journals] == NULL) {
		fprintf(log_file, "not enough memory\n");
		non_chargee = TRUE;
		return;
		}
	memcpy(j_libel[tot_journals], journal, l); 
	j_libel[tot_journals][l] = 0;
	rank = tot_journals;
	sprintf(aux, "02%s", j_code[tot_journals]);
	number = find_key(aux, smj_binary_tree, TRUE, &extra_journal);
	if(number == 0) {
		fprintf(log_file, "cannot create %s\n", aux);
		non_chargee = TRUE;
		return;
		}
	j_number[tot_journals] = number;
	tot_journals++;
	tot_txt = read_first_rec(ktxt, NULL);
	tot_txt++;
	l = strlen(name);
	memset(ptxt, ' ', lrtxt);
	memcpy(ptxt, name, (l > lrtxt ? lrtxt : l) );
	writetxt(tot_txt);
	write_first_rec(ktxt, tot_txt, 0);
	readsmj(number);
	psmj->libel = tot_txt;
	writesmj(number);
	}
/* recherche annee */
y_number = process_year_bib(annee, &extra_year);
/* creation de la reference */
special_format = ( strcmp(volume, "0") == 0 || strcmp(first_page, "0") == 0 );
/* cas JOURNAL 0:0-0(annee). */
if( special_format ) {
/* ici fabriquer JOURNAL/ANNEE/AUTEUR1 */
	if( addr_ra == 0 ) return; /* refer sans RA */
	saved_addr = ftello(in_flat);
	fseeko(in_flat, addr_ra, SEEK_SET);
	fgets(name, sizeof(name), in_flat);
	fseeko(in_flat, saved_addr, SEEK_SET);
	next_word(name + 5, first_author, ";,", sizeof(first_author), NULL);
	if( *first_author == 0 ) return; /* refer sans auteur */
	remove_initials(first_author);
	sprintf(name, "%s/%s/%s", j_code[rank], annee, first_author);
	}
else
	sprintf(name, "%s/%s/%s", j_code[rank], volume, first_page);
compact(name); majuscules(name);
name[40] = 0;
num_bib = find_key(name, bib_binary_tree, TRUE, &extra_bib);
if(num_bib == 0) {
	fprintf(log_file, "cannot create %s\n", name);
	non_chargee = TRUE;
	return;
	}
/* si refer a au moins une seq associee, ne pas refaire les auteurs */
readbib(num_bib);
compute_authors = (pbib->plsub == 0);
connect_locus_bib(loc_num, seq_num, num_bib, extra_bib);
readbib(num_bib);
pbib->j = j_number[rank];
/* ne pas changer l'annee d'une reference deja avec des sequences */
if( ( pbib->y == 0 ) || compute_authors) {
	pbib->y = y_number;
	}
else	{
	if(y_number != pbib->y) {/*deja avec des seqs mais d'une autre annee! */
		extra_year = NULL;
		y_number = pbib->y;
		}
	}
writebib(num_bib);
add_to_smjlng(j_number[rank], seq_num, extra_journal);
if(y_number != 0) add_to_smjlng(y_number, seq_num, extra_year);
/* traitement auteurs */
if(compute_authors) process_ra(addr_ra, num_bib, special_format);
}


void process_ra(off_t addr_ra, int num_bib, int first_aut_only)
{
int num_aut, num, tot_old;
char *p, *q, *old_q;
static char auteur[100], lu[100];
static int *old_auteurs = NULL, max_old = 0;
off_t saved_addr;

/* traitement auteurs */
if(addr_ra == 0) return; /* reference sans ligne RA */
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
saved_addr = ftello(in_flat);
fseeko(in_flat, addr_ra, SEEK_SET);
fgets(lu, sizeof(lu), in_flat);
do	{
	p = lu + 5;
	do	{
		while( *p == ' ') p++;
		p = next_word(p, auteur, ",;\n", sizeof(auteur), NULL);
		if( auteur[0] == 0) break;
		remove_initials(auteur);
		majuscules(auteur);
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
		if(first_aut_only) break;
		}
	while( p != NULL);
	if( !first_aut_only ) fgets(lu, sizeof(lu), in_flat);
	}
while( (!first_aut_only) && strncmp(lu, "RA ", 3) == 0);
/* enlever les anciens auteurs qui ne sont plus cites */
for(num = 0; num < tot_old; num++) {
	if( (num_aut = old_auteurs[num]) == 0) continue;
	mdshrt(kaut, num_aut, -1, num_bib, NULL);
	mdshrt(kbib, num_bib, -2, num_aut, NULL);
	}
fseeko(in_flat, saved_addr, SEEK_SET);
}


void remove_initials(char *auteur)
{
char *q;

/* suppression des initiales */
q = strchr(auteur, '.');
if(q == NULL) q = auteur + strlen(auteur) - 1;
else	q -= 2;
/* suppression des espaces terminaux */
while( *q == ' ') q--;
*(q + 1) = 0;
if( q >= auteur + 20) { /* tronquer a 20 sans espaces terminaux */
	q = auteur + 20;
	*q = 0;
	while( q > auteur && *(--q) == ' ') *q = 0;
	}
}


void process_book(char *ligne, size_t lligne, off_t addr_ra)
{
char *p;
static char annee[5], auteur[60], reference[100];
int num_bib, num_book, num_year, *extra_journal, *extra_year,
	deja_connu, *extra_bib;
off_t saved_addr;

if( addr_ra == 0 ) return; /* refer sans RA */
do	{ /* dans reference la derniere ligne RL du groupe */
	strcpy(reference, ligne);
	fgets(ligne, lligne, in_flat);
	}
while(strncmp(ligne, "RL ", 3) == 0);
p = reference + strlen(reference) - 1; *p = 0;
if( *(p - 1) == '.' ) *(--p) = 0; /* enlever le . terminal */
if( *(p - 1) == ')' ) p -= 5;
else p -= 4;
memcpy(annee, p, 4); annee[4] = 0;
if(strncmp(annee, "19", 2) != 0 && strncmp(annee, "20", 2) != 0) 
				strcpy(annee, "0");
saved_addr = ftello(in_flat);
fseeko(in_flat, addr_ra, SEEK_SET);
fgets(reference, sizeof(reference), in_flat);
fseeko(in_flat, saved_addr, SEEK_SET);
next_word(reference + 5, auteur, ";,", sizeof(auteur), NULL);
if( *auteur == 0 ) return; /* refer sans auteur */
remove_initials(auteur);
sprintf(reference, "BOOK/%s/%s", annee, auteur);
compact(reference);
reference[40] = 0; majuscules(reference);
num_bib = find_key(reference, bib_binary_tree, TRUE, &extra_bib); 
if(num_bib == 0) {
	fprintf(log_file, "cannot create %s\n", reference);
	non_chargee = TRUE;
	return;
	}
num_book = find_key("02BOOK", smj_binary_tree,  TRUE, &extra_journal); 
if(num_book == 0) {
	fprintf(log_file, "cannot create 02BOOK\n");
	non_chargee = TRUE;
	return;
	}
num_year = process_year_bib(annee, &extra_year);
readbib(num_bib);
deja_connu = ( pbib->plsub != 0 );
connect_locus_bib(loc_num, seq_num, num_bib, extra_bib);
readbib(num_bib);
pbib->j = num_book;
pbib->y = num_year;
writebib(num_bib);
add_to_smjlng(num_book, seq_num, extra_journal);
if(num_year != 0) add_to_smjlng(num_year, seq_num, extra_year);
/* traitement auteurs */
if( ! deja_connu) process_ra(addr_ra, num_bib, TRUE);
}


void process_thesis(char *ligne, size_t lligne, off_t addr_ra)
{
char *p;
static char annee[5], auteur[60], reference[100];
int num_bib, num_book, num_year, *extra_journal, *extra_year,
	deja_connu, *extra_bib;
off_t saved_addr;

if( addr_ra == 0 ) return; /* refer sans RA */
p = strchr(ligne, '(');
if( p == NULL) {
	fprintf(log_file, "cannot find year in %s", ligne);
	non_chargee = TRUE;
	return;
	}	
memcpy(annee, p + 1, 4); annee[4] = 0;
saved_addr = ftello(in_flat);
fseeko(in_flat, addr_ra, SEEK_SET);
fgets(reference, sizeof(reference), in_flat);
fseeko(in_flat, saved_addr, SEEK_SET);
next_word(reference + 5, auteur, ";,", sizeof(auteur), NULL);
if( *auteur == 0 ) return; /* refer sans auteur */
remove_initials(auteur);
sprintf(reference, "THESIS/%s/%s", annee, auteur);
compact(reference); majuscules(reference);
num_bib = find_key(reference, bib_binary_tree, TRUE, &extra_bib); 
if(num_bib == 0) {
	fprintf(log_file, "cannot create %s\n", reference);
	non_chargee = TRUE;
	return;
	}
num_book = find_key("02THESIS", smj_binary_tree, TRUE, &extra_journal); 
if(num_book == 0) {
	fprintf(log_file, "cannot create 02THESIS\n");
	non_chargee = TRUE;
	return;
	}
num_year = process_year_bib(annee, &extra_year);
readbib(num_bib);
deja_connu = ( pbib->plsub != 0 );
connect_locus_bib(loc_num, seq_num, num_bib, extra_bib);
readbib(num_bib);
pbib->j = num_book;
pbib->y = num_year;
writebib(num_bib);
add_to_smjlng(num_book, seq_num, extra_journal);
add_to_smjlng(num_year, seq_num, extra_year);
/* traitement auteurs */
if( ! deja_connu ) process_ra(addr_ra, num_bib, TRUE);
return;
}


void process_submitted(char *ligne, size_t lligne, off_t addr_ra)
{
char *p;
static char annee[5], auteur[60], reference[100];
int num_bib, num_book, num_year, *extra_journal, *extra_year,
	deja_connu, *extra_bib;
off_t saved_addr;

if(addr_ra == 0) return; /* reference sans RA */
p = strchr(ligne, ')');
if( p == NULL) {
	fprintf(log_file, "cannot find year in %s", ligne);
	non_chargee = TRUE;
	return;
	}	
memcpy(annee, p - 4, 4); annee[4] = 0;
if(strncmp(annee, "19", 2) != 0 && strncmp(annee, "20", 2) != 0) 
				strcpy(annee, "0");
saved_addr = ftello(in_flat);
fseeko(in_flat, addr_ra, SEEK_SET);
fgets(reference, sizeof(reference), in_flat);
fseeko(in_flat, saved_addr, SEEK_SET);
next_word(reference + 5, auteur, ";,", sizeof(auteur), NULL);
if( *auteur == 0 ) return; /* refer sans auteur */
remove_initials(auteur);
sprintf(reference, "SUBMITTED/%s/%s", annee, auteur);
compact(reference); majuscules(reference);
num_bib = find_key(reference, bib_binary_tree, TRUE, &extra_bib); 
if(num_bib == 0) {
	fprintf(log_file, "cannot create %s\n", reference);
	non_chargee = TRUE;
	return;
	}
num_book = find_key("02SUBMITTED", smj_binary_tree, TRUE, &extra_journal); 
if(num_book == 0) {
	fprintf(log_file, "cannot create 02SUBMITTED\n");
	non_chargee = TRUE;
	return;
	}
num_year = process_year_bib(annee, &extra_year);
readbib(num_bib);
deja_connu = ( pbib->plsub != 0 );
connect_locus_bib(loc_num, seq_num, num_bib, extra_bib);
readbib(num_bib);
pbib->j = num_book;
pbib->y = num_year;
writebib(num_bib);
add_to_smjlng(num_book, seq_num, extra_journal);
if(num_year != 0) add_to_smjlng(num_year, seq_num, extra_year);
/* traitement auteurs */
if( ! deja_connu ) process_ra(addr_ra, num_bib, TRUE);
return;
}



void process_patent(char *ligne, size_t lligne, off_t addr_ra)
{
char *p, *q;
static char annee[5], reference[200], number[100];
int num_bib, num_book, num_year, *extra_journal, *extra_year, deja_connu,
	*extra_bib;

/* recherche du numero */
p = ligne + 5 + 13;
while( *p == ' ') p++;
q = strchr(p, ',');
if(q == NULL) q = strchr(p, '.');
if( q == NULL) {
	fprintf(log_file, "cannot find patent number in %s", ligne);
	non_chargee = TRUE;
	return;
	}	
memcpy(number, p, q - p ); number[ q - p ] = 0;
/* recherche de l'annee (peut etre sur ligne suivante) */
p = NULL;
if( *q == ',' ) {
	p = strchr(q, '.');
	if(p == NULL) {
		fgets(ligne, lligne, in_flat);
		p = strchr(ligne, '.');
		}
	}
if(p != NULL) { memcpy(annee, p - 4, 4); annee[4] = 0; }
else strcpy(annee, "0");
sprintf(reference, "PATENT/%s", number);
compact(reference); majuscules(reference);
num_bib = find_key(reference, bib_binary_tree, TRUE, &extra_bib); 
if(num_bib == 0) {
	fprintf(log_file, "cannot create %s\n", reference);
	non_chargee = TRUE;
	return;
	}
num_book = find_key("02PATENT", smj_binary_tree, TRUE, &extra_journal); 
if(num_book == 0) {
	fprintf(log_file, "cannot create 02PATENT\n");
	non_chargee = TRUE;
	return;
	}
num_year = process_year_bib(annee, &extra_year);
readbib(num_bib);
deja_connu = ( pbib->plsub != 0 );
connect_locus_bib(loc_num, seq_num, num_bib, extra_bib);
readbib(num_bib);
pbib->j = num_book;
if( !deja_connu )
	pbib->y = num_year;
else	{
	num_year = pbib->y;
	extra_year = NULL;
	}
writebib(num_bib);
add_to_smjlng(num_book, seq_num, extra_journal);
if(num_year != 0) add_to_smjlng(num_year, seq_num, extra_year);
/* traitement auteurs */
if(! deja_connu) process_ra(addr_ra, num_bib, FALSE);
return;
}


int load_j_code_libel(void)
{
int tot_j, tot_smj, num, l;
char *p, aux[lrtxt + 1];

tot_smj = read_first_rec(ksmj, NULL);
tot_j = 0;
for(num = 2; num <= tot_smj; num++) {
	readsmj(num);
	if(strncmp(psmj->name, "02", 2) != 0) continue;
	if(     strncmp(psmj->name, "02BOOK ", 7) == 0 ||
		strncmp(psmj->name, "02THESIS ", 9) == 0 ||
		strncmp(psmj->name, "02SUBMITTED ", 12) == 0 ||
		strncmp(psmj->name, "02PATENT ", 9) == 0 ) continue;
	if( tot_j >= MAX_JOURNALS ) goto too_many;
	j_code[tot_j] = (char *)malloc(19);
	if( j_code[tot_j] == NULL ) goto memoire;
	memcpy( j_code[tot_j], psmj->name + 2, 18);
	j_code[tot_j][18] = 0;
	p = j_code[tot_j] + 18;
	while( *(--p) == ' ') *p = 0;
	readtxt(psmj->libel);
	memcpy(aux, ptxt, lrtxt); aux[lrtxt] = 0;
	p = aux - 1;
	while( *(++p) != 0 ) {
		*p = toupper(*p);
		if(*p == '.') *p = ' ';
		}
	compact(aux);
	l = strlen(aux);
	j_libel[tot_j] = (char *)malloc(l + 1);
	if( j_libel[tot_j] == NULL ) goto memoire;
	memcpy( j_libel[tot_j], aux, l);
	j_libel[tot_j][l] = 0;
	
	j_number[tot_j] = num;
	tot_j++;
	}
return tot_j;

memoire:
fprintf(log_file, "ERROR: not enough memory for load_j_code_libel\n");
exit(ERREUR);

too_many:
fprintf(log_file, "ERROR: increase MAX_JOURNALS\n");
exit(ERREUR);
}


int remove_current_seq(int use_address, char *err_file_name, int current_div)
{
int point, valeur, num_smj, *extra, fullstop, journal, annee;
char seq_name[L_MNEMO + 1];

fullstop = FALSE;
/* traitement LOCUS */
readloc(loc_num);
point = ploc->spec; /* traitement especes */
while( point != 0 ) {
	readshrt(point);
	valeur = pshrt->val;
	point = pshrt->next;
	mdlng(kspec, valeur, -2, seq_num, NULL);
	extra = get_spec_extra(valeur);
	if(extra != NULL) *extra = 0;
	}
point = ploc->plref; /* traitement references */
while( point != 0 ) {
	readshrt(point);
	valeur = pshrt->val;
	point = pshrt->next;
	mdshrt(kbib, valeur, -1, seq_num, NULL);
	extra = get_bib_extra(valeur);
	if( extra == NULL ) {
		fprintf(log_file, "logical error removing biblio\n");
		fullstop = TRUE;
		}
	else	*extra = 0;
	readbib(valeur);
	journal = pbib->j;
	annee = pbib->y;
	if(journal != 0) {
		mdlng(ksmj, journal, -1, seq_num, NULL);
		extra = get_smj_extra(journal);
		if( extra == NULL ) {
			fprintf(log_file, "logical error removing journal\n");
			fullstop = TRUE;
			}
		else	*extra = 0;
		}
	if(annee != 0) {
		mdlng(ksmj, annee, -1, seq_num, NULL);
		extra = get_smj_extra(annee);
		if( extra == NULL ) {
			fprintf(log_file, "logical error removing year\n");
			fullstop = TRUE;
			}
		else	*extra = 0;
		}
	}
point = ploc->placc; /* traitement accession */
while( point != 0 ) {
	readshrt(point);
	point = pshrt->next;
	mdshrt(kacc, pshrt->val, -1, seq_num, NULL);
	}
/* traitement status */
if(ploc->stat != 0) {
	mdlng(ksmj, ploc->stat, -1, seq_num, NULL);
	extra = get_smj_extra(ploc->stat);
	if( extra == NULL ) {
		fprintf(log_file, "logical error removing status\n");
		fullstop = TRUE;
		}
	else	*extra = 0;
	}
/* traitement organelle */
if(ploc->org != 0) {
	mdlng(ksmj, ploc->org, -1, seq_num, NULL);
	extra = get_smj_extra(ploc->org);
	if( extra == NULL ) {
		fprintf(log_file, "logical error removing organelle\n");
		fullstop = TRUE;
		}
	else	*extra = 0;
	}

/* traitement SUBSEQ */
readsub(seq_num);
memcpy(seq_name, psub->name, L_MNEMO); seq_name[L_MNEMO] = 0;
point = psub->plkey; /* traitement mots-cles */
while( point != 0 ) {
	readshrt(point);
	valeur = pshrt->val;
	point = pshrt->next;
	mdlng(kkey, valeur, -2, seq_num, NULL);
	extra = get_keyw_extra(valeur);
	if(extra != NULL) *extra = 0;
	}

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
return fullstop;
}




