#include "dir_acnuc.h"
#include <ctype.h>
#include <time.h>
#ifdef unix
#ifdef __APPLE__
#include <termios.h>  /* for Darwin unix on apple */
#else
#include <termio.h>
#endif
#include <unistd.h>
#include <fcntl.h>
#endif

#ifdef __INTEL__
void attente(void);
void attente(void)
{
char rep[10];
printf("Hit <RETURN>"); gets(rep);
}
char *mygets(char *ligne);
char *mygets(char *ligne)
{
char *p;
if( gets(ligne) == NULL) return NULL;
p = ligne + strlen(ligne) - 1;
while(*p == '\r' || *p == '\n') *(p--) = 0;
return ligne;
}
#define gets(p) mygets(p)
#endif


/* typedefs, globals */
typedef struct noeud {
	char *nom;
	char *libel;
	char *libel_upcase;
	int numero;
	int count;
	struct noeud *parent;
	struct paire *liste_desc;
	struct noeud *syno;
	} noeud;
struct paire {
	noeud *valeur;
	struct paire *next;
	};
typedef void (*line_out_function)(char *);
typedef char Boolean;
#define TOT_COMMS 19
static char command_list[TOT_COMMS][10] = {"", "stop", "help", "terminal", "select",
	"find", "general", "species", "keywords", "info", "codes", "names",
	"save", "delete", "lists", "short", "bases", "extract", "modify"};
typedef enum {unknown, stop, help, terminal, com_name_select,
	find, general, species, keywords, info, codes, names,
	save, com_name_delete, lists, com_name_short, bases, extract, modify
	} command_kind;
/* type des formats d'extraction possibles */
typedef enum { acnuc=1, gcg, fasta, analseq, flat } out_option;
out_option extract_format_choice = acnuc;
char extract_format_name[6][30] = {"", "acnuc", "gcg", "fasta", "analseq", "GenBank/EMBL/Swissprot/PIR" };
/* type des genres d'extraction possibles */
typedef enum { simple=1, translate, fragment, feature, region } extract_option;
extract_option extract_options_choice;

noeud **tab_noeud;
char nom_racine[] = "root";
#define LOG_FILE_NAME "query.out"   /* name of log file */
int tlist = 30; /* nbre de listes de bits en memoire */
int last_hidden_list;
static int der_num_liste=0;
#define MAXLISTES 100
char nom_liste[MAXLISTES][120]; /* pointeurs vers noms des listes */
char *deflnames[MAXLISTES]; /* noms listes pour def */
char deflisdef[MAXLISTES][80]; /* definition listes pour def */
int defoccup[MAXLISTES];
int deflocus[MAXLISTES];
char defgenre[MAXLISTES];
int defllen[MAXLISTES];
int *defbitlist; /* listes de bits pour interf avec routine def */
char default_list[20];
int nblpp; /* lignes par page */
int numlig; /* nbre courant de lignes ecrites */
int line_interrupt;
FILE *line_file;
/* pour les annotations */
#define TOTREC 30 /* max # of annot options */
int tot_options, wid_key;
char *annot_options;
#define ANNOT_OPTIONS_GBK  "\
ALL        ACCession  NID        VERsion    KEYwords   SEGment    SOUrce     \
  ORGanism REFerence    AUThors    TITle    \
  JOUrnal    MEDline    REMark   COMment    FEAtures   \
BASe count ORIgin     SEQuence   "
char annot_options_gbk[] = ANNOT_OPTIONS_GBK;
char annot_options_embl[]=
"ALLAC SV NI DT KW OS OC OG RN RC RP RX RA RT RL DR CC FH FT SQ SEQ";
char annot_options_swissprot[]=
"ALLAC DT GN OS OG OC OX RN RP RC RX RA RL CC DR KW FT SQ SEQ";
#define ANNOT_OPTIONS_NBRF  "\
ALL       \
ALTernate_CONtains  ORGanism  DATe      ACCessions\
REFerence \
COMment   GENetics  CLAssificaKEYwords  FEAture   \
SUMmary   SEQuence  "
char annot_options_nbrf[]= ANNOT_OPTIONS_NBRF;
Boolean info_options[TOTREC];
char *poptions[TOTREC], *poptions_maj[TOTREC];
int multi_spec;


/* prototypes of included functions */
command_kind identify_command(char *command);
char *proc_l_option(char *command);
void process_command(command_kind command_id, char *lname, char *command);
void command_lists(int use_lpt);
void command_delete(char *lname, char *command);
void command_short(char *lname, int use_lpt);
FILE *prep_out_file(int *use_lpt, char *command_name, char *list_name);
char *get_f_option(char *command);
void command_names(char *lname, int use_lpt);
void command_bases(char *lname, int use_lpt);
void command_extract(char *lname);
extract_option extract_dialog(char *fname, char *nom_feature, char *bornes, 
		int *use_min_bornes, char *min_bornes);
void line_to_stdout(char *line);
void line_to_file(char *line);
void command_select(char *lname);
void command_info(char *lname, int use_lpt);
void command_codes(int use_lpt);
void command_save(char *lname, char *given_fname, int save_accnums, int save_gcg);
static char *calc_access_name(void);
void command_keywords(int use_lpt);
int count_spec_key_seqs(int num, int is_spec);
int command_key_desc(int num, FILE *out, int prof, int *total, int *blist);
void command_find(char *lname);
void command_species(int use_lpt);
static void next_der_num_liste(void);
static int prep_new_list(char *nomliste);
void init(void);
extern char *prefixacnuc(char *basefile, char *pathname);
static void cre_new_list(char *nom, int numlist);
void ghelp(char *topic);
static void prep_names_line(int num, char **addr_debut);
static int listseqdef(int num, FILE *out);
int findliste(char *lname, char *topic, int direct, int seq_list_only);
int secrit(int todo);
void lire(int todo);
void desc_arbre(int num);
void charge_arbre(void);
noeud *cre_noeud(int numero);
void ajout_branche(int pere, int fils);
void ajout_synonyme(int secondaire, int principal);
void rech_sp_floue(char *chaine, int maxprof, FILE *out);
int calc_arbre(noeud *pracine, int maxprof, FILE *out, int *deja_vu);
int calc_sous_arbre(noeud *pracine, int currx, int prof, int maxprof, 
	FILE *out, int *deja_vu);
int getcount(noeud *node);
int calc_chemin_vers_racine(noeud *pracine, int prof, int old_ret, 
	FILE *out, char *chemin_vers_racine);
void command_modify(char *lname, char *command);
char *get_nl_option(char *command);
int modif_by_length(int old_list, int newnumlist);
int modif_by_date(int numlist, int newnumlist);
static int getdateseq(void);
static int dateseuil(char *ligne);
int modif_subloc(int numlist, int newnumlist);
int modif_locsub(int numlist, int newnumlist);
static void mod_liste_descript(int numlist, char *prefix);
char *nom_genre(char genre);
int modif_by_edit(int numlist, int newnumlist);
int modif_by_scan(int oldnum, int numlist);
void show_scan_options(void);
void set_scan_options(void);
void command_terminal(void);
int get_nblpp(void);
void process_args(int argc, char *argv[]);
void process_mmap(int argc, char *argv[]);
int prepare_interrupt(void);
int check_interrupt(int int_fd);
void finish_interrupt(int int_fd);
char **calc_feat_list(int *p_total);
int compare_2_chaines(const void *c1, const void *c2);
int show_feat_list(void);


/* prototypes of external functions */
int proc_requete(char *requete, char *message, char *nomliste, int *numliste);
void *mycalloc(int nbr, size_t taille);
void pretty_seq(int numseq, int basesperline, char translate, 
			void (*outline)(char *line));
int accseq(char *nom, int *list);
void showannots(int iseq, Boolean *ctabchoix, char **record, int totrec, 
		int l_key, Boolean *nouveau_choix, void (*outoneline)(char *) );
int prep_extract(out_option option, char *fname, extract_option choix,
	char *feature_name, char *bornes, char *min_bornes, char **message);
int extract_1_seq(int numseq, char *bornes);
void fin_extract(void);
int prepch(char *chaine, char **posmot);
int compch(char *cible, int lcible, char **posmot, int nbrmots);
int isenum_flou(char *nom, int *blist, int *locus);
void subloc(int *source, int *destin);
void locsub(int *source, int *destin);
void descen(DIR_FILE *kan, int recnum, int *blist);
char *check_acnuc_gcgacnuc(void);
void sel_seqs_1_node(DIR_FILE *kan, int num, int *seq_list, int hote);


/* external variables */
extern int acc_length_on_disk;


int main(int argc, char *argv[])
{
char command[200], *lname, *p;
command_kind command_id;

#ifdef __INTEL__
atexit(attente);
#endif

#ifdef unix
process_args(argc, argv);
#endif
if( (p = check_acnuc_gcgacnuc()) != NULL ) {
	fprintf(stderr, "Error: %s!\n", p);
	exit(ERREUR);
	}
dir_acnucopen("RO");
#ifdef unix
process_mmap(argc, argv);
#endif
init();

while(TRUE)	{
	if(nblpp != -1) nblpp = get_nblpp();
	lire(3);
	printf("\nCommand? (or H for help)\n");
	gets(command);
	if(strlen(command) == 0) continue;
	command_id = identify_command(command);
	if(command_id == unknown) continue;
	if(command_id == stop) break;
	lname = proc_l_option(command);
	process_command(command_id, lname, command);
	}
dir_acnucclose();
puts("End of ACNUC retrieval program");
return 0;
}


char *proc_l_option(char *command)
{
static char lname[40], *p, *q;

p = strstr(command, "/l=");
if(p == NULL) p = strstr(command, "/L=");
if(p == NULL) return NULL;
p += 2; while(*p == ' ') p++; q = lname;
while( *(++p) != 0 && *p != '/') *(q++) = toupper(*p);
*q = 0;
return lname;
}


command_kind identify_command(char *command)
{
char *p;
command_kind num, retval;
int l, trouve;

p = command - 1; l = trouve = 0;
while(*(++p) != 0 && *p != '/' && *p != ' ') {*p = tolower(*p); l++;}
for(num = 1; num < TOT_COMMS; num++)
	if(strncmp(command, command_list[num], l ) == 0) {
		trouve++;
		retval = num;
		}
if(trouve == 1)
	return retval;
else if (trouve == 0) 
	printf("Unknown command: %s\n", command);
else 
	printf("Ambiguous command: %s\n", command);
return unknown;
}


void process_command(command_kind command_id, char *lname, char *command)
{
int use_lpt;
FILE *out;

use_lpt = ( strstr(command, "/lpt") != NULL || 
		strstr(command, "/LPT") != NULL );

if(command_id == com_name_select) {
	command_select(lname);
	}
else if(command_id == lists) {
	command_lists(use_lpt);
	}
else if(command_id == help) {
	ghelp("MENU");
	}
else if(command_id == terminal) {
	command_terminal();
	}
else if(command_id == general) {
	ghelp("CONT");
	if( secrit(1) ) return;
	putchar('\n');
	ghelp("GENERAL");
	}
else if(command_id == com_name_short) {
	command_short(lname, use_lpt);
	}
else if(command_id == names ) {
	command_names(lname, use_lpt);
	}
else if(command_id == bases ) {
	command_bases(lname, use_lpt);
	}
else if(command_id == com_name_delete ) {
	command_delete(lname, command);
	}
else if(command_id == info ) {
	command_info(lname, use_lpt);
	}
else if(command_id == extract ) {
	command_extract(lname);
	}
else if(command_id == codes ) {
	command_codes(use_lpt);
	}
else if(command_id == save ) {
	char *fname;
	fname = get_f_option(command);
	majuscules(command);
	command_save(lname, fname, strstr(command, "/AC") != NULL, 
		 strstr(command, "/GCG") != NULL );
	}
else if(command_id == keywords) {
	command_keywords(use_lpt);
	}
else if(command_id == find) {
	command_find(lname);
	}
else if(command_id == species ) {
	command_species(use_lpt);
	}
else if(command_id == modify ) {
	command_modify(lname, command);
	}
}


void command_lists(int use_lpt)
{
int num;
FILE *out;

out = prep_out_file( &use_lpt, command_list[lists], NULL );
if( (!use_lpt) && secrit(1) ) return;
for(num = last_hidden_list + 1; num < tlist; num++) 
	if( defoccup[num]) break;
if(num >= tlist) {
	if( (!use_lpt) && secrit(1) ) return;
	fputs("No list currently defined\n", out);
	}
else	{
	if( (!use_lpt) && secrit(1) ) return;
	fputs("Name,  kind,   size, and origin of currently defined lists\n", 
		out);
	for(num = last_hidden_list + 1; num < tlist; num++) {
		if( !defoccup[num]) continue;
		if( (!use_lpt) && secrit(1) ) break;
		fprintf(out, "%.80s\n", nom_liste[num]);
		}
	}
if(use_lpt) fclose(out);
}


void command_delete(char *lname, char *command)
{
int liste, all;
char reponse[10];

all = ( strstr(command, "/all") != NULL || strstr(command, "/ALL") != NULL );
if(all) {
	lire(1);
	printf("Confirm deletion of all lists? (y/[n]) ");
	gets(reponse); majuscules(reponse);
	if(reponse[0] != 'Y' || strlen(reponse) != 1) return;
	for(liste = last_hidden_list + 1; liste < tlist; liste++) {
		if(defoccup[liste]) {
			free(deflnames[liste]);
			defoccup[liste] = FALSE;
			}
		}
	default_list[0] = 0;
	der_num_liste = 0;
	return;
	}
liste = findliste(lname, "DELETE", FALSE, FALSE);
if(liste == 0) return;
defoccup[liste] = FALSE;
if(strcmp(deflnames[liste], default_list) == 0) default_list[0] = 0;
free(deflnames[liste]);
}


void command_short(char *lname, int use_lpt)
{
int liste, seq, *blist;
FILE *out;

liste = findliste(lname, "SHORT", TRUE, TRUE);
if(liste == 0) return;
out = prep_out_file( &use_lpt, command_list[com_name_short], deflnames[liste] );
seq = 1; blist = defbitlist + liste * lenw;
while( ( seq = irbit(blist, seq, nseq) ) != 0) 
	if( listseqdef(seq, out) ) break;
if(use_lpt) fclose(out);
}


FILE *prep_out_file(int *use_lpt, char *command_name, char *list_name)
{
FILE *log_file;
time_t heure;
char *p;

if( !*use_lpt ) return stdout;
log_file = fopen(LOG_FILE_NAME, "a");
if(log_file == NULL) {
	secrit(1);
	printf("Cannot open output file: %s\n", LOG_FILE_NAME);
	*use_lpt = FALSE;
	return stdout;
	}
secrit(1);
printf("Command output sent to file %s\n", LOG_FILE_NAME);
fprintf(log_file, 
	"\n\n========== output of command %s", command_name);
if(list_name != NULL) 
	fprintf(log_file," on list %s", list_name);
time(&heure);
p = asctime(localtime(&heure));
p[strlen(p) - 1] = 0;
fprintf(log_file," at %s ==========\n", p);
return log_file;
}


void command_select(char *lname)
{
	char lname2[11], requete[200], errmessage[200];
	int numlist, l, erreur, numfilles;
	
	lire(2);
     	printf("Enter your selection criteria, or H(elp) (EX: sp=homo sapiens et k=globin@)\n");
     	requete[0] = 0; 
     	gets(requete);
     	if(strlen(requete) == 0) return;
     	if(strcmp(requete, "H") == 0 || strcmp(requete, "h") == 0) {
     		lire(1);
     		printf("Do you want help for simple(S) or elaborate(E) usage? ");
     		gets(requete); majuscules(requete);
     		if(requete[0] =='E')
     			ghelp("SELE");
     		else
     			ghelp("SELS");
     		return;
     		}
	if(lname == NULL) {
		next_der_num_liste();
		sprintf(lname2, "LIST%d", der_num_liste);
		lname = lname2;
		}
     	majuscules(lname);
	erreur = proc_requete(requete,errmessage,lname,&numlist);
	if(erreur == 0) {
		l = strlen(requete); if(l > 79) l = 79;
		memcpy(deflisdef[numlist],requete,l); deflisdef[numlist][l]=0;
		strcpy(default_list, lname);
		secrit(1);
		printf("List %s contains %d %s ", 
			lname, defllen[numlist], nom_genre(defgenre[numlist]) );
		if(defgenre[numlist] == 'S' && !deflocus[numlist]) {
			/* liste des filles */
			non(defbitlist + 2 * lenw, defbitlist + lenw, lenw);
			/* filles de la liste */
			et(defbitlist + 2 * lenw, defbitlist + numlist * lenw,
				defbitlist + 2 * lenw, lenw);
			numfilles = bcount(defbitlist + 2 * lenw, lenbit);
			if(numfilles == 0) deflocus[numlist] = TRUE;
			else 	printf(" (among which %d subsequences)", 
					numfilles);
			}
		putchar('\n');
		cre_new_list(lname, numlist);
		}
	else	{
		der_num_liste--;
		secrit(1);
		puts(errmessage);
		}
}


void command_info(char *lname, int use_lpt)
{
int liste, seq, largeur, l, *blist;
Boolean nouveau_choix;
char reponse[30], *p;
FILE *out;
line_out_function out_function;

liste = findliste(lname, "INFO", TRUE, TRUE);
if(liste == 0) return;
if(embl || swissprot) largeur = 2;
else	largeur = 3;
l = 0;
puts("Choose one or several of the following topics, finish with two <RETURN>s");
for(seq = 0; seq < tot_options; seq++) {
	l += strlen(poptions[seq]);
	if(l > 80) { putchar('\n'); l = strlen(poptions[seq]); }
	fputs(poptions[seq], stdout);
	info_options[seq] = FALSE;
	}
putchar('\n');
while(TRUE) {
	lire(1);
	fputs("? ", stdout); fflush(stdout);
	gets(reponse);
	if(strlen(reponse) == 0) break;
	majuscules(reponse);
	for(seq = 0; seq < tot_options; seq++) {
		p = poptions_maj[seq]; while(*p == ' ') p++;
		if(strncmp(reponse, p, largeur) == 0) {
			info_options[seq] = TRUE;
			break;
			}
		}
	if(seq >= tot_options) {
		secrit(1);
		printf("Unknown topic: %s\n", reponse);
		}
	}
out = prep_out_file( &use_lpt, command_list[info], deflnames[liste] );
seq = 1; blist = defbitlist + liste * lenw;
nouveau_choix = TRUE;
line_interrupt = FALSE;
if(use_lpt) {
	line_file = out;
	out_function = line_to_file;
	}
else	{
	out_function = line_to_stdout;
	}
while( ( seq = irbit(blist, seq, nseq) ) != 0) {
	showannots(seq, info_options, poptions_maj, tot_options, wid_key,
		&nouveau_choix, out_function);
	if(line_interrupt) break;
	}	
if(use_lpt) fclose(out);
}


char *get_f_option(char *command)
{
char *reponse, *p;
static char fname[100];

reponse = strstr(command, "/F=");
if(reponse == NULL) reponse = strstr(command, "/f=");
if(reponse == NULL) return NULL;
reponse += 2; p = fname;
while( *(++reponse) != 0 && *reponse != '/') *(p++) = *reponse;
*p = 0;
return fname;
}


void command_names(char *lname, int use_lpt)
{
int liste, seq, count, *blist;
double total;
char ligne[90], *p;
FILE *out;

liste = findliste(lname, "NAMES", TRUE, FALSE);
if(liste == 0) return;
out = prep_out_file( &use_lpt, command_list[names], deflnames[liste] );
blist = defbitlist + liste * lenw;
total = 0;
if(defgenre[liste] == 'S') {
	if( out == stdout && secrit(3) ) return;
	fprintf(out, "\nSequence names and lengths in list %s (%d seqs)\n\n",
		deflnames[liste], defllen[liste] );
	count=0; p = ligne;
	seq = 1;
	while( ( seq = irbit(blist, seq, nseq) ) != 0) {
		readsub(seq);
		total += psub->length;
		count++;
		if(count==4 ) { 
			if(out == stdout && secrit(1)) return;
			fputs(ligne, out); fputc('\n', out);
			p=ligne;
			count=1;
			}
		prep_names_line(seq, &p);
		}
	if( out == stdout && secrit(2) ) return;
	fprintf(out, "%s\nCumulated sequence length: %.0f %s\n", ligne, total,
		( nbrf || swissprot ? "aminoacids" : "bases") );
	}
else	{
	DIR_FILE *kan;
	if(defgenre[liste] == 'K') kan = kkey;
	else kan = kspec;
	seq = 2;
	while( ( seq = irbit(blist, seq, maxa) ) != 0) {
		dir_read(kan, seq, 1, pspec);
		if( out == stdout && secrit(1) ) break;
		fprintf(out, "%.*s\n", WIDTH_KS, pspec->name);
		}
	}
if(use_lpt) fclose(out);
}


void command_bases(char *lname, int use_lpt)
{
int liste, seq, first, *blist, nuc_db;
char reponse[5];
FILE *out;
	
liste = findliste(lname, "BASES", TRUE, TRUE);
if(liste == 0) return; 
nuc_db = !(nbrf || swissprot);
out = prep_out_file( &use_lpt, command_list[bases], deflnames[liste] );
blist = defbitlist + liste * lenw;
if(use_lpt) {
	seq = 1; line_file = out;
	while( ( seq = irbit(blist, seq, nseq) ) != 0) {
		readsub(seq);
		if(secrit(1)) break;
		printf("Listing %.*s\n", L_MNEMO, psub->name);
		pretty_seq(seq, 100, nuc_db, line_to_file);
		}
	fclose(out);
	}
else	{
	seq = 1; first = TRUE;
	while( ( seq = irbit(blist, seq, nseq) ) != 0) {
		if(first) first = FALSE;
		else 	{
			readsub(seq);
			lire(1);
			printf("Continue with %.*s ? ([y]/n) ", 
				L_MNEMO, psub->name);
			gets(reponse);
			if(reponse[0] == 'n' || reponse[0] == 'N') break;
			}
		line_interrupt = FALSE;
		pretty_seq(seq, 50, nuc_db, line_to_stdout);
		}
	}
}


void command_extract(char *lname)
{
char fname[100], nom_feature[20], bornes[50], min_bornes[50], *message;
int liste, seq, erreur, total, extracted, use_min_bornes, cds_type, *blist;

liste = findliste(lname, "EXTRACT", TRUE, TRUE);
if(liste == 0) return;
extract_options_choice = extract_dialog(fname, nom_feature, bornes, 
	&use_min_bornes, min_bornes);
if(extract_options_choice == 0) return;
erreur = prep_extract(extract_format_choice, fname, extract_options_choice, nom_feature,
	bornes, (use_min_bornes ? min_bornes : NULL), &message);
if(erreur) {
	secrit(1);
	puts(message);
	return;
	}
if(extract_options_choice == translate) cds_type = fcode(ksmj, "04CDS ", 6);
seq = 1; total = 0; blist = defbitlist + liste * lenw;
while( ( seq = irbit(blist, seq, nseq) ) != 0) {
	readsub(seq);
	if( secrit(1) ) break;
	if(extract_options_choice == translate) {
		if(psub->type != cds_type) {
			printf("Sequence %.*s is not protein coding. "
				"It will not be extracted\n", 
				L_MNEMO, psub->name);
   			continue;
     			}
		printf("Translating and ");
		}
	printf("Extracting %.*s\n", L_MNEMO, psub->name);
	extracted = extract_1_seq(seq, bornes);
	if(extracted == -1) {
		secrit(1);
		printf("Problem while writing to file %s\n", fname);
		break;
		}
	total += extracted;
	}
fin_extract();
secrit(1);
printf("%d extracted sequences.\n", total);
}


extract_option extract_dialog(char *fname, char *nom_feature, char *bornes, 
		int *use_min_bornes, char *min_bornes)
{
char reponse[5];
int choice;

lire(1);
printf("Current output format is: %s. Do you want to change it? (y/[n]) ",
	extract_format_name[extract_format_choice] );
gets(reponse); majuscules(reponse);
if(reponse[0] =='Y') {
	lire( (flat - acnuc) + 3 );
	puts("Choose format among");
	for(choice = acnuc; choice <= flat; choice++)
		printf("%d:   %s\n", choice, extract_format_name[choice]);
	fputs("? ", stdout); gets(reponse);
	choice = acnuc -1;
	sscanf(reponse , "%d", &choice);
	if(choice < acnuc || choice > flat) {
		secrit(1);
		puts("Bad choice. Please repeat command.");
		return 0;
		}
	extract_format_choice = choice;
	}
if(extract_format_choice != gcg) {
	lire(1);
	printf("Name of output file? "); gets(fname);
	if(strlen(fname) == 0) return 0;
	if(extract_format_choice != analseq) {
		FILE *in;
		in = fopen(fname, "r");
		if(in != NULL) {
			fclose(in);
			lire(1);
			puts("This file already exists. Do you want to append to it? (y/[n]) ");
			gets(reponse); majuscules(reponse);
			if(reponse[0] != 'Y') {
				secrit(1);
				puts("Please repeat command.");
				return 0;
				}
			}
		}
	}
if( !nbrf ) {
	lire( (swissprot ? 6 : 7) );
	puts("Do you want:\n"
	"  (1) Simple extraction");
	if(!swissprot) puts("  (2) Translate into protein and extract");
	puts("  (3) Fragments or adjacent sequences\n"
	"  (4) Regions defined by sequence FEATURES\n"
	"  (5) Regions adjacent to sequence FEATURES");
	fputs("? ", stdout); gets(reponse); choice = simple - 1;
	sscanf(reponse, "%d", &choice);
	if(choice < simple || choice > region || 
		(swissprot && (choice == translate) ) ) {
		secrit(1);
		puts("Bad choice. Please repeat command.");
		return 0;
		}
	}
else	choice = simple;
if(choice == feature || choice == region) {
	lire(1);
	printf("Feature name? (or H for list of feature names)  "); 
	gets(nom_feature); majuscules(nom_feature);
	if(nom_feature[0] == 'H' && strlen(nom_feature) == 1) {
		if( !show_feat_list() ) return 0;
		lire(1);
		printf("Feature name?  "); 
		gets(nom_feature); majuscules(nom_feature);
		}
	if(strlen(nom_feature) == 0) return 0;
	if(iknum(nom_feature, kkey) == 0) {
		secrit(1);
		printf("%s is not a valid feature name\n", nom_feature);
		return 0;
		}
	}
if(choice == fragment || choice == region) {
	lire(2);
	puts("Describe fragment: (ex: 132,1600  -10,10  e-20,e+10  -20,e+5) [default= 1,e]");
	gets(bornes); 
	if(bornes[0] == 0) strcpy(bornes, "1,e");
	if(strchr(bornes,'-') != NULL || strstr(bornes,"e+") != NULL || 
			choice == region) {
		*use_min_bornes = TRUE;
		lire(2);
		printf("Describe the minimum fragment size: "
			"(same syntax as before) [default= %s]\n", bornes);
		gets(min_bornes); majuscules(min_bornes);
		if(min_bornes[0] == 0) strcpy(min_bornes, bornes);
		}
	else *use_min_bornes = FALSE;
	majuscules(bornes);
	}
return choice;
}


void line_to_stdout(char *line)
{
if(line_interrupt) return;
line_interrupt = secrit(1);
if( !line_interrupt ) puts(line);
}


void line_to_file(char *line)
{
fputs(line, line_file); fputc('\n', line_file);
}


void command_codes(int use_lpt)
{
char code[3], *p, ligne[81], reponse[50], choix_dyn[10];
int option, last, num, l_start_j;
FILE *out;
static char choix[] = "SMJYTO";
static char choix_name[6][10] = {"status", "molecule", "journal", "year", 
	"type", "organelle"};
static char nuclear[2][70] = { "NUCLEAR             Nuclear genome",
                               "NON-ORGANELLAR      Non-organellar expression"};

choix_dyn[0] = 0;
lire(2);
puts("Do you want to see codes of: (or H for help)");
fputs("Journal(j)  Year(y)  ", stdout);
strcat(choix_dyn, "JY");
if(genbank || embl) {
	fputs("Type(t)  ", stdout);
	strcat(choix_dyn, "T");
	}
fputs("Organelle(o)  ", stdout);
strcat(choix_dyn, "O");
if(genbank || embl) {
	fputs("Molecule(m)  ", stdout);
	strcat(choix_dyn, "M");
	}
if(swissprot) {
	fputs("Status(s)  ", stdout);
	strcat(choix_dyn, "S");
	}
fputs("?  ", stdout);
gets(reponse); majuscules(reponse);
if(reponse[0] == 0) return;
if(reponse[0] == 'H') {
	ghelp("CODES");
	return;
	}
p = strchr(choix_dyn, reponse[0]);
if(p == NULL) {
	secrit(1);
	puts("Bad choice. Please repeat command.");
	return;
	}
p = strchr(choix, reponse[0]);
option = p - choix;
if(option == 2) {
	lire(1);
	fputs("Start of journal name? (or <RETURN> for all journals) ", stdout);
	gets(reponse); 
	majuscules(reponse); compact(reponse);
	l_start_j = strlen(reponse);
	}
sprintf(code,"%2.2d",option);
#ifdef vms
/* bug dans printf de vms */
if(*code==' ') *code='0';
#endif
last = read_first_rec(ksmj,NULL);
out = prep_out_file( &use_lpt, command_list[codes], NULL );
if( out == stdout && secrit(3) ) return;
fprintf(out, "List of codes to be used for selection by: %s\n\n"
"Code                Description\n", choix_name[option]);
for(num=2; num<=last; num++) {
	readsmj(num);
	if( strncmp(code,psmj->name,2) )continue;
	if(option == 2 && l_start_j > 0 && 
		strncmp(psmj->name + 2, reponse, l_start_j) < 0) continue;
	memset(ligne,' ',80);
	ligne[80]=0;
	memcpy(ligne,psmj->name+2,18);
	if(psmj->libel) {
		readtxt(psmj->libel);
		memcpy(ligne+20,ptxt,60);
		}
	if(out == stdout && secrit(1)) return;
	fputs(ligne, out); fputc('\n', out);
	}
if(option == 5) {
	if(swissprot || nbrf) 
		p = nuclear[1];
	else
		p = nuclear[0];
	if(out == stdout && secrit(1)) return;
	fputs(p, out); fputc('\n', out);
	}
if(use_lpt) fclose(out);
}


void command_save(char *lname, char *given_fname, int save_accnums, int save_gcg)
{
int liste, num, *blist;
char fname[90], *p, reponse[90], gcg_prefix[10];
FILE *fich;
time_t heure;
DIR_FILE *kan;

liste = findliste(lname, "SAVE", FALSE, FALSE);
if(liste == 0) return;
if(given_fname != NULL)
	strcpy(fname, given_fname);
else	{
	strcpy(fname, deflnames[liste]);
	p = fname - 1; while(*(++p) != 0) *p = tolower(*p);
	if(defgenre[liste] == 'S')	{
		if(save_accnums)
			strcat(fname,".acc");
		else
			strcat(fname,".mne");
		}
	else if(defgenre[liste] == 'K')
		strcat(fname,".key");
	else if(defgenre[liste] == 'E')
		strcat(fname,".spe");
	lire(1);
	printf("Name of file to write list content? [default= %s] ", fname);
	gets(reponse);
	if(reponse[0] != 0) strcpy(fname, reponse);
	}
fich = fopen(fname,"w");
if(fich == NULL) {
	secrit(1);
	printf("File %s cannot be opened.\n",fname);
	return;
	}
if(defgenre[liste] =='S' && save_gcg) {
	time(&heure);
	fprintf(fich, "Saved from acnuc on %s..\n",
		asctime(localtime(&heure)) );
	p = getenv("acnucgcgprefix");
	if(p != NULL) {
		strcpy(gcg_prefix, p);
		strcat(gcg_prefix, ":");
		}
	else if(nbrf)
		strcpy(gcg_prefix, "p:");
	else if(swissprot)
		strcpy(gcg_prefix, "sw:");
	else if(embl)
		strcpy(gcg_prefix, "embl:");
	else
		strcpy(gcg_prefix, "gb:");
	}
num=1; blist = defbitlist + liste * lenw;
if(defgenre[liste] == 'S') {
	while( (num=irbit(blist,num,nseq)) ) {
		readsub(num);
		if(save_accnums) {
			p = calc_access_name();
			fputs(p,fich); putc('\n',fich);
			if(ferror(fich)) break;
			}
		else 	{
			if(save_gcg) 
				fputs(gcg_prefix, fich);
			fprintf(fich,"%.*s\n", L_MNEMO, psub->name);
			if(ferror(fich)) break;
			}
		}
	}
else	{
	if(defgenre[liste] == 'K') kan=kkey;
	else kan=kspec;
	while( (num=irbit(blist,num,maxa)) ) {
		dir_read(kan,num,1,pspec);
		fprintf(fich,"%.*s\n", WIDTH_KS, pspec->name);
		if(ferror(fich)) break;
		}
	}
if(ferror(fich)) {
	secrit(1);
	printf("Error while writing to save file %s\n", fname);
	}
fclose(fich);
}


static char *calc_access_name(void)
/* contexte: psub charge avec seq utilisee
rend acc # de la seq si mere + extension si fille
*/
{
static char nom[20];
char *ext;
int lext=0;

if( psub->pext > 0) {
	char *p;
	ext=strchr(psub->name,'.');
	p=ext;
	while( ( (p - psub->name) < 16 ) && ( *p != ' ' ) ) p++;
	lext=(p-ext);
	memcpy(nom+ACC_LENGTH,ext,lext);
	readext(psub->pext);
	readsub(pext->mere);
	}
readloc(psub->plinf);
if(ploc->placc) {
	readshrt(ploc->placc);
	readacc(pshrt->val);
	memcpy(nom,pacc->name,ACC_LENGTH);
	}
else
	memset(nom,'?',ACC_LENGTH);
nom[ACC_LENGTH+lext]=0;
compact(nom);
return nom;
}


void command_keywords(int use_lpt)
{
char chaine[50], cible[50], *p, *q, *posmot[40];
int nbrmots, num, v_num, total, totk, *blist;
FILE *out;

lire(3);
puts("(Give part(s) of word(s) you want be present in keywords, ex: RIB TEIN)\n"
"keyword searched? (or ALL for all keywords or H for help)");
gets(chaine);
majuscules(chaine);
if(chaine[0] =='H' && strlen(chaine) == 1) {
	ghelp("KEYW");
	return;
	}
if( (int)strlen(chaine) >= sizeof(cible) ) return;
blist = (int *)calloc(longa , sizeof(int));
if(blist == NULL) {
	secrit(1);
	puts("Not enough memory.");
	return;
	}
out = prep_out_file( &use_lpt, command_list[keywords], NULL );
if(use_lpt) fprintf(out, "Keyword template used: %s\n\n", chaine);
if(strcmp(chaine, "ALL") == 0) strcpy(chaine, " ");
p=chaine; cible[0]='@'; q=cible+1;
do	{
	if(*p==' ')
		*q='@';
	else
		*q = *p;
	p++; q++;
	}
while(*p!=0);
*q=0; strcat(cible, "@");
nbrmots = prepch(cible,posmot);
totk = read_first_rec(kkey,NULL);
total = 0;
for(num = 3; num <= totk; num++) {
	if( testbit(blist, num) ) continue;
	readkey(num);
	if( !compch(pkey->name, 40, posmot, nbrmots) ) continue;
	if(pkey->name[0] == 'x') continue;
	v_num = num;
	while(pkey->syno > 0) {
		v_num = pkey->syno;
		readkey(v_num);
		}
	if( command_key_desc(v_num, out, 1, &total, blist) ) break;
	}
if(total == 0) {
	if(out == stdout) secrit(1);
	fputs("No name match the specified word(s)\n", out);
	}
free(blist);
if( use_lpt ) fclose(out);
}


int count_spec_key_seqs(int num, int is_spec)
{
static int first = TRUE, *blist;

if(first) {
	blist = malloc(lenw * sizeof(int));
	if(blist == NULL) return 0;
	first = FALSE;
	}
memset(blist, 0, lenw * sizeof(int));
if(is_spec)
	sel_seqs_1_node(kspec, num, blist, 0);
else
	sel_seqs_1_node(kkey, num, blist, 0);
return bcount(blist, nseq);
}


int command_key_desc(int num, FILE *out, int prof, int *total, int *blist)
{
int desc, i;
char *p;
static char name[sizeof(pkey->name) + 1];

bit1(blist, num);
readkey(num);
(*total)++;
p = pkey->name + sizeof(pkey->name); *p = 0;
while( p > pkey->name && *(--p) == ' ' ) *p = 0;
strcpy(name, pkey->name);
if(out == stdout && secrit( 1 ) ) return TRUE;
for(i = 0; i < 2 * (prof - 1); i++) fputc(' ', out);
fprintf(out, "%s ==> %d seqs.\n", name, count_spec_key_seqs(num, FALSE) );
readkey(num); /* count may have changed pkey */
if(pkey->syno < 0) { /* traitement des synonymes */
	int syno;
	while( (syno = abs(pkey->syno) ) != num ) {
		readkey(syno);
		p = pkey->name + sizeof(pkey->name); *p = 0;
		while( p > pkey->name && *(--p) == ' ' ) *p = 0;
		if(out == stdout && secrit( 1 ) ) return TRUE;
		for(i = 0; i < 2 * (prof - 1); i++) fputc(' ', out);
		fprintf(out, "synonym: %s\n", pkey->name);
		bit1(blist, syno);
		}
	readkey(num);
	}
if(pkey->libel != 0) {
	if(out == stdout && secrit( 1 ) ) return TRUE;
	readtxt(pkey->libel);
	for(i = 0; i < 2 * (prof - 1); i++) fputc(' ', out);
	p = ptxt + lrtxt;
	while( p > ptxt && *(--p) == ' ' ) *p = 0;
	fprintf(out, "{%.60s}\n", ptxt);
	}
readshrt(pkey->desc);
desc = pshrt->next;
while(desc != 0) {
	readshrt(desc); desc = pshrt->next;
	readshrt(pshrt->val);
	if( command_key_desc( abs(pshrt->val), out, prof + 1, total, blist) ) 
		return TRUE;
	}
return FALSE;
}


void command_find(char *given_lname)
{
int nlist, *blist, num, find_command_mode;
static char fs[]="%.*s", fek[]="%.*s";
char nom[81], *p, lname[11];

/* initialisation de la liste */
nlist = prep_new_list(lname);
if(nlist == -1) {
	secrit(1);
	puts("Not enough free memory. Please delete some lists.");
	return;
	}
secrit(2);
puts("Enter sequence or species or keyword names or accession #s, one per line.\n"
"Finish with two <RETURN>s.   Enter H for help.");
blist = defbitlist + nlist*lenw;
memset(blist, 0, lenw * sizeof(int));
defgenre[nlist] = 'S';
deflocus[nlist] = TRUE;
find_command_mode = 2;
while(TRUE) {
lire(1);
fputs("? ", stdout); gets(nom); majuscules(nom);
if(strlen(nom) == 0) {
	char *format; int width;
	defllen[nlist] = bcount(blist,lenbit);
	if(defllen[nlist] == 0) {
		defoccup[nlist] = FALSE;
		der_num_liste--;
		return;
		}
	num=irbit(blist,1,lenbit);
	if(defgenre[nlist] == 'S') {
		readsub(num);
		p = psub->name; format = fs; width = L_MNEMO;
		}
	else if(defgenre[nlist] == 'E') {
		readspec(num);
		p = pspec->name; format = fek; width = WIDTH_KS;
		}
	else 	{
		readkey(num);
		p = pkey->name; format = fek; width = WIDTH_KS;
		}
	sprintf(deflisdef[nlist], format, width, p);
	p = deflisdef[nlist] + strlen(deflisdef[nlist]) - 1;
	while(*p==' ') p--;
	strcpy(p+1," etc...");
	cre_new_list(lname, nlist);
	if(given_lname != NULL) {
		majuscules(given_lname);
		free(deflnames[nlist]);
		deflnames[nlist] = (char *)mycalloc(strlen(given_lname)+1,1);
		strcpy(deflnames[nlist], given_lname);
		}
	cre_new_list(deflnames[nlist], nlist);
	strcpy(default_list, deflnames[nlist]);
	secrit(1);
	printf("Number of %s selected in list %s: %d\n", 
		nom_genre(defgenre[nlist]), deflnames[nlist], defllen[nlist]);
	return;
	}
if(find_command_mode == 2) { 
/* 1ere lecture d'un elet, identif du genre de la liste */
	if(nom[0] == 'H' && strlen(nom) == 1) {
		defoccup[nlist] = FALSE;
		der_num_liste--;
		ghelp("FIND");
		return;
		}
	num=isenum_flou(nom,blist,&deflocus[nlist]);
	if(num) {
		find_command_mode=3;  /* genre liste de sequences */
		}
	else {
	   num=fcode(kacc,nom,ACC_LENGTH);
	   if(num) {
		find_command_mode=3;
		num=pacc->plsub;
		while(num){
			readshrt(num); num=pshrt->next;
			bit1(blist,pshrt->val);
			}
		num= -1;
		}
	   else {
		num=iknum(nom,kspec);
		if(num) {
		   find_command_mode=4; /* genre liste d'especes */
		   defgenre[nlist]='E';
		   }
		else {
		   num=iknum(nom,kkey);
		   if(num) {
			find_command_mode=5; /* genre liste de mots-cles */
			defgenre[nlist]='K';
			}
		   }
		}
	   }
	}
else if(find_command_mode == 3) { /* seq list */
	num=isenum_flou(nom,blist,&deflocus[nlist]);
	if(num==0) {
		num=fcode(kacc,nom,ACC_LENGTH);
		if(num)	{
			num=pacc->plsub;
			while(num){
				readshrt(num); num=pshrt->next;
				bit1(blist,pshrt->val);
				}
			num= -1;
			}
		}
	}
else if(find_command_mode == 4) { /* species list */
	num=iknum(nom,kspec);
	}
else if(find_command_mode == 5) { /* keyword list */
	num=iknum(nom,kkey);
	}
if(find_command_mode == 3) { /* cas sequence */
/* a-t-on trouve une(des) sequence(s) vraie(s) ou supprimee(s)? */
	if(num < 0) {
		et(blist, blist, defbitlist, lenw);
		if(irbit(blist, 1, nseq) == 0) num = 0;
		}
	else if(num > 0) {
		if( !testbit(defbitlist, num) ) num = 0;
		}
	}
if(num == 0) { /* name not found */
	secrit(1);
	fputs("This ", stdout);
	if(find_command_mode == 4)
		fputs("species name", stdout);
	else if(find_command_mode == 5)
		fputs("keyword name", stdout);
	else
		fputs("sequence or accession number", stdout);
	printf(" does not exist: %s\n", nom);
	}
else	{
	if(num>0) bit1(blist,num);
	}
} /* end of while */
}


static void next_der_num_liste(void)
{
char nomliste[11];
int num;

do	{
	der_num_liste++;
	sprintf(nomliste, "LIST%d", der_num_liste);
	for(num = 0; num < tlist; num++) {
		if(defoccup[num] && strcmp(deflnames[num], nomliste) == 0) 
			break;
		}
	}
while(num < tlist);
}


static int prep_new_list(char *nomliste)
{
int num;

for(num = 0; num < tlist; num++)
	if(!defoccup[num]) break;
if(num >= tlist) return -1; /* pas trouvee */
next_der_num_liste();
sprintf(nomliste, "LIST%d", der_num_liste);
defoccup[num] = TRUE;
deflnames[num] = (char *)mycalloc(strlen(nomliste) + 1,1);
strcpy(deflnames[num], nomliste);
majuscules(deflnames[num]);
return num;
}


void init(void)
{
int i;

#ifdef unix
fprintf(stderr,"Opening a %s database in %d divisions\n",
			(flat_format ? "flat" : "gcg" ) , divisions+1);
#endif
nblpp = get_nblpp(); /* nbre par defaut de lignes par ecran */
numlig = 0;
/* banniere de reconnaissance de banque */
ghelp("CONT");
defbitlist = (int *)mycalloc(tlist*lenw,sizeof(int));
for(i=0; i<MAXLISTES; i++) {
	defoccup[i]=FALSE;
	deflnames[i]=NULL;
	}
lngbit(3,defbitlist);  /* liste des valides */
for(i=0; i< lenw; i++) defbitlist[i]= ~defbitlist[i];
bit0(defbitlist,1);
quick_list_meres(defbitlist+lenw); /* liste des meres */
defoccup[0] = defoccup[1] = defoccup[2] = TRUE;
deflnames[0]=(char *)mycalloc(12,1);
strcpy(deflnames[0],"!VALIDSEQS!");
deflnames[1]=(char *)mycalloc(12,1);
strcpy(deflnames[1],"!MERES!");
deflnames[2]=(char *)mycalloc(50,1); /* liste pour usage direct par nom ou acc # */
strcpy(deflnames[2],"!DIRECT!");
last_hidden_list = 2;
default_list[0] = 0;
der_num_liste = 0;
if(embl) {
	wid_key = 3;
	annot_options = annot_options_embl;
	}
else if(swissprot) {
	wid_key = 3;
	annot_options = annot_options_swissprot;
	}
else if(nbrf) {
	wid_key = 10;
	annot_options = annot_options_nbrf;
	}
else  	{ /* GenBank */
	wid_key = 11;
	annot_options = annot_options_gbk;
	}
tot_options = strlen(annot_options) / wid_key;
for(i=0; i<tot_options; i++) {
	poptions[i] = malloc(wid_key+1);
	poptions_maj[i] = malloc(wid_key+1);
	memcpy(poptions[i],annot_options+i*wid_key,wid_key);
	poptions[i][wid_key]=0;
	strcpy(poptions_maj[i], poptions[i]);
	majuscules(poptions_maj[i]);
	}
printf("[%d free lists available]\n", tlist - last_hidden_list - 1);
multi_spec = (swissprot || nbrf);
}


extern char *prefixacnuc(char *basefile, char *pathname)
{
char *p;

p = prepare_env_var(
#ifndef vms
"acnuc"
#else
"BANK"
#endif
);
if(p==NULL) {
#if defined(unix) || defined(vms)
#ifdef vms
	fprintf(stderr,"Error logical name BANK is not defined or accessible\n");
#else
	fprintf(stderr,"Error environment variable acnuc is not defined\n");
#endif
	exit(ERREUR);
#endif
	pathname[0] = 0;
	}
else
	strcpy(pathname,p);
strcat(pathname,basefile);
return pathname;
}


static void cre_new_list(char *nom, int numlist)
{
char indic[3];

strcpy(indic,"  ");
if(defgenre[numlist]=='E')
	strcpy(indic,"Sp");
else if(defgenre[numlist]=='K')
	strcpy(indic,"Kw");
else if( ! deflocus[numlist])
	strcpy(indic,"* ");
sprintf(nom_liste[numlist],"%-10.10s%2.2s%7d ",nom,indic,defllen[numlist]);
memcpy(nom_liste[numlist]+20, deflisdef[numlist], 80);
}


void ghelp(char *topic)
/* GIVES THE HELP RELATED TO THE COMMAND IN ARGUMENT topic */
{
char line[100], *p;
FILE *fhelp;
int nl, num, trouve;

trouve = FALSE;
prefixacnuc("HELP",line);
fhelp = fopen(line,"r");
if(fhelp != NULL) {
	do	{
		p = fgets(line, sizeof(line), fhelp);
		if( p == NULL) break;
		p = strchr(line, '=');
		sscanf(p + 1, "%d", &nl);
		if(strncmp(line, topic, 4) != 0) {
			for(num = 0; num < nl; num++)
				fgets(line, sizeof(line), fhelp);
			continue;
			}
		for(num = 0; num < nl; num++) {
			trouve = TRUE;
			fgets(line, sizeof(line), fhelp);
			if( secrit(1) ) break;
			fputs(line, stdout);
			}
		}
	while(!trouve);
	fclose(fhelp);
	}
if(! trouve) {
	secrit(1);
	printf("Sorry, no help available for this command: %s\n", topic);
	}
}


static void prep_names_line(int num, char **addr_debut)
/* chargement d'une ligne avec nom_seq.....long_seq */
{
char *p, aux[10];
int lname, lsize, i;

p= *addr_debut;
readsub(num);
sprintf(p,"%.*s", L_MNEMO, psub->name);
i=16; while(p[--i]==' ') p[i]=0;
sprintf(aux,"%d",psub->length);
lname=strlen(p);
lsize=strlen(aux);
for(i=0; i<23-(lname+lsize); i++) p[lname+i]='.';
memcpy(p+23-lsize,aux,lsize);
strcpy(p+23,"  ");
p+=25;
*addr_debut=p;
}


static int listseqdef(int num, FILE *out)
/* charge lignes avec la breve descript de la seq ou sous-seq num.
rend TRUE si interrompu par secrit, FALSE sinon */
{
long pinf;
int div, l, lmnemo;
char *p, longseq[15], ligne[100], l80[81];
int doit;

seq_to_annots(num, &pinf, &div);
readsub(num);
if(psub->pext > 0) goto seqfille;
lmnemo = 0; /* calcul long du nom */
while(psub->name[lmnemo] != ' ' && lmnemo < L_MNEMO) lmnemo++;
if( read_annots(pinf, div) == NULL) return FALSE;
/* pour mere, recherche de la definition */
if(embl || swissprot) 
	while(strncmp(pinfo->line,"DE",2) != 0) next_annots(NULL);
else /* genbank ou nbrf */
	next_annots(NULL);
if(pinfo->line[0] != ' ') /* sauter 1er mot si non espace */
	p=strchr(pinfo->line,' ');
else
	p=pinfo->line;
while(*p == ' ') p++;
l=strlen(p);
while(p[l-1]==' ') l--;
p[l]=0;

/* nom + debut definition */
memcpy(ligne,psub->name,lmnemo);
ligne[lmnemo]=' ';
memcpy(ligne+lmnemo+1,p,l+1);
while (TRUE) { /* lignes suivantes de definition */
	next_annots(NULL);
	if(strncmp(pinfo->line,"  ",2) && strncmp(pinfo->line,"DE",2) ) break;
	if(embl || swissprot) memset(pinfo->line, ' ', 2);
	if( out == stdout && secrit(1) )return TRUE;
	fputs(ligne, out); fputc('\n', out);
	strcpy(ligne, pinfo->line);
	}
/* "=long bp" en fin de derniere ligne */
memset(l80,' ',80);
sprintf(longseq,"=%d ",psub->length);
if ( !( nbrf || swissprot) )
	strcat(longseq,"bp");
else
	strcat(longseq,"aa");
l=strlen(ligne); if(l>80) l=80;
memcpy(l80,ligne,l);
l=strlen(longseq);
memcpy(l80+80-l,longseq,l); l80[80]=0;
if( out == stdout && secrit(1) )return TRUE;;
fputs(l80, out); fputc('\n', out);
return FALSE;

seqfille: /* pour les sous-seqs */
/* nom + longueur */
sprintf(ligne,"%16.16s (%d ", psub->name, psub->length);
if(nbrf || swissprot)
	strcat(ligne, "aas)");
else
	strcat(ligne, "bases)");
if( out == stdout && secrit(1) )return TRUE;
fputs(ligne, out); fputc('\n', out);
/* traitement 1ere ligne de features */
if( read_annots(pinf, div) == NULL) return FALSE;
if( out == stdout && secrit(1) )return TRUE;
fputs(pinfo->line, out); fputc('\n', out);
doit = TRUE;
while (TRUE) { /* les autres lignes */
	p = next_annots(NULL);
	if(p == NULL) break;
	if(*p == 0) continue;
	if(p[5] != ' ' || ( strncmp(p,"  ",2) &&
		strncmp(p,"FT",2) ) ) break;
	if(strstr(p,"/translation") != NULL) /* sauter /translation */
		doit=FALSE;
	else if(!doit && strchr(p,'/') != NULL)
		doit=TRUE;
	if (doit) {
		if( out == stdout && secrit(1) )return TRUE;
		fputs(p, out); fputc('\n', out);
		}
	}
return FALSE;
}


int findliste(char *lname, char *topic, int direct, int seq_list_only)
/* rend numero de la liste dont le nom est dans lname
si pas trouve rend 0
si lname == NULL a l'appel, demande le nom de la liste
si direct == TRUE, le nom demande peut aussi etre une nom de seq ou un acc #
si seq_list_only == TRUE, refuse les listes d'esp ou de keyw
*/
{
int num, search, trouve;
char reponse[30];

if(lname == NULL) {
	lire(1);
	printf("List name, "); 
	if(direct)printf("sequence or accession #");
	printf(" or H(elp) ?");
	if(default_list[0] != 0)
		printf("  [%s] ", default_list);
	gets(reponse);
	if(strcmp(reponse, "H") == 0 || strcmp(reponse, "h") == 0) {
		ghelp(topic);
		return 0;
		}
	else if(strlen(reponse) == 0) {
		if(default_list[0] == 0) return 0;
		strcpy(reponse, default_list);
		}
	lname = reponse;
	}
majuscules(lname);
trouve = FALSE;
for(num = last_hidden_list + 1; num < tlist; num++) { /* chercher si c'est une liste */
	if(!defoccup[num]) continue;
	if(strcmp(lname, deflnames[num]) != 0) continue;
	break;
	}
if(num < tlist) {
	if(seq_list_only && defgenre[num] != 'S') {
		secrit(1);
		puts("This command can only be applied to a list of sequences.");
		return 0;
		}
	strcpy(default_list, lname);
	return num;
	}
if(direct) { /* chercher si c'est une sequence ou un acc # */
	num = isenum(lname);
	if(num != 0 && !testbit(defbitlist, num) ) num = 0;
	trouve = (num != 0);
	if(trouve) { /* trouve' une sequence */
		memset(defbitlist + 2 * lenw, 0, lenw * sizeof(int));
		bit1(defbitlist + 2 * lenw, num);
		defgenre[2] = 'S';
		defllen[2] = 1;
		deflocus[2] = (psub->pext <= 0);
		}
	search = !trouve; /* chercher un acc # */
	if(search) {
		compact(lname);
		if( (int)strlen(lname) > acc_length_on_disk) search = FALSE;
		}
	if(search) {
		num = -1;
		sscanf(lname + 2, "%d", &num);
		if(num == -1) search = FALSE;
		}
	if(search) {
		num = accseq(lname, defbitlist + 2 * lenw);
		if(num == 1) { /* enforce valid seqs only */
			et(defbitlist + 2 * lenw, defbitlist + 2 * lenw, defbitlist, lenw);
			if(irbit(defbitlist + 2 * lenw, 1, nseq) != 0) trouve = TRUE;
			}
		if(trouve) {
			defgenre[2] = 'S';
			defllen[2] = bcount(defbitlist + 2 * lenw, lenbit);
			deflocus[2] = TRUE;
			}	
		}
	}
if(trouve) {
	strcpy(default_list, lname);
	strcpy(deflnames[2], lname);
	return 2;
	}
secrit(1);
printf("Unknown list");
if(direct) printf(", sequence or accession #");
printf(": %s\n", lname);
return 0;
}


int secrit(int todo)
{
char reponse[5];

if(nblpp == -1) return FALSE;
if(numlig + todo + 2 > nblpp) {
	puts("<RETURN> key to continue or ST to stop ");
	gets(reponse); majuscules(reponse);
	nblpp = get_nblpp();
	numlig = 0;
	if(strncmp(reponse, "ST", 2) == 0) return TRUE;
	}
numlig += todo;
return FALSE;
}


void lire(int todo)
{
numlig = todo;
}


void desc_arbre(int num)
{
int l, numdesc, point, next, count;

next = num;
do	{ /* lecture et memorisation du nom et de ses synonmymes */
	char *p;
	readspec(next);
	l = sizeof(pspec->name) - 1;
	while(pspec->name[l] == ' ') l--;
	l++;
	p = (char *)mycalloc(l+1,sizeof(char));
	memcpy(p, pspec->name, l); p[l] = 0;
	if(tab_noeud[next] == NULL) tab_noeud[next] = cre_noeud(next);
	tab_noeud[next]->nom = p;
	if(next != num) ajout_synonyme(next,num);
	next = abs(pspec->syno);
	}
while(next != 0 && next != num);
if(next != 0) {
	ajout_synonyme(num, num); /* fermer la chaine des synonymes */
	readspec(num);
	}
if( ! multi_spec) count = 0;
else count = -1;
if( (! multi_spec) && pspec->plsub != 0) {
	point = pspec->plsub;
	while(point != 0) {
		readlng(point);
		for(l = 0; l < SUBINLNG; l++)
			if(plng->sub[l] != 0) count++;
		point = plng->next;
		}
	}
if(pspec->libel!=0) {
	readtxt(pspec->libel);
	l = lrtxt - 1;
	while(l >= 0 && ptxt[l] == ' ')l--;
	if(l >= 0 && ptxt[l] != '|') { /* sauter les libels avec gc seulement */
		l++;
		tab_noeud[num]->libel = (char *)mycalloc(l+1, sizeof(char));
		tab_noeud[num]->libel_upcase = 
			(char *)mycalloc(l+1, sizeof(char));
		memcpy(tab_noeud[num]->libel,ptxt,l); 
		tab_noeud[num]->libel[l] = 0;
		memcpy(tab_noeud[num]->libel_upcase, tab_noeud[num]->libel, 
			l+1);
		majuscules(tab_noeud[num]->libel_upcase);
		}
	}
point = pspec->desc;
readshrt(point);
point = pshrt->next;
while(point != 0) {
	readshrt(point); next = pshrt->next;
	readshrt(pshrt->val);
	numdesc = abs(pshrt->val);
/* ne pas brancher un noeud deja branche ailleurs auparavant */
	if(tab_noeud[numdesc] == NULL || tab_noeud[numdesc]->parent == NULL) {
		ajout_branche(num, numdesc);
		desc_arbre(numdesc);
		}
	if( ! multi_spec) count += tab_noeud[numdesc]->count;
	point = next;
	}
tab_noeud[num]->count = count;
}


void charge_arbre(void)
{
int totspec;

secrit(2);
fputs("Loading species tree.", stdout); fflush(stdout);
totspec = read_first_rec(kspec, NULL);
tab_noeud = mycalloc(totspec + 1, sizeof(noeud *));
desc_arbre(2);
tab_noeud[2]->nom = nom_racine;
puts("done");
}


noeud *cre_noeud(int numero)
{
noeud *pnoeud;
pnoeud = (noeud *)mycalloc(1, sizeof(noeud));
pnoeud->numero = numero;
return pnoeud;
}


void ajout_branche(int pere, int fils)
{
noeud *pnoeud;
struct paire *point, *nouveau;
static int pre_pere=0;
static struct paire *pre_point;
static int compteur = 0;

if(tab_noeud[fils] == NULL) tab_noeud[fils] = cre_noeud(fils);
pnoeud = tab_noeud[fils];
pnoeud->parent = tab_noeud[pere];
nouveau = (struct paire *)mycalloc(1,sizeof(struct paire));
nouveau->valeur = pnoeud;
if( (point = tab_noeud[pere]->liste_desc) == NULL) {
	tab_noeud[pere]->liste_desc = nouveau;
	}
else	{
	if(pere == pre_pere)point = pre_point;
	while(point->next != NULL) point = point->next;
	point->next = nouveau;
	}
pre_pere = pere;
pre_point = nouveau;
if( (++compteur) % 1000 == 0) {
	putchar('.'); fflush(stdout);
	}
}

void ajout_synonyme(int secondaire, int principal)
{
noeud *point = tab_noeud[principal];
while(point->syno != NULL) {
	point = point->syno;
/* eviter repetition et bouclage infini si boucle deja fermee */
	if(point == tab_noeud[secondaire]  || point == tab_noeud[principal]) return;
	}
point->syno = tab_noeud[secondaire];
}


void command_species(int use_lpt)
{
static int first = TRUE;
FILE *out;
char reponse[50], rep_prof[10];
int maxprof;

lire(3);
puts("(Give part(s) of word(s) you want be present in species names ex: MAMM)\n"
	"species searched? (or ALL for all species  or H for help)");
gets(reponse); majuscules(reponse);
if(strlen(reponse) == 0) return;
if(reponse[0] == 'H' && strlen(reponse) == 1) {
	ghelp("SPECIES");
	return;
	}
if(strcmp(reponse, "ALL") == 0) strcpy(reponse, nom_racine);
lire(1);
printf("Max depth for tree display? (2 to get a general overview) ");
gets(rep_prof);
maxprof = 2;
sscanf(rep_prof, "%d", &maxprof);
if(maxprof <= 0) maxprof = 2;
if(first) {
	first = FALSE;
	charge_arbre();
	}
out = prep_out_file( &use_lpt, command_list[species], NULL );
if(use_lpt) fprintf(out, "Species template: %s  Display depth: %d\n",
			reponse, maxprof);
rech_sp_floue(reponse, maxprof, out);
if(use_lpt) fclose(out);
}


void rech_sp_floue(char *chaine, int maxprof, FILE *out)
{
char cible[50], *p, *q, *posmot[20], *libel_spec;
int nbrmots, num, trouve, match, totspec, *deja_vu;
noeud *pnoeud;

if(strcmp(chaine, nom_racine) == 0) {
	calc_arbre(tab_noeud[2], maxprof, out, NULL);
	return;
	}

if( (int)strlen(chaine) >= 50) return;
p = chaine; cible[0] = '@'; q = cible+1;
do	{
	if(*p == ' ')
		*q = '@';
	else
		*q = *p;
	p++; q++;
	}
while(*p!=0);
*q=0; strcat(cible,"@");
majuscules(cible);
nbrmots = prepch(cible,posmot);
for(num = last_hidden_list + 1; num < tlist; num++)
	if(!defoccup[num]) break;
if(num >= tlist) {
	fputs("Not enough memory available. Please delete some lists.\n", out);
	return;
	}
deja_vu = defbitlist + num * lenw;
memset(deja_vu, 0, lenw * sizeof(int));
trouve = 0;
totspec = read_first_rec(kspec, NULL);
for(num = 3; num <= totspec; num++) {
	if( testbit(deja_vu, num) ) continue;
	if(tab_noeud[num] == NULL) continue;
	p = tab_noeud[num]->nom;
	if(p == NULL) continue;
	libel_spec = tab_noeud[num]->libel_upcase;
	match = compch(p, strlen(p), posmot, nbrmots);
	if(!match && libel_spec != NULL) 
		match = compch(libel_spec, strlen(libel_spec), posmot, nbrmots);
	if( !match ) 
		continue;
	trouve++;
	pnoeud = tab_noeud[num];
	while( pnoeud->syno != NULL && pnoeud->parent == NULL ) 
		pnoeud = pnoeud->syno;
	if( out == stdout && secrit(1) ) break;
	fputc('\n', out);
	if( calc_arbre(tab_noeud[pnoeud->numero], maxprof, out, deja_vu) ) 
		break;
	}
if(trouve == 0) {
	if(out == stdout) secrit(1);
	fputs("No matching sequence name found\n", out);
	}
}


int calc_arbre(noeud *pracine, int maxprof, FILE *out, int *deja_vu)
{
int depart;
char *debut;
static char chemin_vers_racine[81];
const int browser_width = 80;

if( pracine->numero == 2 ) { 
	/* exploration depuis la racine <a cacher> */
	struct paire *desc;
	desc = pracine->liste_desc;
	while(desc != NULL) {
		if( calc_sous_arbre(desc->valeur, 0, 1, maxprof, out, NULL) )
			return TRUE;
		desc = desc->next;
		}
	return FALSE;	
	}
else if(pracine->parent->numero == 2) depart = 1;
else	{
	strcpy(chemin_vers_racine,"[Place in tree: ");
	depart = calc_chemin_vers_racine(pracine->parent, 0, 0, 
					out, chemin_vers_racine) + 1;
	strcat(chemin_vers_racine,"]");
	if( out == stdout && secrit(1) ) return TRUE;
	fputs(chemin_vers_racine, out); fputc('\n', out);
	}
return calc_sous_arbre(pracine, 0, depart, maxprof + depart, out, deja_vu);
}


int calc_sous_arbre(noeud *pracine, int currx, int prof, int maxprof, 
	FILE *out, int *deja_vu)
{
static int l;
struct paire *suivant;
static noeud *next;
static char *ligne = NULL;
static int lligne = 0;
const int deltax = 2;

l = strlen(pracine->nom);
if(l + 4 + currx + 1 + 25 > lligne) {
	if(lligne != 0) free(ligne);
	lligne = l + 4 + currx + 1 + 25 + 80;
	ligne = (char *)malloc(lligne);
	}
memset(ligne,' ',currx);
sprintf(ligne + currx, "%2d/ ", prof);
memcpy(ligne + currx + 4, pracine->nom, l + 1);
sprintf(ligne + currx + 4 + l, " ==> %d seqs.", getcount(pracine));
if(deja_vu != NULL) bit1(deja_vu, pracine->numero);
if( out == stdout && secrit(1) ) return TRUE;
fputs(ligne, out); fputc('\n', out);
if( (next=pracine->syno) != NULL) { /* ecrire les synonymes */
	while(next!=pracine) {
		l=strlen(next->nom);
		if(l + currx + 1 + 9 > lligne) {
			if(lligne != 0) free(ligne);
			lligne = l + currx + 1 + 9;
			ligne = (char *)malloc(lligne);
			}
		memset(ligne,' ',currx);
		memcpy(ligne + currx, "Synonym: ", 9);
		memcpy(ligne + currx + 9, next->nom, l + 1);
		if(deja_vu != NULL) bit1(deja_vu, next->numero);
		if( out == stdout && secrit(1) ) return TRUE;
		fputs(ligne, out); fputc('\n', out);
		next=next->syno;
		}
	}
if( pracine->libel!=NULL) {
	l=strlen(pracine->libel);
	if(l + currx + 1 + 2 > lligne) {
		if(lligne != 0) free(ligne);
		lligne = l + currx + 1 + 2;
		ligne = (char *)malloc(lligne);
		}
	memset(ligne, ' ', currx);
	ligne[currx] = '{';
	memcpy(ligne+currx+1, pracine->libel, l);
	ligne[currx +1 + l] = '}';
	ligne[currx +1 + l + 1] = 0;
	if( out == stdout && secrit(1) ) return TRUE;
	fputs(ligne, out); fputc('\n', out);
	}
suivant = pracine->liste_desc;
if(suivant == NULL || prof >= maxprof) return FALSE;
do	{
	if( calc_sous_arbre(suivant->valeur, currx+deltax, prof+1, maxprof, 
		out, deja_vu) ) return TRUE;
	suivant = suivant->next;
	}
while(suivant != NULL);
return FALSE;
}


int getcount(noeud *node)
{
if( node->count == -1 ) {
	node->count = count_spec_key_seqs(node->numero, TRUE);
	}
return node->count;
}


int calc_chemin_vers_racine(noeud *pracine, int prof, int old_ret, 
	FILE *out, char *chemin_vers_racine)
{
int l;
static int m_prof = 0;

if(pracine == NULL) {
	m_prof = 0;
	return old_ret;
	}
old_ret = calc_chemin_vers_racine(pracine->parent, prof + 1, old_ret, 
					out, chemin_vers_racine);
if(prof > m_prof) m_prof = prof;
if(strcmp(pracine->nom, nom_racine) != 0) {
	l = strlen(chemin_vers_racine);
	if(l + 5 + strlen(pracine->nom) > 80) {
		if(out == stdout && secrit(1)) return 0;
		fputs(chemin_vers_racine, out); fputc('\n', out);
		l = 0;
		}
	sprintf(chemin_vers_racine + l, " %d/ ", m_prof - prof);
	strcat(chemin_vers_racine, pracine->nom);
	if( old_ret <  m_prof - prof) old_ret = m_prof - prof;
	}
return old_ret;
}


void command_modify(char *lname, char *command)
{
char reponse[20], *fixed_name;
int liste, *blist, choix, new_list, done, old_der,
	protein_db, deja;

liste = findliste(lname, "MODIFY", FALSE, FALSE);
if(liste == 0) return;
old_der = der_num_liste;
new_list = prep_new_list(reponse);
if(new_list == -1) {
	secrit(1);
	puts("Cannot create new list. Please delete some lists.");
	return;
	}
done = TRUE; protein_db = (nbrf || swissprot);
if(defgenre[liste] != 'S') {
	choix = 1;
	}
else	{
	lire(4);
	puts("You can modify this sequence list according to:\n"
	     "1. Confirmation/Suppression of sequences from list\n"
	     "2. Sequence length\n"
	     "3. Sequence insertion date");
	if( !protein_db ) {
		lire(2);
		puts("4. Replace subsequences by seq containing them\n"
	     	     "5. Add subsequences of seq in list");
		}
	lire(2);
	puts("6. Scan annotations for a string");
	puts("7. Set scanning options + item 6");
	lire(1);
	fputs("? ", stdout); gets(reponse);
	choix = 0;
	sscanf(reponse, "%d", &choix);
	}
deja = tlist;
if( ( fixed_name = get_nl_option(command) )  != NULL && 
		strcmp(fixed_name, deflnames[new_list]) != 0 ) {
	for(deja = last_hidden_list + 1; deja < tlist; deja++) { 
		if( defoccup[deja] && 
			(strcmp(fixed_name, deflnames[deja]) == 0) ) break;
		}
	free(deflnames[new_list]);
	deflnames[new_list] = malloc(strlen(fixed_name) + 1);
	strcpy(deflnames[new_list], fixed_name);
	der_num_liste = old_der;
	}
if(choix < 1 || choix > 7 || (protein_db && (choix == 4 || choix == 5))) {
	puts("Bad choice. Please repeat command.");
	done = FALSE;
	}
else if(choix == 1) done = modif_by_edit(liste, new_list);
else if(choix == 2) done = modif_by_length(liste, new_list);
else if(choix == 3) done = modif_by_date(liste, new_list);
else if(choix == 4) done = modif_subloc(liste, new_list);
else if(choix == 5) done = modif_locsub(liste, new_list);
else if(choix == 6) {
	show_scan_options(); 
	done = modif_by_scan(liste, new_list);
	}
else if(choix == 7) { 
	set_scan_options(); 
	done = modif_by_scan(liste, new_list);
	}
if(done) {
	if(fixed_name != NULL && deja < tlist) {
		secrit(1);
		printf("Warning previously defined list %s was replaced\n", 
			fixed_name);
		defoccup[deja] = FALSE;
		free(deflnames[deja]);
		}
	secrit(1);
	printf("New list created %s contains %d %s\n", deflnames[new_list], 
		defllen[new_list], nom_genre(defgenre[new_list]) );
	strcpy(default_list, deflnames[new_list]);
	}
else 	{
	defoccup[new_list] = FALSE; 
	free(deflnames[new_list]);
	der_num_liste = old_der;
	}
}


char *get_nl_option(char *command)
{
char *reponse, *p;
static char fname[100];

majuscules(command);
reponse = strstr(command, "/NL=");
if(reponse == NULL) return NULL;
reponse += 3; p = fname;
while( *(++reponse) != 0 && *reponse != '/') *(p++) = *reponse;
*p = 0;
return fname;
}




int modif_by_length(int numlist, int newnumlist)
/* returns TRUE if new list was created, FALSE if not
*/
{
int *blist, *newblist, seuil, plus_grand, num;
char reponse[40], *p;

lire(1);
fputs("Give your length threshold: (ex:  L>200  or   L<1000) ", stdout);
gets(reponse);
p = strchr(reponse, '>');
if(p == NULL) p = strchr(reponse, '<');
plus_grand = (strchr(reponse, '>') != NULL);
seuil = 0;
if(p != NULL) sscanf(p + 1, "%d", &seuil);
if(p == NULL || seuil == 0) {
	secrit(1);
	puts("Bad reply. Please repeat command.");
	return FALSE;
	}
blist = defbitlist + numlist * lenw;
newblist = defbitlist + newnumlist * lenw;
defgenre[newnumlist] = 'S';
deflocus[newnumlist] = deflocus[numlist];
memcpy(newblist, blist, lenw * sizeof(int));
num = 1;
while( (num = irbit(newblist, num, nseq)) != 0 ) {
	readsub(num);
	if(plus_grand && psub->length < seuil )
		bit0(newblist, num);
	else if(!plus_grand && psub->length > seuil )
		bit0(newblist, num);
	}
memcpy(deflisdef[newnumlist], deflisdef[numlist], 80);
cre_new_list(deflnames[newnumlist], newnumlist);
sprintf(reponse, "Modify L%c%d of ", (plus_grand ? '>' : '<') , seuil);
mod_liste_descript(newnumlist, reponse);
return TRUE;
}


int modif_by_date(int numlist, int newnumlist)
/* returns TRUE if new list was created, FALSE if not
*/
{
int *blist, *newblist, seuil, dateseq, ok, num;
char reponse[50], *p, nom_seuil[30];
enum { inferieur = 0, egal, superieur } cas;
char symbols[] = "<=>";

ok = FALSE;
lire(2);
printf("Enter your date threshold: \n"
"(ex: > 1/jun/93 or < 31/dec/89 or < 1/jun/2000 or = 15/may/01)  ");
gets(reponse);
if( (p = strchr(reponse, '>')) != NULL )
	{ cas = superieur; ok = TRUE; }
else if( ( p = strchr(reponse, '<')) != NULL )
	{ cas = inferieur; ok = TRUE; }
else if( ( p = strchr(reponse, '=')) != NULL )
	{ cas = egal; ok = TRUE; }
if(ok) 	{
	compact(p + 1);
	strcpy(nom_seuil, p + 1);
	ok = ( (seuil = dateseuil(p+1)) != 0);
	}
if(!ok)	{
	secrit(1);
	puts("Bad reply. Please repeat command.");
	return FALSE;
	}
blist = defbitlist + numlist * lenw;
newblist = defbitlist + newnumlist * lenw;
defgenre[newnumlist] = 'S';
deflocus[newnumlist] = deflocus[numlist];
memcpy(newblist, blist, lenw * sizeof(int));
num = 1;
while( (num = irbit(newblist, num, nseq)) != 0 ) {
	readsub(num);
	if(psub->pext>0) {
		readext(psub->pext);
		readsub(pext->mere);
		}
	dateseq = getdateseq();
	if( cas == superieur && dateseq < seuil )
		bit0(newblist,num);
	else if( cas == inferieur  && dateseq > seuil )
		bit0(newblist,num);
	else if( cas == egal  && dateseq != seuil )
		bit0(newblist,num);
	}
memcpy(deflisdef[newnumlist], deflisdef[numlist], 80);
cre_new_list(deflnames[newnumlist], newnumlist);
sprintf(reponse, "Modify D%c%s of ", symbols[cas] , nom_seuil);
mod_liste_descript(newnumlist, reponse);
return TRUE;
}


static int getdateseq(void)
{
char *date;
int mois, jour, an;

readloc(psub->plinf);
date=ploc->date;
sscanf(date,"%2d",&mois);
sscanf(date+3,"%2d",&jour);
sscanf(date+6,"%2d",&an);
if(an < 50) an += 100;
return an*10000 + mois*100 + jour;
}


static int dateseuil(char *ligne)
{
static char noms_mois[]="JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC";
char *p, *q;
int jour, mois, an;

p=ligne;
while(*p !=0) {
	*p=toupper(*p);
	p++;
	}
jour = -1;
sscanf(ligne,"%d",&jour);
if( jour == -1 )return 0;
p=strchr(ligne,'/');
if(p==NULL) return 0;
ligne= ++p;
p=strchr(ligne,'/');
if(p==NULL) return 0;
*p=0;
q=strstr(noms_mois,ligne);
if(q==NULL) return 0;
mois=(q-noms_mois)/3+1;
p++;
an = -1;
sscanf(p,"%d",&an);
if( an == -1 )return 0;
if(an >= 100) an = an % 100; /* 19xx ou 20xx donnent xx */
if(an < 50) an += 100;
return an*10000 + mois*100 + jour;
}


int modif_subloc(int numlist, int newnumlist)
/* returns TRUE if new list was created, FALSE if not
*/
{
defgenre[newnumlist] = 'S';
deflocus[newnumlist] = TRUE;
subloc(defbitlist + numlist * lenw, defbitlist + newnumlist * lenw);
memcpy(deflisdef[newnumlist], deflisdef[numlist], 80);
cre_new_list(deflnames[newnumlist], newnumlist);
mod_liste_descript(newnumlist, "Replace subseqs. by seqs. of ");
return TRUE;
}


int modif_locsub(int numlist, int newnumlist)
/* returns TRUE if new list was created, FALSE if not
*/
{
defgenre[newnumlist] = 'S';
deflocus[newnumlist] = FALSE;
locsub(defbitlist + numlist * lenw, defbitlist + newnumlist * lenw);
memcpy(deflisdef[newnumlist], deflisdef[numlist], 80);
cre_new_list(deflnames[newnumlist], newnumlist);
mod_liste_descript(newnumlist, "Add subseqs. to seqs. of ");
return TRUE;
}


static void mod_liste_descript(int numlist, char *prefix)
{
int *blist, total, lprefix;
char reponse[10], *p;

blist = defbitlist + numlist * lenw;
total = bcount(blist, lenbit);
sprintf(reponse, "%7d", total);
memcpy(nom_liste[numlist] + 12, reponse, 7);
defllen[numlist] = total;
lprefix = strlen(prefix);
if(strncmp(nom_liste[numlist]+20, prefix, lprefix) != 0) {
	memcpy(nom_liste[numlist] + 20, prefix, lprefix);
	deflisdef[numlist][79] = 0;
	strcpy(nom_liste[numlist] + 20 + lprefix, deflisdef[numlist]);
	p = nom_liste[numlist] + 20;
	total = strlen(p); if(total > 80) total = 80;
	memset(deflisdef[numlist], ' ', 80);
	memcpy(deflisdef[numlist], p, total);
	}
}


char *nom_genre(char genre)
{
static char noms[3][10] = {"sequences", "species", "keywords"};
static char genres[] = "SEK";
int rank_genre;

rank_genre = strchr(genres, genre) - genres;
return noms[rank_genre];
}


int modif_by_edit(int numlist, int newnumlist)
/* returns TRUE if new list was created, FALSE if not
*/
{
int *blist, *newblist, num, fin, wid, ask;
char reponse[50], *p;
static char explain[] = 
	"C confirm, S suppress, F keep all others, Q suppress all others";

blist = defbitlist + numlist * lenw;
newblist = defbitlist + newnumlist * lenw;
defgenre[newnumlist] = defgenre[numlist];
deflocus[newnumlist] = deflocus[numlist];
memcpy(newblist, blist, lenw * sizeof(int));
secrit(2);
printf("Edition of list: %s\n", deflnames[numlist]);
puts(explain);	
if(defgenre[numlist] == 'S') { 
	num = 1; fin = nseq; 
	}
else 	{ 
	num = 2; fin = maxa; wid = sizeof(pspec->name); 
	}
ask = TRUE;
while( (num = irbit(newblist, num, fin)) != 0 ) {
	if(ask) {
		if(defgenre[numlist] == 'S') { 
			listseqdef(num, stdout);
			}
		else 	{
			if(defgenre[numlist] == 'E') { 
				readspec(num); p = pspec->name; 
				}
			else	{ readkey(num); p = pkey->name; }
			memcpy(reponse, p, wid); reponse[wid] = 0;
			p = reponse + wid; 
			while (p > reponse && *(--p) == ' ') *p = 0;
			fputs(reponse, stdout); 
			fputs(": ", stdout);
			}
		lire(1);
		fputs("[c,s,f,q]? ", stdout);
		gets(reponse); majuscules(reponse);
		}
	if(reponse[0] == 'F') break;
	else if(reponse[0] == 'C') continue;
	else if(reponse[0] == 'S') bit0(newblist,num);
	else if(reponse[0] == 'Q') {
		bit0(newblist,num);
		ask = FALSE;
		}
	else	{
		secrit(1);
		puts(explain);
		num--;
		}
	}
memcpy(deflisdef[newnumlist], deflisdef[numlist], 80);
cre_new_list(deflnames[newnumlist], newnumlist);
mod_liste_descript(newnumlist, "Edition of ");
return TRUE;
}


static int scan_all_annots = TRUE;
static Boolean *scan_options = NULL;
static char *scan_mot;
static int scan_found;
static int scan_definition;
static Boolean scan_nouveau_choix = TRUE;

void scan_one_line(char *line)
{
static int in_def = FALSE;
if(scan_found) return;
if(strncmp(line, "ID ", 3) == 0 || strncmp(line, "LOCUS", 5) == 0 ||
	strncmp(line, "ENTRY", 5) == 0 ) return;
if( !scan_definition ) { /* ignorer DEFINITION/DE/TITLE */
	if(embl) {
		if(strncmp(line, "DE", 2) == 0) return;
		}
	else if(genbank) {
		if(strncmp(line, "DEFINITION", 10) == 0) {
			in_def = TRUE;
			return;
			}
		if(in_def && *line == ' ' && *(line+3) == ' ') return;
		in_def = FALSE;
		}
	else if(nbrf) {
		if(strncmp(line, "TITLE", 5) == 0) {
			in_def = TRUE;
			return;
			}
		if(in_def && *line == ' ') return;
		in_def = FALSE;
		}
	}
majuscules(line);
if(strstr(line, scan_mot) != NULL) scan_found = TRUE;
}


static int scan_one_seq(int numseq, char *mot, char *marq_fin, int l_marq)
{
int fille, div;
long pannots;
char *p;

readsub(numseq); fille = (psub->pext > 0);
if( (!fille) && !scan_all_annots ) { /* scan selected fields */
	scan_found = FALSE;
	scan_mot = mot;
	showannots(numseq, scan_options, poptions_maj, tot_options, wid_key,
		&scan_nouveau_choix, scan_one_line);
	return scan_found;
	}
seq_to_annots(numseq, &pannots, &div);
if( read_annots(pannots, div) == NULL) return FALSE;
do	{ /* scan all annots, mere ou fille */
	p=pinfo->line-1; while(*(++p)) *p=toupper(*p);
	if(strstr(pinfo->line, mot) != NULL) return TRUE;
	next_annots(NULL);
	}
while( ( !fille && ( strncmp(pinfo->line, marq_fin, l_marq) != 0 ) ) ||
	(fille && ( strcmptrail(pinfo->line+2,10,NULL,0)==0 ) ) );
return FALSE;
}


void show_scan_options(void)
{
char def_rec[15];
int i;

secrit(2);
puts("Current annotation scanning options:");
if(scan_all_annots) {
	puts("all sequence annotations scanned");
	}
else	{
	if(nbrf) strcpy(def_rec, "TITLE");
	else if(embl || swissprot) strcpy(def_rec, "DE");
	else strcpy(def_rec, "DEFINITION");
	if(scan_options[0]) fputs(def_rec, stdout);
	for(i = 1; i < tot_options; i++) {
		if(scan_options[i]) printf(" %s", poptions[i]);
		}
	puts("");
	}
}


void set_scan_options(void)
{
int i, l, largeur;
char *p, reponse[20];
static char def_rec[15];

if(scan_options == NULL) {
	scan_options = (Boolean *)calloc(tot_options, sizeof(Boolean));
	if(scan_options == NULL) {
		scan_all_annots = TRUE;
		return;
		}
	if(nbrf) strcpy(def_rec, "TITLE");
	else if(embl || swissprot) strcpy(def_rec, "DE");
	else strcpy(def_rec, "DEFINITION");
	}
encore :
show_scan_options();
lire(1);
printf("Do you want to change options? (y/[n])  ");
gets(reponse);
if(toupper(reponse[0]) != 'Y') return;
if(embl || swissprot) largeur = 2;
else	largeur = 3;
puts("Choose one or several of the following topics, finish with two <RETURN>s");
printf("ALL  %s ", def_rec); l = strlen(def_rec) + 6;
for(i = 1; i < tot_options - 1; i++) {
	l += strlen(poptions[i]);
	if(l > 80) { putchar('\n'); l = strlen(poptions[i]); }
	fputs(poptions[i], stdout);
	scan_options[i] = FALSE;
	}
putchar('\n');
scan_all_annots = FALSE; scan_nouveau_choix = TRUE;
for(i = 0; i < tot_options - 1; i++) scan_options[i] = FALSE;
while(TRUE) {
	lire(1);
	fputs("? ", stdout); fflush(stdout);
	gets(reponse);
	if(strlen(reponse) == 0) break;
	majuscules(reponse);
	if(strcmp(reponse, "ALL") == 0) {
		scan_all_annots = TRUE;
		}
	else if(strncmp(reponse, def_rec, largeur) == 0) {
		scan_options[0] = TRUE;
		}
	else	{
		for(i = 1; i < tot_options - 1; i++) {
			p = poptions_maj[i]; while(*p == ' ') p++;
			if(strncmp(reponse, p, largeur) == 0) {
				scan_options[i] = TRUE;
				break;
				}
			}
		if(i >= tot_options - 1) {
			secrit(1);
			printf("Unknown topic: %s\n", reponse);
			}
		}
	}
goto encore;
}


int modif_by_scan(int oldnum, int numlist)
/* returns TRUE if new list was created, FALSE if not
*/
{
char mot[40], *p;
int err, numseq, point, *new_list, lmot, todo, pourcent,
	processed, fille, div, *current_liste, previous, int_fd;
long pannots;
static int first=TRUE, l_marq;
static char marq_fin[12];

if(first) {
	first=FALSE;
	if(embl || swissprot) strcpy(marq_fin, "SQ ");
	else if(nbrf) strcpy(marq_fin, "SUMMARY");
	else strcpy(marq_fin, "ORIGIN");
	l_marq = strlen(marq_fin);
	}
lire(1);
printf("Enter string to be searched: ");
gets(mot); majuscules(mot);
p = mot + strlen(mot); while( p > mot && *(--p) == ' ') *p = 0;
if(strlen(mot) == 0) return FALSE;
current_liste = defbitlist + oldnum * lenw;
new_list = defbitlist + numlist * lenw;
memset(new_list, 0, lenw * sizeof(int));
defgenre[numlist] = 'S';
deflocus[numlist] = TRUE;
todo = bcount(current_liste, nseq);
numseq = 1; processed = 0; previous = -1;
if(!scan_all_annots) {
	scan_definition = scan_options[0];
	scan_options[0] = FALSE;
	}
printf("Scanning sequence annotations for string %s\n", mot);
int_fd = prepare_interrupt();
while( (numseq = irbit(current_liste, numseq, nseq)) != 0 ) {
	if( check_interrupt(int_fd) ) break;
	processed++;
	if( scan_one_seq(numseq, mot, marq_fin, l_marq) ) {
		bit1(new_list, numseq);
		if ( !testbit(defbitlist+lenw, numseq) ) 
			deflocus[numlist] = FALSE;
		}
	pourcent = (processed * 100) / todo; 
	if(pourcent != previous && int_fd != -1) {
		printf("\r%d%%", pourcent); fflush(stdout);
		previous = pourcent;
		}
	}
finish_interrupt(int_fd);
putchar('\n');
if(!scan_all_annots) scan_options[0] = scan_definition;
sprintf(deflisdef[numlist], "scan of %s for \"", deflnames[oldnum]);
point = strlen(deflisdef[numlist]);
strcat(mot, "\""); lmot = strlen(mot);
lmot= ( lmot > 79-point ? 79-point : lmot ); mot[lmot] = 0;
sprintf(deflisdef[numlist] + point, "%s", mot);
defllen[numlist] = bcount(new_list, nseq);
defoccup[numlist] = TRUE;
cre_new_list(deflnames[numlist], numlist);
return TRUE;
}


void command_terminal(void)
{
char reponse[10];
int current;

current = get_nblpp();
if(current == -1) return; /* output is not to a terminal */
lire(1);
printf("Do you want output to stop when screen is full? [y, n, H(elp)] ");
gets(reponse); majuscules(reponse);
if(reponse[0] == 'H') {
	ghelp("TERMINAL");
	return;
	}
if(reponse[0] == 'N') {
	nblpp = -1;
	printf("\nOutput will not stop when screen is full.\n");
	}
else	{
	nblpp = current;
	printf("\nOutput will stop and program will ask\n"
	"<RETURN> key to continue or ST to stop\n"
	"when screen is full.\n");
	}
lire(0);
}


int get_nblpp(void)
{
int num, retval;
#ifdef unix
struct winsize ws_struct;

if( isatty(fileno(stdout)) == 0) return -1;
#endif

num = 24;
#ifdef unix
retval = ioctl(fileno(stdout), TIOCGWINSZ, &ws_struct);
if(retval != -1) num = ws_struct.ws_row;
#endif
return num;
}


#ifdef unix
void process_args(int argc, char *argv[])
{
int num;
char *p, *q;
static char env[2][100];

for(num = 1; num < argc; num++) {
	if(strncmp(argv[num], "-h", 2) == 0) {
		puts("Usage:\nquery  [db name]  [-mmap kname]...");
		exit(0);
		}
	if(strcmp(argv[num], "-mmap") == 0) {
		num++;
		continue;
		}
	p = getenv(argv[num]);
	if(p == NULL) {
		fprintf(stderr, "Bad argument: %s\n", argv[num]);
		exit(ERREUR);
		}
	q = strtok(p, " ");
	sprintf(env[0], "acnuc=%s", q);
	putenv(env[0]);
	q = strtok(NULL, " ");
	if(q != NULL) {
		sprintf(env[1], "gcgacnuc=%s", q);
		putenv(env[1]);
		}
	}
}


void process_mmap(int argc, char *argv[])
{
int num;

for(num = 1; num < argc; num++) {
	if(strcmp(argv[num], "-mmap") == 0) {
		num++;
		if(strcmp(argv[num], "kshrt") == 0) dir_set_mmap(kshrt);
		else if(strcmp(argv[num], "klng") == 0) dir_set_mmap(klng);
		else if(strcmp(argv[num], "ksub") == 0) dir_set_mmap(ksub);
		else if(strcmp(argv[num], "kloc") == 0) dir_set_mmap(kloc);
		}
	}
}
#endif


int prepare_interrupt(void)
{
char *tty;
int int_fd = -1;

#ifdef unix
tty = ttyname(fileno(stdout));
if( tty == NULL ) return -1;
int_fd = open(tty, O_RDONLY | O_NDELAY);
#endif
if( int_fd == -1 ) return -1;
puts("hit <return> to stop");
return int_fd;
}


int check_interrupt(int int_fd)
{
static char buf[10];
int lu = 0;

if(int_fd == -1) return 0;
#ifdef unix
lu = read(int_fd, buf, sizeof(buf));
#endif
return (lu > 0);
}


void finish_interrupt(int int_fd)
{
#ifdef unix
if(int_fd != -1) close(int_fd);
#endif
}


char **calc_feat_list(int *p_total)
{
int *features_blist = NULL, racine_features, num, l;
char *p;
static char **liste = NULL;
static int first = TRUE, total = 0;

if( !first ) {
	*p_total = total;
	return liste;
	}
first = FALSE;
racine_features = iknum("MISC_FEATURE", kkey);
if(racine_features == 0) {
	secrit(1);
	puts("No features in this data base.");
	return NULL;
	}
features_blist = (int *)malloc(longa * sizeof(int));
if(features_blist == NULL) goto no_mem;
descen(kkey, racine_features, features_blist);
num = bcount(features_blist, maxa);
liste = (char **)malloc(num * sizeof(char *));
if(liste == NULL) goto no_mem;
num = 2;
while( (num = irbit(features_blist, num, maxa) ) != 0) {
	readkey(num);
	l = sizeof(pkey->name);
	while(l > 0 && pkey->name[l - 1] == ' ') l--;
	p = malloc(l + 1);
	if(p == NULL) goto no_mem;
	memcpy(p, pkey->name, l); p[l] = 0;
	liste[total++] = p;
	}
free(features_blist);
qsort(liste, total, sizeof(char *), compare_2_chaines);
*p_total = total;
return liste;
no_mem:
	secrit(1);
	puts("Not enough memory available.");
	num = total - 1;
	while(num >= 0) free(liste[num--]);
	if(liste != NULL) free(liste);
	if(features_blist != NULL) free(features_blist);
	first = TRUE; liste = NULL; total = 0;
	return NULL;
}


int compare_2_chaines(const void *c1, const void *c2)
{
return strcmp( *(char **)c1, *(char **)c2 );
}


int show_feat_list(void)
{
int total, num, l, deja;
char **list;


list = calc_feat_list(&total);
if(list == NULL) {
	return FALSE;
	}
secrit(1);
puts("Sorted list of feature names:");
deja = 0;
for(num = 0; num < total; num++) {
	l = strlen(list[num]);
	if(deja + l + 3 > 80) { putchar('\n'); deja = 0; }
	fputs(list[num], stdout); fputs(";  ", stdout);
	deja += l + 3;
	}
putchar('\n');
return TRUE;
}

