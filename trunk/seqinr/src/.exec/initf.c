#include "dir_acnuc.h"
#include <ctype.h>

/* global variables */
int allow_punct = FALSE;
char div_name[20], b_kind[20];
int tot_smj, tot_txt;
struct rsub s_sub;
struct rloc s_loc;
struct rkey s_key;
struct rspec s_spec;
struct rbib s_bib;
struct racc s_acc;
struct rsmj s_smj;
struct rext s_ext;
struct raut s_aut;
struct rshrt s_shrt;
struct rlng s_lng;


/* extern variables (allocated in use_acnuc.c) */
extern int lracc_on_disk; /*taille records du fichier ACCESSION sur le disque */


void add_smj(char *name, char *libel);
void usage(void);
void *dialogue(void);
void proc_args(int argc, char **argv);


void *dialogue(void)
{
char rep[20];

printf("Use `initf -h` for description of the command line interface\n\n");

printf("GenBank(g), EMBL(e), SWISS-PROT(s), or NBRF/PIR(n) ?"); gets(b_kind);
*b_kind = tolower(*b_kind);
printf("Format flat(f) or GCG(g) ?"); gets(rep);
*rep = tolower(*rep);
flat_format = (*rep == 'f');
printf("Name of division? (without extension) "); gets(div_name);
printf("HSUB value [%d]? ", hsub); gets(rep); 
sscanf(rep, "%d", &hsub);
printf("HKWSP value [%d]? ", hkwsp); gets(rep); 
sscanf(rep, "%d", &hkwsp);
return;
}


void proc_args(int argc, char **argv)
{
int num, found = FALSE;
char *p;

if(strncmp(argv[1], "-h", 2) == 0) {
	usage();
	exit(0);
	}
if(strcmp(argv[1], "genbank") != 0 && strcmp(argv[1], "embl") != 0 && 
	strcmp(argv[1], "swissprot") != 0 && strcmp(argv[1], "nbrf") != 0 ) {
	fprintf(stderr, "Incorrect argument %s\n", argv[1]);
	usage();
	exit(ERREUR);
	}
strcpy(b_kind, argv[1]);

num = 1;
while(++num < argc) {
	if(strcmp(argv[num], "punctuation") == 0) allow_punct = TRUE;
	else if(strcmp(argv[num], "gcg") == 0) {
		flat_format = FALSE;
		}
	else if( (p = strchr(argv[num], '=') ) != NULL) {
		*p = 0;
		if(strcmp(argv[num], "hsub") != 0 && 
				strcmp(argv[num], "hkwsp") != 0 && 
				strcmp(argv[num], "div") != 0 ) {
			usage();
			exit(ERREUR);
			}
		if(strcmp(argv[num], "hsub") == 0) 
			sscanf(p + 1, "%d", &hsub);
		else if(strcmp(argv[num], "hkwsp") == 0) 
			sscanf(p + 1, "%d", &hkwsp);
		else	{
			strcpy(div_name, p + 1);
			found = TRUE;
			}
		}
	else	{
		usage();
		exit(ERREUR);
		}
	}
if(!found) {
	usage();
	exit(ERREUR);
	}
}


void usage(void)
{
fprintf(stderr, "Usage:\n"
"initf { genbank | embl | swissprot | nbrf | -h } [gcg] [punctuation]\n"
"	[hsub=x] [hkwsp=x] div=name\n"
);
}



int main(int argc, char **argv)
{
char aux[40];
int num, fin, tot_shrt, big_annots;

hsub = 9973; hkwsp = 1999;
flat_format = TRUE;
if(argc >= 2) proc_args(argc, argv);
else dialogue();

psub = &s_sub; ploc = &s_loc; pkey = &s_key; pspec = &s_spec;
pbib = &s_bib; pacc = &s_acc; psmj = &s_smj; pext = &s_ext; paut = &s_aut;
pshrt = &s_shrt; plng = &s_lng;

if( hsub % 2 == 0 ) hsub++;
if( hkwsp % 2 == 0 ) hkwsp++;

genbank = embl = swissprot = nbrf = FALSE;
big_annots = TRUE;
if(*b_kind == 'e')
	embl = TRUE;
else if(*b_kind == 's')
	swissprot = TRUE;
else if(*b_kind == 'n') {
	nbrf = TRUE; big_annots = FALSE;
	}
else
	genbank = TRUE;


ksub = dir_open("SUBSEQ", NULL, "w+", lrsub, 1000);
write_first_rec(ksub, 1, 0);
kloc = dir_open("LOCUS", NULL, "w+", lrloc, 1000);
write_first_rec(kloc, 1, 0);

kshrt = dir_open("SHORTL", NULL, "w+", lrshrt, 1000);
pshrt->val = - hsub; pshrt->next = - hkwsp;
writeshrt(2);
memset(pshrt, 0, lrshrt);
fin = 2 + (hsub+1)/2 + (hkwsp + 1);
for(num = 3; num <= fin; num++) writeshrt(num);
tot_shrt = fin;

klng = dir_open("LONGL", NULL, "w+", lrlng, 1000);
write_first_rec(klng, 3, 0);
memset(plng, 0, lrlng);
writelng(2);
plng->sub[0] = 1;
writelng(3);

if( !(swissprot || nbrf) ) {
	kext = dir_open("EXTRACT", NULL, "w+", lrext, 1000);
	write_first_rec(kext, 1, 0);
	}
ksmj = dir_open("SMJYT", NULL, "w+", lrsmj, 1000);
tot_smj = 1;
kaut = dir_open("AUTHOR", NULL, "w+", lraut, 1000);
write_first_rec(kaut, 1, 0);
kbib = dir_open("BIBLIO", NULL, "w+", lrbib, 1000);
write_first_rec(kbib, 1, 0);
kkey = dir_open("KEYWORDS", NULL, "w+", lrkey, 1000);
write_first_rec(kkey, 2, 0);
tot_shrt++;
pshrt->val = 2; pshrt->next = 0;
writeshrt(tot_shrt);
memset(pkey, 0, lrkey);
padtosize(pkey->name, "RACINE", 40);
pkey->desc = tot_shrt;
writekey(2);

kspec = dir_open("SPECIES", NULL, "w+", lrspec, 1000);
write_first_rec(kspec, 2, 0);
tot_shrt++;
pshrt->val = 2; pshrt->next = 0;
writeshrt(tot_shrt);
memset(pspec, 0, lrspec);
padtosize(pspec->name, "RACINE", 40);
pspec->desc = tot_shrt;
writespec(2);
write_first_rec(kshrt, tot_shrt, 0);

ktxt = dir_open("TEXT", NULL, "w+", lrtxt, 1000);
tot_txt = 1;

lracc_on_disk = 6 + 2*sizeof(int); /* total + SORTED + fin_sort */
if( lracc_on_disk < sizeof(struct racc) )
	lracc_on_disk = sizeof(struct racc);
#if defined(vms) || defined(__alpha) 
/* machines avec records multiples de 4 octets en F77 */
lracc_on_disk = ( (lracc_on_disk + 3)/4 ) * 4;
#endif
kacc = dir_open("ACCESS", NULL, "w+", lracc_on_disk, 1000);
write_first_rec(kacc, 1, 0);

/* identifiant dans SMJYT de chaque variante de banque */
if(swissprot) add_smj("01PRT", NULL);
else if(embl) add_smj("04ID", "EMBL sequence data library entry");
else if(nbrf) add_smj("01PROTEIN", NULL);
else add_smj("04LOCUS", "sequenced DNA fragment");

/* partie en 00 */
if(swissprot || embl) {
	add_smj("00STANDARD", "Fully annotated sequence");
	add_smj("00PRELIMINARY", "Partially annotated sequence");
	}
else if(genbank) {
	add_smj("00FULL", "Fully annotated sequence");
	add_smj("00SIMPLE", "Partially annotated sequence");
	add_smj("00UNANNOTATED", "unannotated sequence");
	}
if( embl ) {
	add_smj("00UNANNOTATED", "unannotated sequence");
	add_smj("00UNREVIEWED", "unreviewed sequence");
	}


/* partie en 01 */
if(genbank || embl) {
	add_smj("01DNA",      "Sequenced molecule is DNA");
	add_smj("01RNA",      "Sequenced molecule is RNA");
	if(genbank) {
		add_smj("01MRNA",     "Sequenced molecule is mRNA");
		add_smj("01RRNA",     "Sequenced molecule is rRNA");
		add_smj("01TRNA",     "Sequenced molecule is tRNA");
		add_smj("01URNA",     "Sequenced molecule is snRNA");
		}
	}

/* partie en 02 */
add_smj("02BOOK", "published in a book");
if(nbrf) {
	add_smj("02CITATION", "citation");
	}
else	{
	add_smj("02THESIS", "published in a thesis");
	add_smj("02PATENT", "patented sequence");
	}
if(genbank || embl) add_smj("02UNPUBLISHED", "unpublished");
add_smj("02NAR",       "Nucleic Acids Res.");
add_smj("02PNAS",      "Proc. Natl. Acad. Sci. U.S.A.");
add_smj("02JMB",       "J. Mol. Biol.");
add_smj("02MBE",       "Mol. Biol. Evol.");
add_smj("02JBC",       "J. Biol. Chem.");
add_smj("02JBACT",     "J. Bacteriol.");
add_smj("02JME",       "J. Mol. Evol.");
add_smj("02MPE",       "Mol. Phylogenet. Evol.");
add_smj("02BBA",       "Biochim. Biophys. Acta");
add_smj("02FEBS",      "FEBS lett.");
add_smj("02EMBOJ",     "EMBO J.");
add_smj("02JCB",       "J. Cell Biol.");
add_smj("02JGV",       "J. Gen. Virol.");
add_smj("02MCB",       "Mol. Cell. Biol.");
add_smj("02MGG",       "Mol. Gen. Genet.");

/* partie en 04 */
if(genbank || embl) {
	add_smj("04CDS",      ".PE protein coding region");
	add_smj("04TRNA",     ".TR transfer RNA coding region");
	add_smj("04RRNA",     ".RR ribosomal RNA coding region");
	add_smj("04MISC_RNA", ".RN other structural RNA coding region");
	add_smj("04SCRNA",    ".SC small cytoplasmic RNA coding region");
	add_smj("04SNRNA",    ".SN small nuclear RNA coding region");
	}

/* partie en 05 */
add_smj("05CHLOROPLAST",    "Chloroplast genome");
add_smj("05MITOCHONDRION",  "Mitochondrial genome");
if(genbank || embl) add_smj("05KINETOPLAST",    "Kinetoplast genome");
else if (swissprot)  add_smj("05CYANELLE",    "Cyanelle genome");

/* partie en 06 */
if(flat_format)
	strcpy(aux, "06FLT");
else
	strcpy(aux, "06GCG");
strcat(aux, div_name);
if(big_annots)
	add_smj(aux, "rank:0");
else
	add_smj(aux, "rank:0  size:0");

/* partie en 07 */
add_smj("07ACCESSION", "Length of accession numbers = 8");
if(big_annots) add_smj("07BIG_ANNOTS", "Annotation/Sequence addresses on 2*32 bits");
if(allow_punct) 
	add_smj("07ALLOW_PUNCTUATION", "Punctuation in sequence data");

dir_acnucflush();
printf("Normal end.\n");
return 0;
}


void add_smj(char *name, char *libel)
{
static char name2[21], libel2[61];

padtosize(name2, name, 20);
if(libel != NULL) padtosize(libel2, libel, 60);
tot_smj++;
memset(psmj, 0, lrsmj);
memcpy(psmj->name, name2, 20);
if(libel != NULL) {
	tot_txt++;
	memcpy(ptxt, libel2, 60);
	psmj->libel = tot_txt;
	writetxt(tot_txt);
	write_first_rec(ktxt, tot_txt, 0);
	}
writesmj(tot_smj);
write_first_rec(ksmj, tot_smj, 0);
}
