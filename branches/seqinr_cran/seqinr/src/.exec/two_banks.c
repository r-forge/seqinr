#include "dir_acnuc.h"

/* declaration de la structure de memorisation d'une banque acnuc */
struct environ {
	char *acnuc, *gcgacnuc;
	};

struct indexes {
	DIR_FILE *ksub, *kloc, *kkey, *kspec, *kbib, *kacc, *ktxt,
		*ksmj, *kext, *kaut, *kshrt, *klng;
	};

struct db_struct {
	int genbank, nbrf, swissprot, embl;
	int hsub, hkwsp, flat_format, big_annots, use_div_sizes;
	int ACC_LENGTH, hashing_algorithm;
	};

struct divisions {
	int divisions;
	FILE **divannot, **divseq;
	int *annotopened, *seqopened, *div_offset;
	char **gcgname;
	};

struct db_size {
	int nseq, lenbit, lenw, maxa, longa;
	};

struct acnuc_status {
	struct environ environ;
	struct indexes indexes;
	struct db_struct db_struct;
	struct divisions divisions;
	struct db_size db_size;
	int fcode_nsortd[5];
	};


/* included functions */
struct acnuc_status *store_acnuc_status(void);
void set_current_acnuc_db(struct acnuc_status *);
int sizeof_acnuc_status(void);

/* globals */
static struct acnuc_status *current_acnuc_status = NULL;

/* called global variables */
extern int use_div_sizes, hashing_algorithm;
extern FILE **divannot, **divseq;
extern int *annotopened, *seqopened;
extern unsigned *div_offset;
extern char **gcgname;
extern int gfrag_premier;
extern int fcode_nsortd[5];
extern int chg_acnuc(char *acnucvar, char *gcgacnucvar);


struct acnuc_status *store_acnuc_status(void)
/* memorizes the currently opened acnuc db
returns a pointer to all these data
or NULL if not enough memory
*/
{
struct acnuc_status *status;
char *p;
int l, i;

status = (struct acnuc_status *) malloc(sizeof(struct acnuc_status));
if(status == NULL) return NULL;

/* environ */
p = prepare_env_var("acnuc");	
l = strlen(p);
status->environ.acnuc = malloc(l + 1);
if(status->environ.acnuc == NULL) return NULL;
memcpy(status->environ.acnuc, p, l + 1);

p = prepare_env_var("gcgacnuc");	
l = strlen(p);
status->environ.gcgacnuc = malloc(l + 1);
if(status->environ.gcgacnuc == NULL) return NULL;
memcpy(status->environ.gcgacnuc, p, l + 1);

/* indexes */
status->indexes.ksub = ksub;
status->indexes.kloc = kloc;
status->indexes.kkey = kkey;
status->indexes.kspec = kspec;
status->indexes.kbib = kbib;
status->indexes.kacc = kacc;
status->indexes.ktxt = ktxt;
status->indexes.ksmj = ksmj;
status->indexes.kext = kext;
status->indexes.kaut = kaut;
status->indexes.kshrt = kshrt;
status->indexes.klng = klng;

/* db_struct */
status->db_struct.genbank = genbank;
status->db_struct.nbrf = nbrf;
status->db_struct.swissprot = swissprot;
status->db_struct.embl = embl;
status->db_struct.hsub = hsub;
status->db_struct.hkwsp = hkwsp;
status->db_struct.flat_format = flat_format;
status->db_struct.big_annots = big_annots;
status->db_struct.use_div_sizes = use_div_sizes;
status->db_struct.ACC_LENGTH = ACC_LENGTH;
status->db_struct.hashing_algorithm = hashing_algorithm;

/* divisions */
status->divisions.divisions = divisions;
status->divisions.divannot = (FILE **) malloc( (divisions+1) * sizeof(FILE *) );
if( status->divisions.divannot == NULL ) return NULL;
for(i = 0; i <= divisions; i++ )
	status->divisions.divannot[i] = divannot[i];
if( ! flat_format ) {
	status->divisions.divseq = 
		(FILE **) malloc( (divisions+1) * sizeof(FILE *) );
	if( status->divisions.divseq == NULL ) return NULL;
	for(i = 0; i <= divisions; i++ )
		status->divisions.divseq[i] = divseq[i];
	}
status->divisions.annotopened = 
		(int *) malloc( (divisions+1) * sizeof(int) );
if( status->divisions.annotopened == NULL ) return NULL;
for(i = 0; i <= divisions; i++ )
	status->divisions.annotopened[i] = annotopened[i];
if( ! flat_format ) {
	status->divisions.seqopened = 
		(int *) malloc( (divisions+1) * sizeof(int) );
	if( status->divisions.seqopened == NULL ) return NULL;
	for(i = 0; i <= divisions; i++ )
		status->divisions.seqopened[i] = seqopened[i];
	}
if( ! big_annots ) {
	status->divisions.div_offset = 
		(int *) malloc( (divisions+1) * sizeof(int) );
	if( status->divisions.div_offset == NULL ) return NULL;
	for(i = 0; i <= divisions; i++ )
		status->divisions.div_offset[i] = div_offset[i];
	}
status->divisions.gcgname = (char **) malloc( (divisions+1) * sizeof(char *) );
if( status->divisions.gcgname == NULL ) return NULL;
for(i = 0; i <= divisions; i++ )
	status->divisions.gcgname[i] = gcgname[i];

/* db_size */
status->db_size.nseq = nseq;
status->db_size.lenbit = lenbit;
status->db_size.lenw = lenw;
status->db_size.maxa = maxa;
status->db_size.longa = longa;

/* fcode_nsortd */
for(i = 0; i < 5 ; i++)
	status->fcode_nsortd[i] = fcode_nsortd[i];

current_acnuc_status = status;

return status;
} /* end of store_acnuc_status */


void set_current_acnuc_db(struct acnuc_status *status)
/* orienter la banque acnuc courante vers ce qui a ete memorise dans status
*/
{
int i;

/* sauver modifs faites dans la banque courante */
if(current_acnuc_status != NULL) {
	if( divisions > current_acnuc_status->divisions.divisions ) {
		dir_acnucclose();
		fprintf(stderr, 
			"Cannot handle changes in number of divisions!\n");
		exit(ERREUR);
		}
	for(i=0; i <= divisions; i++) {
		current_acnuc_status->divisions.divannot[i] = divannot[i];
		current_acnuc_status->divisions.annotopened[i] = annotopened[i];
		}
	if( ! flat_format ) {
		for(i=0; i <= divisions; i++) {
			current_acnuc_status->divisions.divseq[i] = divseq[i];
			current_acnuc_status->divisions.seqopened[i] = 
								seqopened[i];
			}
		}
	for(i=0; i < 5; i++) 
		current_acnuc_status->fcode_nsortd[i] = fcode_nsortd[i];	
	}

/* environ */
chg_acnuc(status->environ.acnuc, status->environ.gcgacnuc);

/* indexes */
ksub = status->indexes.ksub;
kloc = status->indexes.kloc;
kkey = status->indexes.kkey;
kspec = status->indexes.kspec;
kbib = status->indexes.kbib;
kacc = status->indexes.kacc;
ktxt = status->indexes.ktxt;
ksmj = status->indexes.ksmj;
kext = status->indexes.kext;
kaut = status->indexes.kaut;
kshrt = status->indexes.kshrt;
klng = status->indexes.klng;

/* db_struct */
genbank = status->db_struct.genbank;
nbrf = status->db_struct.nbrf;
swissprot = status->db_struct.swissprot;
embl = status->db_struct.embl;
hsub = status->db_struct.hsub;
hkwsp = status->db_struct.hkwsp;
flat_format = status->db_struct.flat_format;
big_annots = status->db_struct.big_annots;
use_div_sizes = status->db_struct.use_div_sizes;
ACC_LENGTH = status->db_struct.ACC_LENGTH;
hashing_algorithm = status->db_struct.hashing_algorithm;

/* divisions */
divisions = status->divisions.divisions;
divannot = status->divisions.divannot;
annotopened = status->divisions.annotopened;
for(i=0; i <= divisions; i++) {
	divannot[i] = status->divisions.divannot[i];
	annotopened[i] = status->divisions.annotopened[i];
	}
if( ! flat_format ) {
	divseq = status->divisions.divseq;
	seqopened = status->divisions.seqopened;
	for(i=0; i <= divisions; i++) {
		divseq[i] = status->divisions.divseq[i];
		seqopened[i] = status->divisions.seqopened[i];
		}
	}
if( ! big_annots ) {
	div_offset = status->divisions.div_offset;
	for(i=0; i <= divisions; i++) 
		div_offset[i] = status->divisions.div_offset[i];
	}
gcgname = status->divisions.gcgname;
for(i=0; i <= divisions; i++) 
	gcgname[i] = status->divisions.gcgname[i];

/* db_size */
nseq = status->db_size.nseq;
lenbit = status->db_size.lenbit;
lenw = status->db_size.lenw;
maxa = status->db_size.maxa;
longa = status->db_size.longa;

/* fcode_nsortd */
for(i = 0; i < 5 ; i++)
	fcode_nsortd[i] = status->fcode_nsortd[i];

current_acnuc_status = status;

gfrag_premier = TRUE;
} /* end of set_current_acnuc_db */



int sizeof_acnuc_status(void)
{
return sizeof(struct acnuc_status);
}


/* POUR LE TESTER  

main()
{
void *db1, *db2;
int totloc, seqnum;
char seq[100];

chg_acnuc("/banques0/genbank/index", "/banques0/genbank/flat_files");
acnucopen();
db1 = (void *)store_acnuc_status();

chg_acnuc("/banques0/swissprot/index", "/banques0/swissprot/flat_files");
acnucopen();
db2 = (void *)store_acnuc_status();

set_current_acnuc_db(db1);
printf("hsub=%d hkwsp=%d %s nseq=%d\n", hsub,hkwsp,ksub->filename,nseq);
totloc = read_first_rec(kloc, NULL);
readloc(totloc); seqnum = ploc->sub;
gfrag(seqnum, 1, 60, seq);
readsub(seqnum);
printf("%.16s %s\n", psub->name, seq);
seqnum = fcode(kacc, "M00001", ACC_LENGTH);
printf("M00001 = %d\n", seqnum);

set_current_acnuc_db(db2);
printf("hsub=%d hkwsp=%d %s nseq=%d\n", hsub,hkwsp,ksub->filename,nseq);
totloc = read_first_rec(kloc, NULL);
readloc(totloc); seqnum = ploc->sub;
gfrag(seqnum, 1, 60, seq);
readsub(seqnum);
printf("%.16s %s\n", psub->name, seq);
seqnum = fcode(kacc, "Q28799", ACC_LENGTH);
printf("Q28799 = %d\n", seqnum);

set_current_acnuc_db(db1);
printf("hsub=%d hkwsp=%d %s nseq=%d\n", hsub,hkwsp,ksub->filename,nseq);
gfrag(2, 1, 60, seq);
readsub(2);
printf("%.16s %s\n", psub->name, seq);
seqnum = fcode(kacc, "M00001", ACC_LENGTH);
printf("M00001 = %d\n", seqnum);

set_current_acnuc_db(db2);
printf("hsub=%d hkwsp=%d %s nseq=%d\n", hsub,hkwsp,ksub->filename,nseq);
gfrag(2, 1, 60, seq);
readsub(2);
printf("%.16s %s\n", psub->name, seq);
seqnum = fcode(kacc, "Q28799", ACC_LENGTH);
printf("Q28799 = %d\n", seqnum);
}

 FIN DU TEST */
