#ifndef DIR_ACNUC_H
#define DIR_ACNUC_H

#include "dir_io.h"
#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif

/*			STRUCTURE OF RECORDS OF ACNUC INDEX FILES 	*/

#define L_MNEMO 16
extern struct rsub {         /* SUBSEQ   one record for each parent or sub- sequence */
	char name[L_MNEMO]; /* seq name padded by blanks to L_MNEMO chars */
	int length, /* seq length; or 0 if record was deleted */
	   type, /* to SMJYT, for seq type */
	   pext, /* if > 0 this is a subsequence, pext points to EXTRACT for list of exons;
	   	    if <= 0 this is a parent sequence, -pext points to LONGL for list of subseqs */
	   plkey, /* to SHORTL for list of keywords */
	   plinf, /* if parent sequence, plkey points to LOCUS for corresponding record;
	   	     if subsequence, points to SHORTL for list of address of start of annotations; 
	   	     	this list contains only one element to be combined with the division rank
	   	     	for access to annotations */
	   phase, /* 100 * code_number + reading_frame_0_1_2 */
	   h; /* to SUBSEQ for next record with same hashing value or 0  */
	} *psub; 
#define lrsub sizeof(struct rsub)
#define readsub(x) if(dir_read(ksub,x,1,psub)!=1) dir_readerr(ksub,x)
#define writesub(x) if(dir_write(ksub,x,1,psub)) dir_writeerr(ksub,x)

extern struct rloc {         /* LOCUS    one record for each parent sequence */
	int sub, /* to SUBSEQ for corresponding record; or 0 if record was deleted */
	    pnuc, /* offset within flat file of beginning of sequence lines */
	    pinf, /* offset within flat file of beginning of annotation lines */
	    bef,next, /* to LOCUS for previous and next member of segment relationship */
	    spec, /* if genbank of embl, to SPECIES for corresponding species
	    	     if swissprot or nbrf, to SHORTL for list of corresponding species */
	    host, /* to SPECIES, for host species of virus or plasmid, or 0 */
	    plref, /* to SHORTL for list of references */
	    molec, /* to SMJYT for molecule */
	    placc, /* to SHORTL for list of accession numbers */
	    stat, /* to SMJYT for status */
	    org, /* to SMJYT for organelle, or 0 if none */
	    div; /* gives rank of corresponding division, useful for div argument of read_annots */
	char date[16]; /* seq date in format MM/DD/YYMM/DD/YY */
	} *ploc;
#define lrloc sizeof(struct rloc)
#define readloc(x) if(dir_read(kloc,x,1,ploc)!=1) dir_readerr(kloc,x)
#define writeloc(x) if(dir_write(kloc,x,1,ploc)) dir_writeerr(kloc,x)

#define WIDTH_KS 40
extern struct rkey {    /* KEYWORDS one record for each keyword, number 2 is RACINE (tree root) */
	char name[WIDTH_KS];  /* name padded by spaces or truncated to WIDTH_KS chars */
	int libel, /* to TEXT for corresponding label, or 0 if none */
	    plsub, /* to LONGL for list of associated seqs, plsub = 0 for synonyms */
	    desc, /* to SHORTL for list for descendants in tree structure;
	          the absolute value of the first elt of this list is the rank of corresponding
	          record in KEYWORDS; the sign of this number is negative iff there are
	          sequences associated to this record; other elements of list are "desc" values 
	          of records of descendants in tree; desc = 0 for synonyms */
	    syno, /* to KEYWORDS for next synonymous keyword, or 0 if none (major has value < 0) */
	    h; /* to KEYWORDS for next record with same hashing value, or 0 */
	} *pkey;
#define lrkey sizeof(struct rkey)
#define readkey(x) if(dir_read(kkey,x,1,pkey)!=1) dir_readerr(kkey,x)
#define writekey(x) if(dir_write(kkey,x,1,pkey)) dir_writeerr(kkey,x)

extern struct rspec {   /* SPECIES one record for each taxon, number 2 is RACINE (tree root) */
	char name[WIDTH_KS];  /* taxon name padded by spaces or truncated to WIDTH_KS chars */
	int libel, /* to TEXT for corresponding label, or 0 if none */
	    plsub, /* to LONGL for list of associated seqs, or 0 for synonyms */
	    desc, /* to SHORTL for list for descendants in tree structure;
	          the absolute value of the first elt of this list is the rank of corresponding
	          record in SPECIES; the sign of this number is negative iff there are
	          sequences associated to this record; other elements of list are "desc" values 
	          of records of descendants in tree; desc = 0 for synonyms */
	    syno, /* to SPECIES for next synonymous species, or 0 if none (major has value < 0) */
	    h, /* to SPECIES for next record with same hashing value, or 0 */
	    plhost;  /* to LONGL for list of seqs for which this species is a host, 0 for synos */
	} *pspec;
#define lrspec sizeof(struct rspec)
#define readspec(x) if(dir_read(kspec,x,1,pspec)!=1) dir_readerr(kspec,x)
#define writespec(x) if(dir_write(kspec,x,1,pspec)) dir_writeerr(kspec,x)

extern struct rshrt {        /* SHORTL mostly series of linked records containing various values;
		record #2 contains values of the hashing parameters as -hsub , -hkwsp;
		hashing data for seq names, keywords and species are stored from record #3 
		*/
	int val, /* value of an element of the short list */
	    next; /* to SHORTL for next element of the short list, or 0 when list is finished */
	} *pshrt;
#define lrshrt sizeof(struct rshrt)
#define readshrt(x) if(dir_read(kshrt,x,1,pshrt)!=1) dir_readerr(kshrt,x)
#define writeshrt(x) if(dir_write(kshrt,x,1,pshrt)) dir_writeerr(kshrt,x)

#define SUBINLNG 63
extern struct rlng {         /* LONGL series of linked records containing lists of SUBSEQ ranks */
	int sub[SUBINLNG]; /* array of ranks of SUBSEQ records or of 0s */
	int next; /* to LONGL for next element of the long list, or 0 when list is finished */
	} *plng;
#define lrlng sizeof(struct rlng)
#define readlng(x) if(dir_read(klng,x,1,plng)!=1) dir_readerr(klng,x)
#define writelng(x) if(dir_write(klng,x,1,plng)) dir_writeerr(klng,x)

extern struct rext {     /* EXTRACT  series of linked records, one for each exon of each subseq */
	int mere, /* to SUBSEQ for rank of parent seq containing this exon */
	    deb,fin, /* start and end positions of exon within parent sequence */
	    pnuc, /* offset within flat file of beginning of sequence lines of parent seq */
	    next; /* to EXTRACT for next exon of same subseq , or 0 */
	} *pext;
#define lrext sizeof(struct rext)
#define readext(x) if(dir_read(kext,x,1,pext)!=1) dir_readerr(kext,x)
#define writeext(x) if(dir_write(kext,x,1,pext)) dir_writeerr(kext,x)

extern struct rsmj {/* SMJYT one record for each Status, Molec, Journal, Year, Type, organ, div */
	char name[20]; /* first 2 chars 00, 01,... , 07 give nature of record, others give name */
	int plong, /* to LONGL for list of associated seqs */
	    libel; /* to TEXT for corresponding label, or 0 if none */
	} *psmj;
#define lrsmj sizeof(struct rsmj)
#define readsmj(x) if(dir_read(ksmj,x,1,psmj)!=1) dir_readerr(ksmj,x)
#define writesmj(x) if(dir_write(ksmj,x,1,psmj)) dir_writeerr(ksmj,x)

extern struct raut {       /* AUTHOR  one record for each author name (initials ignored) */
	char name[20]; /* author name padded with spaces; if "xxx...xxx" record was deleted */
	int plref, /* to SHORTL for list of references this author belongs to */
	    fut; /* unused */
	} *paut;
#define lraut sizeof(struct raut)
#define readaut(x) if(dir_read(kaut,x,1,paut)!=1) dir_readerr(kaut,x)
#define writeaut(x) if(dir_write(kaut,x,1,paut)) dir_writeerr(kaut,x)

extern struct rbib {       /* BIBLIO one record for each reference, book, thesis, or unpublished */
	char name[40]; /* reference name padded by spaces or truncated to 40 chars 
		journal citations appear as JOURNALCODE/vol/first_page
		book citations as BOOK/year/first_author
		theses citations as THESIS/year/first_author
		patent citations as PATENT/number
		other citations as UNPUBL/year/first_author
		*/
	int plsub, /* to SHORTL for list of associated parent sequences */
	    plaut, /* to SHORTL for list of associated authors (only 1 for book,thesis,unpubl) */
	    j, /* to SMJYT for rank of corresponding journal, or generic values BOOK, THESIS */
	    y; /* to SMJYT for rank of publication year or for rank of year 0 if unknown */
	} *pbib;
#define lrbib sizeof(struct rbib)
#define readbib(x) if(dir_read(kbib,x,1,pbib)!=1) dir_readerr(kbib,x)
#define writebib(x) if(dir_write(kbib,x,1,pbib)) dir_writeerr(kbib,x)

#define lrtxt 60         /* TEXT  taxa and keywords can have labels; elts of SMJYT have labels */
extern char ptxt[]; /* a label padded to lrtxt chars WITHOUT \n or \0 at end */
#define readtxt(x) if(dir_read(ktxt,x,1,ptxt)!=1) dir_readerr(ktxt,x)
#define writetxt(x) if(dir_write(ktxt,x,1,ptxt)) dir_writeerr(ktxt,x)

extern struct racc {     /* ACCESS  one record for each accession number */
	int plsub; /* to SHORTL for list of associated parent seqs */
	char name[1]; /* real size is ACC_LENGTH (global variable) padded by spaces */
	} *pacc;
void readacc(int recnum);  /* !!!!! always use readacc() and writeacc() */
void writeacc(int recnum);

extern struct rinfo{	/* to access seq annotations, always use read_annots/next_annots */
        char line[256]; /* to hold one line of annotations */
        } *pinfo;
#define lrinfo sizeof(struct rinfo)


extern int lmot,hoffst,hsub,hkwsp,nseq,nbrf,lenbit,lenw,maxa,longa,ACC_LENGTH;
extern DIR_FILE *ksub,*kloc,*kkey,*kspec,*kbib,*kacc,*ktxt,*ksmj,*kext,
*kaut,*kshrt,*klng;
extern int nbmrfa,flat_format,gcgcod,unixos,genbank,embl,divisions,swissprot,
	big_annots;
extern char nucbuf[]; /* to hold in memory the last part of sequence read */
extern char **gcgname; /* names of division files */
extern int *annotopened; /* says whether any division is currently opened */
extern FILE **divannot; /* streams associated to opened divisions */


/* prototypes for acnuc access*/
void acnucopen(void);  /* open the full acnuc data base in readonly mode */
void quick_list_meres(int *blist); /* puts in blist of length lenw list of parent sequences */
void simpleopen(void); /* open acnuc using only sequences and annotations */
void dir_acnucclose(void);  /* close it */
	/* get seq number, length, reading frame and genetic code from name */
int gsnuml(char *name,  int *length, int *frame, int *gencode);
int gfrag(int nsub,int first,int lfrag,char *dseq); /*read part of a sequence*/
void seq_to_annots(int numseq, long *faddr, int *div);
char *read_annots(long faddr, int div);
char *next_annots(long *pfaddr);
char *short_descr(int seqnum, char *text, int maxlen);/* get short description */
char *short_descr_p(int seqnum, char *text, int maxlen);/*get parent's descrip*/
char codaa(char *codon,int code);/*translate a codon with correct genetic code*/
char *translate_cds(int seqnum); /*complete translation of a seq in dyn memory*/
char *get_code_descr(int code);/*get short description of variant genetic code*/
int fcode(DIR_FILE *fp, char *search, int lcompar); /*search key in index file*/
int isenum(char *name); /*get sequence number from its name;upper or lowercase*/
int iknum(char *name, DIR_FILE *fp); /*get species or keyword number from name */
int hashmn(char *nom);  /* used by isenum */
int hasnum(char *name);  /* used by iknum */
void lngbit(int point, int *blist); /* read a long list in bit list in memory */
void dir_readerr(DIR_FILE *fich, int recnum);  /*write message after read err*/
int decode_locus_rec(char *line, char **pname, char **pmolec, 
	int *circular, char **pdivision, char **pdate);
/* for any acnuc index file fp, returns the rank of last valid record;
if endsort != NULL, stores in *endsort the rank of the last alphabetically sorted record  */
int read_first_rec(DIR_FILE *fp, int *endsort); 


/* prototypes of utility functions */
char complementer_base(char nucl); /* returns Watson-Crick complement of nucl */
void complementer_seq(char *seq, int len); /*in place complement strand of seq*/
void padtosize(/*out*/ char *paddedname, /*in*/ char *name, int size); /*add trailing spaces */
/* alphabetically compares l1 first chars of s1 and l2 of s2 ignoring trailing spaces
and returns 0 iff strings are equal */
int strcmptrail(char *s1, int l1, char *s2, int l2); 
void compact(char *chaine); /* removes all spaces in place */
void majuscules(char *name); /* converts to uppercase in place */
int trim_key(char *name); /* removes trailing spaces */


/* prototypes for acnuc management*/
void dir_acnucopen(char *); /* use RO or WP or WA to control file access */
void dir_acnucflush(void);  /* flush all data base files */
void dir_writeerr(DIR_FILE *fich, int recnum);/*write message after write err*/
void write_first_rec(DIR_FILE *fp, int total, int endsort);/*update file size*/
         /* modify short list */
int mdshrt(DIR_FILE *kan, int nrec, int offset, int val, int *newplist);
	/* modify long list */
int mdlng(DIR_FILE *kan, int nrec, int offset, int val, int *newplist);
int supshrt(int point, int val); 
int addshrt(int point, int val); 
int addlng(int point, int val); 
int suplng(int point, int val); 
void delseq(int nsub);  /* remove sequence from data base */
void suphsh(int numrec, DIR_FILE *kan); /* remove key from hash chains */
void addhsh(int recnum, DIR_FILE *kan); /* add key to hash chains */
/* call when kacc, kaut, kbib or ksmj were updated */
void write_sorted_part(DIR_FILE *fp); 
/* to create a species and/or get its number, use ascend=NULL for a new root*/
int crespecies(char *ascend, char *nom);
/* to create a keyword and/or get its number, use ascend=NULL for a new root*/
int crekeyword(char *ascend, char *nom);

/* prototypes of bit-manipulation functions */
int irbit(int *pdeblist, int deb, int fin);
void bit1(int *plist, int num);
void bit0(int *plist, int num);
int testbit(int *plist, int num);
int bcount(int *plist, int last);
void et(int *listet, int *list1, int *list2, int len);
void ou(int *listou, int *list1, int *list2, int len);
void non(int *listnon, int *list, int len);

#ifdef vms
#define MINRECSIZE 52
#define MAXDOCSTRING 256
#define MAXRECLENGTH 20000 /* the max rec size under VMS of .seq files for gcg*/
#endif

#endif /* DIR_ACNUC_H */
