#define SIMEXT_NFMAX 500
typedef struct {
	int nfmax;  /* maximum # of fragments in a virtual subseq */
	int newsim; /* TRUE when a new virtual subseq has just been created*/
	int pinfme; /* pointer to annots of parent of virtual subseq */
	int siminf; /* pointer to annots of virtual subseq */
	int div;    /* div of these annots */
	int nfrags; /* # of fragments in virtual subseq (from 1)*/
	int valext[SIMEXT_NFMAX][4];
/* valext(.,0)=beginning of fragment */
/* valext(.,1)=end of fragment	*/
/* valext(.,2)=pointer to nucleot for that fragment */
/* valext(.,3)=# of parent seq of that fragment	*/
	} simext_struct;
