/* convertir une banque acnuc vers le format big_annots
usage unique
*/

#include "dir_acnuc.h"

/* included functions */
void add_smj(char *name, char *libel);


int main(void)
{
int nloc, totloc, newpinf, newpnuc, div1, div2, nsub, point, doit, num_mere;
char mnemo[L_MNEMO + 1], mere[L_MNEMO + 1], *p;
int newpnuc2, newpinf2;

dir_acnucopen("WP");
if(big_annots) {
	fprintf(stderr, "Bank is already in the bigannots format!\n");
	exit(0);
	}

/* conversion de LOCUS */
printf("Processing LOCUS...\n");
totloc = read_first_rec(kloc, NULL);
for(nloc = 2; nloc <= totloc; nloc++) {
	readloc(nloc);
	if(ploc->sub == 0) continue;
	goffset(ploc->pinf, &div1, &newpinf);
	goffset(ploc->pnuc, &div2, &newpnuc);
	if(div1 != div2 || div1 < 0 || div1 > divisions ) {
		readsub(ploc->sub);
		fprintf(stderr, "Pb div1 != div2  %.16s\n", psub->name);
		continue;
		}
	if(ploc->div != 0) {
		readsub(ploc->sub);
		fprintf(stderr, "Pb div != 0  %.16s\n", psub->name);
		continue;
		}
	ploc->pinf = newpinf;
	ploc->pnuc = newpnuc;
	ploc->div = div1;
	writeloc(nloc);
	}
dir_acnucflush();

if( ! (nbrf || swissprot) ) {
/* conversion de EXTRACT & SHORTL */
   printf("Processing EXTRACT & SHORTL...\n");
   for(nsub = 2; nsub <= nseq; nsub++) {
	readsub(nsub);
	if(psub->pext <= 0) continue;
	doit = TRUE;
	memcpy(mnemo, psub->name, L_MNEMO); mnemo[L_MNEMO] = 0;
	memcpy(mere, mnemo, L_MNEMO + 1);
	p = strchr(mere, '.');
	if(p == NULL) doit = FALSE;
	if(doit) {
		*p = 0;
		num_mere = isenum(mere);
		if(num_mere == 0) doit = FALSE;
		}
	if(doit) {
		readsub(num_mere);
		readloc(psub->plinf);
		readsub(nsub);
		readshrt(psub->plinf);
		}
	if(pshrt->next != 0) doit = FALSE;
	if(doit) {
		goffset(pshrt->val, &div1, &newpinf);
		if(div1 != ploc->div || newpinf < ploc->pinf) doit = FALSE;
		}
	if(doit) {
		pshrt->val = newpinf;
		writeshrt(psub->plinf);
		}
	if( !doit ) fprintf(stderr, "Pb annots de fille %s\n", mnemo);
	point =  psub->pext;
	do	{
		readext(point);
		readsub(pext->mere);
		readloc(psub->plinf);
		goffset(pext->pnuc, &div1, &newpnuc);
		if(newpnuc != ploc->pnuc || div1 != ploc->div) {
			fprintf(stderr, "Pb extract %s\n", mnemo);
			break;
			}
		pext->pnuc = newpnuc;
		writeext(point);
		point = pext->next;
		}
	while(point != 0);
	}
   }

add_smj("07BIG_ANNOTS", "Annotation/Sequence addresses on 2*32 bits");

dir_acnucclose();
return 0;
}


void add_smj(char *name, char *libel)
{
static char name2[21], libel2[61];
int tot_smj, tot_txt, totsort;

padtosize(name2, name, 20);
if(libel != NULL) padtosize(libel2, libel, 60);
tot_smj = read_first_rec(ksmj, &totsort) + 1;
tot_txt = read_first_rec(ktxt, NULL) + 1;
memset(psmj, 0, lrsmj);
memcpy(psmj->name, name2, 20);
if(libel != NULL) {
	memcpy(ptxt, libel2, 60);
	psmj->libel = tot_txt;
	writetxt(tot_txt);
	write_first_rec(ktxt, tot_txt, 0);
	}
writesmj(tot_smj);
write_first_rec(ksmj, tot_smj, totsort);
}
