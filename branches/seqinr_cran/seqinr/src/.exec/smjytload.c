#include "dir_acnuc.h"

void job(int option);
int mod_j_name(void);


int main(int argc, char **argv)
{
int option;
char reponse[10];

dir_acnucopen("WP");

while(TRUE) {
	do	{
		printf(
		"\n\nStatus(0)  Molec(1)  Journal(2)  Year(3)  Type(4)  "
	   "Organelle(5)\nsuppr. unused Journals+Years(6) Misc(7) Stop(8)?  ");
		gets(reponse);
		option = -1; sscanf(reponse, "%d", &option);
		}
	while(option < 0 || option > 8);
	
	if(option == 8) {
		dir_acnucclose();
		exit(0);
		}
	else 	job(option);
	}
}


void job(int option)
{
int last, num, sup, ltext;
char reponse[50], code[21], operation, libel[100];

if(option == 6) { /* SUPPRESSION TOUS JOURNAUX ET ANNEES INUTILES */
	printf("Confirmez la suppression de tous les journaux et "
	"annees inutiles: (y/[n])?  ");
	gets(reponse);
	majuscules(reponse);
	if(reponse[0] != 'Y') return;
	last = read_first_rec(ksmj, NULL); sup = 0;
	for(num = 2; num <= last; num++) {
		readsmj(num);
		if(memcmp(psmj->name, "02", 2) != 0 && 
			memcmp(psmj->name, "03", 2) != 0) continue;
		if(psmj->plong != 0) continue;
		printf("%.20s   suppressed\n", psmj->name); sup++;
		memset(psmj->name, 'x', 20);
		psmj->plong = psmj->libel = 0;
		writesmj(num);
		}
	if(sup != 0) {
		write_first_rec(ksmj, last, 0);
		dir_acnucflush();
		}
	return;
	}

do	{
	printf(
 "\nDo you want to add(A) or suppress(S) elements or change a text line(T)?  ");
	if(option == 2) printf("\nor modify(M) a journal name?  ");
	gets(reponse); majuscules(reponse); 
	operation = reponse[0];
	}
while(strchr("ASTM", operation) == NULL);
sprintf(code, "%2.2d", option);
if(operation == 'M' && option == 2) {
	do	num = mod_j_name();
	while(num);
	return;
	}
while(TRUE) {
	printf("%-20s*\n", "Code? or <return>"); gets(reponse);
	if(strcmptrail(reponse, strlen(reponse), NULL, 0) == 0) break;
	padtosize(code + 2, reponse, 18);
	majuscules(code);
	num = fcode(ksmj, code, 20);
	if(operation == 'A') {  /* AJOUT */
		if(num != 0) {
			printf("Error: this code already exists: %s\n", code);
			continue;
			}
		printf("%-60s*\n", "Libelle?");
		gets(libel);
		ltext = read_first_rec(ktxt, NULL);
		last = read_first_rec(ksmj, NULL);
		memcpy(psmj->name, code, 20);
		psmj->plong = 0;
		psmj->libel = ++ltext;
		writesmj(++last);
		padtosize(libel, libel, 60);
		memcpy(ptxt, libel, 60);
		writetxt(ltext);
		write_first_rec(ksmj, last, 0);
		write_first_rec(ktxt, ltext, 0);
		dir_acnucflush();
		printf("%s created\n", code);
		}
	else if(operation == 'S') { /* SUPPRESSION */
		if(num == 0) {
			printf("Error: this code does not exist: %s\n", code);
			continue;
			}
		readsmj(num);
		if(psmj->plong != 0) {
			printf("Error: this code code cannot be suppressed, "
				"it has associated sequences: %s\n", code);
			continue;
			}
		memset(psmj, 0, lrsmj);
		memset(psmj->name, 'x', 20);
		writesmj(num);
		write_first_rec(ksmj, read_first_rec(ksmj, NULL), 0);
		dir_acnucflush();
		printf("%s suppressed\n", code);
		}
	else if(operation == 'T') { /* changement du texte */
		if(num == 0) {
			printf("Error: this code does not exist: %s\n", code);
			continue;
			}
		readsmj(num);
		strcpy(ptxt, "*** empty ***");
		if(psmj->libel != 0) readtxt(psmj->libel);
		printf("Present libel value:\n%-60s*\n%-60s*\n", ptxt,
			"New libel value ?");
		gets(libel);
		padtosize(libel, libel, 60);
		printf("Confirm new libel value\n%s\n? (y/[n])   ", libel);
		gets(reponse); majuscules(reponse);
		if(reponse[0] != 'Y') continue;
		if(psmj->libel == 0) {
			psmj->libel = read_first_rec(ktxt, NULL) + 1;
			write_first_rec(ktxt, psmj->libel, 0);
			writesmj(num);
			}
		memcpy(ptxt, libel, 60);
		writetxt(psmj->libel);
		dir_acnucflush();
		printf("Libel has been changed\n");
		}
	}
}


int mod_j_name(void)
{
char reponse[50], rep[10], code_o[21] = "02", code_n[21] = "02";
int old, num, numseq, ps, i, j, k, num2, nr, ns, point2, pln, tn, pa, pl, y, y2;
char bibname[41], aux[100];
static int *blist = NULL;

if(blist == NULL) {
	blist = calloc(lenw, sizeof(int));
	write_first_rec(kbib, read_first_rec(kbib, NULL), 0);
#ifdef unix
	dir_set_mmap(kbib);
#endif
	}
if(blist == NULL ) {
	printf("Not enough memory\n");
	return FALSE;
	}
printf("\nPresent journal name? (or <RETURN>)  ");
gets(reponse);
if(strlen(reponse) == 0) return FALSE;
majuscules(reponse);
padtosize(code_o+2, reponse, 18);
old = fcode(ksmj, code_o, 20);
if(old == 0) {
	printf("Unknown journal name: %s\n", reponse);
	return TRUE;
	}
printf("%-18s*\n", "New name?");
gets(reponse);
majuscules(reponse); compact(reponse);
padtosize(code_n+2, reponse, 18);
num = fcode(ksmj, code_n, 20);
if(num == 0) {
	printf("Old: %s New: %s  ? (y/[n])   ", code_o+2, code_n+2);
	gets(rep); majuscules(rep);
	if(rep[0] != 'Y') return TRUE;
	num  = old;
	}
else	{
	readsmj(num); pln = psmj->plong; tn = psmj->libel;
	if(tn != 0) { readtxt(tn); }
	else strcpy(ptxt, "*** no libel ***");
	printf("Old: %s New: %s\n"
		"The new name is that of another yet known journal\n%.60s\n"
		"Do you really want to merge the 2 names? (y/[n])   )",
		code_o+2, code_n+2, ptxt);
	gets(rep); majuscules(rep);
	if(rep[0] != 'Y') return TRUE;
	}
readsmj(old);
if(psmj->plong == 0) {
	if(num == old) {
		memcpy(psmj->name, code_n, 20);
		writesmj(old);
		}
	else	{
		memset(psmj, 0, lrsmj);
		memset(psmj->name, 'x', 20);
		writesmj(old);
		}
	write_first_rec(ksmj, read_first_rec(ksmj, NULL), 0);
	dir_acnucflush();
	printf("Change has been done\n");
	return TRUE;
	}
lngbit(psmj->plong, blist);
numseq = 1;
while((numseq = irbit(blist, numseq, nseq)) != 0) {
	readsub(numseq);
	readloc(psub->plinf);
	ps = ploc->plref;
	compact(code_n);
	j = strlen(code_n);
	while(ps != 0) {
		readshrt(ps); ps = pshrt->next;
		nr = pshrt->val;
		readbib(nr); pa = pbib->plaut; pl = pbib->plsub; y = pbib->y;
		if(pbib->j == old) {
			memcpy(bibname, pbib->name, 40); bibname[40] = 0;
			trim_key(bibname);
			k = strlen(bibname);
			i = strchr(bibname, '/') - bibname + 1;
			memcpy(aux, code_n+2, j-2); aux[j-2] = 0;
			strcat(aux, bibname+i-1);
			padtosize(aux, aux, 40);
			num2 = 0;
			if(num != old) num2 = fcode(kbib, aux, 40);
			if(num2 != 0) {
				readbib(num2); y2 = pbib->y;
				memcpy(aux, pbib->name, 40); aux[40] = 0;
				point2 = pl;
				while(point2 != 0) {
					readshrt(point2); point2 = pshrt->next;
					ns = pshrt->val;
					readsub(ns);
					mdshrt(kloc, psub->plinf, -8, nr, NULL);
					mdshrt(kloc, psub->plinf, 8, num2, &ps);
					mdshrt(kbib, num2, 1, ns, NULL);
					if(y2 != y) {
/* DETERMINER SI L'ANNEE Y EST NECESSITEE A CAUSE D'UNE AUTRE REF DE LA SEQ */
						int tmp, need = FALSE;
						tmp = ps;
						while(tmp != 0) {
							readshrt(tmp);
							tmp = pshrt->next;
							readbib(pshrt->val);
							if(pbib->y == y) {
								need = TRUE;
								break;
								}
							}
						if(!need) 
						    mdlng(ksmj,y,-1,ns,NULL);
						mdlng(ksmj,y2,1,ns,NULL);
						}
					}
				point2 = pa;
				while(point2 != 0) {
					readshrt(point2); point2 = pshrt->next;
					mdshrt(kaut, pshrt->val, -1, nr, NULL);
					}
				memset(pbib, 0, lrbib);
				memset(pbib->name, 'x', 40);
				writebib(nr);
				nr = num2;
				y = y2;
				readbib(nr); pl = pbib->plsub; pa = pbib->plaut;
				memcpy(aux, pbib->name, 40); aux[40] = 0;
				}
			memcpy(pbib->name, aux, 40);
			pbib->plsub = pl; pbib->plaut = pa;
			pbib->j = num; pbib->y = y;
			writebib(nr);
			printf("%s  %s\n", bibname, aux);
			if(k-i+j-1 > 40) printf(
				"Warning: reference has been truncated\n");
			}
		}
	if(num != old && pln != 0) addlng(pln, numseq);
	}
padtosize(code_n, code_n, 20);
if(num != old && pln == 0) {
	readsmj(old);
	memcpy(psmj->name, code_n, 20);
	psmj->libel = tn;
	writesmj(num);
	}
if(num == old) {
	readsmj(old);
	memcpy(psmj->name, code_n, 20);
	writesmj(num);
	}
else	{
	memset(psmj, 0, lrsmj);
	memset(psmj->name, 'x', 20);
	writesmj(old);
	}
write_first_rec(ksmj, read_first_rec(ksmj, NULL), 0);
dir_acnucflush();
printf("Change has been done\n");
return TRUE;
}


int fbib(char *name, int nr)
{
static int inmem = 1;
int num, debut;
static int totbib;
static char **allbib;

if(name == NULL) {
	char *tmp;
	totbib = read_first_rec(kbib, NULL);
	tmp = (char *)malloc( totbib * 40);
	if(tmp != NULL) allbib = (char **)malloc((totbib+1) * sizeof(char *));
	if(tmp == NULL || allbib == NULL) return TRUE;
	for(num = 2; num <= totbib; num++)
		allbib[num] = tmp + (num-1)*40;
	return FALSE;
	}

if(nr != 0) {
	memcpy(allbib[nr], name, 40);
	return 0;
	}

for(num = 2; num <= inmem; num++) {
	if(memcmp(name, allbib[num], 40) == 0) return num;
	}
debut = inmem + 1;
for(num = debut; num <= totbib; num++) {
	readbib(num);
	memcpy(allbib[num], pbib->name, 40);
	inmem = num;
	if(memcmp(name, allbib[num], 40) == 0) return num;
	}
return 0;
}

