#include "dir_acnuc.h"
#define WID sizeof(psmj->name)


int main(int argc, char **argv)
{
char *oldname, *newname, smjname[21];
int num1, num2, totsmj;

if(argc != 3) {
	fprintf(stderr, "Usage:     %s old new\n", argv[0]);
	exit(ERREUR);
	}
dir_acnucopen("RW");
oldname = argv[1];
newname = argv[2];
for(num1 = 0; num1 <= divisions; num1++)
	if(strcmp(gcgname[num1], oldname) == 0) break;
if(num1 > divisions) {
	fprintf(stderr, "%s is not the name of an existing division\n",
		oldname);
	exit(ERREUR);
	}
for(num2 = 0; num2 <= divisions; num2++)
	if(strcmp(gcgname[num2], newname) == 0) break;
if(num2 <= divisions) {
	fprintf(stderr, "%s already exists as a division name\n",
		newname);
	exit(ERREUR);
	}
if(strlen(newname) + 5 > WID) {
	fprintf(stderr, "%s is too long for a division name\n",
		newname);
	exit(ERREUR);
	}
if(flat_format)
	strcpy(smjname, "06FLT");
else
	strcpy(smjname, "06GCG");
memset(smjname + 5, ' ', WID-5);
memcpy(smjname + 5, oldname, strlen(oldname));
totsmj = read_first_rec(ksmj, NULL);
for(num2 = 2; num2 <= totsmj; num2++) {
	readsmj(num2);
	if(strncmp(smjname, psmj->name, WID) == 0) break;
	}
if(num2 > totsmj) {
	fprintf(stderr, "change impossible\n");
	exit(ERREUR);
	}
readsmj(num2);
memset(psmj->name + 5, ' ', WID - 5);
memcpy(psmj->name + 5, newname, strlen(newname));
writesmj(num2);
dir_acnucclose();
dir_acnucopen("RO");
gcgors("inf", num1, FALSE);
if( ! (annotopened[num1]) ) {
	printf("Warning: the file corresponding to the new name"
		" cannot be opened\n");
	}
dir_acnucclose();
return 0;
}

