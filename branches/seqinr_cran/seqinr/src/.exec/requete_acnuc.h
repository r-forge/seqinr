#ifndef REQUETE_ACNUC_H
#define REQUETE_ACNUC_H

extern int tlist;
extern char *deflnames[]; /* noms listes pour def */
extern int defoccup[];
extern int deflocus[];
extern char defgenre[];
extern int defllen[];
extern int *defbitlist; /* listes de bits pour interf avec routine def */

/* prototypes */
void prep_acnuc_requete(void);
int proc_requete(char *requete, char *message, char *nomliste, int *numliste);
void free_list(int num);
#endif /*  REQUETE_ACNUC_H  */
