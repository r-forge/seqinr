/* public interface used with generic pointers:
void *init_dynlist(void);
char *find_dynlist(char *name, void *tree, int create_if_not_found);
int sizeof_dynlist(void *arbre);
int remove_dynlist(char *name, void *tree);
char *next_dynlist(void *arbre, void **p_noeud);
int optimize_dynlist(void *tree); 
void free_dynlist(void *tree);

En plus pour dynlist avec pointeur pour extra_data:
char *find_g_dynlist(char *name, void *tree, int create_if_not_found, void **node);
void *dynlist_get_extra(void *node);
void dynlist_set_extra(void *node, void *extra);
*/

#include <stdlib.h>
#include <string.h>
#ifndef FALSE
#define FALSE 0
#define TRUE (!FALSE)
#endif
#define MAX_PROF 300

typedef struct _bt_node {
	char *name;
	struct _bt_node *before;
	struct _bt_node *after;
	struct _bt_node *pere;
	} bt_node;

typedef struct _bt_tree {
	bt_node *racine;
	int total;
	} bt_tree;

/* true prototypes */
bt_tree *init_dynlist(void);

char *find_dynlist(char *name, bt_tree *tree, int create_if_not_found);
/* trouver noeud de nom name dans dynlist tree 
si trouve, rend pointeur vers ce nom, 
sinon, si create_if_not_found==TRUE celui-ci est cree
	mais rend NULL si pas assez de memoire
sinon rend NULL
*/

int sizeof_dynlist(bt_tree *arbre);

int remove_dynlist(char *name, bt_tree *tree);
/* returns TRUE iff name was in list and was removed */

char *next_dynlist(bt_tree *arbre, bt_node **p_noeud);
/* to loop around a dynlist by:
void *loop = NULL;
char *p;
while( (p = next_dynlist(dynlist, &loop)) != NULL) {
	process p coming out sorted
	}
pour parcours avec suppressions:
char *q;
p = next_dynlist(dynlist, &loop);
while(p != NULL) {
	if(ne pas supprimer) {
		p = next_dynlist(dynlist, &loop);
		continue;
		}
	q = next_dynlist(dynlist, &loop);
	remove_dynlist(p, dynlist);
	p = q;
	}
*/

int optimize_dynlist(bt_tree *tree); 
/* done internally, so not much necesary */

void free_dynlist(bt_tree *tree);

/* private stuff */
static bt_node *find_node_internal(char *name, bt_node *racine, 
	int create_if_not_found, int *total, int *p_maxprof);
static void free_bt_internal(bt_node *racine);
static void remove_node_internal(bt_tree *tree, bt_node *noeud);
static bt_node *add_mid_node(int debut, int fin, bt_node **list);
static bt_node **calc_list(bt_node *noeud, bt_node **list);


bt_tree *init_dynlist(void)
{
return (bt_tree *)calloc(1, sizeof(bt_tree));
}


int sizeof_dynlist(bt_tree *arbre)
{
return arbre->total;
}


char *find_dynlist(char *name, bt_tree *tree, int create_if_not_found)
{
bt_node *noeud;
int maxprof = 0;

noeud = find_node_internal(name, tree->racine, create_if_not_found, 
	&(tree->total), &maxprof );
if(tree->racine == NULL) tree->racine = noeud;
if(noeud == NULL) return NULL;
if(create_if_not_found && maxprof >= MAX_PROF) {
	if( optimize_dynlist(tree) == 2 ) return NULL;
	}
return noeud->name;
}


char *next_dynlist(bt_tree *arbre, bt_node **p_noeud)
{
bt_node *anc, *noeud;

if(arbre == NULL) return NULL;
noeud = *p_noeud;
if(noeud == NULL) {
	noeud = arbre->racine;
	if(noeud == NULL) return NULL;
	while(noeud->before != NULL) noeud = noeud->before;
	*p_noeud = noeud;
	return noeud->name;
	}
else if(noeud->after != NULL) {
	noeud = noeud->after;
	while(noeud->before != NULL) noeud = noeud->before;
	*p_noeud = noeud;
	return noeud->name;
	}
else	{
	while(TRUE) {
		anc = noeud->pere;
		if(anc == NULL) return NULL;
		if(anc->before == noeud) {
			*p_noeud = anc;
			return anc->name;
			}
		noeud = anc;
		}
	}
}


void free_dynlist(bt_tree *tree)
{
free_bt_internal(tree->racine);
free(tree);
}


int optimize_dynlist(bt_tree *tree)
{
bt_node **list;
int num = 0;

if(tree->total <= 2) return 0;
list = (bt_node **)malloc(tree->total * sizeof(bt_node *));
if(list == NULL) return 1;
if( calc_list(tree->racine, list) != list + tree->total) { 
	/* dynlist incoherente */
	free(list);
	return 2;
	}
tree->racine = add_mid_node(0, tree->total - 1, list);
tree->racine->pere = NULL;
free(list);
return 0;
}


static bt_node **calc_list(bt_node *noeud, bt_node **list)
{
if(noeud == NULL) return list;
list = calc_list(noeud->before, list);
*(list++) = noeud;
list = calc_list(noeud->after, list);
return list;
}


static bt_node *find_node_internal(char *name, bt_node *noeud, 
	int create_if_not_found, int *total, int *p_maxprof)
{
int test;

if(noeud == NULL) {
	int l;
	if( !create_if_not_found ) return NULL;
	noeud = (bt_node *)calloc(1, sizeof(bt_node));
	if(noeud == NULL) return NULL;
	l = strlen(name) + 1;
	noeud->name = (char *)malloc(l);
	if(noeud->name == NULL) return NULL;
	memcpy(noeud->name, name, l);
	(*total)++;
	return noeud;
	}

test = strcmp(name, noeud->name);
if(test == 0) return noeud;
(*p_maxprof)++;
if(test > 0) {
	if(noeud->after == NULL) {
		if(create_if_not_found) {
			noeud->after = find_node_internal(name, NULL, TRUE, 
				total, p_maxprof);
			if(noeud->after != NULL) noeud->after->pere = noeud;
			}
		return noeud->after;
		}
	else	return find_node_internal(name, noeud->after,
					create_if_not_found, total, p_maxprof);
	}
else	{
	if(noeud->before == NULL) {
		if(create_if_not_found) {
			noeud->before =find_node_internal(name, NULL, TRUE, 
				total, p_maxprof);
			if(noeud->before != NULL) noeud->before->pere = noeud;
			}
		return noeud->before;
		}
	else	return find_node_internal(name, noeud->before, 
					create_if_not_found, total, p_maxprof);
	}
}


int remove_dynlist(char *name, bt_tree *arbre)
{
bt_node *anc, *noeud;
int maxprof = 0;

anc = NULL;
noeud = find_node_internal(name, arbre->racine, FALSE, NULL, &maxprof);
if(noeud != NULL) {
	remove_node_internal(arbre, noeud);
	return TRUE;
	}
else return FALSE;
}


static void remove_node_internal(bt_tree *tree, bt_node *noeud)
{
bt_node *desc, *anc, *n_noeud;

desc = noeud->before;
anc = noeud->pere;
if(desc == NULL) {
	if(anc != NULL) {
		if(anc->before == noeud) anc->before = noeud->after;
		else	anc->after = noeud->after;
		}
	else	tree->racine = noeud->after;
	if(noeud->after != NULL) noeud->after->pere = anc;
	}
else	{
	n_noeud = desc;
	while(n_noeud->after != NULL) {
		n_noeud = n_noeud->after;
		}
	if(n_noeud != desc) {
		n_noeud->pere->after = n_noeud->before;
		if(n_noeud->before != NULL) 
			n_noeud->before->pere = n_noeud->pere;
		}
	if(anc != NULL) {
		if(anc->before == noeud) anc->before = n_noeud;
		else	anc->after = n_noeud;
		}
	else	tree->racine = n_noeud;
	n_noeud->pere = anc;
	if(n_noeud != desc) {
		n_noeud->before = desc;
		desc->pere = n_noeud;
		}
	n_noeud->after = noeud->after;
	if(n_noeud->after != NULL) n_noeud->after->pere = n_noeud;
	}
free(noeud->name);
free(noeud);
(tree->total)--;
}


static void free_bt_internal(bt_node *racine)
{
if(racine == NULL) return;
free_bt_internal(racine->before);
free_bt_internal(racine->after);
free(racine->name);
free(racine);
}


static bt_node *add_mid_node(int debut, int fin, bt_node **list)
{
int milieu = (debut + fin)/2;

if(debut > fin) return NULL;
list[milieu]->before = add_mid_node(debut, milieu - 1, list);
if(list[milieu]->before != NULL) list[milieu]->before->pere = list[milieu];
list[milieu]->after = add_mid_node(milieu + 1, fin, list);
if(list[milieu]->after != NULL) list[milieu]->after->pere = list[milieu];
return list[milieu];
}

/* gdynlist: liste dynamique pour objet generique: nom + extra_data
*/
typedef struct _g_bt_node {
	char *name;
	struct _g_bt_node *before;
	struct _g_bt_node *after;
	struct _g_bt_node *pere;
	void *extra_data;
	} g_bt_node;

typedef struct _g_bt_tree {
	g_bt_node *racine;
	int total;
	} g_bt_tree;

/* defined functions */
char *find_g_dynlist(char *name, g_bt_tree *tree, int create_if_not_found, void **node);
void *dynlist_get_extra(void *noeud);
void dynlist_set_extra(void *noeud, void *extra);


/* private stuff */
static g_bt_node *g_find_node_internal(char *name, g_bt_node *racine,  
	int create_if_not_found, int *total, int *p_maxprof);


char *find_g_dynlist(char *name, g_bt_tree *tree, int create_if_not_found, void **node)
{
g_bt_node *noeud;
int maxprof = 0;

noeud = g_find_node_internal(name, tree->racine, 
	create_if_not_found, &(tree->total), &maxprof );
if(tree->racine == NULL) tree->racine = noeud;
if(noeud == NULL) return NULL;
if(create_if_not_found && maxprof >= MAX_PROF) {
	if( optimize_dynlist( (bt_tree *)tree ) == 2 ) return NULL;
	}
if(node != NULL) *node = noeud;
return noeud->name;
}


void *dynlist_get_extra(void *noeud)
{
return ((g_bt_node *)noeud)->extra_data;
}


void dynlist_set_extra(void *noeud, void *extra)
{
((g_bt_node *)noeud)->extra_data = extra;
}


static g_bt_node *g_find_node_internal(char *name, g_bt_node *noeud,  
	int create_if_not_found, int *total, int *p_maxprof)
{
int test;


if(noeud == NULL) {
	int l;
	if( !create_if_not_found ) return NULL;
	noeud = (g_bt_node *)calloc(1, sizeof(g_bt_node));
	if(noeud == NULL) return NULL;
	l = strlen(name) + 1;
	noeud->name = (char *)malloc(l);
	if(noeud->name == NULL) return NULL;
	memcpy(noeud->name, name, l);
	(*total)++;
	return noeud;
	}

test = strcmp(name, noeud->name);
if(test == 0) return noeud;
(*p_maxprof)++;
if(test > 0) {
	if(noeud->after == NULL) {
		if(create_if_not_found) {
			noeud->after = g_find_node_internal(name, NULL, TRUE, 
				total, p_maxprof);
			if(noeud->after != NULL) noeud->after->pere = noeud;
			}
		return noeud->after;
		}
	else	return g_find_node_internal(name, noeud->after, 
					create_if_not_found, total, p_maxprof);
	}
else	{
	if(noeud->before == NULL) {
		if(create_if_not_found) {
			noeud->before = g_find_node_internal(name, NULL, TRUE, 
				total, p_maxprof);
			if(noeud->before != NULL) noeud->before->pere = noeud;
			}
		return noeud->before;
		}
	else	return g_find_node_internal(name, noeud->before,  
					create_if_not_found, total, p_maxprof);
	}
}

