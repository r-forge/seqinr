#include "dir_acnuc.h"

/* prototype */
void bit1(int *, int);

void lngbit(int point, int *blist)
{
int i, last;
memset(blist,0,lenw*sizeof(int));
last = read_first_rec(klng, NULL);
while(point != 0 && point <= last) {
	readlng(point); point=plng->next;
	for(i=0; i<SUBINLNG; i++) 
		if(plng->sub[i]) bit1(blist,plng->sub[i]);
	}
} 
