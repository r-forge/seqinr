#define lmot (8*sizeof(int))

/* prototypes of included functions */
int irbit(int *pdeblist, int deb, int fin);
void bit1(int *plist, int num);
void bit0(int *plist, int num);
int testbit(int *plist, int num);
int bcount(int *plist, int last);
void et(int *listet, int *list1, int *list2, int len);
void ou(int *listou, int *list1, int *list2, int len);
void non(int *listnon, int *list, int len);


/*
******* Pour tester toutes les fonctions *********
main()
{
int tab[100],i,j;
for(i=0;i<100;i++) tab[i]=0;
for(i=1;i<=64;i++) bit1(tab,i);
printf("%d %d %d %d\n",tab[0],tab[1],tab[2],bcount(tab,64));
for(i=1;i<=64;i++) {
	j=testbit(tab,i);
	if(!j)printf("pb %d\n",i);
	}
if(testbit(tab,65))printf("pb 65\n");
printf("irbit=%d\n",irbit(tab,0,3200));
for(i=64;i>=1;i--) {
	if(bcount(tab,i)!=i)printf("Pb bcount\n");
	}
for(i=1;i<=64;i++) {
	bit0(tab,i); 
	j=irbit(tab,i,64);
	if(j!=i+1 && i<64)printf("pb irbit %d\n",i);
	j=irbit(tab,0,i);
	if(j!=0 )printf("pb irbit %d\n",i);
	if(bcount(tab,128)!=64-i)printf("pb bcount\n");
	}
printf("%d %d %d\n",tab[0],tab[1],tab[2]);
printf("irbit=%d\n",irbit(tab,0,3200));
printf(" ffs was NOT used\n");
}
*/


int irbit(int *pdeblist, int deb, int fin)
{
unsigned int mot;
int *plist, *debw, *finw, retval;

if(deb>=fin) return 0;
finw = pdeblist+(fin-1)/lmot;
debw = pdeblist+deb/lmot;
plist = debw-1;
do	{
	mot= *(++plist);
	if(plist==debw) {
		mot &= ( (~0)<<(deb%lmot) );
		}
	if( plist==finw && (fin%lmot)!=0 ) {
		mot &= ( ~((~0)<<(fin%lmot)) );
		}
	}
while (plist<finw && mot==0);
if(!mot)return 0;
retval = (plist-pdeblist)*lmot;
retval++;
while((mot&1)==0) {
	mot>>=1;
	retval++;
	}
return retval;
}

void bit1(int *plist, int num)
{
num--;
plist+=(num/lmot);
*plist |= (1<<(num%lmot));
}

void bit0(int *plist, int num)
{
num--;
plist+=(num/lmot);
*plist &=  ~(1<<(num%lmot));
}

int testbit(int *plist, int num)
{
num--;
plist += (num/lmot);
return (*plist) & (1<<(num%lmot));
}

int bcount(int *plist, int last)
{
unsigned int mot;
int *lasta,retval=0,i;
lasta=plist+(last-1)/lmot;
do	{
	mot= *plist;
	if( plist==lasta && (last%lmot)!=0 )mot &=  ~((~0)<<(last%lmot));
	while(mot) {
		if(mot&1)retval++;
		mot >>= 1;
		}
	}
while((plist++)<lasta);
return retval;
}

void et(int *listet, int *list1, int *list2, int len)
{
int i;
for(i=0; i< len; i++)
	listet[i]= list1[i] & list2[i];
}

void ou(int *listou, int *list1, int *list2, int len)
{
int i;
for(i=0; i< len; i++)
	listou[i]= list1[i] | list2[i];
}

void non(int *listnon, int *list, int len)
{

int i;
for(i=0; i< len; i++)
	listnon[i]= ~list[i];
}
