#include <sys/times.h>
#include <limits.h>

/* usage:
#include <sys/times.h>
declarer:
clock_t elapsed=0, cpu=0;
avant mesure faire:
debut_timing();
apres mesure faire:
cumul_timing(&elapsed, &cpu); 
ce qui augmente les 2 arguments du temps depuis appel a debut_timing
pour obtenir des secondes faire:
conv_timing_to_secnds(&elapsed, &cpu); 
*/

static struct tms buffer;
static clock_t depart_elapsed, depart_cpu;


void debut_timing(void)
{
depart_elapsed = times(&buffer);

depart_cpu = buffer.tms_utime+buffer.tms_stime;
}

void cumul_timing(clock_t *elapsed, clock_t *cpu)
{
clock_t current_cpu, current_elapsed;

current_elapsed = times(&buffer);

current_cpu = buffer.tms_utime+buffer.tms_stime;
*cpu += (current_cpu - depart_cpu);

*elapsed += (current_elapsed - depart_elapsed);
}

void conv_timing_to_secnds(clock_t *elapsed, clock_t *cpu)
{
*elapsed /= CLK_TCK;
*cpu /= CLK_TCK;
}


/* pour tester 

main()
{
int tab[2250], i;
clock_t elap, cpu;
elap=cpu=0;

for(i=0; i<10000; i++){
	debut_timing();
	non(tab,tab,2250);
	cumul_timing(&elap,&cpu);
	}
conv_timing_to_secnds(&elap,&cpu);
printf("elapsed=%d cpu=%d\n",elap,cpu);
}
*/
