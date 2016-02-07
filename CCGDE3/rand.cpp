
#include "rand.h"


//	rand
double seed;
double oldrand[55];
int jrand;


/* Fisher¨CYates shuffle algorithm */
void shuffle(int *x,int size){
	int i,aux,k=0;
	for(i = size-1; i>0; i--){
		/* get a value between cero and i  */
		k = rnd(0,i); 
		/* exchange of values */
		aux = x[i]; 
		x[i] = x[k];
		x[k] = aux;
	}
}

int flip(float prob){
	if(randomperc() <= prob){
		return(1);
	}else{
		return(0);
	}
}

void randomize(int* n)
{
	int j1;
	for(j1=0; j1<=54; j1++)
	{
		oldrand[j1] = 0.0;
	}
	jrand=0;
	warmup_random (seed, n);
	return;
}

void warmup_random (double seed, int* n)
{
	int j1, ii;
	double new_random, prev_random;
	oldrand[54] = seed;
	new_random = 0.000000001;
	prev_random = seed;
	for(j1=1; j1<=54; j1++)
	{
		ii = (21*j1)%54;
		oldrand[ii] = new_random;
		new_random = prev_random-new_random;
		if(new_random<0.0)
		{
			new_random += 1.0;
		}
		prev_random = oldrand[ii];
	}
	if(n)
		for(int i=0;i<*n+3;i++)
			advance_random();
	advance_random ();
	advance_random ();
	advance_random ();
	jrand = 0;
	return;
}

void advance_random()
{
	int j1;
	double new_random;
	for(j1=0; j1<24; j1++)
	{
		new_random = oldrand[j1]-oldrand[j1+31];
		if(new_random<0.0)
		{
			new_random = new_random+1.0;
		}
		oldrand[j1] = new_random;
	}
	for(j1=24; j1<55; j1++)
	{
		new_random = oldrand[j1]-oldrand[j1-24];
		if(new_random<0.0)
		{
			new_random = new_random+1.0;
		}
		oldrand[j1] = new_random;
	}
}


double randomperc()
{
	jrand++;
	if(jrand>=55)
	{
		jrand = 1;
		advance_random();
	}
	return((double)oldrand[jrand]);
}

int rnd(int low, int high)
{
	int res;
	if (low >= high)
	{
		res = low;
	}
	else
	{
		res = low + (randomperc()*(high-low+1));
		if (res > high)
		{
			res = high;
		}
	}
	return (res);
}

double rndreal (double low, double high)
{
	return (low + (high-low)*randomperc());
}
