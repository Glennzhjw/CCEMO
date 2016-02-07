/*
  rand.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  

#include "CCGDE3.h"

/* Fisher–Yates shuffle algorithm */
void CCGDE3::shuffle(int *x,int size){
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

int CCGDE3::flip(float prob){
	if(randomperc() <= prob){
		return(1);
	}else{
		return(0);
	}
}

void CCGDE3::randomize()
{
    int j1;
    for(j1=0; j1<=54; j1++)
    {
        oldrand[j1] = 0.0;
    }
    jrand=0;
    warmup_random (seed);
    return;
}

void CCGDE3::warmup_random (double seed)
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
    advance_random ();
    advance_random ();
    advance_random ();
    jrand = 0;
    return;
}

void CCGDE3::advance_random()
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


double CCGDE3::randomperc()
{
    jrand++;
    if(jrand>=55)
    {
        jrand = 1;
        advance_random();
    }
    return((double)oldrand[jrand]);
}

int CCGDE3::rnd(int low, int high)
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

double CCGDE3::rndreal (double low, double high)
{
    return (low + (high-low)*randomperc());
}
