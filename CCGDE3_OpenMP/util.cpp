/*
  util.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
  
#include "CCGDE3.h"

void CCGDE3::copyIndividual(Individual *original, Individual *copy){
	int i;
	copy->rank = original->rank;
	copy->fitness = original->fitness;
	copy->restViol = original->restViol;
	copy->crowdDist = original->crowdDist;
	for(i=0; i<DSp;i++){
		copy->x[i] = original->x[i];
	}
	for(i=0; i<Nobj;i++){
		copy->fx[i] = original->fx[i];
	}
	for(i=0; i<D;i++){
		copy->cx[i] = original->cx[i];
	}
}


// double CCGDE3::timeval_diff(struct timeval *a, struct timeval *b){
//     return
//     (double)(a->tv_sec + (double)a->tv_usec/1000000) -
//     (double)(b->tv_sec + (double)b->tv_usec/1000000);
// }
