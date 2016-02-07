/*
  crossover.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
#include "CCGDE3.h"

/* ED/rand/1/bin */
void CCGDE3::ED_rand_1_bin(Species *species){
	int r1,r2,r3;
	int i,j,k;
#pragma omp parallel for private(j,k,r1,r2,r3)
	for(i = 0; i < NP; i++){		
		do{
			r1 = rnd(0,NP - 1);
		}while(r1 == i);
		do{
			r2 = rnd(0,NP - 1);
		}while(r2 == i || r2 == r1);
		do{
			r3 = rnd(0,NP - 1);
		}while(r3 == i || r3 == r1 || r3 == r2);

		j = rnd(0,DSp - 1);
		
		for(k = 0; k < DSp; k++){
			if(flip(CR) || k == j){
				species->uG->individual[i].x[k] = species->xG->individual[r1].x[k] + F*(species->xG->individual[r2].x[k] - species->xG->individual[r3].x[k]);
				
				if(species->uG->individual[i].x[k]<minLimit[species->indexes[k]]){
					species->uG->individual[i].x[k] = minLimit[species->indexes[k]];
				}
				if(species->uG->individual[i].x[k]>maxLimit[species->indexes[k]]){
					species->uG->individual[i].x[k] = maxLimit[species->indexes[k]];
				}
			}else{
				species->uG->individual[i].x[k] = species->xG->individual[i].x[k];
			}
		}				
	}
}