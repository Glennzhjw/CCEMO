/*
  dominanceComparator.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
  
#include "CCGDE3.h"

int CCGDE3::dominanceComparator(Individual *individual1, Individual *individual2){
	int i;

	/* Indicates if some objetive in solution 1 dominates the objetive in solution 2 */
	int dominates1 = 0;
    /* Indicates if some objetive in solution 2 dominates the objetive in solution 1 */
	int dominates2 = 0;
	int result;
	double value1, value2;    
	
	for(i = 0; i < Nobj; i++){
		value1 = individual1->fx[i];
		value2 = individual2->fx[i];
		if(value1 < value2){
            result = -1;
		}else if(value1 > value2){
            result = 1;
		}else{
            result = 0;
		}
		if (result == -1) {
            dominates1 = 1;
		}
		if (result == 1) {
            dominates2 = 1;
		}
	}
    
	if (dominates1 == dominates2){ /* non-dominated solutions */
		return 0;
	}
	if(dominates1 == 1){ /* solution1 dominates */
		return 1;
	}
	return -1;  /* solucion2 dominates */
    
}
