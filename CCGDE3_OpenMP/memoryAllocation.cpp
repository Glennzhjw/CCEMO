/*
  memoryAllocation.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
  
#include "CCGDE3.h"

void CCGDE3::memoryAllocation(){
	int i;
	
	if((minLimit = (double *)calloc(D,sizeof(double))) == NULL){
			printf("ERROR!! --> calloc: no memory for vector\n");
			exit(1);
	}
    if((maxLimit = (double *)calloc(D,sizeof(double))) == NULL){
			printf("ERROR!! --> calloc: no memory for vector\n");
			exit(1);
	}
	
	if((species = (Species *)calloc(numSpecies,sizeof(Species))) == NULL){
		printf("ERROR!! --> calloc: no memory for species\n");
		exit(1);
	}
	
	for(i=0;i<numSpecies;i++){
		
		if((species[i].indexes = (int *)calloc(DSp,sizeof(int))) == NULL){
				printf("ERROR!! --> calloc: no memory for population\n");
				exit(1);
		}
		
		if((species[i].xG = (Population *)calloc(1,sizeof(Population))) == NULL){
				printf("ERROR!! --> calloc: no memory for population\n");
				exit(1);
		}
		if((species[i].uG = (Population *)calloc(1,sizeof(Population))) == NULL){
				printf("ERROR!! --> calloc: no memory for population\n");
				exit(1);
		}
		if((species[i].xGm = (Population *)calloc(1,sizeof(Population))) == NULL){
				printf("ERROR!! --> calloc: no memory for population\n");
				exit(1);
		}
		populationMemoryAllocation(species[i].xG,NP);
		populationMemoryAllocation(species[i].uG,NP);
		populationMemoryAllocation(species[i].xGm,2*NP);
	}

    return;
}


void CCGDE3::populationMemoryAllocation(Population *population, int size){
	int i;
    if((population->individual = (Individual *)calloc(size,sizeof(Individual))) == NULL){
		printf("ERROR!! --> calloc: no memomry for new population\n");
		exit(1);
	}
    for (i=0; i<size; i++){
        individualMemoryAllocation(&(population->individual[i]));
    }
    return;
}

void CCGDE3::individualMemoryAllocation(Individual *individual){
    
	if((individual->x = (double *)calloc(DSp,sizeof(double))) == NULL){
		printf("ERROR!! --> calloc: no memory for vector x\n");
		exit(1);
	}
    
    if((individual->fx = (double *)calloc(Nobj,sizeof(double))) == NULL){
			printf("ERROR!! --> calloc: no memory for vector fx\n");
			exit(1);
	}
	
	if((individual->cx = (double *)calloc(D,sizeof(double))) == NULL){
		printf("ERROR!! --> calloc: no memory for vector x\n");
		exit(1);
	}
	
    return;
}


void CCGDE3::freeMemory(){
	int i;
	
	free(minLimit);
    free(maxLimit);

    for(i=0;i<numSpecies;i++){
		freePopulationMemory(species[i].xG,NP);
		freePopulationMemory(species[i].uG,NP);
		freePopulationMemory(species[i].xGm,2*NP);
		free(species[i].xG);
		free(species[i].uG);
		free(species[i].xGm);
	}
	free(species);
	if(sizeSol!=0){
		freePopulationMemory(Solution,sizeSol);
		free(Solution);
	}

    return;
}


void CCGDE3::freePopulationMemory(Population *population, int size){
    int i;
    for (i=0; i<size; i++){
        freeIndividualMemory(&(population->individual[i]));
    }
    free(population->individual);
    return;
}

void CCGDE3::freeIndividualMemory(Individual *individual){
    
    free(individual->x);    
    free(individual->fx);
    free(individual->cx);    
    return;
}

