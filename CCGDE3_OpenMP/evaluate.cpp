/*
  evaluate.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
  
#include "CCGDE3.h"

void CCGDE3::evaluateSpeciesInitial(){
    int i;
    for (i=0; i<numSpecies; i++){        
        /* Populations evaluation */
        evaluatePopulationRandom(species[i].xG,i);
        /* Crowding distance assignation */
        assignCrowdingDistanceIndexes(species[i].xG,0,NP-1);        
        /* nonDominates Sorting */
        nonDominatedSorting(species[i].xG,NP,NP);
    }
    return;
}

void CCGDE3::evaluatePopulationRandom(Population *population, int indP){
    int i;
    for (i=0; i<NP; i++){
        evaluateIndividualRandom(&(population->individual[i]),indP);
    }
    return;
}

void CCGDE3::evaluateIndividualRandom(Individual *individual,int indP){        
    joinRepresentativesRandom(individual->cx,individual->x,indP);    
    EMO_TEST_SUITE::evaluate_problems(problemInstance,individual->cx, individual->fx, D,1,Nobj);   
    return;
}

void CCGDE3::evaluatePopulation(Population *population, int size, int indP){
	//	size==NP
    int i;
#pragma omp parallel for 
    for (i=0; i<size; i++){
        evaluateIndividual(&(population->individual[i]),indP);
	}
    return;
}

void CCGDE3::evaluateIndividual(Individual *individual, int indP){        
    joinRepresentativesBest(individual->cx,individual->x,indP);    
//     problem(individual->cx, individual->fx, D);    
	EMO_TEST_SUITE::evaluate_problems(problemInstance,individual->cx, individual->fx, D,1,Nobj);   
    return;
}

