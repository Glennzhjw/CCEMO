/*
  species.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
#include "CCGDE3.h"

/* Crear species in a random way */
void CCGDE3::createSpecies(){
	int *indexes;
	int i,j,k;
	
	indexes = (int*)calloc(D,sizeof(int));	
	
	/* Create indexes */
	for(i=0;i<D;i++){
		indexes[i]= i;
	}
	/* shuffle indexes */
	shuffle(indexes,D);
	/* Indexes assignation to each species */
	for(i=0,k=0;i<numSpecies;i++){
		for(j=0;j<DSp;j++,k++){
			species[i].indexes[j] = indexes[k];
		}
	}
	free(indexes);
}


void CCGDE3::joinRepresentativesRandom(double *x,double *indX,int indP){
	int i,j,k;
	/* Sort particles values in a vector of objetives in order to be evaluated */
	for(i=0;i<numSpecies;i++){
        if(i == indP){
            for(j=0;j<DSp;j++){
                x[species[i].indexes[j]] = indX[j];
            }
        }else{
            k = rnd(0,NP - 1);
            for(j=0;j<DSp;j++){
                x[species[i].indexes[j]] = species[i].xG->individual[k].x[j];
            }
        }
		
	}
}

void CCGDE3::joinRepresentativesBest(double *x,double *indX,int indP){
	int i,j,k;
	/* Sort particles values in a vector of objetives in order to be evaluated */
	for(i=0;i<numSpecies;i++){
        if(i == indP){
            for(j=0;j<DSp;j++){
                x[species[i].indexes[j]] = indX[j];
            }
        }else{                        
            do{
                k = rnd(0,NP - 1);
            }while(species[i].xG->individual[k].rank != 1);
            for(j=0;j<DSp;j++){
                x[species[i].indexes[j]] = species[i].xG->individual[k].x[j];
            }
            
        }
		
	}
}

void CCGDE3::generateSolutions(){
	int i,j,k;
	sizeSol=0;
	solutionCounter = 0;
	for(i=0;i<numSpecies;i++){
		sizeSol += species[i].frontSize;
	}

	if((Solution = (Population *)calloc(1,sizeof(Population))) == NULL){
		printf("ERROR!! --> calloc: no memory for population\n");
		exit(1);
	}
    
	populationMemoryAllocation(Solution,sizeSol);
    
	for(i=0,k=0;i<numSpecies;i++){
		for(j=0;j<NP;j++){
			if(species[i].xG->individual[j].rank == 1){
				copyIndividual(&species[i].xG->individual[j],&Solution->individual[k]);
				k++;
			}
		}
	}
       
    nonDominatedSorting(Solution,sizeSol,sizeSol);
}


