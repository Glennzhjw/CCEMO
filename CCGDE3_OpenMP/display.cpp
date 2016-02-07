/*
  display.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/

# include "CCGDE3.h"

void CCGDE3::displayPopulation(Population *population, int size){
	int i, j;
    
    printf("x: \n");
    for(i=0; i<size; i++){
        for(j=0; j<D; j++){
            printf("%d) %f ",i+1,population->individual[i].x[j]);
        }
        printf("\n");
		
    }
    
    printf("f(x):\n");
    for(i=0; i<size; i++){
        for(j=0; j<Nobj; j++){
            printf("%d) %f ",i+1,population->individual[i].fx[j]);
        }
        printf("\n");
		
    }
    
    printf("\nCD:\n");
    for(i=0; i<size; i++){        
        printf("%d) %f ",i+1,population->individual[i].crowdDist);        
        printf("\n");
		
    }
    
    return;
}

void CCGDE3::displayIndividual(Individual *individual){
	int j;
	printf("\nx: \n");
	for(j=0; j<D; j++){
		printf(" %f ",individual->x[j]);
    }    
    printf("\n");
    printf("f(x):\n");
    for(j=0; j<Nobj; j++){
		printf(" %f ",individual->fx[j]);
    }
    printf("\n");
	
}

void CCGDE3::exportPopulation(Population *population, int size, int gener, FILE *fpt){
	int i, j;
	fprintf(fptm,"\nGen:%d\n",gener);
    fprintf(fpt,"------- f(x) --------------\t--------- x ------------\n");
    for(i=0; i<size; i++){        
        for(j=0; j<Nobj; j++){
            fprintf(fpt,"%e\t",population->individual[i].fx[j]);
        }
        if (D != 0){
          for (j=0; j<D; j++){
              fprintf(fpt,"%e\t",population->individual[i].x[j]);
          }
        }
    }
    return;
}


void CCGDE3::exportPopulationObjetivesValues(Population *population, int size, FILE *fpt){
	int i, j;    
    for(i=0; i<size; i++){  
		fprintf(fpt,"\n");      
        for(j=0; j<Nobj; j++){
            fprintf(fpt,"%e\t",population->individual[i].fx[j]);
        }
    }
    return;
}

void CCGDE3::exportNonDominatedPopulationObjetivesValues(Population *population, int size, FILE *fpt){
	int i, j;    
	sizeSolND = 0;
    for(i=0; i<size; i++){
		if(population->individual[i].rank==1){ 
			sizeSolND++;
			for(j=0; j<Nobj; j++){
				fprintf(fpt,"%e\t",population->individual[i].fx[j]);
			}
			fprintf(fpt,"\n");
		}
    }
    return;
}

void CCGDE3::exportPopulationSolutionValues(Population *population, int size, FILE *fpt){
	int i, j;    
    for(i=0; i<size; i++){   
        for(j=0; j<D; j++){
            fprintf(fpt,"%e\t",population->individual[i].x[j]);
        }
        fprintf(fpt,"\n");
    }
    return;
}

void CCGDE3::exportNonDominatedPopulationSolutionValues(Population *population, int size, FILE *fpt){
	int i, j;    
    for(i=0; i<size; i++){  
		if(population->individual[i].rank==1){ 			
			for(j=0; j<D; j++){
				fprintf(fpt,"%e\t",population->individual[i].x[j]);
			}
			fprintf(fpt,"\n");
		}
    }
    return;
}
