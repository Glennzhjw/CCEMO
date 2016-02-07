/*
  selection.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
  
#include "CCGDE3.h"

void CCGDE3::addSolution(Population *xGm, Individual *individual){
	copyIndividual(individual,&(xGm->individual[counter]));
	counter++;
}

void CCGDE3::selection(Species *species){
	
	int result;	
	int i;
	counter = 0;
	for(i=0;i<NP;i++){
		result = dominanceComparator(&(species->xG->individual[i]), &(species->uG->individual[i])) ;
		if(result == 1){ /*  xG[i] dominates uG[i] */
			addSolution(species->xGm,&(species->xG->individual[i]));
		}else if(result == -1){ /* uG[i] dominates xG[i] */
			addSolution(species->xGm,&(species->uG->individual[i]));
		}else{ /* non comparable solutions */
			addSolution(species->xGm,&(species->xG->individual[i]));
			addSolution(species->xGm,&(species->uG->individual[i]));
		}
	}
}

void CCGDE3::nextGeneration(Species *species){
    int result;
    int i,j;
    int end;
    int frontSize;
    int popSize;
    int rank = 1;
    list *pool;
    list *elite;
    list *temp1, *temp2;
    pool = createList(-1);
    elite = createList(-1);
    frontSize = 0;
    popSize = 0;
    
    temp1 = pool;
    for (i=0; i<counter; i++){
        insert(temp1,i);
        temp1 = temp1->child;
    }
    i=0;
    do{
        temp1 = pool->child;
        insert(elite, temp1->index);
        frontSize = 1;
        temp2 = elite->child;
        temp1 = deleteNode(temp1);
        temp1 = temp1->child;
        do{
            temp2 = elite->child;
            if (temp1==NULL){
                break;
            }
            do{
                end = 0;
                result = dominanceComparator(&(species->xGm->individual[temp1->index]), &(species->xGm->individual[temp2->index]));
                if (result == 1){
                    insert(pool, temp2->index);
                    temp2 = deleteNode(temp2);
                    frontSize--;
                    temp2 = temp2->child;
                }
                if (result == 0){
                    temp2 = temp2->child;
                }
                if (result == -1){
                    end = 1;
                }
            }while ((end != 1) && (temp2 != NULL));
            
            if(result == 0 || result == 1){
                insert(elite, temp1->index);
                frontSize++;
                temp1 = deleteNode(temp1);
            }
            temp1 = temp1->child;
        }while(temp1 != NULL);
        
        if(rank == 1){
			if(frontSize <= NP){
				species->frontSize = frontSize;
			}else{
				species->frontSize = NP;
			}
		}
        temp2 = elite->child;
        j=i;
        if((popSize + frontSize)<= NP){
            do{
                copyIndividual(&species->xGm->individual[temp2->index], &species->xG->individual[i]);
                species->xG->individual[i].rank = rank;
                popSize+=1;
                temp2 = temp2->child;
                i+=1;
            }while(temp2 != NULL);
//             assignCrowdingDistanceIndexes(species->xG, j, i-1);
            rank+=1;
        }else{
            fillCrowdingDistance(species->xGm, species->xG, i, frontSize, elite);
            popSize = NP;
            for (j=i; j<NP; j++){
                species->xG->individual[j].rank = rank;
            }
        }
        temp2 = elite->child;
        do{
            temp2 = deleteNode(temp2);
            temp2 = temp2->child;
        }while (elite->child !=NULL);
    }while(popSize < NP);
    
    deleteList(pool);
    deleteList(elite);
    return;
}

void CCGDE3::nonDominatedSorting(Population *xG, int size, int sizeNP){
    int result;
    int i;
    int end;
    int frontSize;
    int popSize;
    int rank = 1;
    list *pool;
    list *elite;
    list *temp1, *temp2;
    pool = createList(-1);
    elite = createList(-1);
    frontSize = 0;
    popSize = 0;
    
    temp1 = pool;
    for (i=0; i<size; i++){
        insert(temp1,i);
        temp1 = temp1->child;
    }
    do{
        temp1 = pool->child;
        insert(elite, temp1->index);
        frontSize = 1;
        temp2 = elite->child;
        temp1 = deleteNode(temp1);
        temp1 = temp1->child;
        do{
            temp2 = elite->child;
            if (temp1==NULL){
                break;
            }
            do{
                end = 0;
                result = dominanceComparator(&(xG->individual[temp1->index]), &(xG->individual[temp2->index]));
                if (result == 1){
                    insert(pool, temp2->index);
                    temp2 = deleteNode(temp2);
                    frontSize--;
                    temp2 = temp2->child;
                }
                if (result == 0){
                    temp2 = temp2->child;
                }
                if (result == -1){
                    end = 1;
                }
            }while ((end != 1) && (temp2 != NULL));
            
            if(result == 0 || result == 1){
                insert(elite, temp1->index);
                frontSize++;
                temp1 = deleteNode(temp1);
            }
            temp1 = temp1->child;
        }while(temp1 != NULL);
        
        
        temp2 = elite->child;
        if((popSize + frontSize)<= sizeNP){
            do{
                xG->individual[temp2->index].rank = rank;
                popSize+=1;
                temp2 = temp2->child;                
            }while(temp2 != NULL);
            rank+=1;
        }
        temp2 = elite->child;
        do{
            temp2 = deleteNode(temp2);
            temp2 = temp2->child;
        }while (elite->child !=NULL);
    }while(popSize < sizeNP);
    
    deleteList(pool);
    deleteList(elite);
    return;
}

void CCGDE3::fillCrowdingDistance(Population *xGm, Population *xG, int count, int frontSize, list *elite){
    int *distance;
    list *temp;
    int i, j;
    int missing = NP - count;
    while(missing<frontSize){
		assignCrowdingDistanceList(xGm, elite->child, frontSize);
		distance = (int *)calloc(frontSize,sizeof(int));
		temp = elite->child;
		for(j=0; j<frontSize; j++){
			distance[j] = temp->index;
			temp = temp->child;
		}
// 		quickSortDistance(xGm, distance, frontSize);
		double min_dist=xGm->individual[distance[0]].crowdDist;
		int min_ind=distance[0];
		for(int k=1;k<frontSize;k++)
		{
			if(xGm->individual[distance[k]].crowdDist<min_dist)
			{
				min_dist=xGm->individual[distance[k]].crowdDist;
				min_ind=distance[k];
			}
		}
// 		deleteInd(elite,distance[0]);
		deleteInd(elite,min_ind);
		frontSize--;
		if(missing<frontSize){
			free(distance);
		}
	}
    for(i=count, j=frontSize-1; i<NP; i++, j--){
        copyIndividual(&xGm->individual[distance[j]], &xG->individual[i]);
    }
    free(distance);
    return;
}
