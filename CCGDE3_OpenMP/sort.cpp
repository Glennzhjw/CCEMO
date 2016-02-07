/*
  sort.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
  
# include "CCGDE3.h"

void CCGDE3::quickSortFrontObj(Population *population, int objcount, int arrayFx[], int sizeArrayFx){
    qSortFrontObj(population, objcount, arrayFx, 0, sizeArrayFx-1);
    return;
}

void CCGDE3::qSortFrontObj(Population *population, int objcount, int arrayFx[], int left, int right){
    int index;
    int temp;
    int i, j;
    double pivot;
    if (left<right){
        index = rnd(left, right);
        temp = arrayFx[right];
        arrayFx[right] = arrayFx[index];
        arrayFx[index] = temp;
        pivot = population->individual[arrayFx[right]].fx[objcount];
        i = left-1;
        for (j=left; j<right; j++){
            if (population->individual[arrayFx[j]].fx[objcount] <= pivot){
                i+=1;
                temp = arrayFx[j];
                arrayFx[j] = arrayFx[i];
                arrayFx[i] = temp;
            }
        }
        index=i+1;
        temp = arrayFx[index];
        arrayFx[index] = arrayFx[right];
        arrayFx[right] = temp;
        qSortFrontObj(population, objcount, arrayFx, left, index-1);
        qSortFrontObj(population, objcount, arrayFx, index+1, right);
    }
    return;
}

void CCGDE3::quickSortDistance(Population *population, int *distance, int frontSize){
    qSortDistance(population, distance, 0, frontSize-1);
    return;
}

void CCGDE3::qSortDistance(Population *population, int *distance, int left, int right){
    int index;
    int temp;
    int i, j;
    double pivot;
    if (left<right){
        index = rnd(left, right);
        temp = distance[right];
        distance[right] = distance[index];
        distance[index] = temp;
        pivot = population->individual[distance[right]].crowdDist;
        i = left-1;
        for (j=left; j<right; j++){
            if (population->individual[distance[j]].crowdDist <= pivot){
                i+=1;
                temp = distance[j];
                distance[j] = distance[i];
                distance[i] = temp;
            }
        }
        index=i+1;
        temp = distance[index];
        distance[index] = distance[right];
        distance[right] = temp;
        qSortDistance(population, distance, left, index-1);
        qSortDistance(population, distance, index+1, right);
    }
    return;
}

int CCGDE3::indexWorstCrowdingDistance(Population *population,int size){
	int i,index;
	double aux;
	aux = population->individual[0].crowdDist;
	index = 0;
	for(i=0;i<size;i++){
		if(population->individual[i].crowdDist < aux){
			aux = population->individual[i].crowdDist;
			index = i;
		}
	}	
	return index;
}
