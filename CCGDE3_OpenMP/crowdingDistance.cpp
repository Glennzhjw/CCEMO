/*
  crowdingDistance.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
#include "CCGDE3.h"

void CCGDE3::assignCrowdingDistanceList(Population *population, list *lst, int frontSize){
    int **arrayFx;
    int *distance;
    int i, j;
    list *temp;
    temp = lst;
    if(frontSize==1){
        population->individual[lst->index].crowdDist = INF;
        return;
    }
    if(frontSize==2){
        population->individual[lst->index].crowdDist = INF;
        population->individual[lst->child->index].crowdDist = INF;
        return;
    }
    arrayFx = (int **)calloc(Nobj,sizeof(int*));
    distance = (int *)calloc(frontSize,sizeof(int));
    for(i=0; i<Nobj; i++){
        arrayFx[i] = (int *)calloc(frontSize,sizeof(int));
    }
    for(j=0; j<frontSize; j++){
        distance[j] = temp->index;
        temp = temp->child;
    }
    assignCrowdingDistance(population, distance, arrayFx, frontSize);
    free(distance);
    for(i=0; i<Nobj; i++){
        free(arrayFx[i]);
    }
    free(arrayFx);
    return;
}

void CCGDE3::assignCrowdingDistanceIndexes(Population *population, int c1, int c2){
    int **arrayFx;
    int *distance;
    int i, j;
    int frontSize;
    frontSize = c2-c1+1;
    if(frontSize==1){
        population->individual[c1].crowdDist = INF;
        return;
    }
    if(frontSize==2){
        population->individual[c1].crowdDist = INF;
        population->individual[c2].crowdDist = INF;
        return;
    }
    arrayFx = (int **)calloc(Nobj,sizeof(int*));
    distance = (int *)calloc(frontSize,sizeof(int));
    for(i=0; i<Nobj; i++){
        arrayFx[i] = (int *)calloc(frontSize,sizeof(int));
    }
    for(j=0; j<frontSize; j++){
        distance[j] = c1++;
    }
    assignCrowdingDistance(population, distance, arrayFx, frontSize);
    free(distance);
    for(i=0; i<Nobj; i++){
        free (arrayFx[i]);
    }
    free(arrayFx);
    return;
}

void CCGDE3::assignCrowdingDistance(Population *population, int *distance, int **arrayFx, int frontSize){
    int i, j;
    for(i=0; i<Nobj; i++){
        for(j=0; j<frontSize; j++){
            arrayFx[i][j] = distance[j];
        }
        quickSortFrontObj(population, i, arrayFx[i], frontSize);
    }
    for(j=0; j<frontSize; j++){
        population->individual[distance[j]].crowdDist = 0.0;
    }
    for(i=0; i<Nobj; i++){
        population->individual[arrayFx[i][0]].crowdDist = INF;
    }
    for(i=0; i<Nobj; i++){
        for(j=1; j<frontSize-1; j++){
            if(population->individual[arrayFx[i][j]].crowdDist != INF){
                if(population->individual[arrayFx[i][frontSize-1]].fx[i] == population->individual[arrayFx[i][0]].fx[i]){
                    population->individual[arrayFx[i][j]].crowdDist += 0.0;
                }else{
                    population->individual[arrayFx[i][j]].crowdDist += (population->individual[arrayFx[i][j+1]].fx[i] - population->individual[arrayFx[i][j-1]].fx[i])/(population->individual[arrayFx[i][frontSize-1]].fx[i] - population->individual[arrayFx[i][0]].fx[i]);
                }
            }
        }
    }

    for(j=0; j<frontSize; j++){
        if (population->individual[distance[j]].crowdDist != INF){
            population->individual[distance[j]].crowdDist = (population->individual[distance[j]].crowdDist)/Nobj;
        }
    }
    return;
}
