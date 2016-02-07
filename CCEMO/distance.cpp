# include "global.h"

void fillCrowdingDistance(int count, int frontSize, list *elite){
	int *distance;
	list *temp;
	int i, j;
	int missing = nArch - count;
	while(missing<frontSize){
		assignCrowdingDistanceList(elite->child, frontSize);
		distance = (int *)calloc(frontSize,sizeof(int));
		temp = elite->child;
		for(j=0; j<frontSize; j++){
			distance[j] = temp->index;
			temp = temp->child;
		}
// 		quickSortDistance(distance, frontSize);
		double min_dist=repertoryDensity[distance[0]];
		int min_ind=distance[0];
		for(int k=1;k<frontSize;k++)
		{
			if(repertoryDensity[distance[k]]<min_dist)
			{
				min_dist=repertoryDensity[distance[k]];
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
	for(i=count, j=frontSize-1; i<nArch; i++, j--){
		copyToArchiveFromRepertory(i,distance[j]);
	}
	free(distance);
	return;
}

void assignCrowdingDistanceList(list *lst, int frontSize){
	int **arrayFx;
	int *distance;
	int i, j;
	list *temp;
	temp = lst;
	if(frontSize==1){
		repertoryDensity[lst->index] = INF;
		return;
	}
	if(frontSize==2){
		repertoryDensity[lst->index] = INF;
		repertoryDensity[lst->child->index] = INF;
		return;
	}
	arrayFx = (int **)calloc(nObj,sizeof(int*));
	distance = (int *)calloc(frontSize,sizeof(int));
	for(i=0; i<nObj; i++){
		arrayFx[i] = (int *)calloc(frontSize,sizeof(int));
	}
	for(j=0; j<frontSize; j++){
		distance[j] = temp->index;
		temp = temp->child;
	}
	assignCrowdingDistance(distance, arrayFx, frontSize);
	free(distance);
	for(i=0; i<nObj; i++){
		free(arrayFx[i]);
	}
	free(arrayFx);
	return;
}

void assignCrowdingDistanceIndexes(int c1, int c2){
	int **arrayFx;
	int *distance;
	int i, j;
	int frontSize;
	frontSize = c2-c1+1;
	if(frontSize==1){
		repertoryDensity[c1] = INF;
		return;
	}
	if(frontSize==2){
		repertoryDensity[c1] = INF;
		repertoryDensity[c2] = INF;
		return;
	}
	arrayFx = (int **)calloc(nObj,sizeof(int*));
	distance = (int *)calloc(frontSize,sizeof(int));
	for(i=0; i<nObj; i++){
		arrayFx[i] = (int *)calloc(frontSize,sizeof(int));
	}
	for(j=0; j<frontSize; j++){
		distance[j] = c1++;
	}
	assignCrowdingDistance(distance, arrayFx, frontSize);
	free(distance);
	for(i=0; i<nObj; i++){
		free (arrayFx[i]);
	}
	free(arrayFx);
	return;
}

void assignCrowdingDistance(int *distance, int **arrayFx, int frontSize){
	int i, j;
	for(i=0; i<nObj; i++){
		for(j=0; j<frontSize; j++){
			arrayFx[i][j] = distance[j];
		}
		quickSortFrontObj(i, arrayFx[i], frontSize);
	}
	for(j=0; j<frontSize; j++){
		repertoryDensity[distance[j]] = 0.0;
	}
	for(i=0; i<nObj; i++){
		repertoryDensity[arrayFx[i][0]] = INF;
	}
	for(i=0; i<nObj; i++){
		for(j=1; j<frontSize-1; j++){
			if(repertoryDensity[arrayFx[i][j]] != INF){
				if(repertFit[(arrayFx[i][frontSize-1])*nObj+i] == repertFit[(arrayFx[i][0])*nObj+i]){
					repertoryDensity[arrayFx[i][j]] += 0.0;
				}else{
					repertoryDensity[arrayFx[i][j]] += (repertFit[(arrayFx[i][j+1])*nObj+i] - repertFit[(arrayFx[i][j-1])*nObj+i])/(repertFit[(arrayFx[i][frontSize-1])*nObj+i] - repertFit[(arrayFx[i][0])*nObj+i]);
				}
			}
		}
	}

	for(j=0; j<frontSize; j++){
		if (repertoryDensity[distance[j]] != INF){
			repertoryDensity[distance[j]] = (repertoryDensity[distance[j]])/nObj;
		}
	}
	return;
}
