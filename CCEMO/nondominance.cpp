# include "global.h"
#include <math.h>

void refineRepertory_generateArchive()
{
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
	for (i=0; i<nRep; i++){
		insert(temp1,i);
		temp1 = temp1->child;
	}
	i=0;
	do{
		temp1 = pool->child;
		if (temp1==NULL){
			break;
		}
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
				result = dominanceComparator(&(repertFit[(temp1->index)*nObj]), &(repertFit[(temp2->index)*nObj]));
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

// 		if(rank == 1){
// 			if(frontSize <= nArch){
// 				cnArch = frontSize;
// 			}else{
// 				cnArch = nArch;
// 			}
// 		}
		temp2 = elite->child;
		j=i;
		if((popSize + frontSize) <= nArch){
			do{
				copyToArchiveFromRepertory(i,temp2->index);
				archiveClass[i]=rank;
				popSize+=1;
				temp2 = temp2->child;
				i+=1;
			}while(temp2 != NULL);
// 			assignCrowdingDistanceIndexes(j,i-1);
			rank+=1;
		}else{
			fillCrowdingDistance(i, frontSize, elite);
			popSize = nArch;
			for (j=i; j<popSize; j++){
				archiveClass[j] = rank;
			}
		}
		temp2 = elite->child;
		do{
			temp2 = deleteNode(temp2);
			temp2 = temp2->child;
		}while (elite->child !=NULL);
	}while(popSize < nArch);

	if (nRep<nArch)
	{
		cnArch=nRep;
	} 
	else
	{
		cnArch=nArch;
	}

	deleteList(pool);
	deleteList(elite);
	return;
}

void refineRepertory_generateArchive_SDE()
{
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
	for (i=0; i<nRep; i++){
		insert(temp1,i);
		temp1 = temp1->child;
	}
	i=0;
	do{
		temp1 = pool->child;
		if (temp1==NULL){
			break;
		}
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
				result = dominanceComparator(&(repertFit[(temp1->index)*nObj]), &(repertFit[(temp2->index)*nObj]));
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

// 		if(rank == 1){
// 			if(frontSize <= nArch){
// 				cnArch = frontSize;
// 			}else{
// 				cnArch = nArch;
// 			}
// 		}
		temp2 = elite->child;
		j=i;
		if((popSize + frontSize) <= nArch){
			do{
				copyToArchiveFromRepertory(i,temp2->index);
				archiveClass[i]=rank;
				popSize+=1;
				temp2 = temp2->child;
				i+=1;
			}while(temp2 != NULL);
// 			assignCrowdingDistanceIndexes(j,i-1);
			rank+=1;
		}else{
			K_Neighbor_Nearest_SDE(i, frontSize, elite);
			popSize = nArch;
			for (j=i; j<popSize; j++){
				archiveClass[j] = rank;
			}
		}
		temp2 = elite->child;
		do{
			temp2 = deleteNode(temp2);
			temp2 = temp2->child;
		}while (elite->child !=NULL);
	}while(popSize < nArch);

	if (nRep<nArch)
	{
		cnArch=nRep;
	} 
	else
	{
		cnArch=nArch;
	}

	deleteList(pool);
	deleteList(elite);
	return;
}

void K_Neighbor_Nearest_SDE(int count, int frontSize, list *elite)
{
	int non_size = frontSize;
	int* non_indexes;

	non_indexes=(int*)malloc(non_size*sizeof(int));
	list* tmp=elite->child;
	for(int i=0;i<non_size;i++)
	{
		non_indexes[i]=tmp->index;
		tmp=tmp->child;
	}

	//���Ҳ���¼ÿ��Ŀ���ϵ����ֵ����Сֵ
	for (int i = 0; i<nObj; ++i)
	{
		fun_min[i] = INF;
		fun_max[i] = -100;
	}
	for (int i = 0; i<non_size; ++i)
	{
		for (int j = 0; j<nObj; ++j)
		{
			if (repertFit[non_indexes[i]*nObj+j] < fun_min[j])
				fun_min[j] = repertFit[non_indexes[i]*nObj+j];
			if (repertFit[non_indexes[i]*nObj+j] > fun_max[j])
				fun_max[j] = repertFit[non_indexes[i]*nObj+j];
		}
	}	

	//���������ռ�Ž⣬����shift-based��ŷ����þ������  
	double **dist;				//��¼������� 
	int **distIndex;			//��¼�����ľ�������Ӧ�ķ�ռ�Ž��±� 
	int *flag;					//����ѱ�ɾ���Ľ�Ϊ0��δ��ɾ���Ľ�Ϊ1 
	dist = (double**)malloc(non_size*sizeof(double*));
	distIndex = (int**)malloc(non_size*sizeof(int*));
	flag = (int*)malloc(non_size*sizeof(int));

	double* distMatrix=generateDistMatrix(repertFit,non_size,non_indexes);

	for (int i = 0; i<non_size; ++i)
	{
		flag[i] = 1;
		dist[i] = (double*)malloc(non_size*sizeof(double));
		distIndex[i] = (int*)malloc(non_size*sizeof(int));
		for (int j = 0; j<non_size; ++j)
		{
			if (i == j)
				dist[i][j] = INF;
			else
				dist[i][j] = distMatrix[i*non_size+j];
			distIndex[i][j] = j;
		}
	}

	//�Ծ��������д��±������ 
	for (int i = 0; i<non_size; ++i)
		sort_dist_index(dist[i], distIndex[i], 0, non_size - 1);

	//�����޳�non_size-NA��ӵ����ռ�Ž� 
	int cur_size = non_size;
	while (cur_size>nArch-count)
	{
		//ɾ���ڽ�������С�Ĵ浵�������ͨ���Ըý�����Ѵ�������ʵ�� 
		int crowd_index;						//��¼�ڽ�������С�Ĵ浵���±� 
		for (int i = 0; i<non_size; ++i)
		{
			if (flag[i])
			{
				crowd_index = i;
				break;
			}
		}
		for (int i = crowd_index + 1; i < non_size; ++i)
		{
			if (flag[i] && cmp_crowd(dist[i], dist[crowd_index]))
				crowd_index = i;
		}
		flag[crowd_index] = 0;

		//�޸��ڽ���������ҵ�ɾ�����±꣬Ȼ���ڽ�����һ��ǰ�� 
		for (int i = 0; i<non_size; ++i)
		{
			if (flag[i] == 0) continue;
			for (int j = 0; j<cur_size; ++j)
			{
				if (distIndex[i][j] == crowd_index)
				{
					while (j<cur_size - 1)
					{
						dist[i][j] = dist[i][j + 1];
						distIndex[i][j] = distIndex[i][j + 1];
						j++;
					}
					dist[i][j] = INF;				//��ʱj=cur_size-1 
					distIndex[i][j] = crowd_index;
					break;
				}
			}
		}
		cur_size -= 1;
	}
	//�����޳����� 

	//���澭���޳�������Ĵ浵�⵽A�� 
	int p_arch=count;
	for (int i = 0; i < non_size; ++i)
	{
		if (flag[i])
		{
			copyToArchiveFromRepertory(p_arch,non_indexes[i]);
			p_arch++;
		}
	}
	for (int i = 0; i<non_size; ++i)
	{
		free(dist[i]);
		free(distIndex[i]);
	}
	free(dist);
	free(distIndex);
	free(flag);
	free(non_indexes);
	free(distMatrix);
}

double* generateDistMatrix(double* f,int non_size,int* non_indexes)
{
	double* distMatrix=(double*)malloc(non_size*non_size*sizeof(double));

	int i,j;
	for(i=0;i<non_size;i++)
	{
		distMatrix[i*non_size+i]=0.0;
		for(j=i+1;j<non_size;j++)
		{
			distMatrix[i*non_size+j]=
				distMatrix[j*non_size+i]=
				Euclid_Dist(&f[non_indexes[i]*nObj],&f[non_indexes[j]*nObj]);
		}
	}
	return distMatrix;
}

// �����һ����ŷ����þ��룬vec1 ������vec2�ľ��� 
double Euclid_Dist(double* vec1, double* vec2)
{
	double sum = 0;
	for(int n=0; n<nObj; n++)
	{
		if (vec1[n] < vec2[n])			//����SPEA2�е�SDE�������ŷ����þ���
			sum += pow(vec1[n] - vec2[n], 2) / pow(fun_max[n] - fun_min[n], 2);
	}
	return sqrt(sum);
}

// record index orser after sorting
void sort_dist_index(double *a, int *b, int left, int right)
{
	int pivot;
	int i, j;
	double temp;
	double temp_a;
	int temp_b;
	if(left<right)
	{
		temp=a[right];
		i=left-1;
		for(j=left;j<right;++j)
		{
			if(a[j]<=temp)
			{
				i+=1;
				temp_a=a[i];
				a[i]=a[j];
				a[j]=temp_a;
				temp_b=b[i];
				b[i]=b[j];
				b[j]=temp_b;
			}
		}
		pivot=i+1;
		temp_a=a[pivot];
		a[pivot]=a[right];
		a[right]=temp_a;
		temp_b=b[pivot];
		b[pivot]=b[right];
		b[right]=temp_b;
		sort_dist_index(a,b,left,pivot-1);
		sort_dist_index(a,b,pivot+1,right);
	}
}

//K���ڱȽϺ���
bool cmp_crowd(double *a, double *b)
{
	for(int i=0;i<10;++i)				//�������Ƚ�ǰNA���ڽ����� 
		if(a[i]!=b[i])	
			return a[i]<b[i];
	return  false;
}