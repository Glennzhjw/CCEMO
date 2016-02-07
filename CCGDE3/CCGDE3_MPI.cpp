#include "CCGDE3_MPI.h"

//	
double* x_variable;	//current population stored in one mpi
double* x_variable_sub;	//
double* x_fitness;
double* u_variable;	//for evaluation, join new sub and best subs
double* u_variable_sub;	//for evolution, new sub
double* u_fitness;	//new fitness
int*    x_rank;

double* collection_variable;	//for selection
double* collection_variable_sub;
double* collection_fitness;	//for selection
double* collection_dist;	//for selection
int     collection_size;

double* x_sub_all;	//for evolution
int*    x_rank_all;

int my_frontSize;
int* frontSize_all;	//rank==1

double* collection_nonDom_x;
double* collection_nonDom_fit;

double* finalSolutions;
double* finalFitness;
int* finalRank;

//	limits
double *minLimit;
double *maxLimit;

//
int* indexes;

//	problem name
char problemInstance[256];

int    	gen;				/* Counter for the number of generations 	*/
int   	Gmax;				/* Maximum number of generations 			*/
int     D;					/* Number of decision variables				*/
int 	DSp;				/* Number of decision variables	per species	*/
int 	Nobj;				/* Number of objetive functions				*/
int 	numSpecies;		
int 	NP;					/* Population size							*/
int 	Cycles;
int 	cycle;
int 	numRest;

//	siae of nondominated solution
int 	sizeSol;
int 	sizeSolND;
int		solutionCounter;

//	coefficient for DE
double 	F;			
double  CR;					

//	output file pointer
FILE *fptx;
FILE *fptu;
FILE *fptv;
FILE *fptm;
FILE *fptt;
FILE *fpts;

/* Problem definition */
void (*problem)(double *x, double *fx,int numVar);


//	MPI
int mpi_size;
int mpi_rank;
int* recv_size;
int* disp_size;



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/*	crossover	*/
void ED_rand_1_bin()
{
	int r1,r2,r3;
	int i,j,k;
	for(i = 0; i < NP; i++){		
		do{
			r1 = rnd(0,NP - 1);
		}while(r1 == i);
		do{
			r2 = rnd(0,NP - 1);
		}while(r2 == i || r2 == r1);
		do{
			r3 = rnd(0,NP - 1);
		}while(r3 == i || r3 == r1 || r3 == r2);

		j = rnd(0,DSp - 1);

		for(k = 0; k < DSp; k++){
			if(flip(CR) || k == j){
				u_variable_sub[i*DSp+k] = x_variable_sub[r1*DSp+k] + 
					F*(x_variable_sub[r2*DSp+k] - x_variable_sub[r3*DSp+k]);

				if(u_variable_sub[i*DSp+k]<minLimit[indexes[mpi_rank*DSp+k]]){
					u_variable_sub[i*DSp+k] = minLimit[indexes[mpi_rank*DSp+k]];
				}
				if(u_variable_sub[i*DSp+k]>maxLimit[indexes[mpi_rank*DSp+k]]){
					u_variable_sub[i*DSp+k] = maxLimit[indexes[mpi_rank*DSp+k]];
				}
			}else{
				u_variable_sub[i*DSp+k] = x_variable_sub[i*DSp+k];
			}
		}				
	}
	return;
}

/*	initialization	*/
void set_parameters(int _Cycles, int _Gmax, int _numSpecies, int _NP, 
	int _Nobj, int _D, double _F, double _CR, char *pro)
{
	Cycles=_Cycles;
	Gmax=_Gmax;
	numSpecies=_numSpecies;
	NP = _NP;
	Nobj = _Nobj;
	D = _D;
	DSp = (int)(D/numSpecies);
	F = _F;
	CR = _CR;
	strcpy(problemInstance,pro);

	EMO_TEST_SUITE::set_para(D,Nobj,2);
	setMPI();
	memoryAllocation();

	return;
}

void run(int iRun)
{
	double t_ini, t_fin;
	double secs;

	char evalsal[50];//={"results/Evaluations_"};
	char timesal[50];//={"results/Time_"};
	char sizesal[50];//={"results/Size_"};
	char funsal[50];//={"results/FUN_"};
	char varsal[50];//={"results/VAR_"};
	int evaluations;

	numRest = 0;
	cycle = 0;
	gen = 0;

	sprintf(evalsal,"results/Evaluations_%s_%d",problemInstance,D);
	sprintf(timesal,"results/Time_%s_%d",problemInstance,D);
	sprintf(sizesal,"results/Size_%s_%d",problemInstance,D);
	sprintf(funsal,"results/FUN_%s_%d_%d",problemInstance,D,iRun);
	sprintf(varsal,"results/VAR_%s_%d_%d",problemInstance,D,iRun);

	if(mpi_rank==0)
		t_ini=clock();

	randomize();    
	initializeProblem(problemInstance);    
	initializePopulation();
	synchronize_x_sub_all();
	evaluateSpeciesInitial();

	int i;
	while(cycle<Cycles){
		cycle++;
// 		if(mpi_rank==0)
// 			printf("Cycle:\t%d\n",cycle);

		gen = 0;
		while(gen<Gmax){
			gen++;
			ED_rand_1_bin();
			//species parallel, results differ from original
			evaluatePopulation();
			addSolution();
			selection();
			synchronize_x_sub_all();
		}
	}

	evaluations = cycle * Gmax * numSpecies * NP;

	generateSolutions();	

	if(mpi_rank==0)
	{
		fptx = fopen(funsal,"w");     

		exportNonDominatedPopulationObjetivesValues(finalFitness,sizeSol,fptx);    
		fptv = fopen(varsal,"w");
		exportNonDominatedPopulationSolutionValues(finalSolutions,sizeSol,fptv);

		fptu = fopen(evalsal,"a");
		fprintf(fptu, "%d\n",evaluations);

		fpts = fopen(sizesal,"a");
		fprintf(fpts, "%d\n",sizeSolND);

		if(mpi_rank==0)
			t_fin=clock();

		secs = (t_fin-t_ini)/CLOCKS_PER_SEC;

		fptt = fopen(timesal,"a");
		fprintf(fptt,"%.16g \n", secs);

		fclose(fptx);
		fclose(fptu);
		fclose(fptv);
		fclose(fptt);
		fclose(fpts);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return;
}

void selection()
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
	for (i=0; i<collection_size; i++){
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
				result = dominanceComparator(&(collection_fitness[(temp1->index)*Nobj]), 
					&(collection_fitness[(temp2->index)*Nobj]));
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
				my_frontSize = frontSize;
			}else{
				my_frontSize = NP;
			}
		}
		temp2 = elite->child;
		j=i;
		if((popSize + frontSize)<= NP){
			do{
				copyIndividual(temp2->index, i);
				x_rank[i] = rank;
				popSize+=1;
				temp2 = temp2->child;
				i+=1;
			}while(temp2 != NULL);
			assignCrowdingDistanceIndexes(j, i-1);
			rank+=1;
		}else{
			fillCrowdingDistance(i, frontSize, elite);
			popSize = NP;
			for (j=i; j<NP; j++){
				x_rank[j] = rank;
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

//for each mpi, compare x and u, add better ones
void addSolution()
{
	int i;
	int result;
	collection_size=0;
	for(i=0;i<NP;i++)
	{
		result=dominanceComparator(&x_fitness[i*Nobj],&u_fitness[i*Nobj]);
		if(result==1)
		{
			memcpy(&collection_variable[collection_size*D],&x_variable[i*D],D*sizeof(double));
			memcpy(&collection_fitness[collection_size*Nobj],&x_fitness[i*Nobj],Nobj*sizeof(double));
			memcpy(&collection_variable_sub[collection_size*DSp],&x_variable_sub[i*DSp],DSp*sizeof(double));
			collection_size++;
		}
		if(result==-1)
		{
			memcpy(&collection_variable[collection_size*D],&u_variable[i*D],D*sizeof(double));
			memcpy(&collection_fitness[collection_size*Nobj],&u_fitness[i*Nobj],Nobj*sizeof(double));
			memcpy(&collection_variable_sub[collection_size*DSp],&u_variable_sub[i*DSp],DSp*sizeof(double));
			collection_size++;
		}
		if(result==0)
		{
			memcpy(&collection_variable[collection_size*D],&x_variable[i*D],D*sizeof(double));
			memcpy(&collection_fitness[collection_size*Nobj],&x_fitness[i*Nobj],Nobj*sizeof(double));
			memcpy(&collection_variable_sub[collection_size*DSp],&x_variable_sub[i*DSp],DSp*sizeof(double));
			collection_size++;
			memcpy(&collection_variable[collection_size*D],&u_variable[i*D],D*sizeof(double));
			memcpy(&collection_fitness[collection_size*Nobj],&u_fitness[i*Nobj],Nobj*sizeof(double));
			memcpy(&collection_variable_sub[collection_size*DSp],&u_variable_sub[i*DSp],DSp*sizeof(double));
			collection_size++;
		}
	}
}

void evaluatePopulation()
{
	int i,j,k;
	for (i=0; i<NP; i++){
		joinRepresentativesBest(&u_variable[i*D],&u_variable_sub[i*DSp]);    
		EMO_TEST_SUITE::evaluate_problems(problemInstance,&u_variable[i*D],&u_fitness[i*Nobj], D,1,Nobj);   
	}
	return;
}

void joinRepresentativesBest(double *x,double *indX)
{
	int i,j,k;
	int var_address;
	/* Sort particles values in a vector of objetives in order to be evaluated */
	for(i=0;i<numSpecies;i++){
		if(i == mpi_rank){
			for(j=0;j<DSp;j++){
				x[indexes[i*DSp+j]] = indX[j];
			}
		}else{                        
			do{
				k = rnd(0,NP - 1);
			}while(x_rank_all[i*NP+k] != 1);
			var_address=i*NP*DSp+k*DSp;
			for(j=0;j<DSp;j++){
				x[indexes[i*DSp+j]] = 
					x_sub_all[var_address+j];
			}
		}
	}

	return;
}

void evaluateSpeciesInitial()
{
	evaluatePopulationRandom();
	nonDominatedSorting(x_fitness,x_rank,NP,NP);

	return;
}

void nonDominatedSorting(double* fitness_var, int* rank_var, int size, int sizeNP)
{
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
				result = dominanceComparator(&(fitness_var[(temp1->index)*Nobj]), &(fitness_var[(temp2->index)*Nobj]));
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
				rank_var[temp2->index] = rank;
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

int dominanceComparator(double *individual1, double *individual2){
	int i;

	/* Indicates if some objetive in solution 1 dominates the objetive in solution 2 */
	int dominates1 = 0;
	/* Indicates if some objetive in solution 2 dominates the objetive in solution 1 */
	int dominates2 = 0;
	int result;
	double value1, value2;    

	for(i = 0; i < Nobj; i++){
		value1 = individual1[i];
		value2 = individual2[i];
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

void evaluatePopulationRandom()
{
	int i,j,k;
	for (i=0; i<NP; i++){
		joinRepresentativesRandom(&x_variable[i*D],&x_variable_sub[i*DSp]);    
		EMO_TEST_SUITE::evaluate_problems(problemInstance,&x_variable[i*D],&x_fitness[i*Nobj], D,1,Nobj);   
	}
	return;
}

void joinRepresentativesRandom(double *x,double *indX)
{
	int i,j,k;
	int var_address;
	/* Sort particles values in a vector of objetives in order to be evaluated */
	for(i=0;i<numSpecies;i++){
		if(i == mpi_rank){
			for(j=0;j<DSp;j++){
				x[indexes[i*DSp+j]] = indX[j];
			}
		}else{
			k = rnd(0,NP - 1);
			var_address=i*NP*DSp+k*DSp;
			for(j=0;j<DSp;j++){
				x[indexes[i*DSp+j]] = 
					x_sub_all[var_address+j];
			}
		}
	}
	return;
}

void synchronize_x_sub_all()
{
	MPI_Allgather(x_variable_sub,NP*DSp,MPI_DOUBLE,
		x_sub_all,NP*DSp,MPI_DOUBLE,
		MPI_COMM_WORLD);
	MPI_Allgather(x_rank,NP,MPI_INT,
		x_rank_all,NP,MPI_INT,
		MPI_COMM_WORLD);

	return;
}

void initializePopulation()
{
	int i,j;
	/* Create indexes */
	for(i=0;i<D;i++){
		indexes[i]= D-1-i;
	}
	/* shuffle indexes */
	shuffle(indexes,D);
	MPI_Bcast(indexes,D,MPI_INT,0,MPI_COMM_WORLD);

	for(i=0;i<NP;i++)
	{
		for(j=0;j<DSp;j++)
		{
			x_variable_sub[i*DSp+j]=
				rndreal(minLimit[indexes[mpi_rank*DSp+j]],maxLimit[indexes[mpi_rank*DSp+j]]);
		}
	}

	return;
}

void initializeProblem(char *problemName)
{
	EMO_TEST_SUITE::setLimits(problemName, minLimit, maxLimit, D);
	return;
}

void copyIndividual(int source, int dest)
{
	//collection to x
	memcpy(&x_variable[dest*D],&collection_variable[source*D],D*sizeof(double));
	memcpy(&x_fitness[dest*Nobj],&collection_fitness[source*Nobj],Nobj*sizeof(double));
	memcpy(&x_variable_sub[dest*DSp],&collection_variable_sub[source*DSp],DSp*sizeof(double));
}

void setMPI()
{
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	recv_size=(int*)malloc(mpi_size*sizeof(int));
	disp_size=(int*)malloc(mpi_size*sizeof(int));

	return;
}

void update_recv_disp(int* num, int n)
{
	int size=n;
	int i;
	for(i=0;i<mpi_size;i++)
	{
		recv_size[i]=num[i]*size;
	}
	disp_size[0]=0;
	for(i=1;i<mpi_size;i++)
	{
		disp_size[i]=disp_size[i-1]+recv_size[i-1];
	}
}

void memoryAllocation()
{
	x_variable=(double*)malloc(NP*D*sizeof(double));
	x_variable_sub=(double*)malloc(NP*DSp*sizeof(double));
	x_fitness=(double*)malloc(NP*Nobj*sizeof(double));
	u_variable=(double*)malloc(NP*D*sizeof(double));
	u_variable_sub=(double*)malloc(NP*DSp*sizeof(double));
	u_fitness=(double*)malloc(NP*Nobj*sizeof(double));
	x_rank=(int*)malloc(NP*sizeof(int));

	collection_variable=(double*)malloc(2*NP*D*sizeof(double));
	collection_variable_sub=(double*)malloc(2*NP*DSp*sizeof(double));
	collection_fitness=(double*)malloc(2*NP*Nobj*sizeof(double));
	collection_dist=(double*)malloc(2*NP*sizeof(double));

	x_sub_all=(double*)malloc(numSpecies*NP*DSp*sizeof(double));
	x_rank_all=(int*)malloc(numSpecies*NP*sizeof(int));

	frontSize_all=(int*)malloc(numSpecies*sizeof(int));

	collection_nonDom_x=(double*)malloc(NP*D*sizeof(double));
	collection_nonDom_fit=(double*)malloc(NP*Nobj*sizeof(double));

	finalSolutions=(double*)malloc(numSpecies*NP*D*sizeof(double));
	finalFitness=(double*)malloc(numSpecies*NP*Nobj*sizeof(double));
	finalRank=(int*)malloc(numSpecies*NP*sizeof(int));

	minLimit=(double*)malloc(D*sizeof(double));
	maxLimit=(double*)malloc(D*sizeof(double));

	indexes=(int*)malloc(D*sizeof(int));

	return;
}

void freeMemory()
{
	free(x_variable);
	free(x_variable_sub);
	free(x_fitness);
	free(u_variable);
	free(u_variable_sub);
	free(u_fitness);
	free(x_rank);

	free(collection_variable);
	free(collection_variable_sub);
	free(collection_fitness);
	free(collection_dist);

	free(x_sub_all);
	free(x_rank_all);

	free(frontSize_all);

	free(collection_nonDom_x);
	free(collection_nonDom_fit);

	free(finalSolutions);
	free(finalFitness);
	free(finalRank);

	free(minLimit);
	free(maxLimit);

	free(indexes);

	free(recv_size);
	free(disp_size);

	return;
}



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void assignCrowdingDistanceList(list *lst, int frontSize)
{
	int **arrayFx;
	int *distance;
	int i, j;
	list *temp;
	temp = lst;
	if(frontSize==1){
		collection_dist[lst->index] = INF;
		return;
	}
	if(frontSize==2){
		collection_dist[lst->index] = INF;
		collection_dist[lst->child->index] = INF;
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
	assignCrowdingDistance(distance, arrayFx, frontSize);
	free(distance);
	for(i=0; i<Nobj; i++){
		free(arrayFx[i]);
	}
	free(arrayFx);
	return;
}

void assignCrowdingDistanceIndexes(int c1, int c2)
{
	int **arrayFx;
	int *distance;
	int i, j;
	int frontSize;
	frontSize = c2-c1+1;
	if(frontSize==1){
		collection_dist[c1] = INF;
		return;
	}
	if(frontSize==2){
		collection_dist[c1] = INF;
		collection_dist[c2] = INF;
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
	assignCrowdingDistance(distance, arrayFx, frontSize);
	free(distance);
	for(i=0; i<Nobj; i++){
		free (arrayFx[i]);
	}
	free(arrayFx);
	return;
}

void assignCrowdingDistance(int *distance, int **arrayFx, int frontSize)
{
	int i, j;
	for(i=0; i<Nobj; i++){
		for(j=0; j<frontSize; j++){
			arrayFx[i][j] = distance[j];
		}
		quickSortFrontObj(i, arrayFx[i], frontSize);
	}
	for(j=0; j<frontSize; j++){
		collection_dist[distance[j]] = 0.0;
	}
	for(i=0; i<Nobj; i++){
		collection_dist[arrayFx[i][0]] = INF;
	}
	for(i=0; i<Nobj; i++){
		for(j=1; j<frontSize-1; j++){
			if(collection_dist[arrayFx[i][j]] != INF){
				if(collection_fitness[arrayFx[i][frontSize-1]*Nobj+i] == collection_fitness[arrayFx[i][0]*Nobj+i]){
					collection_dist[arrayFx[i][j]] += 0.0;
				}else{
					collection_dist[arrayFx[i][j]] += 
						(collection_fitness[arrayFx[i][j+1]*Nobj+i] - 
						collection_fitness[arrayFx[i][j-1]*Nobj+i])/
						(collection_fitness[arrayFx[i][frontSize-1]*Nobj+i] - 
						collection_fitness[arrayFx[i][0]*Nobj+i]);
				}
			}
		}
	}

	for(j=0; j<frontSize; j++){
		if (collection_dist[distance[j]] != INF){
			collection_dist[distance[j]] = (collection_dist[distance[j]])/Nobj;
		}
	}
	return;
}

void quickSortFrontObj(int objcount, int arrayFx[], int sizeArrayFx){
	qSortFrontObj(objcount, arrayFx, 0, sizeArrayFx-1);
	return;
}

void qSortFrontObj(int objcount, int arrayFx[], int left, int right){
	int index;
	int temp;
	int i, j;
	double pivot;
	if (left<right){
		index = rnd(left, right);
		temp = arrayFx[right];
		arrayFx[right] = arrayFx[index];
		arrayFx[index] = temp;
		pivot = collection_fitness[arrayFx[right]*Nobj+objcount];
		i = left-1;
		for (j=left; j<right; j++){
			if (collection_fitness[arrayFx[j]*Nobj+objcount] <= pivot){
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
		qSortFrontObj(objcount, arrayFx, left, index-1);
		qSortFrontObj(objcount, arrayFx, index+1, right);
	}
	return;
}

void quickSortDistance(int *distance, int frontSize){
	qSortDistance(distance, 0, frontSize-1);
	return;
}

void qSortDistance(int *distance, int left, int right){
	int index;
	int temp;
	int i, j;
	double pivot;
	if (left<right){
		index = rnd(left, right);
		temp = distance[right];
		distance[right] = distance[index];
		distance[index] = temp;
		pivot = collection_dist[distance[right]];
		i = left-1;
		for (j=left; j<right; j++)
		{
			if (collection_dist[distance[j]] <= pivot)
			{
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
		qSortDistance(distance, left, index-1);
		qSortDistance(distance, index+1, right);
	}
	return;
}

int indexWorstCrowdingDistance(int size)
{
	int i,index;
	double aux;
	aux = collection_dist[0];
	index = 0;
	for(i=0;i<size;i++)
	{
		if(collection_dist[i] < aux)
		{
			aux = collection_dist[i];
			index = i;
		}
	}	
	return index;
}

void fillCrowdingDistance(int count, int frontSize, list *elite)
{
	int *distance;
	list *temp;
	int i, j;
	int missing = NP - count;
	while(missing<frontSize)
	{
		assignCrowdingDistanceList(elite->child, frontSize);
		distance = (int *)calloc(frontSize,sizeof(int));
		temp = elite->child;
		for(j=0; j<frontSize; j++)
		{
			distance[j] = temp->index;
			temp = temp->child;
		}
		quickSortDistance(distance, frontSize);
		deleteInd(elite,distance[0]);
		frontSize--;
		if(missing<frontSize)
		{
			free(distance);
		}
	}
	for(i=count, j=frontSize-1; i<NP; i++, j--)
	{
		copyIndividual(distance[j], i);
	}
	free(distance);
	return;
}

void generateSolutions()
{
	int i,j,k;
	sizeSol=0;
	solutionCounter = 0;
	MPI_Allgather(&my_frontSize,1,MPI_INT,
		frontSize_all,1,MPI_INT,
		MPI_COMM_WORLD);
	int count=0;
	for(i=0;i<NP;i++)
	{
		if(x_rank[i]==1)
		{
			memcpy(&collection_nonDom_x[count*D],
				&x_variable[i*D],D*sizeof(double));
			memcpy(&collection_nonDom_fit[count*Nobj],
				&x_fitness[i*Nobj],Nobj*sizeof(double));
			count++;
		}
	}

	update_recv_disp(frontSize_all,D);
	MPI_Gatherv(collection_nonDom_x,my_frontSize*D,MPI_DOUBLE,
		finalSolutions,recv_size,disp_size,MPI_DOUBLE,
		0,MPI_COMM_WORLD);
	update_recv_disp(frontSize_all,Nobj);
	MPI_Gatherv(collection_nonDom_fit,my_frontSize*Nobj,MPI_DOUBLE,
		finalFitness,recv_size,disp_size,MPI_DOUBLE,
		0,MPI_COMM_WORLD);

	for(i=0;i<numSpecies;i++)
	{
		sizeSol += frontSize_all[i];
	}

	if(mpi_rank==0)
		nonDominatedSorting(finalFitness,finalRank,sizeSol,sizeSol);
}

void exportNonDominatedPopulationObjetivesValues(double* fun, int size, FILE *fpt)
{
	int i,j;
	sizeSolND=0;
	for(i=0;i<size;i++)
	{
		if(finalRank[i]==1)
		{
			sizeSolND++;
			for(j=0;j<Nobj;j++)
			{
					fprintf(fpt,"%e\t",fun[i*Nobj+j]);
			}
		}
		fprintf(fpt,"\n");
	}
}

void exportNonDominatedPopulationSolutionValues(double* var, int size, FILE *fpt)
{
	int i,j;
	for(i=0;i<size;i++)
	{
		if(finalRank[i]==1)
		{
			for(j=0;j<D;j++)
			{
				fprintf(fpt,"%lf\t",var[i*D+j]);
			}
			fprintf(fpt,"\n");
		}
	}
}