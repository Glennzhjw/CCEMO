# ifndef _CCGDE3_
# define _CCGDE3_
/*
  CCGDE3.h

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/



# include "EMO_test_suite.h"
# include <cstdio>
# include <cmath>
# include <limits.h>
# include <unistd.h>
# include <ctime>
# include <iostream>
# include <string.h>
# include <cstdlib>
# include "ZDT.h"
# include <time.h>
#include <omp.h>

# define INF 1.0e14
using namespace std;

/* Individual definition */
typedef struct{
    int    rank;
    int    generation;		
    double *x;				
	double *fx;	
	double *cx;			
	double fitness;			
	double restViol;
    double *restrictions;
    double crowdDist;
}Individual;

/* Population definition */
typedef struct{
    Individual *individual;
}Population;

/* Species definition */
typedef struct{
	int    *indexes;		
	int 	frontSize;
	Population *xG;			/* Actual population	 					*/
	Population *uG;			/* Trial population							*/
	Population *xGm;			
}Species;

//	list
typedef struct node{
	int index;
	struct node *parent;
	struct node *child;
}list,node;



class CCGDE3
{
public:
	CCGDE3();
	~CCGDE3();

public:
	char hvexe[100];
	double hvpercent;
	double percent;

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
	int 	counter;
	int 	numRest;

	int 	sizeSol;
	int 	sizeSolND;
	int		solutionCounter;

	double 	F;			
	double  CR;					

	double *minLimit;
	double *maxLimit;

	Species *species;
	Population *Solution;


	FILE *fptx;
	FILE *fptu;
	FILE *fptv;
	FILE *fptm;
	FILE *fptt;
	FILE *fpts;

	/* Problem definition */
	void (*problem)(double *x, double *fx,int numVar);

	//	rand
	double seed;
	double oldrand[55];
	int jrand;

public:
	/* crossover */
	void ED_rand_1_bin(Species *species);

	/* crowdingDistance */
	void assignCrowdingDistanceList(Population *population, list *lst, int frontSize);
	void assignCrowdingDistanceIndexes(Population *population, int c1, int c2);
	void assignCrowdingDistance(Population *population, int *distance, int **arrayFx, int frontSize);

	/* display */
	void displayPopulation(Population *population, int size);
	void displayIndividual(Individual *individual);
	void exportPopulation(Population *population, int size, int gener, FILE *fpt);
	void exportPopulationObjetivesValues(Population *population, int size, FILE *fpt);
	void exportPopulationSolutionValues(Population *population, int size, FILE *fpt);
	void exportNonDominatedPopulationObjetivesValues(Population *population, int size, FILE *fpt);
	void exportNonDominatedPopulationSolutionValues(Population *population, int size, FILE *fpt);

	/* dominanceComparator */
	int dominanceComparator(Individual *individual1, Individual *individual2);

	/* initialize */
	void set_parameters(int _Cycles, int _Gmax, int _numSpecies, int _NP, 
		int _Nobj, int _D, double _F, double _CR, char *pro);
	void run(int iRun);
	void initializeSpecies();
	void initializePopulation(Population *population, int size, int *indexes);
	void initializeIndividual(Individual *individual,int *indexes);
	void initializeProblem(char *problemName);
#ifdef MPI_AVAILABLE
	void setMPI();
#endif

	/* evaluate */
	void evaluateSpeciesInitial();
	void evaluatePopulationRandom(Population *population, int indP);
	void evaluateIndividualRandom(Individual *individual,int indP);
	void evaluatePopulation(Population *population, int size, int indP);
	void evaluateIndividual(Individual *individual, int indP);

	/* memoryAllocation */
	void memoryAllocation();
	void populationMemoryAllocation(Population *population, int size);
	void individualMemoryAllocation(Individual *individual);
	void freeMemory();
	void freePopulationMemory(Population *population, int size);
	void freeIndividualMemory(Individual *individual);

	/* selection */

	void addSolution(Population *xGm, Individual *individual);
	void selection(Species *species);
	void nextGeneration(Species *species);
	void nonDominatedSorting(Population *xG, int size, int sizeNP);
	void fillCrowdingDistance(Population *xGm, Population *xG, int count, int frontSize, list *elite);

	/* sort */
	void quickSortFrontObj(Population *population, int objcount, int arrayFx[], int sizeArrayFx);
	void qSortFrontObj(Population *population, int objcount, int arrayFx[], int left, int right);
	void quickSortDistance(Population *population, int *distance, int frontSize);
	void qSortDistance(Population *population, int *distance, int left, int right);
	int indexWorstCrowdingDistance(Population *population,int size);

	/* species */
	void createSpecies();
	void joinRepresentativesRandom(double *x,double *indX,int indP);
	void joinRepresentativesBest(double *x,double *indX,int indP);
	void generateSolutions();

	/* util */
	void copyIndividual(Individual *original, Individual *copy);
	double timeval_diff(struct timeval *a, struct timeval *b);
	//	list
	void insert(node *n, int x);
	node* deleteNode(node *n);
	node* deleteInd(list *l,int index);
	void deleteList(list *l);
	list* createList(int index);
	//	rand
	void shuffle(int *x,int size);
	int flip(float prob);
	void randomize(void);
	void warmup_random (double seed);
	void advance_random (void);
	double randomperc(void);
	int rnd (int low, int high);
	double rndreal (double low, double high);
};

# endif
