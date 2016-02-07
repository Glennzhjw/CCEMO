#ifndef _CCGDE3_MPI_H_
#define _CCGDE3_MPI_H_

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
# include "list.h"
# include "rand.h"

# define INF 1.0e14

/* crossover */
void ED_rand_1_bin();

/* crowdingDistance */
void assignCrowdingDistanceList(list *lst, int frontSize);
void assignCrowdingDistanceIndexes(int c1, int c2);
void assignCrowdingDistance(int *distance, int **arrayFx, int frontSize);

/* display */
void exportNonDominatedPopulationObjetivesValues(double* fun, int size, FILE *fpt);
void exportNonDominatedPopulationSolutionValues(double* var, int size, FILE *fpt);

/* dominanceComparator */
int dominanceComparator(double *individual1, double *individual2);

/* initialize */
void set_parameters(int _Cycles, int _Gmax, int _numSpecies, int _NP, 
	int _Nobj, int _D, double _F, double _CR, char *pro);
void run(int iRun);
void initializePopulation();
void initializeProblem(char *problemName);
void setMPI();
void update_recv_disp(int* num, int n);

/* evaluate */
void evaluateSpeciesInitial();
void synchronize_x_sub_all();
void evaluatePopulationRandom();
void evaluatePopulation();

/* memoryAllocation */
void memoryAllocation();
void freeMemory();

/* selection */
void addSolution();
void selection();
void nonDominatedSorting(double* fitness_var, int* rank_var, int size, int sizeNP);
void fillCrowdingDistance(int count, int frontSize, list *elite);

/* sort */
void quickSortFrontObj(int objcount, int arrayFx[], int sizeArrayFx);
void qSortFrontObj(int objcount, int arrayFx[], int left, int right);
void quickSortDistance(int *distance, int frontSize);
void qSortDistance(int *distance, int left, int right);
int indexWorstCrowdingDistance(int size);

/* species */
void joinRepresentativesRandom(double *x,double *indX);
void joinRepresentativesBest(double *x,double *indX);
void generateSolutions();

/* util */
void copyIndividual(int source, int dest);
double timeval_diff(struct timeval *a, struct timeval *b);


#endif