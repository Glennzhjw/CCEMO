/*
  initialize.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
# include "CCGDE3.h"

CCGDE3::CCGDE3(){}
CCGDE3::~CCGDE3(){}

void CCGDE3::set_parameters(int _Cycles, int _Gmax, int _numSpecies, int _NP, 
							int _Nobj, int _D,
							double _F, double _CR, char *pro)
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
	memoryAllocation();
}

void CCGDE3::run(int iRun)
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

	t_ini=clock();

	randomize();    
	initializeProblem(problemInstance);    
	createSpecies();
	initializeSpecies();
	evaluateSpeciesInitial();

	int i;
	while(cycle<Cycles){
		cycle++;

		for(i=0;i<numSpecies;i++){
			gen = 0;
			while(gen<Gmax){
				gen++;
				ED_rand_1_bin(&species[i]);
				evaluatePopulation(species[i].uG,NP,i);
				selection(&species[i]);
				nextGeneration(&species[i]);
			}
		}               
	}

	evaluations = cycle * Gmax * numSpecies * NP;

	{
		fptx = fopen(funsal,"w");     

		generateSolutions();	
		exportNonDominatedPopulationObjetivesValues(Solution,sizeSol,fptx);    
		fptv = fopen(varsal,"w");
		exportNonDominatedPopulationSolutionValues(Solution,sizeSol,fptv);

		fptu = fopen(evalsal,"a");
		fprintf(fptu, "%d\n",evaluations);

		fpts = fopen(sizesal,"a");
		fprintf(fpts, "%d\n",sizeSolND);

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
}

void CCGDE3::initializeSpecies(){
    int i;
    for (i=0; i<numSpecies; i++){
        initializePopulation(species[i].xG,NP,species[i].indexes);
    }
    return;
}

void CCGDE3::initializePopulation(Population *population, int size, int *indexes){
    int i;
    for (i=0; i<size; i++){
        initializeIndividual(&(population->individual[i]),indexes);
    }
    return;
}

void CCGDE3::initializeIndividual(Individual *individual,int *indexes){
    int j;
    
	for (j=0; j<DSp; j++){
		individual->x[j] = rndreal(minLimit[indexes[j]], maxLimit[indexes[j]]);
	}
    
    return;
}

void CCGDE3::initializeProblem(char *problemName){
	EMO_TEST_SUITE::setLimits(problemName, minLimit, maxLimit, D);
// 	if(strcmp(problemName, "ZDT1") == 0){
//         limitsZDT1(minLimit,maxLimit,D);
//         problem = ZDT1;
//         strcpy(hvexe, "./hv -r \"1.1 1.1\" ");
//         hvpercent = 0.875646 * percent;
//     }else if(strcmp(problemName, "ZDT2") == 0){
//         limitsZDT2(minLimit,maxLimit,D);
//         problem = ZDT2;
//         strcpy(hvexe, "./hv -r \"1.1 1.1\" ");
//         hvpercent = 0.542332 * percent;
//     }else if(strcmp(problemName, "ZDT3") == 0){
//         limitsZDT3(minLimit,maxLimit,D);
//         problem = ZDT3;
//         strcpy(hvexe, "./hv -r \"0.9 1.1\" ");
//         hvpercent = 0.955208 * percent;
//     }else if(strcmp(problemName, "ZDT4") == 0){
//         limitsZDT4(minLimit,maxLimit,D);
//         problem = ZDT4;
//         strcpy(hvexe, "./hv -r \"1.1 1.1\" ");
//         hvpercent = 0.873844 * percent;
//     }else if(strcmp(problemName, "ZDT6") == 0){
//         limitsZDT6(minLimit,maxLimit,D);
//         problem = ZDT6;
//         strcpy(hvexe, "./hv -r \"1.1 1.1\" ");
//         hvpercent = 0.506946 * percent;
//     }
}