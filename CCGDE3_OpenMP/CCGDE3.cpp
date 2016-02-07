/*
  CCGDE3.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
#include "CCGDE3.h"
#include <fstream>

double seed=0.36;

int main(int argc, char **argv)
{
// 	omp_set_num_threads(20);

	int fun_start_num = 16; // the start number of test function, you can find the explanation in testInstance.txt
	int fun_end_num = 16;    // the end number of test function
	int num_run = 30;        //the number of running CCGDE3 independently

	char _problemName[256];        
	int _Nobj;//number of objective
	int _D;//dim of variable
	int _Cycles;
	int _Gmax;
	int _numSpecies;
	int _NP;
	int _DSp;
	double _F=0.5;
	double _CR=0.5;

	int seq;
	std::ifstream readf("testInstance.txt");
	int iPro;
	for(iPro=1;iPro<fun_start_num;iPro++)
	{
		readf>>seq;
		readf>>_problemName;
		readf>>_D;
		readf>>_Nobj;
	}
	for(iPro=fun_start_num;iPro<=fun_end_num;iPro++)
	{
		readf>>seq;
		readf>>_problemName;
		readf>>_D;
		readf>>_Nobj;

		_Gmax=1;
		_numSpecies=2;
		_NP=40;
		_DSp = (int)(_D/_numSpecies);
		_Cycles=(int)(_D*1e4/_Gmax/_numSpecies/_NP);

		{
			printf("\nCCGDE3\n");
			printf("\nCycles:%d ",_Cycles);
			printf("\nGens:%d ",_Gmax);
			printf("\nSpecies:%d ",_numSpecies);
			printf("\nNP:%d ",_NP);
			printf("\nProblem:%s ",_problemName);
			printf("\nNobj:%d ",_Nobj);
			printf("\nD:%d ",_D);   
			printf("\nF:%f ",_F);
			printf("\nCR:%f ",_CR);
			printf("\nSeed:%f \n",seed);
		}

		CCGDE3 CCGDE3_handle;
		CCGDE3_handle.set_parameters(_Cycles,_Gmax,_numSpecies,_NP,_Nobj,_D,_F,_CR,_problemName);
		for(int iRun=1;iRun<=num_run;iRun++)
		{
			CCGDE3_handle.seed=seed;
			printf("--  RUN %d  --\n",iRun);
			CCGDE3_handle.run(iRun);
			seed=CCGDE3_handle.randomperc();
		}
		CCGDE3_handle.freeMemory();
	}

	return 0;
}

