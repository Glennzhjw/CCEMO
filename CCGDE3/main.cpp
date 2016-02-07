#include "CCGDE3_MPI.h"
#include <fstream>
#include <assert.h>

int main(int argc, char **argv)
{
#ifndef MPI_AVAILABLE
	printf("MPI is not used, exiting...\n");
	exit(1);
#endif
	MPI_Init(&argc,&argv);
	int my_rank,group_size;
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&group_size);

	int fun_start_num = 1; // the start number of test function, you can find the explanation in testInstance.txt
	int fun_end_num = 31;    // the end number of test function
	int num_run = 1;        //the number of running CCGDE3 independently

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

	seed=0.36;

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
		_numSpecies=group_size;
		assert(group_size<=_D);
		assert((_D%group_size==0));
		_NP=100;
		_DSp = (int)(_D/_numSpecies);
		_Cycles=(int)(_D*1e4/_Gmax/_numSpecies/_NP);

		if(my_rank==0)
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

		set_parameters(_Cycles,_Gmax,_numSpecies,_NP,_Nobj,_D,_F,_CR,_problemName);
		for(int iRun=1;iRun<=num_run;iRun++)
		{
			MPI_Bcast(&seed,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			if(my_rank==0)
				printf("--  RUN %d  --\n",iRun);
			run(iRun);
			seed=randomperc();
		}
		freeMemory();
	}

	MPI_Finalize();

	return 0;
}