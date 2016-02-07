#include "global.h"

int  seed = 36;		// seed for random number generation
int  NP=40;
int  N_arch=100;

#include <fstream>

int main(int argc, char **argv)
{
#ifndef MPI_AVAILABLE
	printf("\nMPI is not used, exiting...\n");
	exit(-1);
#endif

	MPI_Init(&argc,&argv);
	int my_rank;
	int group_size;
// 	char my_name[256];
// 	int name_len;
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&group_size);

	int fun_start_num = 1; // the start number of test function, you can find the explanation in testInstance.txt
	int fun_end_num = 31;    // the end number of test function
	int num_run = 30;        //the number of running test independently
	int iRun;

	int seq;
	char prob[256];
	int ndim;
	int nobj;

	std::ifstream readf("testInstance.txt");
	int iPro;
	for(iPro=1;iPro<fun_start_num;iPro++)
	{
		readf>>seq;
		readf>>prob;
		readf>>ndim;
		readf>>nobj;
	}
	for(iPro=fun_start_num;iPro<=fun_end_num;iPro++)
	{
		readf>>seq;
		readf>>prob;
		readf>>ndim;
		readf>>nobj;

		if((nobj+1)>group_size)
		{
			printf("\nERROR:\tThe number of mpi processes is \n\
				not enough, exiting...\n");
			exit(-2);
		}
		if(NP*nobj+N_arch<group_size)
		{
			printf("\nERROR:\tThe number of mpi processes is \n\
				too much, exiting...\n");
			exit(-3);
		}
		set_parameter(NP,ndim,nobj,N_arch);

		if(my_rank==0)
		{
			printf("\n\n-- PROBLEM %s\n",prob);
			printf("--  variables: %d\n",ndim);
			printf("--  objectives: %d\n\n",nobj);
		}

		for(iRun=1;iRun<=num_run;iRun++)
		{
			seed = (seed + 111 + my_rank * 10) % 1235;
			rnd_uni_init = -(long)seed;
			if(my_rank==0)
				printf("--   run %d   --\n",iRun);
			run(prob,iRun);
		}
		freeMemory();
	}

	MPI_Finalize();

	return 0;
}