#include "global.h"
#include "time.h"

void set_parameter(int npop, int ndim, int nobj, int narch)
{
	F=0.5;
	CR=0.5;

	nObj=nobj;
	nDim=ndim;
	nPop=npop;
	nArch=narch;
	nSwm=nObj;
	maxIteration=(int)(nDim*1e4);
	//	30
	S_SIZE=3;
	{
		S[0]=2;
		S[1]=3;
		S[2]=6;
		S[3]=10;
		S[4]=15;
		S[5]=30;
	}
	//	50
// 	S_SIZE=3;
// 	{
// 		S[0]=2;
// 		S[1]=5;
// 		S[2]=10;
// 		S[3]=25;
// 		S[4]=50;
// 		S[5]=30;
// 	}
//	archive
	arch_S_SIZE=2;
	{
		arch_S[0]=2;
		arch_S[1]=5;
		arch_S[2]=10;
	}
	allocateMemory();

	setMPI();
#ifdef SUITE_PARALLEL
	EMO_TEST_SUITE::mpi_initialization(mpi_rank,mpi_size,nPop,nObj);
#endif
}

void run(char* func_name, int iRun)
{
	iter=0;
	iter_each=0;
	cnArch=0;

	char objsal[256];
	char varsal[256];
	char timesal[256];
	strcpy(testInstance,func_name);

	double start_time,end_time;
	double duration;

	sprintf(objsal,"results/FUN_%s_%d_%d",testInstance,nDim,iRun);
	sprintf(varsal,"results/VAR_%s_%d_%d",testInstance,nDim,iRun);
	sprintf(timesal,"results/TIME_%s_%d",testInstance,nDim);

	if(mpi_rank==0)
		start_time=clock();

	initializeProblem();
	initializePopulation();
// 	showPopulation();
	initializeIndex();
	if(mpi_color)
		evaluatePopInitial();
	MPI_Barrier(MPI_COMM_WORLD);
	if(master_flag)
		collection_initial();
	if(!mpi_color&&mpi_rank_master==mpi_rank_master_archive)
		refineRepertory_generateArchive();
	MPI_Barrier(MPI_COMM_WORLD);
	if(master_flag)
		SynchronizeArchive();
	SynchronizeArchive_one();
// 	refineRepertory_generateArchive_SDE();
// 	showArchive();
// 	showLimits();
	if(mpi_color)
	{
		update_xBest_initial();
		update_xBest_archive();
	}
	MPI_Barrier(MPI_COMM_WORLD);
	iter=nPop*nObj;

// 	for(int i=0;i<mpi_size;i++)
// 	{
// 		if(i==0&&mpi_rank==0)
// 			showArchive();
// 		if(i!=0&&mpi_rank==i)
// 		{
// 			showGlobalBest();
// 			showPopulation();
// 		}
// 		MPI_Barrier(MPI_COMM_WORLD);
// 	}

	while(iter<maxIteration)
	{
		iter_each=0;
		nRep=0;
		if(mpi_color)
		{
			update_objective();
			permIndexes();
		}
		else
		{
			update_archive();
			perm_archIndex();
		}
		MPI_Barrier(MPI_COMM_WORLD);
// 		showGlobalBest();
		if(master_flag)
		{
			collect2master_archive();
			if(mpi_rank_master==mpi_rank_master_archive)
				refineRepertory_generateArchive();
			SynchronizeArchive();
		}
		SynchronizeArchive_one();
		if(mpi_color)
			update_xBest_archive();
// 		showArchive();
		MPI_Barrier(MPI_COMM_WORLD);
		update_iteration();
	}

	if(mpi_rank==0)
	{
		end_time=clock();
		duration=(end_time-start_time)/CLOCKS_PER_SEC;
	}

	if(mpi_rank==0)
	{
		fpttime=fopen(timesal,"a");
		fprintf(fpttime, "%e\n", duration);

		get_nonDominateSize();

		fptobj=fopen(objsal,"w");
		save_obj(fptobj);

		fptvar=fopen(varsal,"w");
		save_var(fptvar);

		fclose(fpttime);
		fclose(fptobj);
		fclose(fptvar);
	}
}

void update_objective()
{
	cooperativeCoevolution();
	refineRepertory_generateArchive();
	SynchronizeBest();
	collect2repertory0();
	if(master_flag)
		refineRepertory_generateArchive();
}

void update_archive()
{
	cooperativeCoevolution_archive();
	saveArchiveOld();
	refineRepertory_generateArchive();
	collect2repertory0();
	if(master_flag)
	{
		addArchiveOld2repertory();
		refineRepertory_generateArchive();
// 		refineRepertory_generateArchive_SDE();
	}
}
