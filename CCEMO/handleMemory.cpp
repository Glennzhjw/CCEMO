# include "global.h"

void allocateMemory()
{
	minLimit=allocDouble(nDim);
	maxLimit=allocDouble(nDim);
	Indexes=allocInt(nDim);
	archIndex=allocInt(nDim);
	xCurrent=allocDouble(nPop*nDim);
	uTrail=allocDouble(nPop*nDim);
// 	cTrail=allocDouble(nSwm*nPop*nDim);
	xFitness=allocDouble(nPop*nObj);
// 	cFitness=allocDouble(nSwm*nPop*nObj);
	xBest=allocDouble(nDim);
	xBFitness=allocDouble(nObj);
	archive=allocDouble(nArch*nDim);
	archFit=allocDouble(nArch*nObj);
	archiveOld=allocDouble(nArch*nDim);
	archFitOld=allocDouble(nArch*nObj);
	archiveClass=allocInt(nArch);
	uArchive=allocDouble(nArch*nDim);
	repertory=allocDouble(MAX_SIZE*nDim);
	repertFit=allocDouble(MAX_SIZE*nObj);
	repertoryDensity=allocDouble(MAX_SIZE);
	repertoryF=allocDouble(MAX_SIZE);
	repertoryCR=allocDouble(MAX_SIZE);
	fun_max=allocDouble(nObj);
	fun_min=allocDouble(nObj);
	S_F=allocDouble(nPop);
	S_CR=allocDouble(nPop);
	Sflag=allocInt(nPop);
	arch_F=allocDouble(nArch);
	arch_CR=allocDouble(nArch);
	dimFlag=allocInt(nDim);
}

double* allocDouble(int size)
{
	double* tmp;
	if((tmp = (double *)calloc(size,sizeof(double))) == NULL){
		printf("ERROR!! --> calloc: no memory for vector\n");
		exit(1);
	}
	return tmp;
}

int* allocInt(int size)
{
	int* tmp;
	if((tmp = (int *)calloc(size,sizeof(int))) == NULL){
		printf("ERROR!! --> calloc: no memory for vector\n");
		exit(1);
	}
	return tmp;
}

void freeMemory()
{
	free(minLimit);
	free(maxLimit);
	free(Indexes);
	free(archIndex);
	free(xCurrent);
	free(uTrail);
// 	free(cTrail);
	free(xFitness);
// 	free(cFitness);
	free(xBest);
	free(xBFitness);
	free(archive);
	free(archFit);
	free(archiveOld);
	free(archFitOld);
	free(archiveClass);
	free(uArchive);
	free(repertory);
	free(repertFit);
	free(repertoryDensity);
	free(repertoryF);
	free(repertoryCR);
	free(fun_max);
	free(fun_min);
	free(S_F);
	free(S_CR);
	free(Sflag);
	free(arch_F);
	free(arch_CR);
	free(dimFlag);
	freeMPI();
#ifdef SUITE_PARALLEL
	EMO_TEST_SUITE::mpi_finalization();
#endif
}