# include "global.h"

void showRepertory()
{
	int i,j;
	printf("\n Repertory...\n");
	for(i=0;i<nRep;i++)
	{
		printf("ID:%d\n",i+1);
		for(j=0;j<nDim;j++)
		{
			printf("%lf\t",repertory[INDEX(0,i,j)]);
		}
		printf("\n");
		printf("Fitness:\n");
		for(j=0;j<nObj;j++)
		{
			printf("%.16e\t",repertFit[INDEX_f(0,i,j)]);
		}
		printf("\n");
	}
}

void showLimits()
{
	int i,j;
	printf("\n Limits...\n");
	for(i=0;i<nObj;i++)
	{
		printf("ID:%d\n",i+1);
		for(j=0;j<nDim;j++)
		{
			printf("%lf\t",minLimit[j]);
		}
		printf("\n");
		for(j=0;j<nDim;j++)
		{
			printf("%lf\t",maxLimit[j]);
		}
		printf("\n");
	}
}

void showArchive()
{
	int i,j;
	printf("\n Archive...%d\n",cnArch);
	for(i=0;i<cnArch;i++)
	{
		printf("ID:%d\n",i+1);
		for(j=0;j<nDim;j++)
		{
			printf("%lf\t",archive[INDEX(0,i,j)]);
		}
		printf("\n");
		printf("Fitness:\n");
		for(j=0;j<nObj;j++)
		{
			printf("%.16e\t",archFit[INDEX_f(0,i,j)]);
		}
		printf("\n");
		printf("rank:\t%d\n",archiveClass[i]);
	}
}

void showGlobalBest()
{
	int i,j;
	printf("xBest\t%d\n",mpi_rank);
	for(i=0;i<nDim;i++)
	{
		printf("%lf\t",xBest[i]);
	}
	printf("\n");
	for(j=0;j<nObj;j++)
	{
		printf("%.16e\t",xBFitness[j]);
	}
	printf("\n");
}

void showPopulation()
{
	int i,j;

	printf("Swarm:\t%d\n",mpi_rank);
	for(i=0;i<nPop;i++)
	{
		printf("ID:\t%d\n",i+1);
		for(j=0;j<nDim;j++)
		{
			printf("%lf\t",xCurrent[INDEX(0,i,j)]);
		}
		printf("\n");
		for(j=0;j<nObj;j++)
		{
			printf("%.16e\t",xFitness[INDEX_f(0,i,j)]);
		}
		printf("\n");
	}
}