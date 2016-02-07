#include "global.h"

void initializeProblem()
{
	EMO_TEST_SUITE::set_para(nDim,nObj,2);
	nSwm=nObj;
	ParameterSet();
	int i;
	F_mu=0.5;
	CR_mu=0.5;
	c_para=0.1;
	for(i=0;i<nPop;i++)
	{
		Sflag[i]=0;
	}

	for(i=0;i<nArch;i++)
	{
		arch_F[i]=0.5;
		arch_CR[i]=0.9;
	}
	t1=0.1;
	t2=0.1;
}

void initializePopulation()
{
	int i,j;

	for(i=0;i<nPop;i++)
	{
		for(j=0;j<nDim;j++)
		{
			xCurrent[i*nDim+j]=rndreal(minLimit[j], maxLimit[j]);
		}
	}
	for(i=0;i<nObj;i++)
	{
		xBFitness[i]=INF;
	}
}

void initializePopulation_sinusMap()
{
	int i,j;

	for(i=0;i<nPop;i++)
	{
		for(j=0;j<nDim;j++)
		{
			xCurrent[i*nDim+j]=minLimit[j]+sinusMap()*(maxLimit[j]-minLimit[j]);
		}
	}
	for(i=0;i<nObj;i++)
	{
		xBFitness[i]=INF;
	}
}

void initializePopulation_UD()
{
	//	UD for mixture experiments
	int i,j,k;
// 	int a;
	int s=nDim-1;
// 	double sum;
	int* tmp_d=(int*)calloc(nPop,sizeof(int));
	double* tmp_x=(double*)calloc(nDim,sizeof(double));
	for(j=0;j<nPop;j++)
		tmp_d[j]=j+1;
	for(i=0;i<nSwm;i++)
	{
		for(k=0;k<nDim;k++)
		{
			shuffle(tmp_d,nPop);
			shuffle(tmp_d,nPop);
			shuffle(tmp_d,nPop);
			for(j=0;j<nPop;j++)
			{
				xCurrent[INDEX(i,j,k)]=(2*tmp_d[j]-1.0)/(2.0*nPop);
			}
		}
// 		for(j=0;j<nPop;j++)
// 		{
// 			sum=0.0;
// 			for(k=0;k<nDim-1;k++)
// 			{
// 				sum+=xCurrent[INDEX(i,j,k)];
// 			}
// 			for(k=0;k<nDim;k++)
// 			{
// 				xCurrent[INDEX(i,j,k)]/=sum;
// 			}
// 		}
	}
	for(i=0;i<nSwm;i++)
	{
// 		for(j=0;j<nPop;j++)
// 		{
// 			//	transformation
// 			for(k=0;k<s;k++)
// 			{
// 				tmp_x[k]=1.0-pow(xCurrent[INDEX(i,j,k)],1.0/(s-1.0-j));
// 				for(a=0;a<k-1;a++)
// 				{
// 					tmp_x[k]*=pow(xCurrent[INDEX(i,j,a)],1.0/(s-1.0-a));
// 				}
// 			}
// 			tmp_x[s]=1.0;
// 			for(k=0;k<s;k++)
// 			{
// 				tmp_x[s]*=pow(xCurrent[INDEX(i,j,k)],1.0/(s-1.0-j));
// 			}
// 			for(k=0;k<nDim;k++)
// 			{
// 				xCurrent[INDEX(i,j,k)]=/*minLimit[k]+(maxLimit[k]-minLimit[k])*
// 					(2.0**/tmp_x[k]/*-1.0)/(2.0*nPop)*/;
// 			}
// 		}
		for(j=0;j<nObj;j++)
		{
			xBFitness[i*nObj+j]=INF;
		}
	}
	free(tmp_d);
	free(tmp_x);
}

void evaluatePopInitial()
{
	int i;
	for(i=0;i<nPop;i++)
		EMO_TEST_SUITE::evaluate_problems(testInstance,&xCurrent[i*nDim],
			&xFitness[i*nObj],nDim,1,nObj);

	MPI_Bcast(xCurrent,nPop*nDim,MPI_DOUBLE,0,my_comm);
	MPI_Bcast(xFitness,nPop*nObj,MPI_DOUBLE,0,my_comm);
}

void ParameterSet()
{
	EMO_TEST_SUITE::setLimits(testInstance, minLimit, maxLimit, nDim);
}