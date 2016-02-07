# include "global.h"

void copyToArchiveFromRepertory(int iA, int iR)
{
	memcpy(&archive[iA*nDim],&repertory[iR*nDim],nDim*sizeof(double));
	memcpy(&archFit[iA*nObj],&repertFit[iR*nObj],nObj*sizeof(double));
// 	arch_F[iA]=repertoryF[iR];
// 	arch_CR[iA]=repertoryCR[iR];
}

int INDEX(int a, int b, int c)
{
	return ((a)*nPop*nDim+(b)*nDim+(c));
}

int INDEX_f(int a, int b, int c)
{
	return ((a)*nPop*nObj+(b)*nObj+(c));
}

void selectSamples(int pp,int candidate,int *r1,int *r2,int *r3)
{
	if (r1)
	{
		do
		{
			*r1 = rnd(0,pp-1);
		}
		while (*r1 == candidate);
	}

	if (r2)
	{
		do
		{
			*r2 = rnd(0,pp-1);
		}
		while ((*r2 == candidate) || (*r2 == *r1));
	}

	if (r3)
	{
		do
		{
			*r3 = rnd(0,pp-1);
		}
		while ((*r3 == candidate) || (*r3 == *r1)||(*r3 == *r2));
	}
}

void saveArchiveOld()
{
	memcpy(archiveOld,archive,cnArch*nDim*sizeof(double));
	memcpy(archFitOld,archFit,cnArch*nObj*sizeof(double));
	cnArchOld=cnArch;
}

void addArchiveOld2repertory()
{
	if(!cnArchOld)
		return;
	memcpy(&repertory[nRep*nDim],archive,cnArchOld*nDim*sizeof(double));
	memcpy(&repertFit[nRep*nObj],archFit,cnArchOld*nObj*sizeof(double));
// 	memcpy(&repertoryF[nRep],arch_F,cnArch*sizeof(double));
// 	memcpy(&repertoryCR[nRep],arch_CR,cnArch*sizeof(double));
	nRep+=cnArchOld;

	// eleminate duplicate ones
// 	int i=0;
// 	int j;
// 	while(i<nRep)
// 	{
// 		for(j=i+1;j<nRep;j++)
// 		{
// 			if(isDuplicate(repertory,i,j))
// 				break;
// 		}
// 		if(j<nRep-1)
// 		{
// 			memcpy(&repertory[i*nDim],&repertory[(nRep-1)*nDim],nDim*sizeof(double));
// 			memcpy(&repertFit[i*nObj],&repertFit[(nRep-1)*nObj],nObj*sizeof(double));
// 			nRep--;
// 		}
// 		else
// 		{
// 			i++;
// 		}
// 	}
}

bool isDuplicate(double* s, int i, int j)
{
	int k;
	for(k=0;k<nDim;k++)
	{
		if(s[i*nDim+k]!=s[j*nDim+k])
			return false;
	}
	return true;
}

void collect2repertory_select()
{
	memcpy(repertory,xCurrent,nSwm*nPop*nDim*sizeof(double));
	nRep=nSwm*nPop;
	memcpy(repertFit,xFitness,nRep*nObj*sizeof(double));

	memcpy(&repertory[nRep*nDim],archive,cnArch*nDim*sizeof(double));
	memcpy(&repertFit[nRep*nObj],archFit,cnArch*nObj*sizeof(double));
	nRep+=cnArch;
}

void get_nonDominateSize()
{
	nonDominateSize=cnArch;
	while(archiveClass[nonDominateSize-1]>1)
		nonDominateSize--;
}

void save_obj(FILE *fpt)
{
	int i,j;
	for(i=0;i<nonDominateSize;i++)
	{
		for(j=0;j<nObj;j++)
		{
			fprintf(fpt,"%e\t",archFit[INDEX_f(0,i,j)]);
		}
		fprintf(fpt,"\n");
	}
}

void save_var(FILE *fpt)
{
	int i,j;
	for(i=0;i<nonDominateSize;i++)
	{
		for(j=0;j<nDim;j++)
		{
			fprintf(fpt,"%e\t",archive[INDEX(0,i,j)]);
		}
		fprintf(fpt,"\n");
	}
}