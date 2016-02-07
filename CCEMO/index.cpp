# include "global.h"

void initializeIndex()
{
	int j;

	if(mpi_color)
	{
		for(j=0;j<nDim;j++)
		{
			Indexes[j]=nDim-1-j;
		}
		shuffle(Indexes,nDim);
		dimIn1Group=S[rnd(0,S_SIZE-1)];
		nGroup=nDim/dimIn1Group;

		MPI_Bcast(Indexes,nDim,MPI_INT,0,my_comm);
		MPI_Bcast(&dimIn1Group,1,MPI_INT,0,my_comm);
		MPI_Bcast(&nGroup,1,MPI_INT,0,my_comm);
	}
	else
	{
		for(j=0;j<nDim;j++)
		{
			archIndex[j]=nDim-1-j;
		}
		shuffle(archIndex,nDim);
		arch_dimIn1Group=arch_S[rnd(0,arch_S_SIZE-1)];
		arch_nGroup=nDim/arch_dimIn1Group;

		MPI_Bcast(archIndex,nDim,MPI_INT,0,my_comm);
		MPI_Bcast(&arch_dimIn1Group,1,MPI_INT,0,my_comm);
		MPI_Bcast(&arch_nGroup,1,MPI_INT,0,my_comm);
	}
}

void permIndexes()
{
	shuffle(Indexes,nDim);
	dimIn1Group=S[rnd(0,S_SIZE-1)];
	nGroup=nDim/dimIn1Group;

	MPI_Bcast(Indexes,nDim,MPI_INT,0,my_comm);
	MPI_Bcast(&dimIn1Group,1,MPI_INT,0,my_comm);
	MPI_Bcast(&nGroup,1,MPI_INT,0,my_comm);
// 	for(i=0;i<nSwm;i++)
// 	{
// 		for(j=0;j<nDim;j++)
// 		{
// 			printf("%d\t",Indexes[i*nDim+j]);
// 		}
// 	}
}

void perm_archIndex()
{
	shuffle(archIndex,nDim);
	arch_dimIn1Group=arch_S[rnd(0,arch_S_SIZE-1)];
	arch_nGroup=nDim/arch_dimIn1Group;

	MPI_Bcast(archIndex,nDim,MPI_INT,0,my_comm);
	MPI_Bcast(&arch_dimIn1Group,1,MPI_INT,0,my_comm);
	MPI_Bcast(&arch_nGroup,1,MPI_INT,0,my_comm);
}