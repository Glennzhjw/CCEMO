#include "global.h"

void setMPI()
{
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Get_processor_name(my_name,&name_len);

	recv_size=(int*)malloc(MAX_SIZE*sizeof(int));
	disp_size=(int*)malloc(MAX_SIZE*sizeof(int));
	each_size=(int*)malloc(MAX_SIZE*sizeof(int));

	int i;
	for(i=0;i<MAX_SIZE;i++)
		recv_size[i]=disp_size[i]=each_size[i]=0;

	allocateMPI();
	if(mpi_color) allocateTask(nPop);//fixed number
}

void update_recv_disp(int* num, int n, int l)
{
	int size=n;
	int i;
	int len=l;
	for(i=0;i<len;i++)
	{
		recv_size[i]=num[i]*size;
	}
	disp_size[0]=0;
	for(i=1;i<len;i++)
	{
		disp_size[i]=disp_size[i-1]+recv_size[i-1];
	}
}

void freeMPI()
{
	free(recv_size);
	free(disp_size);
	free(each_size);

	MPI_Comm_free(&my_comm);
	MPI_Comm_free(&comm_masters);
}

void collection_initial()
{
	int i;
	int my_send;
	for(i=0;i<=nObj;i++) each_size[i]=nPop;
	each_size[mpi_rank_master_archive]=0;

	update_recv_disp(each_size,nDim,nObj+1);
	my_send=each_size[mpi_rank_master]*nDim;
	MPI_Gatherv(xCurrent,my_send,MPI_DOUBLE,
		repertory,recv_size,disp_size,MPI_DOUBLE,
		mpi_rank_master_archive,comm_masters);
	update_recv_disp(each_size,nObj,nObj+1);
	my_send=each_size[mpi_rank_master]*nObj;
	MPI_Gatherv(xFitness,my_send,MPI_DOUBLE,
		repertFit,recv_size,disp_size,MPI_DOUBLE,
		mpi_rank_master_archive,comm_masters);
	if(mpi_rank_master==mpi_rank_master_archive)
		nRep=nPop*nObj;
}

void SynchronizeBest()
{
	MPI_Allgather(&task_num,1,MPI_INT,
		each_size,1,MPI_INT,my_comm);
	update_recv_disp(each_size,nDim,mpi_size_self);
	MPI_Gatherv(&xCurrent[task_l*nDim],task_num*nDim,MPI_DOUBLE,
		repertory,recv_size,disp_size,MPI_DOUBLE,
		0,my_comm);
	update_recv_disp(each_size,nObj,mpi_size_self);
	MPI_Gatherv(&xFitness[task_l*nObj],task_num*nObj,MPI_DOUBLE,
		repertFit,recv_size,disp_size,MPI_DOUBLE,
		0,my_comm);

	memcpy(xCurrent,repertory,nPop*nDim*sizeof(double));
	memcpy(xFitness,repertFit,nPop*nObj*sizeof(double));

	MPI_Bcast(xCurrent,nPop*nDim,MPI_DOUBLE,0,my_comm);
	MPI_Bcast(xFitness,nPop*nObj,MPI_DOUBLE,0,my_comm);

	double* best_fit=(double*)malloc(mpi_size_self*nObj*sizeof(double));
	int best_ind=-1;

	MPI_Gather(xBFitness,nObj,MPI_DOUBLE,
		best_fit,nObj,MPI_DOUBLE,
		0,my_comm);

	int i;
	if(mpi_color&&master_flag)
	{
		best_ind=0;
		for(i=1;i<mpi_size_self;i++)
		{
			if(check_dominance(&best_fit[best_ind*nObj],&best_fit[i*nObj])!=1)
				best_ind=i;
		}
	}

	MPI_Bcast(&best_ind,1,MPI_INT,0,my_comm);
	MPI_Bcast(xBest,nDim,MPI_DOUBLE,best_ind,my_comm);
	MPI_Bcast(xBFitness,nObj,MPI_DOUBLE,best_ind,my_comm);

	free(best_fit);
}

void collect2repertory0()
{
	MPI_Allgather(&cnArch,1,MPI_INT,
		each_size,1,MPI_INT,my_comm);
	update_recv_disp(each_size,nDim,mpi_size_self);
	MPI_Gatherv(archive,cnArch*nDim,MPI_DOUBLE,
		repertory,recv_size,disp_size,MPI_DOUBLE,
		0,my_comm);
	update_recv_disp(each_size,nObj,mpi_size_self);
	MPI_Gatherv(archFit,cnArch*nObj,MPI_DOUBLE,
		repertFit,recv_size,disp_size,MPI_DOUBLE,
		0,my_comm);
	MPI_Reduce(&cnArch,&nRep,1,MPI_INT,MPI_SUM,0,my_comm);
}

//for all master nodes of my_comm
void SynchronizeArchive()
{
	MPI_Bcast(&cnArch,1,MPI_INT,mpi_rank_master_archive,comm_masters);
	MPI_Bcast(archive,cnArch*nDim,MPI_DOUBLE,mpi_rank_master_archive,comm_masters);
	MPI_Bcast(archFit,cnArch*nObj,MPI_DOUBLE,mpi_rank_master_archive,comm_masters);
	MPI_Bcast(archiveClass,cnArch,MPI_INT,mpi_rank_master_archive,comm_masters);
}

//for each my_comm
void SynchronizeArchive_one()
{
	MPI_Bcast(&cnArch,1,MPI_INT,0,my_comm);
	MPI_Bcast(archive,cnArch*nDim,MPI_DOUBLE,0,my_comm);
	MPI_Bcast(archFit,cnArch*nObj,MPI_DOUBLE,0,my_comm);
	MPI_Bcast(archiveClass,cnArch,MPI_INT,0,my_comm);
}

void collect2master_archive()
{
	int my_send;

	MPI_Allgather(&cnArch,1,MPI_INT,
		each_size,1,MPI_INT,comm_masters);

	update_recv_disp(each_size,nDim,mpi_size_master);
	my_send=each_size[mpi_rank_master]*nDim;
	MPI_Gatherv(archive,my_send,MPI_DOUBLE,
		repertory,recv_size,disp_size,MPI_DOUBLE,
		mpi_rank_master_archive,comm_masters);
	update_recv_disp(each_size,nObj,mpi_size_master);
	my_send=each_size[mpi_rank_master]*nObj;
	MPI_Gatherv(archFit,my_send,MPI_DOUBLE,
		repertFit,recv_size,disp_size,MPI_DOUBLE,
		mpi_rank_master_archive,comm_masters);
	MPI_Reduce(&cnArch,&nRep,1,MPI_INT,MPI_SUM,
		mpi_rank_master_archive,comm_masters);
}

void update_iteration()
{
	MPI_Reduce(&iter_each,&iter_sum,1,MPI_INT,MPI_SUM,0,my_comm);
// 	printf("\nSUM=%d\n",iter_sum);
	if(master_flag)
	{
		iter_each=iter_sum;
		MPI_Allreduce(&iter_each,&iter_sum,1,MPI_INT,MPI_SUM,comm_masters);
		iter+=iter_sum;
	}
	MPI_Bcast(&iter,1,MPI_INT,0,my_comm);
// 	printf("\nITER=%d\n",iter);
}

/************************************************************************/
/* 
Allocate MPI processes for objectives and archive,
number of MPI processes:
	objective problems --- [1, nPop]
	archive            --- [1, nArch]

				archive   obj1   obj2   ...   objn
	mpi_color	   0        1      2            n
There are (nObj+1) communication regions (MPI_Comm)

*/
/************************************************************************/
void allocateMPI()
{
	int quo,rem;
	quo=mpi_size/(nObj+1);
	rem=mpi_size%(nObj+1);

	int i;
	if(quo<nPop)
	{
		for(i=0;i<=nObj;i++)
		{
			each_size[i]=quo;
			if(i<rem) each_size[i]++;
		}
	}
	else
	{
		for(i=1;i<=nObj;i++)
		{
			each_size[i]=nPop;
		}
		each_size[0]=mpi_size-nObj*nPop;
	}
	update_recv_disp(each_size,1,nObj+1);
	mpi_color=0;
	while(mpi_rank>=recv_size[mpi_color]+disp_size[mpi_color]) mpi_color++;
	//MPI processes are split to (nObj+1) MPI_Comm
	MPI_Comm_split(MPI_COMM_WORLD,mpi_color,mpi_rank,&my_comm);
	MPI_Comm_size(my_comm,&mpi_size_self);
	MPI_Comm_rank(my_comm,&mpi_rank_self);
	//select the master of each my_comm
	if(mpi_rank_self==0) master_flag=1;
	else master_flag=0;
	//MPI processes are split to 2 MPI_Comm,
	//if master_flag==1, comm_masters contains the masters of each my_comm
	//else, comm_masters contain the remain MPIs
	MPI_Comm_split(MPI_COMM_WORLD,master_flag,mpi_rank,&comm_masters);
	MPI_Comm_size(comm_masters,&mpi_size_master);
	MPI_Comm_rank(comm_masters,&mpi_rank_master);
	if(master_flag)
	{
		for(i=0;i<mpi_size_master;i++)
		{
			mpi_rank_master_archive=mpi_color;
			MPI_Bcast(&mpi_rank_master_archive,1,MPI_INT,i,comm_masters);
			MPI_Barrier(comm_masters);
			if(mpi_rank_master_archive==0)
			{
				mpi_rank_master_archive=i;
				break;
			}
		}
	}

	//	output mpi info
	for(i=0;i<mpi_size;i++)
	{
		if(mpi_rank==i)
		{
			printf("\n//----------------------------------------------\n");
			printf("My rank = %d in %d MPI processes, name: \t%s.\n",
				mpi_rank,mpi_size,my_name);
			printf("My color = %d, so I'm in ",mpi_color);
			if(!mpi_color)
			{
				printf("archive.\n");
				printf("And the size for archive is %d.\n",mpi_size_self);
			}
			else
			{
				printf("objective %d.\n",mpi_color);
				printf("And the size for the objective is %d.\n",mpi_size_self);
			}
			if(master_flag)
			{
				printf("I am a master.\n");
				printf("My rank in the master group is %d.\n",mpi_rank_master);
				printf("The size of the master group is %d.\n",mpi_size_master);
				if(mpi_rank_master==mpi_rank_master_archive)
					printf("And I am the master of the master group.\n");
			}
			else
			{
				printf("I am not a master.\n");
				printf("My rank in the non-master group is %d.\n",mpi_rank_master);
				printf("The size of the non-master group is %d.\n",mpi_size_master);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

void allocateTask(int n)
{
	int quo,rem;
	int i;
	int size=n;
	{
		quo=size/mpi_size_self;
		rem=size%mpi_size_self;
		for(i=0;i<mpi_size_self;i++)
		{
			each_size[i]=quo;
			if(i<rem) each_size[i]++;
		}
		update_recv_disp(each_size,1,mpi_size_self);
		task_l=disp_size[mpi_rank_self];
		task_r=task_l+recv_size[mpi_rank_self];
		task_num=recv_size[mpi_rank_self];
	}
}