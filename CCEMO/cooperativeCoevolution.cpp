# include "global.h"

void cooperativeCoevolution()
{
	DE_gen_uTrail();
	mainLoop();
// 	if(iter>=maxIteration)
// 		return;
}

void cooperativeCoevolution_ada()
{
	generate_F_CR();
	DE_gen_uTrail_ada();
	mainLoop2();
// 	if(iter>=maxIteration)
// 		return;
	update_F_CR_mu();
}

void cooperativeCoevolution_archive()
{
	DE_gen_archive();
	mainLoop_archive();
}

void update_xBest_initial()
{
	int i;
	int best_idx;
	double* best_fit;
	int probIdx;

	probIdx=mpi_rank-1;
	best_fit=(double*)malloc(nObj*sizeof(double));
	memcpy(best_fit,xBFitness,nObj*sizeof(double));

	best_idx=-1;
	for(i=0;i<nPop;i++)
	{
		if(xFitness[i*nObj+probIdx]<=best_fit[probIdx]
		 &&check_dominance(best_fit,&xFitness[i*nObj])!=1)
		{
			memcpy(best_fit,&xFitness[i*nObj],nObj*sizeof(double));
			best_idx=i;
		}
	}
	if(best_idx!=-1)
	{
		memcpy(xBest,&xCurrent[best_idx*nDim],nDim*sizeof(double));
		memcpy(xBFitness,&xFitness[best_idx*nObj],nObj*sizeof(double));
	}
	free(best_fit);
}

void update_xCurrent()
{
	int i;
	int best_idx;
	double* best_fit=(double*)calloc(nObj,sizeof(double));
	int probIdx;
	best_idx=-1;

	probIdx=mpi_rank-1;
	memcpy(best_fit,xBFitness,nObj*sizeof(double));
	for(i=0;i<nPop;i++)
	{
		if(cFitness[INDEX_f(0,i,probIdx)]<=xFitness[INDEX_f(0,i,probIdx)]
		 &&check_dominance(&xFitness[i*nObj],&cFitness[i*nObj])!=1)
		{
			Sflag[i]=1;
			memcpy(&xFitness[INDEX_f(0,i,0)],&cFitness[INDEX_f(0,i,0)],nObj*sizeof(double));
			memcpy(&xCurrent[INDEX(0,i,0)],&cTrail[INDEX(0,i,0)],nDim*sizeof(double));
		}
		else
			Sflag[i]=0;
		if(cFitness[INDEX_f(0,i,probIdx)]<=best_fit[probIdx]
		 &&check_dominance(best_fit,&cFitness[i*nObj])!=1)
		{
			memcpy(best_fit,&cFitness[i*nObj],nObj*sizeof(double));
			best_idx=i;
		}
	}
	if(best_idx!=-1)
	{
		memcpy(xBest,&cTrail[INDEX(0,best_idx,0)],nDim*sizeof(double));
		memcpy(xBFitness,&cFitness[INDEX_f(0,best_idx,0)],nObj*sizeof(double));
	}

	free(best_fit);
}

void update_xCurrent_one(int iPop)
{
	int i;
	int probIdx;
	double* best_fit=(double*)malloc(nObj*sizeof(double));
	int best_idx;
	double* best_fit_g=(double*)malloc(nObj*sizeof(double));
	int best_idx_g;
	probIdx=mpi_rank-1;

	memcpy(best_fit,&xFitness[iPop*nDim],nObj*sizeof(double));
	best_idx=-1;
	memcpy(best_fit_g,xBFitness,nObj*sizeof(double));
	best_idx_g=-1;

	for(i=0;i<nRep-nRep_tmp;i++)
	{
		if(cFitness[INDEX_f(0,i,probIdx)]<=best_fit[probIdx]
		&&check_dominance(best_fit,&cFitness[i*nObj])!=1)
		{
			memcpy(best_fit,&cFitness[INDEX_f(0,i,0)],nObj*sizeof(double));
			best_idx=i;
		}

		if(cFitness[INDEX_f(0,i,probIdx)]<=best_fit_g[probIdx]
		&&check_dominance(best_fit_g,&cFitness[i*nObj])!=1)
		{
			memcpy(best_fit_g,&cFitness[i*nObj],nObj*sizeof(double));
			best_idx_g=i;
		}
	}

	if(best_idx!=-1)
	{
		memcpy(&xCurrent[iPop*nDim],&cTrail[best_idx*nDim],nDim*sizeof(double));
		memcpy(&xFitness[iPop*nObj],&cFitness[best_idx*nObj],nObj*sizeof(double));
		Sflag[iPop]=1;
	}
	else
	{
		Sflag[iPop]=0;
	}

	if(best_idx_g!=-1)
	{
		memcpy(xBest,&cTrail[best_idx_g*nDim],nDim*sizeof(double));
		memcpy(xBFitness,&cFitness[best_idx_g*nObj],nObj*sizeof(double));
	}

	free(best_fit);
	free(best_fit_g);
}

void update_xBest_archive()
{
	int i;
	int best_idx;
	double* best_fit=(double*)malloc(nObj*sizeof(double));
	int probIdx;

	memcpy(best_fit,xBFitness,nObj*sizeof(double));
	probIdx=mpi_rank-1;

	best_idx=-1;
	for(i=0;i<cnArch;i++)
	{
		if(archFit[i*nObj+probIdx]<=best_fit[probIdx]
		 &&check_dominance(best_fit,&archFit[i*nObj])!=1)
		{
			memcpy(best_fit,&archFit[i*nObj],nObj*sizeof(double));
			best_idx=i;
		}
	}
	if(best_idx!=-1)
	{
		memcpy(xBest,&archive[best_idx*nDim],nDim*sizeof(double));
		memcpy(xBFitness,&archFit[best_idx*nObj],nObj*sizeof(double));
	}
	free(best_fit);
}

void mainLoop()
{
	int i,j;

	nRep=0;
	iter_each=0;
	for(i=task_l;i<task_r;i++)
	{
		nRep_tmp=nRep;
		cTrail=&repertory[nRep*nDim];
		cFitness=&repertFit[nRep*nObj];
// 		cF=&repertoryF[nRep];
// 		cCR=&repertoryCR[nRep];
// 		memcpy(cF,&S_F[i*nPop],nPop*sizeof(double));
// 		memcpy(cCR,&S_CR[i*nPop],nPop*sizeof(double));
		joinSolutions(i);
		for(j=0;j<nRep-nRep_tmp;j++)
			EMO_TEST_SUITE::evaluate_problems(testInstance,&cTrail[j*nDim],
				&cFitness[j*nObj],nDim,1,nObj);
		iter_each+=nRep-nRep_tmp;
// 		if(iter>=maxIteration)
// 			return;
		update_xCurrent_one(i);
	}
}

void joinSolutions(int iPop)
{
	int i;
	int j;
	int num;
	int rx;
	int copy_ind;
	int a;
	int a_idx;
	int a_low,a_high;

	for(i=0;i<cnArch;i++)
	{
		if(archiveClass[i]>1)
			break;
	}
	num=i;

	copy_ind=0;
	for(i=0;i<nGroup;i++)
	{
		rx=rnd(0,num-1);
		memcpy(&cTrail[copy_ind*nDim],&archive[rx*nDim],nDim*sizeof(double));
		a_low=i*dimIn1Group;
		a_high=(i+1)*dimIn1Group;
		for(a=a_low;a<a_high;a++)
		{
			a_idx=Indexes[a];
			cTrail[INDEX(0,copy_ind,a_idx)]=uTrail[INDEX(0,iPop,a_idx)];
		}

// 		for(j=0;j<nGroup;j++)
// 		{
// 			a_low=j*dimIn1Group;
// 			a_high=(j+1)*dimIn1Group;
// 			if(j==i)
// 			{
// 				for(a=a_low;a<a_high;a++)
// 				{
// 					a_idx=Indexes[a];
// 					cTrail[INDEX(0,copy_ind,a_idx)]=uTrail[INDEX(0,iPop,a_idx)];
// 				}
// 			}
// 			else
// 			{
// // 				selectSamples(nPop,iPop,&rx,NULL,NULL);
// 				for(a=a_low;a<a_high;a++)
// 				{
// 					a_idx=Indexes[a];
// 					cTrail[INDEX(0,copy_ind,a_idx)]=archive[INDEX(0,rx,a_idx)];
// 				}
// 			}
// 		}
		copy_ind++;
		nRep++;
	}
}

void mainLoop2()
{
	int i;
	nRep=0;
	joinSolutions2();
	nRep+=nPop;
	for(i=0;i<nRep;i++)
		EMO_TEST_SUITE::evaluate_problems(testInstance,&cTrail[i*nDim],
			&cFitness[i*nObj],nDim,1,nObj);
	iter_each+=nRep;
// 	if(iter>=maxIteration)
// 		return;
	update_xCurrent();
}

void joinSolutions2()
{
	int iPop;
	int iGrp;
	int g_ind;
	int r;
	int a;
	int a_idx;
	int a_low,a_high;

	cTrail=&repertory[nRep*nDim];
	cFitness=&repertFit[nRep*nObj];
// 	cF=&repertoryF[nRep];
// 	cCR=&repertoryCR[nRep];
	for(iPop=0;iPop<nPop;iPop++)
	{
		g_ind=rnd(0,nGroup-1);
		for(iGrp=0;iGrp<nGroup;iGrp++)
		{
			a_low=iGrp*dimIn1Group;
			a_high=(iGrp+1)*dimIn1Group;
			if(flip_r((float)0.5)&&iGrp!=g_ind)
			{
				r=rnd(0,cnArch-1);
				for(a=a_low;a<a_high;a++)
				{
					a_idx=Indexes[a];
					cTrail[INDEX(0,iPop,a_idx)]=archive[INDEX(0,r,a_idx)];
				}
			}
			else
			{
				for(a=a_low;a<a_high;a++)
				{
					a_idx=Indexes[a];
					cTrail[INDEX(0,iPop,a_idx)]=uTrail[INDEX(0,iPop,a_idx)];
				}
			}
		}
	}
// 	memcpy(cF,&S_F[iSwm*nPop],nPop*sizeof(double));
// 	memcpy(cCR,&S_CR[iSwm*nPop],nPop*sizeof(double));
}

void mainLoop_archive()
{
	int i,j;

	nRep=0;
	iter_each=0;
	for(i=task_l;i<task_r;i++)
	{
		nRep_tmp=nRep;
		cArchive=&repertory[nRep*nDim];
		cArchFit=&repertFit[nRep*nObj];
// 		cF=&repertoryF[nRep];
// 		cCR=&repertoryCR[nRep];
// 		memcpy(cF,&S_F[i*nPop],nPop*sizeof(double));
// 		memcpy(cCR,&S_CR[i*nPop],nPop*sizeof(double));
		joinSolutions_archive(i);
		for(j=0;j<nRep-nRep_tmp;j++)
			EMO_TEST_SUITE::evaluate_problems(testInstance,&cArchive[j*nDim],
				&cArchFit[j*nObj],nDim,1,nObj);
		iter_each+=nRep-nRep_tmp;
	}
}

void joinSolutions_archive(int iArch)
{
	int i;
	int j;
	int num;
	int rArch;
	int copy_ind;
	int a;
	int a_idx;
	int a_low,a_high;

	for(i=0;i<cnArch;i++)
	{
		if(archiveClass[i]>1)
			break;
	}
	num=i;

	copy_ind=0;
	for(i=0;i<arch_nGroup;i++)
	{
		selectSamples(num,iArch,&rArch,NULL,NULL);
		memcpy(&cArchive[copy_ind*nDim],&archive[rArch*nDim],nDim*sizeof(double));
		a_low=i*arch_dimIn1Group;
		a_high=(i+1)*arch_dimIn1Group;
		for(a=a_low;a<a_high;a++)
		{
			a_idx=Indexes[a];
			cArchive[INDEX(0,copy_ind,a_idx)]=uArchive[INDEX(0,iArch,a_idx)];
		}

// 		for(j=0;j<arch_nGroup;j++)
// 		{
// 			a_low=j*arch_dimIn1Group;
// 			a_high=(j+1)*arch_dimIn1Group;
// 			if(j==i)
// 			{
// 				for(a=a_low;a<a_high;a++)
// 				{
// 					a_idx=Indexes[a];
// 					cArchive[INDEX(0,copy_ind,a_idx)]=uArchive[INDEX(0,iArch,a_idx)];
// 				}
// 			}
// 			else
// 			{
// 				for(a=a_low;a<a_high;a++)
// 				{
// 					a_idx=Indexes[a];
// 					cArchive[INDEX(0,copy_ind,a_idx)]=archive[INDEX(0,rArch,a_idx)];
// 				}
// 			}
// 		}
		copy_ind++;
		nRep++;
	}
}