# include "global.h"

void DE_gen_uTrail()
{
	int r1,r2,r3;
	int j,k;
	int k_ind;
	for(j=task_l;j<task_r;j++)
	{
		selectSamples(nPop,j,&r1,&r2,NULL);

		do 
		{
			r3 = rnd(0, cnArch - 1);
		} while (archiveClass[r3]!=1);

		k_ind = rnd(0,nDim - 1);

		for(k=0;k<nDim;k++)
		{
			if(flip_r((float)CR) || k == k_ind){
				uTrail[j*nDim+k] = xCurrent[j*nDim+k] + F*(xBest[k] - xCurrent[j*nDim+k]) + 
					F*(xCurrent[r1*nDim+k] - xCurrent[r2*nDim+k]) + 
					F*(archive[r3*nDim+k] - xCurrent[j*nDim+k]);
				if(uTrail[j*nDim+k]<minLimit[k]){
					uTrail[j*nDim+k] = (minLimit[k]+xCurrent[j*nDim+k])/2.0;
				}
				if(uTrail[j*nDim+k]>maxLimit[k]){
					uTrail[j*nDim+k] = (maxLimit[k]+xCurrent[j*nDim+k])/2.0;
				}
			}
			else
			{
				uTrail[j*nDim+k] = xCurrent[j*nDim+k];
			}
		}
	}
}

void DE_gen_uTrail_ada()
{
	int r1,r2,r3;
	int j,k;
	int k_ind;
	int j_idx;
	for(j=0;j<nPop;j++)
	{
		selectSamples(nPop,j,&r1,&r2,NULL);

		do 
		{
			r3 = rnd(0, cnArch - 1);
		} while (archiveClass[r3]!=1);

		k_ind = rnd(0,nDim - 1);

		j_idx = j;

		for(k=0;k<nDim;k++)
		{
			if(flip_r((float)S_CR[j_idx]/*CR*/) || k == k_ind){
				uTrail[j*nDim+k] = xCurrent[j*nDim+k] + S_F[j_idx]/*F*/*(xBest[k] - xCurrent[j*nDim+k]) + 
					S_F[j_idx]/*F*/*(xCurrent[r1*nDim+k] - xCurrent[r2*nDim+k]) + 
					S_F[j_idx]/*F*/*(archive[r3*nDim+k] - xCurrent[j*nDim+k]);
			}
			else
			{
				uTrail[j*nDim+k] = xCurrent[j*nDim+k];
			}
			if(uTrail[j*nDim+k]<minLimit[k]){
				uTrail[j*nDim+k] = (minLimit[k]+xCurrent[j*nDim+k])/2.0;
			}
			if(uTrail[j*nDim+k]>maxLimit[k]){
				uTrail[j*nDim+k] = (maxLimit[k]+xCurrent[j*nDim+k])/2.0;
			}
		}
	}
}

void DE_gen_archive()
{
	if(cnArch<4||archiveClass[3]>1)
		return;
	int i,j;
	int r1,r2,r3;
	int num;
	int j_ind;

	for(i=0;i<cnArch;i++)
	{
		if(archiveClass[i]>1)
			break;
	}
	num=i;

	allocateTask(cnArch);

	for(i=task_l;i<task_r;i++)
	{
		selectSamples(num,i,&r1,&r2,&r3);
		j_ind=rnd(0,nDim - 1);

		for(j=0;j<nDim;j++)
		{
			if(flip_r((float)/*arch_CR[i]*/CR) || j == j_ind)
			{
				uArchive[INDEX(0,i,j)]=archive[INDEX(0,r1,j)]+
					/*arch_F[i]*/F*(archive[INDEX(0,r2,j)]-archive[INDEX(0,r3,j)]);
			}
			else
			{
				uArchive[INDEX(0,i,j)]=archive[INDEX(0,i,j)];
			}
			if(uArchive[INDEX(0,i,j)]<minLimit[j]){
				uArchive[INDEX(0,i,j)] = minLimit[j];
			}
			if(uArchive[INDEX(0,i,j)]>maxLimit[j]){
				uArchive[INDEX(0,i,j)] = maxLimit[j];
			}
		}
	}
}

void DE_gen_archive2()
{
	if(cnArch<4||archiveClass[3]>1)
		return;
	int i;
	int r1,r2,r3;
	int num;
	int j_ind;
	int g_ind;
	int iGrp;
	int a;
	int a_low,a_high;
	int a_idx;
	int r;

	for(i=0;i<cnArch;i++)
	{
		if(archiveClass[i]>1)
			break;
	}
	num=i;

	cArchive=&repertory[nRep*nDim];
	cArchFit=&repertFit[nRep*nObj];
// 	uF=&repertoryF[nRep];
// 	uCR=&repertoryCR[nRep];

// 	generate_arch_F_CR();

	for(i=0;i<cnArch;i++)
	{
		selectSamples(num,i,&r1,&r2,&r3);
		j_ind=rnd(0,nDim - 1);
		g_ind=rnd(0,arch_nGroup-1);

		for(iGrp=0;iGrp<arch_nGroup;iGrp++)
		{
			a_low=iGrp*arch_dimIn1Group;
			a_high=(iGrp+1)*arch_dimIn1Group;
			if(flip_r((float)0.63)&&iGrp!=g_ind)
			{
				r=rnd(0,cnArch-1);
				for(a=a_low;a<a_high;a++)
				{
					a_idx=archIndex[a];
					cArchive[INDEX(0,i,a_idx)]=archive[INDEX(0,r,a_idx)];
				}
			}
			else
			{
				for(a=a_low;a<a_high;a++)
				{
					a_idx=archIndex[a];
					if(flip_r((float)CR)||a_idx==j_ind)
					{
						cArchive[INDEX(0,i,a_idx)]=archive[INDEX(0,r1,a_idx)]+
							/*arch_F[i]*/F*(archive[INDEX(0,r2,a_idx)]-archive[INDEX(0,r3,a_idx)]);
						if(cArchive[INDEX(0,i,a_idx)]<minLimit[a_idx]){
							cArchive[INDEX(0,i,a_idx)] = minLimit[a_idx];
						}
						if(uArchive[INDEX(0,i,a_idx)]>maxLimit[a_idx]){
							cArchive[INDEX(0,i,a_idx)] = maxLimit[a_idx];
						}
					}
					else
					{
						cArchive[INDEX(0,i,a_idx)]=archive[INDEX(0,i,a_idx)];
					}
				}
			}
		}
	}

// 	memcpy(uF,arch_F,cnArch*sizeof(double));
// 	memcpy(uCR,arch_CR,cnArch*sizeof(double));

	nRep+=cnArch;
	for(i=0;i<cnArch;i++)
		EMO_TEST_SUITE::evaluate_problems(testInstance,
			&cArchive[i*nDim],&cArchFit[i*nObj],nDim,1,nObj);
	iter_each+=cnArch;
// 	if(iter>=maxIteration)
// 		return;
}

void DE_gen_archive3()
{
	if(cnArch<4||archiveClass[3]>1)
		return;
	int i,j;
	int num;
	int j_ind;

	for(i=0;i<cnArch;i++)
	{
		if(archiveClass[i]>1)
			break;
	}
	num=i;

	cArchive=&repertory[nRep*nDim];
	cArchFit=&repertFit[nRep*nObj];

	int count=0;

	for(j=0;j<nDim;j++)
	{
		if(isAllTheSame(num,j))
		{
			dimFlag[count++]=j;
		}
	}
	if(!count)
		return;

	for(i=0;i<num;i++)
	{
		j_ind=rnd(0,count - 1);
		j_ind=dimFlag[j_ind];

		memcpy(&cArchive[INDEX(0,i,0)],&archive[INDEX(0,i,0)],nDim*sizeof(double));
		cArchive[INDEX(0,i,j_ind)]+=(maxLimit[j_ind]-minLimit[j_ind])/100000.0*gaussrand(0.0,1.0);
		if(cArchive[INDEX(0,i,j_ind)]<minLimit[j_ind]){
			cArchive[INDEX(0,i,j_ind)] = (minLimit[j_ind]+archive[INDEX(0,i,j_ind)])/2.0;
		}
		if(cArchive[INDEX(0,i,j_ind)]>maxLimit[j_ind]){
			cArchive[INDEX(0,i,j_ind)] = (maxLimit[j_ind]+archive[INDEX(0,i,j_ind)])/2.0;
		}
	}

	nRep+=cnArch;
	for(i=0;i<cnArch;i++)
		EMO_TEST_SUITE::evaluate_problems(testInstance,
			&cArchive[i*nDim],&cArchFit[i*nObj],nDim,1,nObj);
	iter_each+=cnArch;
// 	if(iter>=maxIteration)
// 		return;
}

int isAllTheSame(int size,int iDim)
{
	int i;
	double value=archive[INDEX(0,0,iDim)];
	for(i=1;i<size;i++)
	{
		if(archive[INDEX(0,i,iDim)]!=value)
			return 0;
	}
	return 1;
}

void update_F_CR_mu()
{
	int i;
	double sum1,sum2,sum3;
	int count;
	sum1=0.0;
	sum2=0.0;
	sum3=0.0;
	count=0;
	for(i=0;i<nPop;i++)
	{
		if(Sflag[i])
		{
			sum1+=(S_F[i]*S_F[i]);
			sum2+=S_F[i];
			sum3+=S_CR[i];
			count++;
		}
	}
	if(count)
	{
		F_mu=(1.0-c_para)*F_mu+c_para*sum1/sum2;
		CR_mu=(1.0-c_para)*CR_mu+c_para*sum3/count;
	}
}

void generate_F_CR()
{
	int j;
	int c_idx;
	for(j=0;j<nPop;j++)
	{
		c_idx=j;
		do {
			S_F[c_idx]=cauchyrand(F_mu,0.1);
		} while (S_F[c_idx]<=0.0);
		if(S_F[c_idx]>1.0)
			S_F[c_idx]=1.0;

		S_CR[c_idx]=gaussrand(CR_mu,0.1);
		if(S_CR[c_idx]<0.0)
			S_CR[c_idx]=0.0;
		else if(S_CR[c_idx]>1.0)
			S_CR[c_idx]=1.0;
	}
}

void generate_arch_F_CR()
{
	int i;
	for(i=0;i<cnArch;i++)
	{
		if(rnd_uni(&rnd_uni_init)<t1)
			arch_F[i]=0.1+0.9*rnd_uni(&rnd_uni_init);
		if(rnd_uni(&rnd_uni_init)<t2)
			arch_CR[i]=rnd_uni(&rnd_uni_init);
	}
}