#include "EMO_test_suite.h"

//////////////////////////////////////////////////////////////////////////
#ifdef MPI_AVAILABLE
#include "mpi.h"
#endif
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <vector>
#include "assert.h"


//////////////////////////////////////////////////////////////////////////
using namespace std;

namespace EMO_TEST_SUITE
{

	//////////////////////////////////////////////////////////////////////////

	// #define E  2.7182818284590452353602874713526625
	#define PI 3.1415926535897932384626433832795029
	#define MYSIGN(x) ((x)>0?1.0:-1.0)

	int nvar,  nobj;                    //  the number of variables and objectives
	int position_parameters;

	//////////////////////////////////////////////////////////////////////////
	#if (defined(MPI_AVAILABLE) && defined(SUITE_PARALLEL))
	mpi_info suite_mpi_info;
	double tmp_f[MAX_VL];
	MPI_Comm comm_range;
	void set_MPI_Comm(MPI_Comm cur_comm);
	#endif

	void set_para(int dim, int nObj, int pos_para)
	{
		nvar=dim;
		nobj=nObj;
	// 	position_parameters=pos_para;
		position_parameters=2*(nobj-1);
	}

	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////// 
	//	Implementation
	//////////////////////////////////////////////////////////////////////////
	void setLimits(char *pro, double *minLimit, double *maxLimit, int dim)
	{
		//	DTLZ
		if(strcmp(pro, "DTLZ1") == 0
			||strcmp(pro, "DTLZ2") == 0
			||strcmp(pro, "DTLZ3") == 0
			||strcmp(pro, "DTLZ4") == 0
			||strcmp(pro, "DTLZ5") == 0
			||strcmp(pro, "DTLZ6") == 0
			||strcmp(pro, "DTLZ7") == 0)
		{
			for(int i=0;i<dim;i++)
			{
				minLimit[i]=0.0;
				maxLimit[i]=1.0;
			}
			return;
		}
		//	UF
		if(strcmp(pro, "UF1") == 0
			||strcmp(pro, "UF2") == 0
			||strcmp(pro, "UF5") == 0
			||strcmp(pro, "UF6") == 0
			||strcmp(pro, "UF7") == 0)
		{
			minLimit[0]=0.0;
			maxLimit[0]=1.0;
			for(int i=1;i<dim;i++)
			{
				minLimit[i]=-1.0;
				maxLimit[i]=1.0;
			}
			return;
		}
		if(strcmp(pro, "UF3") == 0)
		{
			for(int i=0;i<dim;i++)
			{
				minLimit[i]=0.0;
				maxLimit[i]=1.0;
			}
			return;
		}
		if(strcmp(pro, "UF4") == 0)
		{
			minLimit[0]=0.0;
			maxLimit[0]=1.0;
			for(int i=1;i<dim;i++)
			{
				minLimit[i]=-2.0;
				maxLimit[i]=2.0;
			}
			return;
		}
		if(strcmp(pro, "UF8") == 0
			||strcmp(pro, "UF9") == 0
			||strcmp(pro, "UF10") == 0)
		{
			minLimit[0]=0.0;
			maxLimit[0]=1.0;
			minLimit[1]=0.0;
			maxLimit[1]=1.0;
			for(int i=2;i<dim;i++)
			{
				minLimit[i]=-2.0;
				maxLimit[i]=2.0;
			}
			return;
		}
		//	CF
		if(strcmp(pro, "CF1") == 0)
		{
			for(int i=0;i<dim;i++)
			{
				minLimit[i]=0.0;
				maxLimit[i]=1.0;
			}
			return;
		}
		if(strcmp(pro, "CF2") == 0)
		{
			minLimit[0]=0.0;
			maxLimit[0]=1.0;
			for(int i=1;i<dim;i++)
			{
				minLimit[i]=-1.0;
				maxLimit[i]=1.0;
			}
			return;
		}
		if(strcmp(pro, "CF3") == 0
			||strcmp(pro, "CF4") == 0
			||strcmp(pro, "CF5") == 0
			||strcmp(pro, "CF6") == 0
			||strcmp(pro, "CF7") == 0)
		{
			minLimit[0]=0.0;
			maxLimit[0]=1.0;
			for(int i=1;i<dim;i++)
			{
				minLimit[i]=-2.0;
				maxLimit[i]=2.0;
			}
			return;
		}
		if(strcmp(pro, "CF8") == 0)
		{
			minLimit[0]=0.0;
			maxLimit[0]=1.0;
			minLimit[1]=0.0;
			maxLimit[1]=1.0;
			for(int i=2;i<dim;i++)
			{
				minLimit[i]=-4.0;
				maxLimit[i]=4.0;
			}
			return;
		}
		if(strcmp(pro, "CF9") == 0
			||strcmp(pro, "CF10") == 0)
		{
			minLimit[0]=0.0;
			maxLimit[0]=1.0;
			minLimit[1]=0.0;
			maxLimit[1]=1.0;
			for(int i=2;i<dim;i++)
			{
				minLimit[i]=-2.0;
				maxLimit[i]=2.0;
			}
			return;
		}
		//	WFG
		if(strcmp(pro, "WFG1") == 0
			||strcmp(pro, "WFG2") == 0
			||strcmp(pro, "WFG3") == 0
			||strcmp(pro, "WFG4") == 0
			||strcmp(pro, "WFG5") == 0
			||strcmp(pro, "WFG6") == 0
			||strcmp(pro, "WFG7") == 0
			||strcmp(pro, "WFG8") == 0
			||strcmp(pro, "WFG9") == 0)
		{
			for(int i=0;i<dim;i++)
			{
				minLimit[i]=0.0;
				maxLimit[i]=2*(i+1.0);
			}
			return;
		}
		if(strcmp(pro, "I1") == 0
			||strcmp(pro, "I2") == 0
			||strcmp(pro, "I3") == 0
			||strcmp(pro, "I4") == 0
			||strcmp(pro, "I5") == 0)
		{
			for(int i=0;i<dim;i++)
			{
				minLimit[i]=0.0;
				maxLimit[i]=1.0;
			}
			return;
		}
		{
			fprintf(stderr,"Unknown problem %s\n", pro );
			assert( false );
			return;
		}
	}

	void evaluate_problems(char *pro, double *xreal, double *obj, int dim, int nx, int _nobj)
	{
		// check the number of particles to calculate
	#if (defined(MPI_AVAILABLE) && defined(SUITE_PARALLEL))
		if((nx)!=suite_mpi_info.num || (_nobj)!=suite_mpi_info.num_obj)
			suite_mpi_info.update_n(nx, _nobj);
	#endif

		nobj=_nobj;
		nvar=dim;

	#if (defined(MPI_AVAILABLE) && defined(SUITE_PARALLEL))
		int count = 0;
		int n_start = suite_mpi_info.getStart();
		int n_end = suite_mpi_info.getEnd();
		int  my_tasks = suite_mpi_info.getMyTaskNum();
		int* tasks = suite_mpi_info.my_task;
		int* disps = suite_mpi_info.displacement;
		int vl_send = suite_mpi_info.vl_send;
		int* vl_task = suite_mpi_info.vl_task;
		int* vl_disp = suite_mpi_info.vl_disp;
		int i_idx;

		if( my_tasks )
		{

			for(int i=n_start; i<n_end; i++)
			{
				i_idx=i-n_start;

				//---DTLZ------------------//
				if (strcmp(pro, "DTLZ1") == 0)
					dtlz1(&xreal[i*dim], &tmp_f[i_idx*nobj], dim, nobj);
				else if (strcmp(pro, "DTLZ2") == 0)
					dtlz2(&xreal[i*dim], &tmp_f[i_idx*nobj], dim, nobj);
				else if (strcmp(pro, "DTLZ3") == 0)
					dtlz3(&xreal[i*dim], &tmp_f[i_idx*nobj], dim, nobj);
				else if (strcmp(pro, "DTLZ4") == 0)
					dtlz4(&xreal[i*dim], &tmp_f[i_idx*nobj], dim, nobj);
				else if (strcmp(pro, "DTLZ5") == 0)
					dtlz5(&xreal[i*dim], &tmp_f[i_idx*nobj], dim, nobj);
				else if (strcmp(pro, "DTLZ6") == 0)
					dtlz6(&xreal[i*dim], &tmp_f[i_idx*nobj], dim, nobj);
				else if (strcmp(pro, "DTLZ7") == 0)
					dtlz7(&xreal[i*dim], &tmp_f[i_idx*nobj], dim, nobj);
				else

					if(strcmp(pro, "UF1") == 0)
						UF1(&xreal[i*dim], &tmp_f[i_idx*nobj], dim);
					else if(strcmp(pro, "UF2") == 0)
						UF2(&xreal[i*dim], &tmp_f[i_idx*nobj], dim);
					else if(strcmp(pro, "UF3") == 0)
						UF3(&xreal[i*dim], &tmp_f[i_idx*nobj], dim);
					else if(strcmp(pro, "UF4") == 0)
						UF4(&xreal[i*dim], &tmp_f[i_idx*nobj], dim);
					else if(strcmp(pro, "UF5") == 0)
						UF5(&xreal[i*dim], &tmp_f[i_idx*nobj], dim);
					else if(strcmp(pro, "UF6") == 0)
						UF6(&xreal[i*dim], &tmp_f[i_idx*nobj], dim);
					else if(strcmp(pro, "UF7") == 0)
						UF7(&xreal[i*dim], &tmp_f[i_idx*nobj], dim);
					else if(strcmp(pro, "UF8") == 0)
						UF8(&xreal[i*dim], &tmp_f[i_idx*nobj], dim);
					else if(strcmp(pro, "UF9") == 0)
						UF9(&xreal[i*dim], &tmp_f[i_idx*nobj], dim);
					else if(strcmp(pro, "UF10") == 0)
						UF10(&xreal[i*dim], &tmp_f[i_idx*nobj], dim);
					else

						// CF is omitted
	// 					if(strcmp(pro, "CF1") == 0)
	// 						CF1(&xreal[i*dim], &tmp_f[i_idx*nobj], dim);

						//---WFG------------------//
						{
							wfg_eval(&xreal[i*dim],dim,nobj,pro,&tmp_f[i_idx*nobj]);
						}
			}
		}

		MPI_Allgatherv(tmp_f, vl_send, MPI_DOUBLE,
			obj, vl_task, vl_disp, MPI_DOUBLE, MPI_COMM_WORLD);

		tasks = NULL;
		disps = NULL;

		vl_task = NULL;
		vl_disp = NULL;

		MPI_Barrier(MPI_COMM_WORLD);

	#else

		for(int i=0; i<nx; i++)
		{
			if (strcmp(pro, "DTLZ1") == 0)
				dtlz1(&xreal[i*dim], &obj[i*nobj], dim, nobj);
			else if (strcmp(pro, "DTLZ2") == 0)
				dtlz2(&xreal[i*dim], &obj[i*nobj], dim, nobj);
			else if (strcmp(pro, "DTLZ3") == 0)
				dtlz3(&xreal[i*dim], &obj[i*nobj], dim, nobj);
			else if (strcmp(pro, "DTLZ4") == 0)
				dtlz4(&xreal[i*dim], &obj[i*nobj], dim, nobj);
			else if (strcmp(pro, "DTLZ5") == 0)
				dtlz5(&xreal[i*dim], &obj[i*nobj], dim, nobj);
			else if (strcmp(pro, "DTLZ6") == 0)
				dtlz6(&xreal[i*dim], &obj[i*nobj], dim, nobj);
			else if (strcmp(pro, "DTLZ7") == 0)
				dtlz7(&xreal[i*dim], &obj[i*nobj], dim, nobj);
			else

				if(strcmp(pro, "UF1") == 0)
					UF1(&xreal[i*dim], &obj[i*nobj], dim);
				else if(strcmp(pro, "UF2") == 0)
					UF2(&xreal[i*dim], &obj[i*nobj], dim);
				else if(strcmp(pro, "UF3") == 0)
					UF3(&xreal[i*dim], &obj[i*nobj], dim);
				else if(strcmp(pro, "UF4") == 0)
					UF4(&xreal[i*dim], &obj[i*nobj], dim);
				else if(strcmp(pro, "UF5") == 0)
					UF5(&xreal[i*dim], &obj[i*nobj], dim);
				else if(strcmp(pro, "UF6") == 0)
					UF6(&xreal[i*dim], &obj[i*nobj], dim);
				else if(strcmp(pro, "UF7") == 0)
					UF7(&xreal[i*dim], &obj[i*nobj], dim);
				else if(strcmp(pro, "UF8") == 0)
					UF8(&xreal[i*dim], &obj[i*nobj], dim);
				else if(strcmp(pro, "UF9") == 0)
					UF9(&xreal[i*dim], &obj[i*nobj], dim);
				else if(strcmp(pro, "UF10") == 0)
					UF10(&xreal[i*dim], &obj[i*nobj], dim);
				else

					//---WFG------------------//
					wfg_eval(&xreal[i*dim],dim,nobj,pro,&obj[i*nobj]);
		}

	#endif	// MPI_AVAILABLE
	}

	//	DTLZ
	void dtlz1 (double *xreal, double *obj, int dim, int nobj)
	{
		double sum=0;
		double gx;
		int i, j;

		for (i=nobj-1; i<dim; i++)
		{
			sum += pow ((xreal[i]-0.5), 2.0) - cos(20*PI*(xreal[i]-0.5));
		}
		gx = 100 * (sum+dim-nobj+1) + 1.0;
		sum = gx;
		for (j=0; j<nobj-1; j++)
		{
			sum = sum * xreal[j];
		}
		obj[0] = 0.5 * sum;

		for (i=1; i<nobj; i++)
		{
			sum = gx;
			for (j=0; j<nobj-1-i; j++)
			{
				sum = sum * xreal[j];
			}
			sum = sum * (1.0-xreal[nobj-1-i]);
			obj[i] = 0.5 * sum;
		}
		return;
	}


	void dtlz2 (double *xreal, double *obj, int dim, int nobj)
	{
		double sum=0;
		double gx;
		int i, j;

		for (i=nobj-1; i<dim; i++)
		{
			sum += pow ((xreal[i]-0.5), 2.0);
		}
		gx = 1.0 + sum;
		sum = gx;
		for (j=0; j<nobj-1; j++)
		{
			sum = sum * cos(xreal[j]*PI/2.0);
		}
		obj[0] = sum;

		for (i=1; i<nobj; i++)
		{
			sum = gx;
			for (j=0; j<nobj-1-i; j++)
			{
				sum = sum * cos(xreal[j]*PI/2.0);
			}
			sum = sum * sin(xreal[nobj-1-i]*PI/2.0);
			obj[i] = sum;
		}
		return;
	}


	void dtlz3 (double *xreal, double *obj, int dim, int nobj)
	{
		double sum=0;
		double gx;
		int i, j;

		for (i=nobj-1; i<dim; i++)
		{
			sum += pow ((xreal[i]-0.5), 2.0) - cos(20*PI*(xreal[i]-0.5));
		}
		gx = 100 * (sum+dim-nobj+1) + 1.0;
		sum = gx;
		for (j=0; j<nobj-1; j++)
		{
			sum = sum * cos(xreal[j]*PI/2.0);
		}
		obj[0] = sum;

		for (i=1; i<nobj; i++)
		{
			sum = gx;
			for (j=0; j<nobj-1-i; j++)
			{
				sum = sum * cos(xreal[j]*PI/2.0);
			}
			sum = sum * sin(xreal[nobj-1-i]*PI/2.0);
			obj[i] = sum;
		}
		return;

	}


	void dtlz4 (double *_xreal, double *obj, int dim, int nobj)
	{
		double sum=0;
		double gx;
		int i, j;
		double *xreal;
		xreal = (double *)malloc(dim * sizeof(double));
		memcpy(xreal, _xreal, dim * sizeof(double));

		for (i=nobj-1; i<dim; i++)
		{
			sum += pow ((xreal[i]-0.5), 2.0);
		}
		for (i=0; i<nobj-1; i++)
		{
			xreal[i]=pow((float)xreal[i],(float)100);
		}
		gx = 1.0 + sum;
		sum = gx;
		for (j=0; j<nobj-1; j++)
		{
			sum = sum * cos(xreal[j]*PI/2.0);
		}
		obj[0] = sum;

		for (i=1; i<nobj; i++)
		{
			sum = gx;
			for (j=0; j<nobj-1-i; j++)
			{
				sum = sum * cos(xreal[j]*PI/2.0);
			}
			sum = sum * sin(xreal[nobj-1-i]*PI/2.0);
			obj[i] = sum;
		}
		free(xreal);
		return;
	}


	void dtlz5 (double *xreal, double *obj, int dim, int nobj)
	{

		double sum=0;
		double gx;
		int i, j;
		double *x;
		x=(double*)malloc((nobj-1)*sizeof(double));

		for (i=nobj-1; i<dim; i++)
		{
			sum += pow ((xreal[i]-0.5), 2.0);
		}
		for (i=1; i<nobj-1; i++)
		{
			x[i] = PI/(4*(1+sum))*(1+2*sum*xreal[i]);
		}
		gx = 1.0 + sum;
		sum = gx;
		for (j=1; j<nobj-1; j++)
		{
			sum = sum * cos(x[j]);
		}
		sum = sum * cos(xreal[0]*PI/2.0);
		obj[0] = sum;

		for (i=1; i<nobj; i++)
		{
			sum = gx;
			for (j=1; j<nobj-1-i; j++)
			{
				sum = sum * cos(x[j]);
			}
			if (i == nobj-1)
			{
				sum = sum * sin(xreal[0]*PI/2.0);
			}
			else 
			{
				sum = sum * sin(x[nobj-1-i]);
				sum = sum * cos(xreal[0]*PI/2.0);
			}		
			obj[i] = sum;
		}
		free (x);
		return;
	}


	void dtlz6 (double *xreal, double *obj, int dim, int nobj)
	{
	
		double sum=0;
		double gx;
		int i, j;
		double *x;
		x=(double*)malloc((nobj-1)*sizeof(double));

		for (i=nobj-1; i<dim; i++)
		{
			sum += pow ((xreal[i]), 0.1);
		}
		for (i=1; i<nobj-1; i++)
		{
			x[i] = PI/(4*(1+sum))*(1+2*sum*xreal[i]);
		}
		gx = 1.0 + sum;
		sum = gx;
		for (j=1; j<nobj-1; j++)
		{
			sum = sum * cos(x[j]);
		}
		sum = sum * cos(xreal[0]*PI/2.0);
		obj[0] = sum;

		for (i=1; i<nobj; i++)
		{
			sum = gx;
			for (j=1; j<nobj-1-i; j++)
			{
				sum = sum * cos(x[j]);
			}
			if (i == nobj-1)
			{
				sum = sum * sin(xreal[0]*PI/2.0);
			}
			else 
			{
				sum = sum * sin(x[nobj-1-i]);
				sum = sum * cos(xreal[0]*PI/2.0);
			}		
			obj[i] = sum;
		}
		free (x);
		return;
	}

 
	void dtlz7 (double *xreal, double *obj, int dim, int nobj)
	{
		double sum=0, temp=0;
		double gx;
		int i;

		for (i=nobj-1; i<dim; i++)
		{
			sum += xreal[i];
		}
		gx = 1.0 + 9.0*sum/(dim-nobj+1.0);
		sum = gx;
		for (i=0; i<nobj-1; i++)
		{
			obj[i] = xreal[i];
		}
		for (i=0; i<nobj-1; i++)
		{
			temp += (obj[i]/(sum+1))*(1+sin(3*PI*obj[i]));
		}
		temp = nobj-temp;
		obj[nobj-1] = (sum+1)*temp;
		return;
	}

	/*
	#ifdef dtlz5_I_M 
	void test_problem (double *xreal, double *obj, double *constr)
	{
	
		double sum=0;
		double gx;
		int i, j;
		double *x;
		x=(double*)malloc((nobj-1)*sizeof(double));

		for (i=nobj-1; i<nreal; i++)
		{
			sum += pow ((xreal[i]-0.5), 2.0);
		}
	
		for (i=0; i<I_number-1; i++)
		{
			x[i] = xreal[i]*PI/2.0;
		}
		for (i=I_number-1; i<nobj-1; i++)
		{
			x[i] = PI/(4*(1+sum))*(1+2*sum*xreal[i]);
		}
		gx = 1.0 + 100 * sum;
	
		sum = gx;
		for (j=0; j<nobj-1; j++)
		{
			sum = sum * cos(x[j]);
		}
		obj[0] = sum;

		for (i=1; i<nobj; i++)
		{
			sum = gx;
			for (j=0; j<nobj-1-i; j++)
			{
				sum = sum * cos(x[j]);
			}
			sum = sum * sin(x[nobj-1-i]);
			obj[i] = sum;
		}
		free (x);
		return;
	}
	#endif
	*/

	//	DTLZ
	//////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////// 



	//////////////////////////////////////////////////////////////////////////
	//	UF	&	CF
	/****************************************************************************/
	// unconstraint test instances
	/****************************************************************************/
	void UF1(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;

		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}
		f[0] = x[0]				+ 2.0 * sum1 / (double)count1;
		f[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double)count2;
	}

	void UF2(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;

		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			if(j % 2 == 0) 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(24.0*PI*x[0]+4.0*j*PI/nx)+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			} 
			else 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(24.0*PI*x[0]+4.0*j*PI/nx)+2.0)*cos(6.0*PI*x[0]+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
		}
		f[0] = x[0]				+ 2.0 * sum1 / (double)count1;
		f[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double)count2;
	}

	void UF3(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;

		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		f[0] = x[0]				+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		f[1] = 1.0 - sqrt(x[0]) + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
	}

	void UF4(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj;

		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		f[0] = x[0]				+ 2.0*sum1 / (double)count1;
		f[1] = 1.0 - x[0]*x[0]	+ 2.0*sum2 / (double)count2;
	}

	void UF5(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;

		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 10.0; E = 0.1;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		hj = (0.5/N + E)*fabs(sin(2.0*N*PI*x[0]));
		f[0] = x[0]	      + hj + 2.0*sum1 / (double)count1;
		f[1] = 1.0 - x[0] + hj + 2.0*sum2 / (double)count2;
	}

	void UF6(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, hj, pj, N, E;
		N = 2.0; E = 0.1;

		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}

		hj = 2.0*(0.5/N + E)*sin(2.0*N*PI*x[0]);
		if(hj<0.0) hj = 0.0;
		f[0] = x[0]	      + hj + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		f[1] = 1.0 - x[0] + hj + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
	}

	void UF7(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;

		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		yj = pow(x[0],0.2);
		f[0] = yj	    + 2.0*sum1 / (double)count1;
		f[1] = 1.0 - yj + 2.0*sum2 / (double)count2;
	}

	void UF8(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj;

		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
			if(j % 3 == 1) 
			{
				sum1  += yj*yj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += yj*yj;
				count2++;
			}
			else
			{
				sum3  += yj*yj;
				count3++;
			}
		}
		f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
		f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
		f[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;
	}

	void UF9(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj, E;

		E = 0.1;
		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
			if(j % 3 == 1) 
			{
				sum1  += yj*yj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += yj*yj;
				count2++;
			}
			else
			{
				sum3  += yj*yj;
				count3++;
			}
		}
		yj = (1.0+E)*(1.0-4.0*(2.0*x[0]-1.0)*(2.0*x[0]-1.0));
		if(yj<0.0) yj = 0.0;
		f[0] = 0.5*(yj + 2*x[0])*x[1]		+ 2.0*sum1 / (double)count1;
		f[1] = 0.5*(yj - 2*x[0] + 2.0)*x[1] + 2.0*sum2 / (double)count2;
		f[2] = 1.0 - x[1]                   + 2.0*sum3 / (double)count3;
	}

	void UF10(double *x, double *f, const unsigned int nx)
	{
		unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj, hj;

		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
			hj = 4.0*yj*yj - cos(8.0*PI*yj) + 1.0;
			if(j % 3 == 1) 
			{
				sum1  += hj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += hj;
				count2++;
			}
			else
			{
				sum3  += hj;
				count3++;
			}
		}
		f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
		f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
		f[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;
	}

	/****************************************************************************/
	// constraint test instances
	/****************************************************************************/
	void CF1(double *x, double *f, double *c, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, N, a;
		N = 10.0; a = 1.0;

		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			if (j % 2 == 1) 
			{
				sum1  += yj*yj;
				count1++;
			} 
			else 
			{
				sum2  += yj*yj;
				count2++;
			}
		}
		f[0] = x[0]		  + 2.0*sum1 / (double)count1;
		f[1] = 1.0 - x[0] + 2.0*sum2 / (double)count2;
		c[0] = f[1] + f[0] - a*fabs(sin(N*PI*(f[0]-f[1]+1.0))) - 1.0; 
	}

	void CF2(double *x, double *f, double *c, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, N, a, t;
		N = 2.0; a = 1.0;

		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			if (j % 2 == 1) 
			{
				yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
				sum1  += yj*yj;
				count1++;
			} 
			else 
			{
				yj = x[j-1] - cos(6.0*PI*x[0] + j*PI/nx);
				sum2  += yj*yj;
				count2++;
			}
		}
		f[0] = x[0]		        + 2.0*sum1 / (double)count1;
		f[1] = 1.0 - sqrt(x[0]) + 2.0*sum2 / (double)count2;
		t	 = f[1] + sqrt(f[0]) - a*sin(N*PI*(sqrt(f[0])-f[1]+1.0)) - 1.0;
		c[0] = MYSIGN(t)*fabs(t)/(1+exp(4.0*fabs(t)));
	}

	void CF3(double *x, double *f, double *c, const unsigned int nx)
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj, N, a;
		N = 2.0; a = 1.0;

		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}

		f[0] = x[0]	           + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		f[1] = 1.0 - x[0]*x[0] + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		c[0] = f[1] + f[0]*f[0] - a*sin(N*PI*(f[0]*f[0]-f[1]+1.0)) - 1.0;
	}

	void CF4(double *x, double *f, double *c, const unsigned int nx)
	{
		unsigned int j;
		double sum1, sum2, yj, t;

		sum1   = sum2   = 0.0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
			if (j % 2 == 1) 
			{
				sum1  += yj*yj;
			} 
			else
			{
				if (j==2)
					sum2 += yj < 1.5-0.75*sqrt(2.0) ? fabs(yj) : (0.125+(yj-1)*(yj-1));
				else
					sum2  += yj*yj;
			}
		}
		f[0] = x[0]		  + sum1;
		f[1] = 1.0 - x[0] + sum2;
		t	 = x[1] - sin(6.0*x[0]*PI+2.0*PI/nx) - 0.5*x[0] + 0.25;
		c[0] = MYSIGN(t)*fabs(t)/(1+exp(4.0*fabs(t)));
	}

	void CF5(double *x, double *f, double *c, const unsigned int nx)
	{
		unsigned int j;
		double sum1, sum2, yj;

		sum1   = sum2   = 0.0;
		for(j = 2; j <= nx; j++) 
		{
			if (j % 2 == 1) 
			{
				yj    = x[j-1] - 0.8*x[0]*cos(6.0*PI*x[0] + j*PI/nx);
				sum1 += 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			} 
			else 
			{
				yj = x[j-1] - 0.8*x[0]*sin(6.0*PI*x[0] + j*PI/nx);
				if (j==2)
					sum2 += yj < 1.5-0.75*sqrt(2.0) ? fabs(yj) : (0.125+(yj-1)*(yj-1));
				else
					sum2 += 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			}
		}
		f[0] = x[0]		  + sum1;
		f[1] = 1.0 - x[0] + sum2;
		c[0] = x[1] - 0.8*x[0]*sin(6.0*x[0]*PI+2.0*PI/nx) - 0.5*x[0] + 0.25;
	}

	void CF6(double *x, double *f, double *c, const unsigned int nx)
	{
		unsigned int j;
		double sum1, sum2, yj;

		sum1   = sum2   = 0.0;
		for(j = 2; j <= nx; j++)
		{
			if (j % 2 == 1) 
			{
				yj     = x[j-1] - 0.8*x[0]*cos(6.0*PI*x[0] + j*PI/nx);
				sum1  += yj*yj;
			} 
			else 
			{
				yj     = x[j-1] - 0.8*x[0]*sin(6.0*PI*x[0] + j*PI/nx);
				sum2  += yj*yj;
			}
		}
		f[0] = x[0]		                 + sum1;
		f[1] = (1.0 - x[0])*(1.0 - x[0]) + sum2;
		c[0] = x[1]-0.8*x[0]*sin(6.0*x[0]*PI+2.0*PI/nx) - MYSIGN((x[0]-0.5)*(1.0-x[0]))*sqrt(fabs((x[0]-0.5)*(1.0-x[0])));
		c[1] = x[3]-0.8*x[0]*sin(6.0*x[0]*PI+4.0*PI/nx) - MYSIGN(0.25*sqrt(1-x[0])-0.5*(1.0-x[0]))*sqrt(fabs(0.25*sqrt(1-x[0])-0.5*(1.0-x[0])));
	}

	void CF7(double *x, double *f, double *c, const unsigned int nx)
	{
		unsigned int j;
		double sum1, sum2, yj;

		sum1   = sum2   = 0.0;
		for(j = 2; j <= nx; j++)
		{
			if (j % 2 == 1) 
			{
				yj     = x[j-1] - cos(6.0*PI*x[0] + j*PI/nx);
				sum1  += 2.0*yj*yj-cos(4.0*PI*yj)+1.0;
			} 
			else 
			{
				yj     = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
				if (j==2 || j==4)
					sum2 += yj*yj;
				else
					sum2  += 2.0*yj*yj-cos(4.0*PI*yj)+1.0;
			}
		}
		f[0] = x[0]		                 + sum1;
		f[1] = (1.0 - x[0])*(1.0 - x[0]) + sum2;
		c[0] = x[1]-sin(6.0*x[0]*PI+2.0*PI/nx) - MYSIGN((x[0]-0.5)*(1.0-x[0]))*sqrt(fabs((x[0]-0.5)*(1.0-x[0])));
		c[1] = x[3]-sin(6.0*x[0]*PI+4.0*PI/nx) - MYSIGN(0.25*sqrt(1-x[0])-0.5*(1.0-x[0]))*sqrt(fabs(0.25*sqrt(1-x[0])-0.5*(1.0-x[0])));
	}

	void CF8(double *x, double *f, double *c, const unsigned int nx)
	{
		unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj, N, a;
		N = 2.0; a = 4.0;

		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
			if(j % 3 == 1) 
			{
				sum1  += yj*yj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += yj*yj;
				count2++;
			}
			else
			{
				sum3  += yj*yj;
				count3++;
			}
		}
		f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
		f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
		f[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;
		c[0] = (f[0]*f[0]+f[1]*f[1])/(1-f[2]*f[2]) - a*fabs(sin(N*PI*((f[0]*f[0]-f[1]*f[1])/(1-f[2]*f[2])+1.0))) - 1.0;
	}

	void CF9(double *x, double *f, double *c, const unsigned int nx)
	{
		unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj, N, a;
		N = 2.0; a = 3.0;

		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
			if(j % 3 == 1) 
			{
				sum1  += yj*yj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += yj*yj;
				count2++;
			}
			else
			{
				sum3  += yj*yj;
				count3++;
			}
		}
		f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
		f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
		f[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;
		c[0] = (f[0]*f[0]+f[1]*f[1])/(1-f[2]*f[2]) - a*sin(N*PI*((f[0]*f[0]-f[1]*f[1])/(1-f[2]*f[2])+1.0)) - 1.0;
	}

	void CF10(double *x, double *f, double *c, const unsigned int nx)
	{
		unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj, hj, N, a;
		N = 2.0; a = 1.0;

		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
			hj = 4.0*yj*yj - cos(8.0*PI*yj) + 1.0;
			if(j % 3 == 1) 
			{
				sum1  += hj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += hj;
				count2++;
			}
			else
			{
				sum3  += hj;
				count3++;
			}
		}
		f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
		f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
		f[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;
		c[0] = (f[0]*f[0]+f[1]*f[1])/(1-f[2]*f[2]) - a*sin(N*PI*((f[0]*f[0]-f[1]*f[1])/(1-f[2]*f[2])+1.0)) - 1.0;
	}
	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////



	//////////////////////////////////////////////////////////////////////////
	//	WFG
	void wfg_eval( double* x, int n, int M, char* problem, double* fit )
	{
		vector< double > z;  // the decison vector
		vector< double > f;  // the fitness vector
		for( int i = 0; i < n; i++ )
		{
			//z.push_back( x[i]*(2*(i+1)) );
			z.push_back( x[i]);
		}	
		f = problem_calc_fitness( z, position_parameters, M, problem);
		for( int i = 0; i < (int)f.size(); i++ )
			fit[i] = f[i]; 
		vector<double>(z).swap(z);
		vector<double>(f).swap(f);
	}

	vector< double > problem_calc_fitness
	(
	 const vector< double >& z,
	 const int k,
	 const int M,
	 char* fn
	 )
	{
		if ( strcmp("WFG1",fn) == 0 )
		{
			return WFG1( z, k, M );
		}	
		else if ( strcmp("WFG2",fn)== 0  )
		{
			return WFG2( z, k, M );
		}
		else if ( strcmp("WFG3",fn)== 0  )
		{
			return WFG3( z, k, M );
		}
		else if ( strcmp("WFG4",fn) == 0 )
		{
			return WFG4( z, k, M );
		}
		else if ( strcmp("WFG5",fn) == 0 )
		{
			return WFG5( z, k, M );
		}
		else if ( strcmp("WFG6",fn) == 0)
		{
			return WFG6( z, k, M );
		}
		else if ( strcmp("WFG7",fn) == 0 )
		{
			return WFG7( z, k, M );
		}
		else if ( strcmp("WFG8",fn) == 0 )
		{
			return WFG8( z, k, M );
		}
		else if ( strcmp("WFG9",fn) == 0 )
		{
			return WFG9( z, k, M );
		}
		else if ( strcmp("I1",fn) == 0 )
		{
			return I1( z, k, M );
		}
		else if ( strcmp("I2",fn) == 0 )
		{
			return I2( z, k, M );
		}
		else if ( strcmp("I3",fn) == 0 )
		{
			return I3( z, k, M );
		}
		else if ( strcmp("I4",fn) == 0 )
		{
			return I4( z, k, M );
		}
		else if ( strcmp("I5",fn) == 0)
		{
			return I5( z, k, M );
		}
		else
		{
			fprintf(stderr,"Unknown problem %s\n", fn );
			assert( false );
			return vector< double >();
		}
	}

	//	Example problems
	//	//////////////////////////////////////////////////////////////////////////
	//** Reduces each paramer in "z" to the domain [0,1]. ***********************
	vector< double > WFG_normalise_z( const vector< double >& z )
	{
	  vector< double > result;

	  for( int i = 0; i < static_cast< int >( z.size() ); i++ )
	  {
		const double bound = 2.0*( i+1 );

		assert( z[i] >= 0.0   );
		assert( z[i] <= bound );

		result.push_back( z[i] / bound );
	  }

	  return result;
	}

	vector< double > WFG1
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );

		vector< double > y = WFG_normalise_z( z );

		y = WFG1_t1( y, k );
		y = WFG1_t2( y, k );
		y = WFG1_t3( y );
		y = WFG1_t4( y, k, M );

		return WFG1_shape( y );
	}

	vector< double > WFG2
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );
		assert( ( static_cast< int >( z.size() )-k ) % 2 == 0 );

		vector< double > y = WFG_normalise_z( z );

		y = WFG1_t1( y, k );
		y = WFG2_t2( y, k );
		y = WFG2_t3( y, k, M );

		return WFG2_shape( y );
	}

	vector< double > WFG3
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );
		assert( ( static_cast< int >( z.size() )-k ) % 2 == 0 );

		vector< double > y = WFG_normalise_z( z );

		y = WFG1_t1( y, k );
		y = WFG2_t2( y, k );
		y = WFG2_t3( y, k, M );

		return WFG3_shape( y );
	}

	vector< double > WFG4
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );

		vector< double > y = WFG_normalise_z( z );

		y = WFG4_t1( y );
		y = WFG2_t3( y, k, M );

		return WFG4_shape( y );
	}

	vector< double > WFG5
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );

		vector< double > y = WFG_normalise_z( z );

		y = WFG5_t1( y );
		y = WFG2_t3( y, k, M );

		return WFG4_shape( y );
	}

	vector< double > WFG6
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );

		vector< double > y = WFG_normalise_z( z );

		y = WFG1_t1( y, k );
		y = WFG6_t2( y, k, M );

		return WFG4_shape( y );
	}

	vector< double > WFG7
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );

		vector< double > y = WFG_normalise_z( z );

		y = WFG7_t1( y, k );
		y = WFG1_t1( y, k );
		y = WFG2_t3( y, k, M );

		return WFG4_shape( y );
	}

	vector< double > WFG8
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );

		vector< double > y = WFG_normalise_z( z );

		y = WFG8_t1( y, k );
		y = WFG1_t1( y, k );
		y = WFG2_t3( y, k, M );

		return WFG4_shape( y );
	}

	vector< double > WFG9
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );

		vector< double > y = WFG_normalise_z( z );

		y = WFG9_t1( y );
		y = WFG9_t2( y, k );
		y = WFG6_t2( y, k, M );

		return WFG4_shape( y );
	}

	vector< double > I1
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );

		vector< double > y = z;

		y = I1_t2( y, k );
		y = I1_t3( y, k, M );

		return I1_shape( y );
	}

	vector< double > I2
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );

		vector< double > y = z;

		y = I2_t1( y );
		y = I1_t2( y, k );
		y = I1_t3( y, k, M );

		return I1_shape( y );
	}

	vector< double > I3
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );

		vector< double > y = z;

		y = I3_t1( y );
		y = I1_t2( y, k );
		y = I1_t3( y, k, M );

		return I1_shape( y );
	}

	vector< double > I4
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );

		vector< double > y = z;

		y = I1_t2( y, k );
		y = I4_t3( y, k, M );

		return I1_shape( y );
	}

	vector< double > I5
	(
	 const vector< double >& z,
	 const int k,
	 const int M
	 )
	{
		assert( ArgsOK( z, k, M ) );

		vector< double > y = z;

		y = I3_t1( y );
		y = I1_t2( y, k );
		y = I4_t3( y, k, M );

		return I1_shape( y );
	}

	//	Example Shapes
	//////////////////////////////////////////////////////////////////////////
	//** Construct a vector of length M-1, with values "1,0,0,..." if ***********
	//** "degenerate" is true, otherwise with values "1,1,1,..." if   ***********
	//** "degenerate" is false.                                       ***********
	vector< short > WFG_create_A( const int M, const bool degenerate )
	{
	  assert( M >= 2 );

	  if ( degenerate )
	  {
		vector< short > A( M-1, 0 );
		A[0] = 1;

		return A;
	  }
	  else
	  {
		return vector< short >( M-1, 1 );
	  }
	}

	//** Given the vector "x" (the last value of which is the sole distance ****
	//** parameter), and the shape function results in "h", calculate the   ****
	//** scaled fitness values for a WFG problem.                           ****
	vector< double > WFG_calculate_f
	(
	  const vector< double >& x,
	  const vector< double >& h
	)
	{
	  assert( vector_in_01( x ) );
	  assert( vector_in_01( h ) );
	  assert( x.size() == h.size() );

	  const int M = static_cast< int >( h.size() );

	  vector< double > S;

	  for( int m = 1; m <= M; m++ )
	  {
		S.push_back( m*2.0 );
	  }

	  return calculate_f( 1.0, x, h, S );
	}

	vector< double > WFG1_shape( const vector< double >& t_p )
	{
		assert( vector_in_01( t_p ) );
		assert( t_p.size() >= 2 );

		const int M = static_cast< int >( t_p.size() );

		const vector< short >&  A = WFG_create_A( M, false );
		const vector< double >& x = calculate_x( t_p, A );

		vector< double > h;

		for( int m = 1; m <= M-1; m++ )
		{
			h.push_back( convex( x, m ) );
		}
		h.push_back( mixed( x, 5, 1.0 ) );

		return WFG_calculate_f( x, h );
	}

	vector< double > WFG2_shape( const vector< double >& t_p )
	{
		assert( vector_in_01( t_p ) );
		assert( t_p.size() >= 2 );

		const int M = static_cast< int >( t_p.size() );

		const vector< short >&  A = WFG_create_A( M, false );
		const vector< double >& x = calculate_x( t_p, A );

		vector< double > h;

		for( int m = 1; m <= M-1; m++ )
		{
			h.push_back( convex( x, m ) );
		}
		h.push_back( disc( x, 5, 1.0, 1.0 ) );

		return WFG_calculate_f( x, h );
	}

	vector< double > WFG3_shape( const vector< double >& t_p )
	{
		assert( vector_in_01( t_p ) );
		assert( t_p.size() >= 2 );

		const int M = static_cast< int >( t_p.size() );

		const vector< short >&  A = WFG_create_A( M, true );
		const vector< double >& x = calculate_x( t_p, A );

		vector< double > h;

		for( int m = 1; m <= M; m++ )
		{
			h.push_back( linear( x, m ) );
		}

		return WFG_calculate_f( x, h );
	}

	vector< double > WFG4_shape( const vector< double >& t_p )
	{
		assert( vector_in_01( t_p ) );
		assert( t_p.size() >= 2 );

		const int M = static_cast< int >( t_p.size() );

		const vector< short >&  A = WFG_create_A( M, false );
		const vector< double >& x = calculate_x( t_p, A );

		vector< double > h;

		for( int m = 1; m <= M; m++ )
		{
			h.push_back( concave( x, m ) );
		}

		return WFG_calculate_f( x, h );
	}

	vector< double > I1_shape( const vector< double >& t_p )
	{
		assert( vector_in_01( t_p ) );
		assert( t_p.size() >= 2 );

		const int M = static_cast< int >( t_p.size() );

		const vector< short >&  A = WFG_create_A( M, false );
		const vector< double >& x = calculate_x( t_p, A );

		vector< double > h;

		for( int m = 1; m <= M; m++ )
		{
			h.push_back( concave( x, m ) );
		}

		return calculate_f( 1.0, x, h, vector< double >( M, 1.0 ) );
	}

	//	Example Transitions
	//	//////////////////////////////////////////////////////////////////////////
	//** Construct a vector with the elements v[head], ..., v[tail-1]. **********
	vector< double > subvector
	(
	  const vector< double >& v,
	  const int head,
	  const int tail
	)
	{
	  assert( head >= 0 );
	  assert( head < tail );
	  assert( tail <= static_cast< int >( v.size() ) );

	  vector< double > result;

	  for( int i = head; i < tail; i++ )
	  {
		result.push_back( v[i] );
	  }

	  return result;
	}

	vector< double > WFG1_t1
	(
	 const vector< double >& y,
	 const int k
	 )
	{
		const int n = static_cast< int >( y.size() );

		assert( vector_in_01( y ) );
		assert( k >= 1 );
		assert( k <  n );

		vector< double > t;

		for( int i = 0; i < k; i++ )
		{
			t.push_back( y[i] );
		}

		for( int i = k; i < n; i++ )
		{
			t.push_back( s_linear( y[i], 0.35 ) );
		}

		return t;
	}

	vector< double > WFG1_t2
	(
	 const vector< double >& y,
	 const int k
	 )
	{
		const int n = static_cast< int >( y.size() );

		assert( vector_in_01( y ) );
		assert( k >= 1 );
		assert( k <  n );

		vector< double > t;

		for( int i = 0; i < k; i++ )
		{
			t.push_back( y[i] );
		}

		for( int i = k; i < n; i++ )
		{
			t.push_back( b_flat( y[i], 0.8, 0.75, 0.85 ) );
		}

		return t;
	}

	vector< double > WFG1_t3( const vector< double >& y )
	{
		const int n = static_cast< int >( y.size() );

		assert( vector_in_01( y ) );

		vector< double > t;

		for( int i = 0; i < n; i++ )
		{
			t.push_back( b_poly( y[i], 0.02 ) );
		}

		return t;
	}

	vector< double > WFG1_t4
	(
	 const vector< double >& y,
	 const int k,
	 const int M
	 )
	{
		const int n = static_cast< int >( y.size() );

		assert( vector_in_01( y ) );
		assert( k >= 1 );
		assert( k <  n );
		assert( M >= 2 );
		assert( k % ( M-1 ) == 0 );

		vector< double > w;

		for( int i = 1; i <= n; i++ )
		{
			w.push_back( 2.0*i );
		}

		vector< double > t;

		for( int i = 1; i <= M-1; i++ )
		{
			const int head = ( i-1 )*k/( M-1 );
			const int tail = i*k/( M-1 );

			const vector< double >& y_sub = subvector( y, head, tail );
			const vector< double >& w_sub = subvector( w, head, tail );

			t.push_back( r_sum( y_sub, w_sub ) );
		}

		const vector< double >& y_sub = subvector( y, k, n );
		const vector< double >& w_sub = subvector( w, k, n );

		t.push_back( r_sum( y_sub, w_sub ) );

		return t;
	}

	vector< double > WFG2_t2
	(
	 const vector< double >& y,
	 const int k
	 )
	{
		const int n = static_cast< int >( y.size() );
		const int l = n-k;

		assert( vector_in_01( y ) );
		assert( k >= 1 );
		assert( k <  n );
		assert( l % 2 == 0 );

		vector< double > t;

		for( int i = 0; i < k; i++ )
		{
			t.push_back( y[i] );
		}

		for( int i = k+1; i <= k+l/2; i++ )
		{
			const int head = k+2*( i-k )-2;
			const int tail = k+2*( i-k );

			t.push_back( r_nonsep( subvector( y, head, tail ), 2 ) );
		}

		return t;
	}

	vector< double > WFG2_t3
	(
	 const vector< double >& y,
	 const int k,
	 const int M
	 )
	{
		const int n = static_cast< int >( y.size() );

		assert( vector_in_01( y ) );
		assert( k >= 1 );
		assert( k <  n );
		assert( M >= 2 );
		assert( k % ( M-1 ) == 0 );

		const vector< double > w( n, 1.0 );

		vector< double > t;

		for( int i = 1; i <= M-1; i++ )
		{
			const int head = ( i-1 )*k/( M-1 );
			const int tail = i*k/( M-1 );

			const vector< double >& y_sub = subvector( y, head, tail );
			const vector< double >& w_sub = subvector( w, head, tail );

			t.push_back( r_sum( y_sub, w_sub ) );
		}

		const vector< double >& y_sub = subvector( y, k, n );
		const vector< double >& w_sub = subvector( w, k, n );

		t.push_back( r_sum( y_sub, w_sub ) );

		return t;
	}

	vector< double > WFG4_t1( const vector< double >& y )
	{
		const int n = static_cast< int >( y.size() );

		assert( vector_in_01( y ) );

		vector< double > t;

		for( int i = 0; i < n; i++ )
		{
			t.push_back( s_multi( y[i], 30, 10, 0.35 ) );
		}

		return t;
	}

	vector< double > WFG5_t1( const vector< double >& y )
	{
		const int n = static_cast< int >( y.size() );

		assert( vector_in_01( y ) );

		vector< double > t;

		for( int i = 0; i < n; i++ )
		{
			t.push_back( s_decept( y[i], 0.35, 0.001, 0.05 ) );
		}

		return t;
	}

	vector< double > WFG6_t2
	(
	 const vector< double >& y,
	 const int k,
	 const int M
	 )
	{
		const int n = static_cast< int >( y.size() );

		assert( vector_in_01( y ) );
		assert( k >= 1 );
		assert( k <  n );
		assert( M >= 2 );
		assert( k % ( M-1 ) == 0 );

		vector< double > t;

		for( int i = 1; i <= M-1; i++ )
		{
			const int head = ( i-1 )*k/( M-1 );
			const int tail = i*k/( M-1 );

			const vector< double >& y_sub = subvector( y, head, tail );

			t.push_back( r_nonsep( y_sub, k/( M-1 ) ) );
		}

		const vector< double >& y_sub = subvector( y, k, n );

		t.push_back( r_nonsep( y_sub, n-k ) );

		return t;
	}

	vector< double > WFG7_t1
	(
	 const vector< double >& y,
	 const int k
	 )
	{
		const int n = static_cast< int >( y.size() );

		assert( vector_in_01( y ) );
		assert( k >= 1 );
		assert( k <  n );

		const vector< double > w( n, 1.0 );

		vector< double > t;

		for( int i = 0; i < k; i++ )
		{
			const vector< double >& y_sub = subvector( y, i+1, n );
			const vector< double >& w_sub = subvector( w, i+1, n );

			const double u = r_sum( y_sub, w_sub );

			t.push_back( b_param( y[i], u, 0.98/49.98, 0.02, 50 ) );
		}

		for( int i = k; i < n; i++ )
		{
			t.push_back( y[i] );
		}

		return t;
	}

	vector< double > WFG8_t1
	(
	 const vector< double >& y,
	 const int k
	 )
	{
		const int n = static_cast< int >( y.size() );

		assert( vector_in_01( y ) );
		assert( k >= 1 );
		assert( k <  n );

		const vector< double > w( n, 1.0 );

		vector< double > t;

		for( int i = 0; i < k; i++ )
		{
			t.push_back( y[i] );
		}

		for( int i = k; i < n; i++ )
		{
			const vector< double >& y_sub = subvector( y, 0, i );
			const vector< double >& w_sub = subvector( w, 0, i );

			const double u = r_sum( y_sub, w_sub );

			t.push_back( b_param( y[i], u, 0.98/49.98, 0.02, 50 ) );
		}

		return t;
	}

	vector< double > WFG9_t1( const vector< double >& y )
	{
		const int n = static_cast< int >( y.size() );

		assert( vector_in_01( y ) );

		const vector< double > w( n, 1.0 );

		vector< double > t;

		for( int i = 0; i < n-1; i++ )
		{
			const vector< double >& y_sub = subvector( y, i+1, n );
			const vector< double >& w_sub = subvector( w, i+1, n );

			const double u = r_sum( y_sub, w_sub );

			t.push_back( b_param( y[i], u, 0.98/49.98, 0.02, 50 ) );
		}

		t.push_back( y.back() );

		return t;
	}

	vector< double > WFG9_t2
	(
	 const vector< double >& y,
	 const int k
	 )
	{
		const int n = static_cast< int >( y.size() );

		assert( vector_in_01( y ) );
		assert( k >= 1 );
		assert( k <  n );

		vector< double > t;

		for( int i = 0; i < k; i++ )
		{
			t.push_back( s_decept( y[i], 0.35, 0.001, 0.05 ) );
		}

		for( int i = k; i < n; i++ )
		{
			t.push_back( s_multi( y[i], 30, 95, 0.35 ) );
		}

		return t;
	}

	vector< double > I1_t2
	(
	 const vector< double >& y,
	 const int k
	 )
	{
		return WFG1_t1( y, k );
	}

	vector< double > I1_t3
	(
	 const vector< double >& y,
	 const int k,
	 const int M
	 )
	{
		return WFG2_t3( y, k, M );
	}

	vector< double > I2_t1( const vector< double >& y )
	{
		return WFG9_t1( y );
	}

	vector< double > I3_t1( const vector< double >& y )
	{
		const int n = static_cast< int >( y.size() );

		assert( vector_in_01( y ) );

		const vector< double > w( n, 1.0 );

		vector< double > t;

		t.push_back( y.front() );

		for( int i = 1; i < n; i++ )
		{
			const vector< double >& y_sub = subvector( y, 0, i );
			const vector< double >& w_sub = subvector( w, 0, i );

			const double u = r_sum( y_sub, w_sub );

			t.push_back( b_param( y[i], u, 0.98/49.98, 0.02, 50 ) );
		}

		return t;
	}

	vector< double > I4_t3
	(
	 const vector< double >& y,
	 const int k,
	 const int M
	 )
	{
		return WFG6_t2( y, k, M );
	}

	//	Framework Functions
	//	//////////////////////////////////////////////////////////////////////////
	vector< double > normalise_z
	(
	 const vector< double >& z,
	 const vector< double >& z_max
	 )
	{
		vector< double > result;

		for( int i = 0; i < static_cast< int >( z.size() ); i++ )
		{
			assert( z[i] >= 0.0 );
			assert( z[i] <= z_max[i] );
			assert( z_max[i] > 0.0 );

			result.push_back( z[i] / z_max[i] );
		}

		return result;
	}

	vector< double >  calculate_x
	(
	 const vector< double >& t_p,
	 const vector< short >& A
	 )
	{
		assert( vector_in_01( t_p ) );
		assert( t_p.size() != 0 );
		assert( A.size() == t_p.size()-1 );

		vector< double > result;

		for( int i = 0; i < static_cast< int >( t_p.size() ) - 1; i++ )
		{
			assert( A[i] == 0 || A[i] == 1 );

			const double tmp1 = (t_p.back() > A[i]) ? t_p.back() : A[i];// std::max< double >(t_p.back(), A[i]);
			result.push_back( tmp1*( t_p[i] - 0.5 ) + 0.5 );
		}

		result.push_back( t_p.back() );

		return result;
	}

	vector< double >  calculate_f
	(
	 const double&           D,
	 const vector< double >& x,
	 const vector< double >& h,
	 const vector< double >& S
	 )
	{
		assert( D > 0.0 );
		assert( vector_in_01( x ) );
		assert( vector_in_01( h ) );
		assert( x.size() == h.size() );
		assert( h.size() == S.size() );

		vector< double > result;

		for( int i = 0; i < static_cast< int >( h.size() ); i++ )
		{
			assert( S[i] > 0.0 );

			result.push_back( D*x.back() + S[i]*h[i] );
		}

		return result;
	}

	//	Misc
	//	//////////////////////////////////////////////////////////////////////////
	double correct_to_01( const double& a, const double& epsilon )
	{
		assert( epsilon >= 0.0 );

		const double min = 0.0;
		const double max = 1.0;

		const double min_epsilon = min - epsilon;
		const double max_epsilon = max + epsilon;

		if ( a <= min && a >= min_epsilon )
		{
			return min;
		}
		else if ( a >= max && a <= max_epsilon )
		{
			return max;
		}
		else
		{
			return a;
		}
	}

	bool vector_in_01( const vector< double >& x )
	{
		for( int i = 0; i < static_cast< int >( x.size() ); i++ )
		{
			if( x[i] < 0.0 || x[i] > 1.0 )
			{
				return false;
			}
		}

		return true;
	}

	//	Shape Functions
	//	//////////////////////////////////////////////////////////////////////////
	//** True if all elements of "x" are in [0,1], and m is in [1, x.size()]. ***
	bool shape_args_ok( const vector< double >& x, const int m )
	{
	  const int M = static_cast< int >( x.size() );

	  return vector_in_01( x ) && m >= 1 && m <= M;
	}

	double linear( const vector< double >& x, const int m )
	{
		assert( shape_args_ok( x, m ) );

		const int M = static_cast< int >( x.size() );
		double result = 1.0;

		for( int i=1; i <= M-m; i++ )
		{
			result *= x[i-1];
		}

		if( m != 1 )
		{
			result *= 1 - x[M-m];
		}

		return correct_to_01( result );
	}

	double convex( const vector< double >& x, const int m )
	{
		assert( shape_args_ok( x, m ) );

		const int M = static_cast< int >( x.size() );
		double result = 1.0;

		for( int i=1; i <= M-m; i++ )
		{
			result *= 1.0 - cos( x[i-1]*PI/2.0 );
		}

		if( m != 1 )
		{
			result *= 1.0 - sin( x[M-m]*PI/2.0 );
		}

		return correct_to_01( result );
	}

	double concave( const vector< double >& x, const int m )
	{
		assert( shape_args_ok( x, m ) );

		const int M = static_cast< int >( x.size() );
		double result = 1.0;

		for( int i=1; i <= M-m; i++ )
		{
			result *= sin( x[i-1]*PI/2.0 );
		}

		if( m != 1 )
		{
			result *= cos( x[M-m]*PI/2.0 );
		}

		return correct_to_01( result );
	}

	double mixed
	(
	 const vector< double >& x,
	 const int A,
	 const double& alpha
	 )
	{
		assert( vector_in_01( x ) );
		assert( x.size() != 0   );
		assert( A        >= 1   );
		assert( alpha    >  0.0 );

		const double tmp = 2.0*A*PI;

		return correct_to_01( pow( 1.0-x[0]-cos( tmp*x[0] + PI/2.0 )/tmp, alpha ) );
	}

	double disc
	(
	 const vector< double >& x,
	 const int A,
	 const double& alpha,
	 const double& beta
	 )
	{
		assert( vector_in_01( x ) );
		assert( x.size() != 0   );
		assert( A        >= 1   );
		assert( alpha    >  0.0 );
		assert( beta     >  0.0 );

		const double tmp1 = A*pow( x[0], beta )*PI;
		return correct_to_01( 1.0 - pow( x[0], alpha )*pow( cos( tmp1 ), 2.0 ) );
	}

	//	TransFunctions
	//	//////////////////////////////////////////////////////////////////////////
	double b_poly( const double& y, const double& alpha )
	{
		assert( y >= 0.0 );
		assert( y <= 1.0 );
		assert( alpha >  0.0 );
		assert( alpha != 1.0 );

		return correct_to_01( pow( y, alpha ) );
	}

	double b_flat
	(
	 const double& y,
	 const double& A,
	 const double& B,
	 const double& C
	 )
	{
		assert( y >= 0.0 );
		assert( y <= 1.0 );
		assert( A >= 0.0 );
		assert( A <= 1.0 );
		assert( B >= 0.0 );
		assert( B <= 1.0 );
		assert( C >= 0.0 );
		assert( C <= 1.0 );
		assert( B < C );
		assert( B != 0.0 || A == 0.0 );
		assert( B != 0.0 || C != 1.0 );
		assert( C != 1.0 || A == 1.0 );
		assert( C != 1.0 || B != 0.0 );

		const double tmp1 = min( 0.0, floor( y-B ) ) * A*( B-y )/B;
		const double tmp2 = min( 0.0, floor( C-y ) ) * ( 1.0-A )*( y-C )/( 1.0-C );

		return correct_to_01( A+tmp1-tmp2 );
	}

	double b_param
	(
	 const double& y,
	 const double& u,
	 const double& A,
	 const double& B,
	 const double& C
	 )
	{
		assert( y >= 0.0 );
		assert( y <= 1.0 );
		assert( u >= 0.0 );
		assert( u <= 1.0 );
		assert( A > 0.0 );
		assert( A < 1.0 );
		assert( B > 0.0 );
		assert( B < C );

		const double v = A - ( 1.0-2.0*u )*fabs( floor( 0.5-u )+A );

		return correct_to_01( pow( y, B + ( C-B )*v ) );
	}

	double s_linear( const double& y, const double& A )
	{
		assert( y >= 0.0 );
		assert( y <= 1.0 );
		assert( A > 0.0 );
		assert( A < 1.0 );

		return correct_to_01( fabs( y-A )/fabs( floor( A-y )+A ) );
	}

	double s_decept
	(
	 const double& y,
	 const double& A,
	 const double& B,
	 const double& C
	 )
	{
		assert( y >= 0.0 );
		assert( y <= 1.0 );
		assert( A > 0.0 );
		assert( A < 1.0 );
		assert( B > 0.0 );
		assert( B < 1.0 );
		assert( C > 0.0 );
		assert( C < 1.0 );
		assert( A - B > 0.0 );
		assert( A + B < 1.0 );

		const double tmp1 = floor( y-A+B )*( 1.0-C+( A-B )/B )/( A-B );
		const double tmp2 = floor( A+B-y )*( 1.0-C+( 1.0-A-B )/B )/( 1.0-A-B );

		return correct_to_01( 1.0 + ( fabs( y-A )-B )*( tmp1 + tmp2 + 1.0/B ) );
	}

	double s_multi
	(
	 const double& y,
	 const int A,
	 const double& B,
	 const double& C
	 )
	{
		assert( y >= 0.0 );
		assert( y <= 1.0 );
		assert( A >= 1 );
		assert( B >= 0.0 );
		assert( ( 4.0*A+2.0 )*PI >= 4.0*B );
		assert( C > 0.0 );
		assert( C < 1.0 );

		const double tmp1 = fabs( y-C )/( 2.0*( floor( C-y )+C ) );
		const double tmp2 = ( 4.0*A+2.0 )*PI*( 0.5-tmp1 );

		return correct_to_01( ( 1.0 + cos( tmp2 ) + 4.0*B*pow( tmp1, 2.0 ) )/( B+2.0 ) );
	}

	double r_sum
	(
	 const vector< double >& y,
	 const vector< double >& w
	 )
	{
		assert( y.size() != 0        );
		assert( w.size() == y.size() );
		assert( vector_in_01( y ) );

		double numerator   = 0.0;
		double denominator = 0.0;

		for( int i = 0; i < static_cast< int >( y.size() ); i++ )
		{
			assert( w[i] > 0.0 );

			numerator   += w[i]*y[i];
			denominator += w[i];
		}

		return correct_to_01( numerator / denominator );
	}

	double r_nonsep( const std::vector< double >& y, const int A )
	{
		const int y_len = static_cast< int >( y.size() );

		assert( y_len != 0 );
		assert( vector_in_01( y ) );
		assert( A >= 1 );
		assert( A <= y_len );
		assert( y.size() % A == 0 );

		double numerator = 0.0;

		for( int j = 0; j < y_len; j++ )
		{
			numerator += y[j];

			for( int k = 0; k <= A-2; k++ )
			{
				numerator += fabs( y[j] - y[( j+k+1 ) % y_len] );
			}
		}

		const double tmp = ceil( A/2.0 );
		const double denominator = y_len*tmp*( 1.0 + 2.0*A - 2.0*tmp )/A;

		return correct_to_01( numerator / denominator );
	}
}
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////// 
////////////////////////////////////////////////////////////////////////// 