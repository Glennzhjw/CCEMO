#ifndef _EMO_TEST_SUITE_
#define _EMO_TEST_SUITE_

//////////////////////////////////////////////////////////////////////////
//	controlling macro define
#define MPI_AVAILABLE
// #define SUITE_PARALLEL

//////////////////////////////////////////////////////////////////////////
#ifdef MPI_AVAILABLE
#include "mpi.h"
#endif
#include <stdio.h>
#include <iostream>
#include <vector>


//////////////////////////////////////////////////////////////////////////
using namespace std;

//////////////////////////////////////////////////////////////////////////

namespace EMO_TEST_SUITE
{

#define MAX_POP (256)
#define MAX_OBJ (64)
#define MAX_VL  (MAX_POP*MAX_OBJ)
#define MAX_DIM (100)
#define MAX_MPI (256)

//////////////////////////////////////////////////////////////////////////
//	MPI setting
#if (defined(MPI_AVAILABLE) && defined(SUITE_PARALLEL))
// Store mpi information
struct mpi_info
{
	int id;// mpi rank, the id_th process
	int size; // group size
	int* my_task; // store num of tasks of each rank
	int* displacement; // 
	int vl_send;
	int* vl_task;
	int* vl_disp;

	int num;//number of particles to calculate
	int num_obj;//number of objectives

	void initialize_mpi(int rank, int group, int n, int nObj)
	{
		id = rank;
		size = group;
		num = n;
		num_obj = nObj;

		my_task = (int*)malloc(size*sizeof(int));
		displacement = (int*)malloc(size*sizeof(int));
		vl_task = (int*)malloc(size*sizeof(int));
		vl_disp = (int*)malloc(size*sizeof(int));

		update_n(n, nObj);
	}

	// Update num of particles to calculate
	void update_n(int n, int nObj)
	{
		num = n;
		num_obj = nObj;
		int quotent, remain;

		if(num<size){
			if(num){
				for(int i=0;i<size;i++){
					if(i<num)
					{
						my_task[i]=1;
						displacement[i]=i;
					}
					else
					{
						my_task[i]=0;
						displacement[i]=num;
					}
				}
			}
			else{
				for(int i=0;i<size;i++){
					my_task[i]=0;
					displacement[i]=0;
				}
			}
		}
		else{
			quotent=num/size;
			remain=num%size;
			for(int i=0;i<size;i++){
				my_task[i]=quotent;
				if(i<remain)
					my_task[i]++;
			}
			displacement[0]=0;
			for(int i=1;i<size;i++){
				displacement[i]=displacement[i-1]+my_task[i-1];
			}
		}

		vl_send = my_task[id]*num_obj;
		for(int i=0;i<size;i++)
		{
			vl_task[i]=my_task[i]*num_obj;
			vl_disp[i]=displacement[i]*num_obj;
		}
	}

	// free memory
	void clean()
	{
		if(my_task) free(my_task);
		if(displacement) free(displacement);

		if(vl_task) free(vl_task);
		if(vl_disp) free(vl_disp);
	}

	// which particle is the first one for rank of id to calculate
	int getStart()
	{
		return displacement[id];
	}

	// the last one + 1
	int getEnd()
	{
		return (displacement[id] + my_task[id]);
	}

	// number of particles for rank of id to calculate
	int getMyTaskNum()
	{
		return my_task[id];
	}
};
#endif


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



	void set_para(int dim, int nObj, int pos_para);

	extern int nvar,  nobj;                    //  the number of variables and objectives
	extern int position_parameters;

	//////////////////////////////////////////////////////////////////////////
#if (defined(MPI_AVAILABLE) && defined(SUITE_PARALLEL))
	extern mpi_info suite_mpi_info;
	extern double tmp_f[MAX_VL];
	extern MPI_Comm comm_range;
	void set_MPI_Comm(MPI_Comm cur_comm)
	{
		MPI_Comm_dup(cur_comm,&comm_range);
	}
	inline void mpi_initialization(int rank, int group, int n, int nObj)
	{
		suite_mpi_info.initialize_mpi(rank, group, n, nObj);
	}
	inline void mpi_finalization()
	{
		suite_mpi_info.clean();
	}
#endif

	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////// 
	void setLimits(char *pro, double *minLimit, double *maxLimit, int dim);

	void evaluate_problems(char *pro, double *xreal, double *obj, int dim, int nx, int nobj);

	//////////////////////////////////////////////////////////////////////////
	//	DTLZ
	void dtlz1(double *xreal, double *obj, int dim, int nobj);
	void dtlz2(double *xreal, double *obj, int dim, int nobj);
	void dtlz3(double *xreal, double *obj, int dim, int nobj);
	void dtlz4(double *xreal, double *obj, int dim, int nobj);
	void dtlz5(double *xreal, double *obj, int dim, int nobj);
	void dtlz6(double *xreal, double *obj, int dim, int nobj);
	void dtlz7(double *xreal, double *obj, int dim, int nobj);

	//////////////////////////////////////////////////////////////////////////
	//	UF
	void UF1(double *x, double *f, const unsigned int nx);
	void UF2(double *x, double *f, const unsigned int nx);
	void UF3(double *x, double *f, const unsigned int nx);
	void UF4(double *x, double *f, const unsigned int nx);
	void UF5(double *x, double *f, const unsigned int nx);
	void UF6(double *x, double *f, const unsigned int nx);
	void UF7(double *x, double *f, const unsigned int nx);
	void UF8(double *x, double *f, const unsigned int nx);
	void UF9(double *x, double *f, const unsigned int nx);
	void UF10(double *x, double *f, const unsigned int nx);

	//////////////////////////////////////////////////////////////////////////
	//	CF
	void CF1(double *x, double *f, double *c, const unsigned int nx);
	void CF2(double *x, double *f, double *c, const unsigned int nx);
	void CF3(double *x, double *f, double *c, const unsigned int nx);
	void CF4(double *x, double *f, double *c, const unsigned int nx);
	void CF5(double *x, double *f, double *c, const unsigned int nx);
	void CF6(double *x, double *f, double *c, const unsigned int nx);
	void CF7(double *x, double *f, double *c, const unsigned int nx);
	void CF8(double *x, double *f, double *c, const unsigned int nx);
	void CF9(double *x, double *f, double *c, const unsigned int nx);
	void CF10(double *x, double *f, double *c, const unsigned int nx);

	//////////////////////////////////////////////////////////////////////////
	//	WFG
	void wfg_eval( double* x, int n, int M, char* problem, double* fit );
	vector< double > problem_calc_fitness( const vector< double >& z, 
		const int k, const int M, char* fn );

	//////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////// 
	////////////////////////////////////////////////////////////////////////// 
	//	Example problems
	//** True if "k" in [1,z.size()), "M" >= 2, and "k" mod ("M"-1) == 0. *******
	inline bool ArgsOK( const vector< double >& z, const int k, const int M )
	{
		const int n = static_cast< int >( z.size() );

		return k >= 1 && k < n && M >= 2 && k % ( M-1 ) == 0;
	}

	//** Reduces each paramer in "z" to the domain [0,1]. ***********************
	vector< double > WFG_normalise_z( const vector< double >& z );

	//** The WFG1 problem. ******************************************************
	vector< double > WFG1(const vector< double >& z,const int k,const int M);
	//** The WFG2 problem. ******************************************************
	vector< double > WFG2(const vector< double >& z,const int k,const int M);
	//** The WFG3 problem. ******************************************************
	vector< double > WFG3(const vector< double >& z,const int k,const int M);
	//** The WFG4 problem. ******************************************************
	vector< double > WFG4(const vector< double >& z,const int k,const int M);
	//** The WFG5 problem. ******************************************************
	vector< double > WFG5(const vector< double >& z,const int k,const int M);
	//** The WFG6 problem. ******************************************************
	vector< double > WFG6(const vector< double >& z,const int k,const int M);
	//** The WFG7 problem. ******************************************************
	vector< double > WFG7(const vector< double >& z,const int k,const int M);
	//** The WFG8 problem. ******************************************************
	vector< double > WFG8(const vector< double >& z,const int k,const int M);
	//** The WFG9 problem. ******************************************************
	vector< double > WFG9(const vector< double >& z,const int k,const int M);

	//** The I1 problem. ********************************************************
	vector< double > I1(const vector< double >& z,const int k,const int M);
	//** The I2 problem. ********************************************************
	vector< double > I2(const vector< double >& z,const int k,const int M);
	//** The I3 problem. ********************************************************
	vector< double > I3(const vector< double >& z,const int k,const int M);
	//** The I4 problem. ********************************************************
	vector< double > I4(const vector< double >& z,const int k,const int M);
	//** The I5 problem. ********************************************************
	vector< double > I5(const vector< double >& z,const int k,const int M);

	//////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////// 
	////////////////////////////////////////////////////////////////////////// 
	//	Example Shapes
	//** Construct a vector of length M-1, with values "1,0,0,..." if ***********
	//** "degenerate" is true, otherwise with values "1,1,1,..." if   ***********
	//** "degenerate" is false.                                       ***********
	vector< short > WFG_create_A( const int M, const bool degenerate );

	//** Given the vector "x" (the last value of which is the sole distance ****
	//** parameter), and the shape function results in "h", calculate the   ****
	//** scaled fitness values for a WFG problem.                           ****
	vector< double > WFG_calculate_f(const vector< double >& x,
		const vector< double >& h);

	//** Given the last transition vector, get the fitness values for WFG1. *****
	vector< double > WFG1_shape( const vector< double >& t_p );
	//** Given the last transition vector, get the fitness values for WFG2. *****
	vector< double > WFG2_shape( const vector< double >& t_p );
	//** Given the last transition vector, get the fitness values for WFG3. *****
	vector< double > WFG3_shape( const vector< double >& t_p );
	//** Given the last transition vector, get the fitness values for WFG4. *****
	vector< double > WFG4_shape( const vector< double >& t_p );
	//** Given the last transition vector, get the fitness values for I1. *******
	vector< double > I1_shape( const vector< double >& t_p );

	//////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////// 
	////////////////////////////////////////////////////////////////////////// 
	//	Example Transitions
	//** Construct a vector with the elements v[head], ..., v[tail-1]. **********
	vector< double > subvector(const vector< double >& v,
	  const int head,const int tail);

	//** t1 from WFG1. **********************************************************
	vector< double > WFG1_t1( const vector< double >& y, const int k );
	//** t2 from WFG1. **********************************************************
	vector< double > WFG1_t2( const vector< double >& y, const int k );
	//** t3 from WFG1. **********************************************************
	vector< double > WFG1_t3( const vector< double >& y );
	//** t4 from WFG1. **********************************************************
	vector< double > WFG1_t4(const vector< double >& y,const int k,const int M);

	//** t2 from WFG2. **********************************************************
	vector< double > WFG2_t2( const vector< double >& y, const int k );
	//** t3 from WFG2. Effectively as per WFG4, t2. *****************************
	vector< double > WFG2_t3(const vector< double >& y,const int k,const int M);

	//** t1 from WFG4. **********************************************************
	vector< double > WFG4_t1( const vector< double >& y );

	//** t1 from WFG5. **********************************************************
	vector< double > WFG5_t1( const vector< double >& y );

	//** t2 from WFG6. **********************************************************
	vector< double > WFG6_t2(const vector< double >& y,const int k,const int M);

	//** t1 from WFG7. **********************************************************
	vector< double > WFG7_t1( const vector< double >& y, const int k );

	//** t1 from WFG8. **********************************************************
	vector< double > WFG8_t1( const vector< double >& y, const int k );

	//** t1 from WFG9. **********************************************************
	vector< double > WFG9_t1( const vector< double >& y );
	//** t2 from WFG9. **********************************************************
	vector< double > WFG9_t2( const vector< double >& y, const int k );

	//** t2 from I1. ************************************************************
	vector< double > I1_t2( const vector< double >& y, const int k );
	//** t3 from I1. ************************************************************
	vector< double > I1_t3(const vector< double >& y,const int k,const int M);

	//** t1 from I2. ************************************************************
	vector< double > I2_t1( const vector< double >& y );

	//** t1 from I3. ************************************************************
	vector< double > I3_t1( const vector< double >& y );

	//** t3 from I4. ************************************************************
	vector< double > I4_t3(const vector< double >& y,const int k,const int M);

	//////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////// 
	////////////////////////////////////////////////////////////////////////// 
	//	Framework Functions
	//** Normalise the elements of "z" to the domain [0,1]. *********************
	vector< double > normalise_z(const vector< double >& z,const vector< double >& z_max);

	//** Degenerate the values of "t_p" based on the degeneracy vector "A". *****
	vector< double > calculate_x(const vector< double >& t_p,const vector< short >& A);

	//** Calculate the fitness vector using the distance scaling constant D, ****
	//** the distance parameter in "x", the shape function values in "h",    ****
	//** and the scaling constants in "S".                                   ****
	vector< double > calculate_f(const double& D,const vector< double >& x,
	  const vector< double >& h,const vector< double >& S);

	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////// 
	//	Misc
	//** Used to correct values in [-epislon,0] to 0, and [1,epsilon] to 1. *****
	double correct_to_01( const double& a, const double& epsilon = 1.0e-10 );

	//** Returns true if all elements of "x" are in [0,1], false otherwise. *****
	bool vector_in_01( const std::vector< double >& x );

	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////// 
	//	Shape Functions
	//** True if all elements of "x" are in [0,1], and m is in [1, x.size()]. ***
	bool shape_args_ok( const vector< double >& x, const int m );

	//** The linear shape function. (m is indexed from 1.) **********************
	double linear ( const vector< double >& x, const int m );

	//** The convex shape function. (m is indexed from 1.) **********************
	double convex ( const vector< double >& x, const int m );

	//** The concave shape function. (m is indexed from 1.) *********************
	double concave( const vector< double >& x, const int m );

	//** The mixed convex/concave shape function. *******************************
	double mixed(const vector< double >& x,const int A,const double& alpha);

	//** The disconnected shape function. ***************************************
	double disc(const vector< double >& x,const int A,const double& alpha,
		const double& beta);

	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////// 
	//	TransFunctions
	//** Calculate the minimum of two doubles. **********************************
	inline double min( const double& a, const double& b )
	{
		double m = a > b ? b : a;
		return m;// std::min< const double >(a, b);
	}

	//** The polynomial bias transformation function. ***************************
	double b_poly( const double& y, const double& alpha );

	//** The flat region bias transformation function. **************************
	double b_flat(const double& y,const double& A,const double& B,const double& C);

	//** The parameter dependent bias transformation function. ******************
	double b_param(const double& y,const double& u,
	  const double& A,const double& B,const double& C);

	//** The linear shift transformation function. ******************************
	double s_linear( const double& y, const double& A );

	//** The deceptive shift transformation function. ***************************
	double s_decept(const double& y,const double& A,const double& B,const double& C);

	//** The multi-modal shift transformation function. *************************
	double s_multi(const double& y,const int A,const double& B,const double& C);

	//** The weighted sum reduction transformation function. ********************
	double r_sum(const vector< double >& y,const vector< double >& w);

	//** The non-separable reduction transformation function. *******************
	double r_nonsep( const vector< double >& y, const int A );
}

#endif	//	_EMO_TEST_SUITE_