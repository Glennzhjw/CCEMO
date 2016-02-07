#include "global.h"

#define MAX_STR 256

int S[MAX_STR];//{2,5,10,50,100,250};
char testInstance[MAX_STR];

int nObj,nSwm,nDim,nPop;
int nArch,cnArch;
int nRep;
int nRep_tmp;
int nonDominateSize;

double* minLimit;
double* maxLimit;

int* Indexes;
int  nGroup;
int  dimIn1Group;
int S_SIZE;
// arch
int* archIndex;
int arch_nGroup;
int arch_dimIn1Group;
int arch_S[MAX_STR];
int arch_S_SIZE;

double* xCurrent;
double* uTrail;
double* cTrail;
double* xFitness;
double* cFitness;
double* xBest;
double* xBFitness;
double* cF;
double* cCR;

double* archive;
double* archFit;
double* archiveOld;
double* archFitOld;
int cnArchOld;
double* uArchive;
double* cArchive;
double* cArchFit;
int* archiveClass;
double* uF;
double* uCR;

double* repertory;
double* repertFit;
double* repertoryDensity;
double* repertoryF;
double* repertoryCR;

double* fun_max;
double* fun_min;

double F,CR;
double *S_F,*S_CR;
int *Sflag;
double F_mu,CR_mu;
double c_para;
double t1,t2;
double *arch_F,*arch_CR;
int iter,maxIteration;
int iter_each;
int iter_sum;

int* dimFlag;

FILE *fptobj;
FILE *fptvar;
FILE *fpttime;

//	MPI
int mpi_size;
int mpi_rank;
char my_name[MAX_STR];
int name_len;
int task_l;
int task_r;
int task_num;
MPI_Comm my_comm;
int mpi_size_self;
int mpi_rank_self;
int mpi_color;
MPI_Comm comm_masters;
int mpi_size_master;
int mpi_rank_master;
int master_flag;
int mpi_rank_master_archive;

int* recv_size;
int* disp_size;
int* each_size;