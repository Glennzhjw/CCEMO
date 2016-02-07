//  [1/26/2016 John]
#include "list.h"
#include "rand.h"
#include "EMO_test_suite.h"
#include <string.h>

# define INF 1.0e14
# define MAX_SIZE 65536
//global variables and functions

//---------------------------------------------------------
//variables

extern char testInstance[];//name of test problem
extern int nObj;//number of objective problems
extern int nSwm;//number of swarm (nSwm==nObj)
extern int nDim;//dimension of variable
extern int nPop;//number of population
extern int nArch;//size of archive
extern int cnArch;//number of solutions in archive
extern int nRep;//number of solutions collected for selection
extern int nRep_tmp;
extern int nonDominateSize;//

extern double* minLimit;//min limit of variable
extern double* maxLimit;//max limit of variable

//dynamic grouping
extern int S[];//{2,5,10,50,100,250};//sizes of group
extern int* Indexes;//permed index
extern int  nGroup;//number of group
extern int  dimIn1Group;//size of 1 group
extern int S_SIZE;//number of elements in S[]

//dynamic grouping for arch
extern int* archIndex;
extern int arch_nGroup;
extern int arch_dimIn1Group;
extern int arch_S[];
extern int arch_S_SIZE;

//evolution
extern double* xCurrent;
extern double* uTrail;
extern double* cTrail;
extern double* xFitness;
extern double* cFitness;
extern double* xBest;
extern double* xBFitness;
extern double* cF;
extern double* cCR;

//archive
extern double* archive;
extern double* archFit;
extern double* archiveOld;
extern double* archFitOld;
extern int cnArchOld;
extern double* uArchive;
extern double* cArchive;
extern double* cArchFit;
extern int* archiveClass;
extern double* uF;
extern double* uCR;

//collection for selection
extern double* repertory;
extern double* repertFit;
extern double* repertoryDensity;
extern double* repertoryF;
extern double* repertoryCR;

extern double* fun_max;
extern double* fun_min;

//adaptive DE
extern double F,CR;
extern double *S_F,*S_CR;
extern int *Sflag;
extern double F_mu,CR_mu;
extern double c_para;
extern double t1,t2;
extern double *arch_F,*arch_CR;
extern int iter,maxIteration;
extern int iter_each;
extern int iter_sum;

extern int* dimFlag;

//output file pointer
extern FILE *fptobj;//objective value
extern FILE *fptvar;//solution variable
extern FILE *fpttime;//operation time

//	MPI
extern int mpi_size;
extern int mpi_rank;
extern char my_name[];
extern int name_len;
extern int task_l;
extern int task_r;
extern int task_num;
extern MPI_Comm my_comm;
extern int mpi_size_self;
extern int mpi_rank_self;
extern int mpi_color;
extern MPI_Comm comm_masters;
extern int mpi_size_master;
extern int mpi_rank_master;
extern int master_flag;
extern int mpi_rank_master_archive;

extern int* recv_size;
extern int* disp_size;
extern int* each_size;


//---------------------------------------------------------
//functions

//	control
void set_parameter(int npop, int ndim, int nobj, int narch);
void run(char* func_name,int iRun);
void update_objective();
void update_archive();

//	differentialEvolution
void DE_gen_uTrail();
void DE_gen_uTrail_ada();
void DE_gen_archive();
void DE_gen_archive2();
void DE_gen_archive3();
int  isAllTheSame(int size,int iDim);
void update_F_CR_mu();
void generate_F_CR();
void generate_arch_F_CR();

//	cooperativeCoevolution
void update_xBest_initial();
void update_xCurrent();
void update_xCurrent_one(int iPop);
void update_xBest_archive();
void cooperativeCoevolution();
void cooperativeCoevolution_ada();
void cooperativeCoevolution_archive();
void mainLoop();
void joinSolutions(int iPop);
void mainLoop2();
void joinSolutions2();
void mainLoop_archive();
void joinSolutions_archive(int iArch);

//	handleMemory
void allocateMemory();
double* allocDouble(int size);
int* allocInt(int size);
void freeMemory();

//	index
void initializeIndex();
void permIndexes();
void perm_archIndex();

//	initialization
void initializeProblem();
void initializePopulation();
void initializePopulation_sinusMap();
void initializePopulation_UD();
void ParameterSet();
void evaluatePopInitial();

//	nondominance
void refineRepertory_generateArchive();
void refineRepertory_generateArchive_SDE();
void K_Neighbor_Nearest_SDE(int count, int frontSize, list *elite);
double* generateDistMatrix(double* f,int non_size,int* non_indexes);
double Euclid_Dist(double* vec1, double* vec2);
void sort_dist_index(double *a, int *b, int left, int right);
bool cmp_crowd(double *a, double *b);

//	dominanceComparator
int dominanceComparator(double* obj1, double* obj2);
int check_dominance(double* obj1, double* obj2);

//	distance
void fillCrowdingDistance(int count, int frontSize, list *elite);
void assignCrowdingDistanceList(list *lst, int frontSize);
void assignCrowdingDistanceIndexes(int c1, int c2);
void assignCrowdingDistance(int *distance, int **arrayFx, int frontSize);

//	sort
void quickSortDistance(int *distance, int frontSize);
void qSortDistance(int *distance, int left, int right);
void quickSortFrontObj(int objcount, int arrayFx[], int sizeArrayFx);
void qSortFrontObj(int objcount, int arrayFx[], int left, int right);

//	display
void showPopulation();
void showArchive();
void showRepertory();
void show_uArchive();
void showLimits();
void showGlobalBest();

//	utility
void copyToArchiveFromRepertory(int iA, int iR);
int INDEX(int a, int b, int c);
int INDEX_f(int a, int b, int c);
void selectSamples(int pp,int candidate,int *r1,int *r2,int *r3);
void saveArchiveOld();
void addArchiveOld2repertory();
bool isDuplicate(double* s, int i, int j);
void collect2repertory_select();
void get_nonDominateSize();
void save_obj(FILE *fpt);
void save_var(FILE *fpt);

//	MPI
void setMPI();
void update_recv_disp(int* num, int n, int l);
void freeMPI();
void collection_initial();
void SynchronizeBest();
void collect2repertory0();
void SynchronizeArchive();
void SynchronizeArchive_one();
void collect2master_archive();
void update_iteration();
void allocateMPI();
void allocateTask(int n);