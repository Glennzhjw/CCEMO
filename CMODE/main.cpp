/*
此程序实现CMDE算法的DE/current-to-best/1变异策略，即在分种群中采用 DE/current-to-best/1变异策略.
与 main(DE current-to-pbest 1).cpp中的DE/current-to-pbest/1变异策略相比，相对较好。
越界处理都采取移动步长的方法效果较好
*/

#include "EMO_test_suite.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

#include <fstream>

using namespace std;

int subpopSize = 20;		// size of every subpop
int archiveSize = 100;	// limited size of archive
// demensionality of variables and objectives
/* for dtlz1 nvar=nobj+4;
for dtlz2-6 nvar=nobj+9;
for dtlz7 nvar=nobj+19;
for uf   nvar = 30
*/
int numObjectives = 2;		// number of objectives
int numVariables  = 30;		// number of decision variables
int runtimes =1;			// run times
int max_fun = 300000;		// the maximal function evaluations

FILE *file;				// write archive solutions
FILE *filet;			// write run times
FILE *filex;
/*
FILE *file1;
FILE *file2;
FILE *file3;
FILE *file4;
FILE *file5;
FILE *file6;
FILE *file7;
FILE *file8;
FILE *file9;
char outfile1[20];
char outfile2[20];
char outfile3[20];
char outfile4[20];
char outfile5[20];
char outfile6[20];
char outfile7[20];
char outfile8[20];
char outfile9[20];
*/

char prob[30];				// test problem name
char outfile[256];		// output arhive name
char timefile[256];		// output time name
char xfile[256];		// output variables

double *lowBound;		// low bound of decision variable
double *uppBound;		// upper bound of decision variable
int  seed = 237;		// seed for random number generation
long rnd_uni_init;		// parameters for random number generation
int	 fes;				// number of function evaluations 

double cc = 1.0 / 10;	//	1/c ={1,2,5,10,20,50,100}
vector<double> mu_cr(numObjectives);							// init average value of normal distribution 
vector<double> mu_f(numObjectives);								// init average value of cauchy distribution 
vector<vector<double> > SF(numObjectives);						// store F of successful mutation 
vector<vector<double> > SCR(numObjectives); 					// store CR of successful crossover

double *fun_min;
double *fun_max;

#include "ran.h"
#include "common.h"
#include "individual.h"
// #include "recombination.h"
#include "rcmode.h"

//主程序
int main(int argc, char *argv[])
{
// 	if (argc < 2)
// 	{
// 		printf("error parameters.");
// 		exit(1);
// 	}
	
	int fun_start_num = 1; // the start number of test function, you can find the explanation in testInstance.txt
	int fun_end_num = 1;    // the end number of test function

	int seq;
	std::ifstream readf("testInstance.txt");
	int iPro;
	for(iPro=1;iPro<fun_start_num;iPro++)
	{
		readf>>seq;
		readf>>prob;
		readf>>numVariables;
		readf>>numObjectives;
	}
	for(iPro=fun_start_num;iPro<=fun_end_num;iPro++)
	{
		readf>>seq;
		readf>>prob;
		readf>>numVariables;
		readf>>numObjectives;

		max_fun=(int)(1e4*numVariables);

// 		strcpy_s(outfile, 20, prob);
		//strcat_s(outfile, ".5D");
// 		fopen_s(&file, outfile, "w");
		/*
		strcpy_s(timefile, 20, outfile);
		strcat_s(timefile, " time");
		fopen_s(&filet, timefile, "w");

		strcpy_s(outfile1, 20, "WFG9.10D 40");
		strcpy_s(outfile2, 20, "WFG9.10D 60");
		strcpy_s(outfile3, 20, "WFG9.10D 80");
		strcpy_s(outfile4, 20, "WFG9.10D 100");
		strcpy_s(outfile5, 20, "WFG3.10D 120");
		strcpy_s(outfile6, 20, "WFG3.10D 140");
		strcpy_s(outfile7, 20, "WFG3.10D 160");
		strcpy_s(outfile8, 20, "WFG3.10D 180");
		strcpy_s(outfile9, 20, "WFG3.10D 200");
		fopen_s(&file1, outfile1, "w");
		fopen_s(&file2, outfile2, "w");
		fopen_s(&file3, outfile3, "w");
		fopen_s(&file4, outfile4, "w");
		fopen_s(&file5, outfile5, "w");
		fopen_s(&file6, outfile6, "w");
		fopen_s(&file7, outfile7, "w");
		fopen_s(&file8, outfile8, "w");
		fopen_s(&file9, outfile9, "w");
		*/
		fun_min = new double[numObjectives];
		fun_max = new double[numObjectives];
		lowBound = new double[numVariables];
		uppBound = new double[numVariables];
		/* low and upper boundary of UF11 and UF12 , while [0, 2i] is boundary of UF13*/	
		//double lb[30] = { -1.773, -1.846, -1.053, -2.370, -1.603, -1.878, -1.677, -0.935, -1.891, -0.964, -0.885, -1.690, -2.235, -1.541, -0.720, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000 };
		//double ub[30] = { 1.403, 1.562, 2.009, 0.976, 1.490, 1.334, 1.074, 2.354, 1.462, 2.372, 2.267, 1.309, 0.842, 1.665, 2.476, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 };
// 		lowBound[0] = 0;
// 		uppBound[0] = 1;
// 		for (int i = 1; i<numVariables; ++i)
// 		{
// 			lowBound[i] = -1;// lb[i];
// 			uppBound[i] = 1;// ub[i];
// 		}
		
		EMO_TEST_SUITE::setLimits(prob,lowBound,uppBound,numVariables);

		double startTime, endTime;
		double avgtime = 0;
		sprintf(timefile,"results/TIME_%s_%d",prob,numVariables);
		filet=fopen(timefile,"a");

		printf("\n\n-- PROBLEM %s\n",prob);
		printf("--  variables %d \n",numVariables);
		printf("--  objectives %d \n\n",numObjectives);

		for (int run = 1; run <= runtimes; ++run)
		{
			sprintf(outfile,"results/FUN_%s_%d_%d",prob,numVariables,run);
			file=fopen(outfile,"w");
			sprintf(xfile,"results/VAR_%s_%d_%d",prob,numVariables,run);
			filex=fopen(xfile,"w");

			// init mu_cr and mu_f in JADE
			for(int i=0; i<mu_cr.size(); ++i)
			{
				mu_cr[i] = 0.5;
				mu_f[i]  = 0.5;
			}
			fes = 0;	// reset 0 every run
			startTime = clock();
			seed = (seed + 111) % 1235;
			rnd_uni_init = -(long)seed;

			printf("--   RUN %d   --\n",run);
			TRCMODE RCMODE;
			RCMODE.run(max_fun, run);

			endTime = clock();
			double aaa = (endTime - startTime) / CLOCKS_PER_SEC;
			avgtime += aaa;
			fprintf(filet, "%e\n", aaa);
			printf("run%d time is %lf s\n", run, aaa);
			fclose(file);
			fclose(filex);
		}
		//fprintf(filet, "average run time is %lfs\n", avgtime / runtimes);
		delete[]fun_min;
		delete[]fun_max;
		delete[]lowBound;
		delete[]uppBound;
		fclose(filet);
		/*
		fclose(file1);
		fclose(file2);
		fclose(file3);
		fclose(file4);
		fclose(file5);	fclose(file6);	fclose(file7);	fclose(file8);	fclose(file9);*/
	}

}
