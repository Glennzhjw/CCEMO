#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

//#include "ZDTfamily.h"
//#include "dtlz.h"
// #include "cec09.h"
//#include "wfg_interface.h"
// #include "EMO_test_suite.h"

// individual
class TIndividual{
public:
	TIndividual();
	virtual ~TIndividual();

	void var_init();							// init decision variables
	void obj_eval();							// evaluate objectives
	bool operator<(const TIndividual &ind2);	// dominate relation
	bool operator==(const TIndividual &ind2);
	void operator=(const TIndividual &ind2);

	vector<double> x_var;	// decision variables
	vector<double> y_obj;	// objective variables

	double fm;				// mutation factor
	double cr;				// crossover probability

	int		rank;			// fast-non-dominated-sorting rank
	double	crowdist;		// crowding distance

	int		subpop;			// store the subpop id for current individual
	bool	isoffspring;	// decide whether current individual is a offspring

};

TIndividual::TIndividual()
{
	for (int i = 0; i<numVariables; ++i)
		x_var.push_back(0.0);
	for (int i = 0; i<numObjectives; ++i)
		y_obj.push_back(0.0);
}

TIndividual::~TIndividual()
{
}

// 初始化个体的决策向量
void TIndividual::var_init()
{
	for (int i = 0; i<numVariables; ++i)
		x_var[i] = lowBound[i] + rnd_uni(&rnd_uni_init)*(uppBound[i] - lowBound[i]);
}

// 适应值评估
void TIndividual::obj_eval()
{
	//objectives(x_var,y_obj);
	double *x = new double[numVariables];
	double *y = new double[numObjectives];
	for (int i = 0; i<numVariables; ++i)
		x[i] = x_var[i];

	//ZDTFamily(prob, x, y, numVariables);
	//dtlz_problems(prob, x,y,numVariables,numObjectives);
// 	uf_problems(prob, x, y, numVariables);
	//wfg_eval(x, numVariables, 8, numObjectives, prob, y);
	EMO_TEST_SUITE::evaluate_problems(prob,x,y,numVariables,1,numObjectives);
	fes++;
	for (int i = 0; i<numObjectives; ++i)
		y_obj[i] = y[i];
	delete[]x;
	delete[]y;
}

bool TIndividual::operator<(const TIndividual &ind2)
{
	for (int i = 0; i<numObjectives; i++)
	if (y_obj[i] > ind2.y_obj[i])
		return false;
	if (y_obj == ind2.y_obj)
		return false;
	return true;
}

bool TIndividual::operator==(const TIndividual &ind2)
{
	if (ind2.x_var == x_var && ind2.y_obj == y_obj)
		return true;
	else
		return false;
}

void TIndividual::operator=(const TIndividual &ind2)
{
	x_var = ind2.x_var;
	y_obj = ind2.y_obj;
	fm	  = ind2.fm;
	cr	  = ind2.cr;

	rank = ind2.rank;
	crowdist = ind2.crowdist;
	/* the isoffspring is set for deciding whether current individual is a successful solution */
	subpop = ind2.subpop;
	isoffspring = ind2.isoffspring;
}



// sub population
class TSPOP{
public:
	TSPOP();
	virtual ~TSPOP();

	void subpop_init();
	void globalbest(int obj);
	void operator=(const TSPOP &sub2);

	vector<TIndividual> indiv;
	TIndividual bestIndiv;			// global best individual in current subpop

	vector<double> utility;
	double psubpop;			// the probability of selecting current subpop
};

TSPOP::TSPOP()
{
	TIndividual t;
	for (int i = 0; i<subpopSize; ++i)
		indiv.push_back(t);
}

TSPOP::~TSPOP()
{
}

// init current subpop
void TSPOP::subpop_init()
{
	for (int i = 0; i<subpopSize; ++i)
	{
		indiv[i].var_init();
		indiv[i].obj_eval();
		indiv[i].fm = 0.5;
		indiv[i].cr = 0.9;
		indiv[i].isoffspring = false;
	}
}

// look for the global best individual in current subpop
void TSPOP::globalbest(int obj)
{
	bestIndiv = indiv[0];
	for (int i = 1; i < subpopSize; ++i)
	{
		if (indiv[i].y_obj[obj] < bestIndiv.y_obj[obj])
			bestIndiv = indiv[i];
	}
}

void TSPOP::operator=(const TSPOP &sub2)
{
	indiv = sub2.indiv;
	bestIndiv = sub2.bestIndiv;
	utility = sub2.utility;
	psubpop = sub2.psubpop;
}

#endif
