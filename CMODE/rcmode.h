#ifndef RCMODE_H
#define RCMODE_H


class TRCMODE{
public:
	TRCMODE();
	virtual ~TRCMODE();

	void init_population();
	void init_archive();
	void cal_utility(int gen);
	void archive_evolution_de(vector <TIndividual> &mix);	// archive evolution in DE/rand/1
	void evolution(int gen);								// mating restriction, recombination, mutation, update
	void run(int mg, int rn);
	void report(vector <TIndividual> archive, FILE *file);
	void report_x(vector <TIndividual> archive, FILE *file);


	int lpdselect();									// select a subpop based on the probability
	int check_dominance(TIndividual a, TIndividual b);	// compare Pareto dominance relation between individual a and b
	void my_sort(vector <TIndividual> &s, int k);		// sort individuals by k'th objective
	void crowding_distance_sorting(vector <TIndividual> arch, int select_size);	// select the last rank solutions based on crowding distance
	void K_Neighbor_Nearest_SDE(vector <TIndividual> arch, int select_size);	// select the last rank solutions based on KNN
	void archive_update();						// fast non-dominated sorting and select the first pops solutions

	vector<TSPOP> population;
	vector<TIndividual> archive;
	vector<TIndividual> child_population;
};

TRCMODE::TRCMODE()
{
	TSPOP tp;
	for (int i = 0; i<numObjectives; ++i)
		population.push_back(tp);
}

TRCMODE::~TRCMODE()
{
}

void TRCMODE::init_population()
{
	if(!child_population.empty())
		child_population.clear();
	for (int i = 0; i<population.size(); ++i)
	{
		population[i].subpop_init();
		for(int j=0; j<population[i].indiv.size(); ++j)
			child_population.push_back(population[i].indiv[j]);
		//population[i].psubpop = 1.0;	// initialize the selecting probability of every subpop
	}
}

void TRCMODE::init_archive()
{
	if(!archive.empty())
		archive.clear();
}

void TRCMODE::report(vector <TIndividual> archive, FILE *file)
{
	fprintf(file, "%d\n", archive.size());
	for (int n = 0; n<archive.size(); n++)
	{
		for (int k = 0; k<numObjectives; k++)
			fprintf(file, "%e\t", archive[n].y_obj[k]);
		fprintf(file, "\n");
	}
	return;
}

void TRCMODE::report_x(vector<TIndividual> archive, FILE *filex)
{
	fprintf(filex, "%d\n", archive.size());
	for (int n = 0; n<archive.size(); n++)
	{
		for (int k = 0; k<numVariables; k++)
			fprintf(filex, "%e\t", archive[n].x_var[k]);
		fprintf(filex, "\n");
	}
	return;
}

// density in objective space

void  TRCMODE::cal_utility(int gen)
{
	// calculate density in objective space
	double *uti = new double[numObjectives];
	for (int i = 0; i<population.size(); ++i)
	{
		double *min_uti = new double[numObjectives];
		double *max_uti = new double[numObjectives];
		for (int k = 0; k<numObjectives; ++k)
		{
			min_uti[k] = population[i].indiv[0].y_obj[k];
			max_uti[k] = population[i].indiv[0].y_obj[k];
		}
		for (int j = 1; j<population[i].indiv.size(); ++j)
		{
			for (int k = 0; k<numObjectives; ++k)
			{
				if (population[i].indiv[j].y_obj[k] < min_uti[k])
					min_uti[k] = population[i].indiv[j].y_obj[k];
				if (population[i].indiv[j].y_obj[k] > max_uti[k])
					max_uti[k] = population[i].indiv[j].y_obj[k];
			}
		}
		uti[i] = 1;
		for (int k = 0; k<numObjectives; ++k)
			uti[i] *= (max_uti[k] - min_uti[k]);
		delete[]min_uti;
		delete[]max_uti;
	}
	// store the utility
	if (gen < 20)
	{
		for (int i = 0; i < population.size(); ++i)
			population[i].utility.push_back(uti[i]);
	}
	// recalculate the selecting probability of every subpop and store the newest utility
	else
	{
		double *delta_uti = new double[numObjectives];
		double max_delta = -100;
		for (int i = 0; i<population.size(); ++i)
		{
			delta_uti[i] = (population[i].utility[0] - uti[i]) / population[i].utility[0];
			if (delta_uti[i] > max_delta)
				max_delta = delta_uti[i];
			population[i].utility.erase(population[i].utility.begin());
			population[i].utility.push_back(uti[i]);
		}
		if (max_delta <= 0)
		{
			for (int i = 0; i < population.size(); ++i)
				population[i].psubpop = 1.0;
		}
		else
		{
			for (int i = 0; i<population.size(); ++i)
				population[i].psubpop = delta_uti[i] / max_delta;
		}
		delete[]delta_uti;
	}
	delete[]uti;
}

// density in decision space
/*
void  TRCMODE::cal_utility(int gen)
{
	// calculate density in objective space
	double *uti = new double[numVariables];
	for (int i = 0; i<population.size(); ++i)
	{
		double *min_uti = new double[numVariables];
		double *max_uti = new double[numVariables];
		for (int k = 0; k<numVariables; ++k)
		{
			min_uti[k] = population[i].indiv[0].x_var[k];
			max_uti[k] = population[i].indiv[0].x_var[k];
		}
		for (int j = 1; j<population[i].indiv.size(); ++j)
		{
			for (int k = 0; k<numVariables; ++k)
			{
				if (population[i].indiv[j].x_var[k] < min_uti[k])
					min_uti[k] = population[i].indiv[j].x_var[k];
				if (population[i].indiv[j].x_var[k] > max_uti[k])
					max_uti[k] = population[i].indiv[j].x_var[k];
			}
		}
		uti[i] = 1;
		for (int k = 0; k<numVariables; ++k)
			uti[i] *= (max_uti[k] - min_uti[k]);
		delete[]min_uti;
		delete[]max_uti;
	}
	// store the utility
	if (gen < 20)
	{
		for (int i = 0; i < population.size(); ++i)
			population[i].utility.push_back(uti[i]);
	}
	// recalculate the selecting probability of every subpop and store the newest utility
	else
	{
		double *delta_uti = new double[numObjectives];
		double max_delta = -1;
		for (int i = 0; i<population.size(); ++i)
		{
			delta_uti[i] = fabs(population[i].utility[0] - uti[i]);
			if (delta_uti[i] > max_delta)
				max_delta = delta_uti[i];
			population[i].utility.erase(population[i].utility.begin());
			population[i].utility.push_back(uti[i]);
		}
		if (max_delta == 0)
		{
			for (int i = 0; i < population.size(); ++i)
				population[i].psubpop = 1.0;
		}
		else
		{
			for (int i = 0; i<population.size(); ++i)
				population[i].psubpop = delta_uti[i] / max_delta;
		}
		delete[]delta_uti;
	}
	delete[]uti;
}
*/

int TRCMODE::check_dominance(TIndividual a, TIndividual b)
{
	int flag1 = 0;
	int flag2 = 0;
	for (int i = 0; i<numObjectives; ++i)
	{
		if (a.y_obj[i] < b.y_obj[i])
			flag1 = 1;
		else if (a.y_obj[i] > b.y_obj[i])
			flag2 = 1;
	}
	if (flag1 == 1 && flag2 == 0)
		return 1;
	else if (flag1 == 0 && flag2 == 1)
		return -1;
	else
		return 0;
}

void TRCMODE::my_sort(vector <TIndividual> &s, int k)
{
	for (int i = 0; i<s.size(); ++i)
	{
		for (int j = i + 1; j<s.size(); ++j)
		{
			if (s[j].y_obj[k] < s[i].y_obj[k])
			{
				TIndividual t = s[i];
				s[i] = s[j];
				s[j] = t;
			}
		}
	}
}

void TRCMODE::crowding_distance_sorting(vector <TIndividual> arch, int select_size)
{
	int archsize = arch.size();
	for (int i = 0; i<archsize; ++i)	// initialize crowding distance
		arch[i].crowdist = 0;
	double fmin;
	double fmax;
	for (int k = 0; k<numObjectives; ++k)
	{
		my_sort(arch, k);
		fmin = arch[0].y_obj[k];
		fmax = arch[archsize - 1].y_obj[k];

		arch[0].crowdist += 1000000.0;
		arch[archsize - 1].crowdist += 1000000.0;
		for (int i = 1; i<archsize - 1; ++i)
			arch[i].crowdist += (arch[i + 1].y_obj[k] - arch[i - 1].y_obj[k]) / (fmax - fmin);
	}
	// select select_size solutions by crowding distance
	for (int i = 0; i<select_size; ++i)
	{
		for (int j = i + 1; j<archsize; ++j)
		{
			if (arch[j].crowdist > arch[i].crowdist)
			{
				TIndividual t = arch[j];
				arch[j] = arch[i];
				arch[i] = t;
			}
		}
		archive.push_back(arch[i]);
	}
}

void TRCMODE::K_Neighbor_Nearest_SDE(vector <TIndividual> arch, int select_size)
{
	int non_size = arch.size();
	
	//查找并记录每个目标上的最大值和最小值
	for (int i = 0; i<numObjectives; ++i)
	{
		fun_min[i] = INF;
		fun_max[i] = -100;
	}
	for (int i = 0; i<non_size; ++i)
	{
		for (int j = 0; j<numObjectives; ++j)
		{
			if (arch[i].y_obj[j] < fun_min[j])
				fun_min[j] = arch[i].y_obj[j];
			if (arch[i].y_obj[j] > fun_max[j])
				fun_max[j] = arch[i].y_obj[j];
		}
	}	

	//计算各个非占优解，经过shift-based的欧几里得距离矩阵  
	double **dist;				//记录距离矩阵 
	int **distIndex;			//记录排序后的距离矩阵对应的非占优解下标 
	int *flag;					//标记已被删除的解为0，未被删除的解为1 
	dist = (double**)malloc(non_size*sizeof(double*));
	distIndex = (int**)malloc(non_size*sizeof(int*));
	flag = (int*)malloc(non_size*sizeof(int));
	for (int i = 0; i<non_size; ++i)
	{
		flag[i] = 1;
		dist[i] = (double*)malloc(non_size*sizeof(double));
		distIndex[i] = (int*)malloc(non_size*sizeof(int));
		for (int j = 0; j<non_size; ++j)
		{
			if (i == j)
				dist[i][j] = INF;
			else
				dist[i][j] = Euclid_Dist(arch[i].y_obj, arch[j].y_obj);
			distIndex[i][j] = j;
		}
	}

	//对距离矩阵进行带下标的排序 
	for (int i = 0; i<non_size; ++i)
		sort_dist_index(dist[i], distIndex[i], 0, non_size - 1);

	//迭代剔除non_size-NA个拥挤非占优解 
	int cur_size = non_size;
	while (cur_size>select_size)
	{
		//删除邻近距离最小的存档解操作，通过对该解进行已处理标记来实现 
		int crowd_index;						//记录邻近距离最小的存档解下标 
		for (int i = 0; i<non_size; ++i)
		{
			if (flag[i])
			{
				crowd_index = i;
				break;
			}
		}
		for (int i = crowd_index + 1; i < non_size; ++i)
		{
			if (flag[i] && cmp_crowd(dist[i], dist[crowd_index]))
				crowd_index = i;
		}
		flag[crowd_index] = 0;

		//修改邻近距离矩阵，找到删除解下标，然后将邻近距离一致前移 
		for (int i = 0; i<non_size; ++i)
		{
			if (flag[i] == 0) continue;
			for (int j = 0; j<cur_size; ++j)
			{
				if (distIndex[i][j] == crowd_index)
				{
					while (j<cur_size - 1)
					{
						dist[i][j] = dist[i][j + 1];
						distIndex[i][j] = distIndex[i][j + 1];
						j++;
					}
					dist[i][j] = INF;				//此时j=cur_size-1 
					distIndex[i][j] = crowd_index;
					break;
				}
			}
		}
		cur_size -= 1;
	}
	//迭代剔除结束 

	//保存经过剔除操作后的存档解到A中 
	for (int i = 0; i < non_size; ++i)
	{
		if (flag[i])
			archive.push_back(arch[i]);
	}
	for (int i = 0; i<non_size; ++i)
	{
		free(dist[i]);
		free(distIndex[i]);
	}
	free(dist);
	free(distIndex);
	free(flag);
}

void TRCMODE::archive_evolution_de(vector <TIndividual> &mix)
{
	double tao = 0.1;	// {0.05, 0.1, 0.15, 0.2, 0.25, 0.3}
	if (archive.size() < 4)
		return;
	for (int r = 0; r<archive.size(); ++r)
	{
		if (rnd_uni(&rnd_uni_init) < tao)
			archive[r].fm = 0.1 + 0.9*rnd_uni(&rnd_uni_init);
		if (rnd_uni(&rnd_uni_init) < tao)
			archive[r].cr = rnd_uni(&rnd_uni_init);
		int r1, r2, r3;
		SelectSamples(archive.size(), r, &r1, &r2, &r3);				// select three mutually exclusive individuals with i
		// mutation and crossover
		TIndividual child;
		int idx_rnd = int(rnd_uni(&rnd_uni_init)*numVariables);	// randomly select a dimension
		for (int j = 0; j<numVariables; j++)
		{
			double rnd = rnd_uni(&rnd_uni_init);
			if (rnd <= archive[r].cr || j == idx_rnd)				//------------- here for setting of crossover rate --------
				child.x_var[j] = archive[r1].x_var[j] + archive[r].fm*(archive[r2].x_var[j] - archive[r3].x_var[j]);
			else
				child.x_var[j] = archive[r].x_var[j];
			if (child.x_var[j] < lowBound[j])
				child.x_var[j] = lowBound[j];
			if (child.x_var[j] > uppBound[j])
				child.x_var[j] = uppBound[j];
			/*
			if(child.x_var[j] < lowBound[j])
				while(child.x_var[j] < lowBound[j])
					child.x_var[j] += uppBound[j] - lowBound[j];
			if(child.x_var[j] > uppBound[j])
				while(child.x_var[j] > uppBound[j])
					child.x_var[j] -= uppBound[j] - lowBound[j];
			*/
		}
		child.subpop = numObjectives;		// indicate it is an archive individual
		child.isoffspring = true;
		// inhert fm and cr from parent
		child.fm = archive[r].fm;
		child.cr = archive[r].cr;
		child.obj_eval();
		mix.push_back(child);
	}
}

void TRCMODE::archive_update()
{
	// if more than archiveSize
	vector <TIndividual> mix;
	if (!mix.empty())
		mix.clear();
	for (int i = 0; i<child_population.size(); ++i)	// put offsprings individuals into mix
		mix.push_back(child_population[i]);
	for (int i = 0; i<archive.size(); ++i)			// put initial archive individuals into mix
		mix.push_back(archive[i]);
	archive_evolution_de(mix);						// put archive evolutinary individuals into mix
	vector <TIndividual> arch;
	if (!arch.empty())
		arch.clear();
	bool ff;
	for (int i = 0; i < mix.size(); ++i)
	{
		ff = true;
		for (int j = 0; j < mix.size(); ++j)
		{
			if (i == j)
				continue;
			if (check_dominance(mix[j], mix[i]) == 1)
			{
				ff = false;
				break;
			}
		}
		if (ff)
			arch.push_back(mix[i]);
	}

	// clear archive, and update archive in the next generation
	if (!archive.empty())
		archive.clear();
	/*--------------------------update archive begin----------------------------*/
	int arch_size = arch.size();
	if (arch_size <= archiveSize)		// put all individuals in mix into archive, no need for fast_rank_sorting
	{
		for (int i = 0; i < arch_size; ++i)
			archive.push_back(arch[i]);
	}
	else
		K_Neighbor_Nearest_SDE(arch, archiveSize); //crowding_distance_sorting(arch, select_size);
	/*--------------------------update archive end----------------------------*/
}

int TRCMODE::lpdselect()
{
	double *sump = new double[numObjectives + 1];
	sump[0] = 0;
	for (int i = 0; i<numObjectives; ++i)
		sump[i + 1] = sump[i] + population[i].psubpop;
	double ra = sump[numObjectives] * rnd_uni(&rnd_uni_init);
	int in = 0;
	while (in <= numObjectives && ra >= sump[in])
		in++;
	if (in > numObjectives)
		in = numObjectives;
	//printf("%f %f %d\n ", sump[pops], ra, in);
	delete[]sump;
	return (in - 1);
}

// recombination, mutation, update in DE/rand/1
void TRCMODE::evolution(int gen)
{
	// new solution generation
	if (!child_population.empty())
		child_population.clear();
	for (int n = 0; n<numObjectives; ++n)
	{
		//if (rnd_uni(&rnd_uni_init) > population[n].psubpop)
			//continue;
		population[n].globalbest(n);
		SF[n].clear();
		SCR[n].clear();
		TSPOP child;
		for (int r = 0; r<subpopSize; r++)
		{
			// adaptation for F and CR
			population[n].indiv[r].cr = gaussrand(mu_cr[n], 0.1);
			while (population[n].indiv[r].cr < 0 || population[n].indiv[r].cr > 1)
				population[n].indiv[r].cr = gaussrand(mu_cr[n], 0.1);

			population[n].indiv[r].fm = cauchyrand(mu_f[n], 0.1);
			while (population[n].indiv[r].fm <= 0)
				population[n].indiv[r].fm = cauchyrand(mu_f[n], 0.1);
			if(population[n].indiv[r].fm > 1)		// must be put after while loop
				population[n].indiv[r].fm = 1;
			/* -------------begin DE/current-to-best/1 evolution-----------------*/
			int r1, r2;
			SelectSamples(subpopSize, r, &r1, &r2, NULL);			// select three mutually exclusive individuals with r from subpop n
			int ar = int(archive.size()*rnd_uni(&rnd_uni_init));	// randomly select an individual from archive
			// mutation and crossover
			int idx_rnd = int(rnd_uni(&rnd_uni_init)*numVariables);	// randomly select a dimension
			for (int j = 0; j<numVariables; j++)
			{
				double rnd = rnd_uni(&rnd_uni_init);
				if (rnd<=population[n].indiv[r].cr || j == idx_rnd)				//------------- here for setting of crossover rate --------
					child.indiv[r].x_var[j] = population[n].indiv[r].x_var[j] + population[n].indiv[r].fm*(population[n].indiv[r1].x_var[j] - population[n].indiv[r2].x_var[j]) + population[n].indiv[r].fm*(population[n].bestIndiv.x_var[j] - population[n].indiv[r].x_var[j]) + population[n].indiv[r].fm*(archive[ar].x_var[j] - population[n].indiv[r].x_var[j]);
				else
					child.indiv[r].x_var[j] = population[n].indiv[r].x_var[j];
				if (child.indiv[r].x_var[j] < lowBound[j])
					child.indiv[r].x_var[j] = lowBound[j];
				if (child.indiv[r].x_var[j] > uppBound[j])
					child.indiv[r].x_var[j] = uppBound[j];
				/*
				if(child[r].x_var[j] < lowBound[j])
					while(child[r].x_var[j] < lowBound[j])
						child[r].x_var[j] += uppBound[j] - lowBound[j];
				if(child[r].x_var[j] > uppBound[j])
					while(child[r].x_var[j] > uppBound[j])
						child[r].x_var[j] -= uppBound[j] - lowBound[j];
				*/
			}
		}
		// update subpop synchronously
		for (int r = 0; r<subpopSize; ++r)
		{
			child.indiv[r].subpop = n;
			child.indiv[r].isoffspring = true;
			child.indiv[r].cr = population[n].indiv[r].cr;	// inherit the f and cr of parent
			child.indiv[r].fm = population[n].indiv[r].fm;
			child.indiv[r].obj_eval();
			child_population.push_back(child.indiv[r]);
			if (child.indiv[r].y_obj[n] < population[n].indiv[r].y_obj[n])		// update parent
			{
				SCR[n].push_back(population[n].indiv[r].cr);
				SF[n].push_back(population[n].indiv[r].fm);
				population[n].indiv[r] = child.indiv[r];
			}
		}
		// uodate adaption parameters 
		if (SCR[n].size() != 0)
		{
			double sum = 0;
			for (int t = 0; t<SCR[n].size(); ++t)
				sum += SCR[n][t];
			mu_cr[n] = (1-cc)*mu_cr[n] + cc*(sum / SCR[n].size());
			double sum1 = 0;
			double sum2 = 0;
			for (int t = 0; t<SF[n].size(); ++t)
			{
				sum1 += SF[n][t] * SF[n][t];
				sum2 += SF[n][t];
			}
			mu_f[n] = (1-cc)*mu_f[n] + cc*(sum1 / sum2);
		}
	}
	// archive update
	archive_update();
	/*
	// contribution of every subpopulation and archive
	int *contribute = new int[numObjectives + 1];
	for (int i = 0; i <= numObjectives; ++i)
		contribute[i] = 0;
	for (int i = 0; i<archive.size(); ++i)
	{
		if (archive[i].isoffspring == true)
		{
			int ind = archive[i].subpop;
			contribute[ind]++;
			archive[i].isoffspring = false;
		}
	}
	fprintf(cont, "%d", gen);
	for (int i = 0; i <= numObjectives; ++i)
		fprintf(cont, "\t%d", contribute[i]);
	fprintf(cont, "\n");
	*/
}

void TRCMODE::run(int mg, int rn)
{
	// mg: maximal number of function evaluations 
	// rn: run number
	init_population();
	init_archive();
	archive_update();
	int realgen = 0;	// real generation
	while (fes <= mg)
	{
		//cout<<"run "<<rn<<": "<<fes<<" evaluations..."<<endl;
		//cal_utility(realgen++);
		evolution(realgen++);
		/*
		if (fes >= 400000 && fes <= 400000 + subpopSize*numObjectives + archive.size())
			report(archive,file1);
		if (fes >= 600000 && fes <= 600000 + subpopSize*numObjectives + archive.size())
			report(archive, file2);
		if (fes >= 800000 && fes <= 800000 + subpopSize*numObjectives + archive.size())
			report(archive, file3);
		if (fes >= 1000000 && fes <= 1000000 + subpopSize*numObjectives + archive.size())
			report(archive, file4);
		if (fes >= 1200000 && fes <= 1200000 + subpopSize*numObjectives + archive.size())
			report(archive, file5);
		if (fes >= 1400000 && fes <= 1400000 + subpopSize*numObjectives + archive.size())
			report(archive, file6);
		if (fes >= 1600000 && fes <= 1600000 + subpopSize*numObjectives + archive.size())
			report(archive, file7);
		if (fes >= 1800000 && fes <= 1800000 + subpopSize*numObjectives + archive.size())
			report(archive, file8);
		if (fes >= 2000000 && fes <= 2000000 + subpopSize*numObjectives + archive.size())
			report(archive, file9);
		*/
	}
	//report the PF
	report(archive, file);
	report_x(archive, filex);
}

#endif
