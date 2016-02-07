#ifndef __COMMON_H_
#define __COMMON_H_


double distanceArray(double vec1[], double vec2[], int dim)
{
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+= (vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}

double distanceVector(vector <double> &vec1, vector <double> &vec2)
{
    double sum = 0;
	for(int n=0; n<vec1.size(); n++)
	    sum+=(vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}


double norm_vector(vector <double> &x)
{
	double sum = 0;
	for(int i=0;i<x.size();i++)
        sum = sum + x[i]*x[i];
    return sqrt(sum);
}

double sum_vector(vector<double>&vec)
{
	double sum = 0;
	for(int i=0;i<vec.size();i++)
        sum = sum + vec[i];
    return sum;
}

double innerproduct(vector <double>&vec1, vector <double>&vec2)
{
    double sum = 0;
	for(int i=0; i<vec1.size(); i++)
		sum+= vec1[i]*vec2[i];
	return sum;
}

void minfastsort(double x[], int idx[], int n, int m)
{
    for(int i=0; i<m; i++)
	{
	    for(int j=i+1; j<n; j++)
			if(x[i]>x[j])
			{
			    double temp = x[i];
				x[i]        = x[j];
				x[j]        = temp;
				int id      = idx[i];
				idx[i]      = idx[j];
				idx[j]      = id;
			}
	}
}

// 计算归一化的欧几里得距离，vec1 到其它vec2的距离 
double Euclid_Dist(vector <double> &vec1, vector <double> &vec2)
{
	double sum = 0;
	for(int n=0; n<vec1.size(); n++)
	{
		//if (vec1[n] < vec2[n])			经过SPEA2中的SDE操作后的欧几里得距离
		sum += pow(vec1[n] - vec2[n], 2) / pow(fun_max[n] - fun_min[n], 2);
	}
	return sqrt(sum);
}

//K近邻比较函数
bool cmp_crowd(double *a, double *b)
{
	for(int i=0;i<10;++i)				//这里最多比较前NA个邻近距离 
		if(a[i]!=b[i])	
  			return a[i]<b[i];
	return  false;
}

// record index orser after sorting
void sort_dist_index(double *a, int *b, int left, int right)
{
	int pivot;
	int i, j;
	double temp;
	double temp_a;
	int temp_b;
	if(left<right)
	{
		temp=a[right];
		i=left-1;
		for(j=left;j<right;++j)
		{
			if(a[j]<=temp)
			{
				i+=1;
				temp_a=a[i];
				a[i]=a[j];
				a[j]=temp_a;
				temp_b=b[i];
				b[i]=b[j];
				b[j]=temp_b;
			}
		}
		pivot=i+1;
		temp_a=a[pivot];
		a[pivot]=a[right];
		a[right]=temp_a;
		temp_b=b[pivot];
		b[pivot]=b[right];
		b[right]=temp_b;
		sort_dist_index(a,b,left,pivot-1);
		sort_dist_index(a,b,pivot+1,right);
	}
}

#endif