#include "rand.h"
#include <stdlib.h>
#include <math.h>

#define PI 3.1415926535897932384626433832795
long rnd_uni_init;
double sinus_init=0.5637112;
double logistic_init=0.3133;
double circle_init=0.43111;
double gauss_init=0.11;
double henon_init1=0.547, henon_init2=0.152;
double sinusoidal_init=0.4323312;
double tent_init=0.27114;
/*------Constants for rnd_uni()---------------------------------------*/
#define INF 1.0e14

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

//generate a uniformly distributed random number in [0,1)
double rnd_uni(long *idum)
{
	long j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0)
	{
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--)
		{
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;

}/*------End of rnd_uni()--------------------------*/

//产生正态分布随机数
double gaussrand(double a, double b)
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if ( phase == 0 ) {
		do {
			double U1 = rnd_uni(&rnd_uni_init);
			double U2 = rnd_uni(&rnd_uni_init);

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	//return X;
	return (X*b+a);
} 

//generate a Cauthy distributed random number,C(a,b) 
double cauchyrand(double a, double b)
{
	double tem=tan(PI*(rnd_uni(&rnd_uni_init)-0.5));
	return (tem*b+a);
}

int rnd(int low, int high)
{
	int res;
	if (low >= high)
	{
		res = low;
	}
	else
	{
		res = low + (int)(rnd_uni(&rnd_uni_init)*(high-low+1));
		if (res > high)
		{
			res = high;
		}
	}
	return (res);
}

/* Fisher鈥揧ates shuffle algorithm */
void shuffle(int *x,int size){
	int i,aux,k=0;
	for(i = size-1; i>0; i--){
		/* get a value between cero and i  */
		k = rnd(0,i); 
		/* exchange of values */
		aux = x[i]; 
		x[i] = x[k];
		x[k] = aux;
	}
}

int flip_r(float prob){
	if(rnd_uni(&rnd_uni_init) <= prob){
		return(1);
	}else{
		return(0);
	}
}

double rndreal (double low, double high)
{
	return (low + (high-low)*rnd_uni(&rnd_uni_init));
}

double HUPRandomExponential(double mu)
{
	double u;
	do
	{   
		u = rnd_uni(&rnd_uni_init);
	}
	while(u == 0);// discard 0
	return -mu * log(u);
}

double HUPRandomLevy(double c, double alpha)
{
	if((alpha <= 0) || (alpha > 2))  
		return 0;

	double u;
	do
	{
		u = rnd_uni(&rnd_uni_init);
	}while(u == 0);        // discard 0
	u = PI * (u - 0.5);

	if(alpha == 1)        // Cauchy case 
		return c * tan(u);

	double v;
	do
	{
		v = HUPRandomExponential(1.0);
	}while(v == 0);

	if(alpha == 2)        // Gaussian case 
		return c * (2 * sin(u) * sqrt(v));

	// general case
	double t = sin(alpha * u) / pow(cos(u), 1 / alpha);
	double s = pow(cos((1 - alpha) * u) / v, (1 - alpha) / alpha);
	return c * t * s;
}

double LevyRand(double c, double alpha)
{
	return HUPRandomLevy(c,alpha);
}

double sinusMap()
{
	sinus_init=2.3*pow(sinus_init,2.0*sin(PI*sinus_init));
	sinus_init=sinus_init-floor(sinus_init);
	return sinus_init;
}

double logisticMap()
{
	double a=4.0;

	logistic_init=a*logistic_init*(1.0-logistic_init);
	return logistic_init;
}

double circleMap()
{
	double a,b;
	a=0.5;
	b=0.2;

	circle_init=circle_init+b-(a/2*PI)*
		sin(2*PI*circle_init);
	circle_init=circle_init-floor(circle_init);
	return circle_init;
}

double gaussMap()
{
	if(gauss_init==0)
	{
		gauss_init=0;
		return gauss_init;
	}
	else
	{
		gauss_init=1.0/gauss_init;
		gauss_init=gauss_init-floor(gauss_init);
		return gauss_init;
	}
}

double henonMap()
{
	double a,b;
	a=1.4;
	b=0.3;
	double newx;

	newx=1.0-a*henon_init2*henon_init2+b*henon_init1;
	henon_init1=henon_init2;
	henon_init2=newx;
	return henon_init2;
}

double sinusoidalMap()
{
	double a;
	a=2.3;

	sinusoidal_init=a*sinusoidal_init*sinusoidal_init*
		sin(PI*sinusoidal_init);
	return sinusoidal_init;
}

double tentMap()
{
	if(tent_init<0.7)
	{
		tent_init=tent_init/0.7;
	}
	else
	{
		tent_init=10.0/3.0*tent_init*(1.0-tent_init);
	}
	return tent_init;
}