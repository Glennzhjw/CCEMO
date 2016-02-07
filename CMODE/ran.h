#ifndef RAN_H
#define RAN_H

//产生正态分布随机数
double gaussrand(double a, double b)
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
     
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
             
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

//generate a Cauthy distributed random number,C(a,b) 
double cauchyrand(double a, double b)
{
	double tem=tan(3.14159*(rnd_uni(&rnd_uni_init)-0.5));
	return (tem*b+a);
}

/*select mutually exclusive integers , which are all different from candidate
 , from the range [1 , NP]	*/
void SelectSamples(int pp,int candidate,int *r1,int *r2,int *r3)
{
	if (r1)
	{
		do
		{
			*r1 = int(pp*rnd_uni(&rnd_uni_init));
		}
		while (*r1 == candidate);
	}

	if (r2)
	{
		do
		{
			*r2 = int(pp*rnd_uni(&rnd_uni_init));
		}
		while ((*r2 == candidate) || (*r2 == *r1));
	}
	
	if (r3)
	{
		do
		{
			*r3 = int(pp*rnd_uni(&rnd_uni_init));
		}
		while ((*r3 == candidate) || (*r3 == *r1)||(*r3 == *r2));
	}
}

#endif
