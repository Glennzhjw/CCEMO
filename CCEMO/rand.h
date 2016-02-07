#ifndef _RAND_H_
#define _RAND_H_

extern long rnd_uni_init;
extern double sinus_init;
extern double logistic_init;
extern double circle_init;
extern double gauss_init;
extern double henon_init1, henon_init2;
extern double sinusoidal_init;
extern double tent_init;

double gaussrand(double a, double b);
double rnd_uni(long *idum);
double cauchyrand(double a, double b);
int rnd(int low, int high);
void shuffle(int *x,int size);
int flip_r(float prob);
double rndreal(double low, double high);

double HUPRandomExponential(double mu);
double HUPRandomLevy(double c, double alpha);
double LevyRand(double c, double alpha);

//CIDE
double sinusMap();
double logisticMap();
double circleMap();
double gaussMap();
double henonMap();
double sinusoidalMap();
double tentMap();

#endif