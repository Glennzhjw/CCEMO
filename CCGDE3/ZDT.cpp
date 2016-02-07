/*
  ZDT.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
  
#include "ZDT.h"

void ZDT1(double *x, double *fx,int numVar){
    double f1, f2, g, h;
    int i;
    f1 = x[0];
    g = 0.0;
    for(i = 1; i<numVar; i++){
        g += x[i];
    }
    g = 9.0*g/((double)numVar-1);
    g += 1.0;
    h = 1.0 - sqrt(f1/g);
    f2 = g*h;
    fx[0] = f1;
    fx[1] = f2;
    return;
}


void ZDT2(double *x, double *fx,int numVar){
    double f1, f2, g, h;
    int i;
    f1 = x[0];
    g = 0.0;
    for(i=1; i<numVar; i++){
        g += x[i];
    }
    g = 9.0*g/((double)numVar-1);
    g += 1.0;
    h = 1.0 - pow((f1/g),2.0);
    f2 = g*h;
    fx[0] = f1;
    fx[1] = f2;
    return;
}

void ZDT3(double *x, double *fx,int numVar){
	double f1, f2, g, h;
    int i;
    f1 = x[0];
    g = 0.0;
    for(i=1; i<numVar; i++){
        g += x[i];
    }
    g = 9.0*g/((double)numVar-1);
    g += 1.0;
    h = 1.0 - sqrt(f1/g) - (f1/g)*sin(10.0*PI*f1);
    f2 = g*h;
    fx[0] = f1;
    fx[1] = f2;
    return;
}

void ZDT4(double *x, double *fx,int numVar){
    double f1, f2, g, h;
    int i;
    f1 = x[0];
    g = 0.0;
    for(i=1; i<numVar; i++){
        g += x[i]*x[i] - 10.0*cos(4.0*PI*x[i]);
    }
    g += 1.0 + (10.0*((double)numVar-1.0));
    h = 1.0 - sqrt(f1/g);
    f2 = g*h;
    fx[0] = f1;
    fx[1] = f2;
    return;
}

void ZDT6(double *x, double *fx,int numVar){
	double f1, f2, g, h;
    int i;
    f1 = 1.0 - ( exp(-4.0*x[0]) ) * pow( (sin(6.0*PI*x[0])),6.0 );
    g = 0.0;
    for (i=1; i<numVar; i++){
        g += x[i];
    }
    g = g/((double)numVar - 1.0);
    g = pow(g,0.25);
    g = 1.0 + 9.0*g;
    h = 1.0 - pow((f1/g),2.0);
    f2 = g*h;
    fx[0] = f1;
    fx[1] = f2;
    return;
}

void limitsZDT1(double *minLimit, double *maxLimit, int numVar){
	int i;
	for(i = 0; i < numVar; i++){
		minLimit[i] = 0.0;
		maxLimit[i] = 1.0;
	}
}

void limitsZDT2(double *minLimit, double *maxLimit, int numVar){
	int i;
	for(i = 0; i < numVar; i++){
		minLimit[i] = 0.0;
		maxLimit[i] = 1.0;
	}
}

void limitsZDT3(double *minLimit, double *maxLimit, int numVar){
	int i;
	for(i = 0; i < numVar; i++){
		minLimit[i] = 0.0;
		maxLimit[i] = 1.0;
	}
}

void limitsZDT4(double *minLimit, double *maxLimit, int numVar){
	int i;
	minLimit[0] = 0.0;
	maxLimit[0] = 1.0;
	for(i = 1; i < numVar; i++){
		minLimit[i] = -5.0;
		maxLimit[i] = 5.0;
	}
}

void limitsZDT6(double *minLimit, double *maxLimit, int numVar){
	int i;
	for(i = 0; i < numVar; i++){
		minLimit[i] = 0.0;
		maxLimit[i] = 1.0;
	}
}
