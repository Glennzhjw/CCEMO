/*
  ZDT.h

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
# ifndef _ZDT_H_
# define _ZDT_H_

# include <math.h>

#ifndef PI
# define PI 3.14159265358979
#endif

void ZDT1(double *x, double *fx,int numVar);
void ZDT2(double *x, double *fx,int numVar);
void ZDT3(double *x, double *fx,int numVar);
void ZDT4(double *x, double *fx,int numVar);
void ZDT6(double *x, double *fx,int numVar);
void limitsZDT1(double *minLimit, double *maxLimit, int numVar);
void limitsZDT2(double *minLimit, double *maxLimit, int numVar);
void limitsZDT3(double *minLimit, double *maxLimit, int numVar);
void limitsZDT4(double *minLimit, double *maxLimit, int numVar);
void limitsZDT6(double *minLimit, double *maxLimit, int numVar);

# endif
