/*
  rand.c

  Author:
       Luis Miguel Antonio <lmiguel@computacion.cs.cinvestav.mx>       

  Copyright (c) 2013 Luis Miguel Antonio		*/
  
# ifndef _RAND_H_
# define _RAND_H_

extern double seed;
extern double oldrand[55];
extern int jrand;


void shuffle(int *x,int size);
int flip(float prob);
void randomize(void);
void warmup_random (double seed);
void advance_random (void);
double randomperc(void);
int rnd (int low, int high);
double rndreal (double low, double high);

# endif
