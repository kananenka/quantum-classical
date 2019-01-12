#include "const.hpp"

#ifndef TY_H
#define TY_H

/*  bond constraints */
typedef struct {
    int N;
    int *inda;
    int *indb;
    double *val;
} CONSTR;

/* interaction matrix */
typedef struct {
  int nLJ;
  int  *LJa;
  int  *LJb;
  double *eps;
  double *sig;
  int nC;
  int  *Ca;
  int  *Cb;
  double *vij;
} INTERM;

#endif
