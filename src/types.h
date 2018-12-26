#include "const.hpp"

#ifndef TY_H
#define TY_H

/*  bond constraints */
typedef struct {
    int N;
    int inda[nbond_max];
    int indb[nbond_max];
    double val[nbond_max];
} const_bonds;

/* angles contraints */
typedef struct {
    int N;
    int inda[nang_max];
    int indb[nang_max];
    int indc[nang_max];
    double val[nang_max];
} const_angles;

#endif
