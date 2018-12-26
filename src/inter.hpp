#include "const.hpp"
#include "util.hpp"

#ifndef INT_H
#define INT_H

void coulomb_m(double* charges, double* vij, int natoms,
               bool* inter);

void force_c(double *forces, double *vij, double* sigma, double* eps,
             double *xyz, double eps_r, double *box, double r_cut, 
             int natoms);


#endif
