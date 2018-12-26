#include "inter.hpp"
#include "const.hpp"
#include "util.hpp"
#include "pes.hpp"
//#include "extern.hpp"
#include "mkl_lapacke.h"

#ifndef QUA_H
#define QUA_H

void quantum(int tag_mol, int tag_atom, double* sigma,
             double* eps, double* vij, double* xyz, int atoms_mol, int nmol,
             int natoms, double* box, double r0, double dr, int Nst, 
             double r_cut, double eps_r);

#endif
