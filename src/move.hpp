#include "qs.hpp"
#include "inter.hpp"
#include "util.hpp"
#include "settle.hpp"

#ifndef MV_H
#define MV_H

void move(Subsystem &QS, Settle &Stl, double* xyz, double* vel, double* mass,
          double* vij, double* sigma, double* eps, double eps_r,
          double* box, double* forces, int natoms, double lj_cut, 
          double dt, int atoms_mol, int nmols, int alpha);

void com_v(double* mass, double* vel, int natoms);

#endif
