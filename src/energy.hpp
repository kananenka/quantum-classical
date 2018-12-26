#include "util.hpp"
#include "const.hpp"

#ifndef EN_H
#define EN_H

void energy(double &ek, double &ep, double &elj,
            double &ec, double &tempK, double* vel,
            double* xyz, double* mass, int natoms,
            int Nconst, double* vij, double* sigma,
            double *eps, double eps_r, double* box,
            double r_cut);

#endif
