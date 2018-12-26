#include "types.h"

#ifndef RW_H
#define RW_H

void read_xyz(std::string xyz_file, double* xyz, double* vel,
              int natoms, double* box, int atoms_mol);

void read_ff(std::string ff_file, double* charges, double* mass, 
             double* sigma, double* eps, int natoms);

void read_const(std::string cst_file, int &nconst, 
                int natoms, bool* inter, const_bonds* cbond,
                const_angles* cang);

#endif
