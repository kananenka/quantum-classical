#include <cmath>
#include <ctime>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "inter.hpp"

/*
   implement particle-mesh Ewald summation here ?
   https://github.com/chemlab/chemlab/blob/180287718380ae635aa6401cbd183d93e3bdc81c/chemlab/md/ewald.py

*/

void coulomb_m(double* charge, double* vij, int natoms,
               bool* inter)
{
/*
    Build Coulomb matrix:
    V_ij = f*q_i*q_j
  
    where f = 1/4*pi*e, and e is the permitivity 

    December 2018
*/

  for(int i=0; i<natoms; ++i)
     for(int j=0; j<natoms; ++j)
        if(inter[natoms*i+j])
           vij[natoms*i+j] = f_electr*charge[i]*charge[j];

  return;
}

