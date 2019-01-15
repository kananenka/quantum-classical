#include <string>
#include <iostream>
#include <cstdlib>
#include "types.h"
#include "const.hpp"
#include "util.hpp"

#ifndef __SYS__
#define __SYS__

/* 
    System class will contain all information about total
    system such as coordinates, velocities, constraints, etc.

    January 2019

*/

class System {

  private:
  
  int nmax, Ndf;
  double tmass;
  FILE *ffile, *sfile;

  public:

  double Ek, Ep, ELJ, Ec, Et;
  double Tk, rho, vol, rcut, dt; 
  std::string data_file;
  int natoms, ncons;
  double *xyz, *vel, *mass, *box, *sig, *eps, *chg, *force, *fc, *flj;   
  double *cst;
  CONSTR *cbond;
  INTERM *imat;

  System(int, std::string, double, double);

  ~System()
  {
     free(xyz), free(vel), free(mass);
     free(box), free(sig), free(eps);
     free(cst), free(chg), free(cbond);
     free(cbond->inda), free(cbond->indb), free(cbond->val);
     free(imat->LJa), free(imat->LJb), free(imat->eps), free(imat->sig);
     free(imat->Ca), free(imat->Cb), free(imat->vij);
     free(force), free(fc), free(flj);
   };

   void read_in();

   void build_interaction();

   void com_v();

   void inbox();

   void kinetic();

   void potential();

   void printe();

   void vv_vel();

   void vv_xyz();

   void save();

};

#endif
