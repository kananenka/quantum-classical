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

  int shake_iter_max;  
  double shake_tol;
  int nmax, Ndf;
  double tmass, dt2;
  FILE *ffile, *sfile;

  public:

  double Ek, Ep, ELJ, Ec, Et;
  double Tk, rho, vol, rcut, dt; 
  std::string data_file;
  int natoms, ncons;
  double *xyz, *xyz_s, *vel, *mass, *box, *sig, *eps, *chg, *force, *fc, *flj;   
  double *cst, *dxt, *dyt, *dzt, *dxx, *dyy, *dzz;
  CONSTR *cbond;
  INTERM *imat;

  System(int, std::string, double, double);

  ~System()
  {
     free(xyz), free(xyz_s), free(vel), free(mass);
     free(box), free(sig), free(eps);
     free(cst), free(chg), free(cbond);
     free(cbond->inda), free(cbond->indb), free(cbond->val);
     free(imat->LJa), free(imat->LJb), free(imat->eps), free(imat->sig);
     free(imat->Ca), free(imat->Cb), free(imat->vij);
     free(force), free(fc), free(flj);
     free(dxt), free(dyt), free(dzt), free(dxx), free(dyy), free(dzz);
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

   void shake();

};

#endif
