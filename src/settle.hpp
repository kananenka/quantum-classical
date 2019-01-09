#include <cmath>
#include <fstream>
#include <iostream>
#include "const.hpp"
#include "util.hpp"

#ifndef __SE__
#define __SE__

class Settle {

private:

public:

  double ra, rb, rc, rc2;
  double mass, wo, wh, invdt;
  double *invmass, *bondsq, *M2, *bond;
  int nmol;
  int max_iter;
  int iatom[3], jatom[3];

  Settle(double dOH, double dHH, int nmol, double dt);

  ~Settle(){
    free (invmass);
    free (bondsq);
    free (M2);
    free (bond);
  };

  void settle1(double*, double*, double*, double*);
  void shake_w(double*, double*, double*);


};


#endif
