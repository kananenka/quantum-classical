#include <cmath>
#include <fstream>
#include <iostream>
#include "const.hpp"

#ifndef __SE__
#define __SE__

class Settle {

private:

public:

  double ra, rb, rc, rc2;
  double mass, wo, wh, invdt;
  int nmol;

  Settle(double dOH, double dHH, int nmol, double dt);
  ~Settle(){};

  void settle1(double* xyzb, double* xyza, double* vel);


};


#endif
