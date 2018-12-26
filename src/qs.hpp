#include <cmath>
#include <ctime>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include "pes.hpp"
#include "const.hpp"
#include "util.hpp"
#include "mkl_lapacke.h"

#ifndef __QS__
#define __QS__

class Subsystem {

private:
  int info;
  char jobz;
  char uplo;
  FILE *Efile;

public:
  int ngrid;
  int tag_mol;
  int tag_atom;
  int indH;
  int cindO, cindH;
  double r0;
  double dr;
  double e0, e1;
  double qav0, qav1;
  double drb, mass;
  double w01, w12, anh;

  double *Es, *Esb, *Ek, *Hqs, *evals, *evecs;

  Subsystem(int Ngrid, int tag_mol, int tag_atom,
            double r0, double dr, int atoms_mol);

  ~Subsystem(){
     free (Ek);
     free (Es);
     free (Esb);
     free (Hqs);
     free (evals);
     free (evecs);
     fclose (Efile);
  };

  void energy(double *, double*, double*, double*,
              int ,int ,int ,double *, double);

  void kinetic();

  void writeEn();

};

#endif
