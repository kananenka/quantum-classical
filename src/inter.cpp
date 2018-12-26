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

void force_c(double *forces, double *vij, double* sigma, double* eps,
             double *xyz, double eps_r, double *box, double r_cut, 
             int natoms)
{
/*
    Calculate classical forces: Coulomb and Lennard-Jones
    Fi = -\sum_j du(r_ij)/dr_ij * vec(r_ij)/r_ij

    Classical force is acting on every degree of freedom
    including quantum subsystem

    December 2018
*/

  double *Fc, *Flj;

  Fc  = (double *)calloc(3*natoms, sizeof(double));
  Flj = (double *)calloc(3*natoms, sizeof(double));

  double tij, elj, ljtemp;
  double rx, ry, rz, rij;
  double urx, ury, urz; 
  double sr, sr2, sr4, sr6, sr12;
  double sr_s, sr2_s, sr4_s, sr6_s, sr12_s;

  for(int i=0; i<natoms; ++i){
     for(int j=0; j<natoms; ++j){
        if(i==j) continue;
        tij = vij[i*natoms+j];
        /* check if two atoms interact via Coulomb
           forces (must belong to different molecules) */
        if(abs(tij) > 1.0e-7){
           rx = xyz[3*i]   - xyz[3*j];
           ry = xyz[3*i+1] - xyz[3*j+1];
           rz = xyz[3*i+2] - xyz[3*j+2];
           rx = minImage(rx, box[0]);
           ry = minImage(ry, box[1]);
           rz = minImage(rz, box[2]);
           rij = sqrt(rx*rx + ry*ry + rz*rz);
           /* unit vector in ij direction */
           urx = rx/rij;
           ury = ry/rij;
           urz = rz/rij;
           Fc[3*i]   = tij*urx/(eps_r*rij*rij);
           Fc[3*i+1] = tij*ury/(eps_r*rij*rij);
           Fc[3*i+2] = tij*urz/(eps_r*rij*rij);
           /* add Lennard-Jones forces */
           elj = eps[i*natoms+j];
           if(abs(elj) > 1.0e-7){
              sr     = sigma[i*natoms+j]/rij;
              sr2    = sr*sr;
              sr4    = sr2*sr2;
              sr6    = sr4*sr2;
              sr12   = sr6*sr6;
              sr_s   = sigma[i*natoms+j]/r_cut;
              sr2_s  = sr_s*sr_s;
              sr4_s  = sr2_s*sr2_s;
              sr6_s  = sr4_s*sr2_s;
              sr12_s = sr6_s*sr6_s;
              ljtemp = 12.0*sr12/rij - 6.0*sr6/rij
                     - 12.0*sr12_s/r_cut + 6.0*sr6_s/r_cut;
              Flj[3*i]   = 4.0*elj*ljtemp*urx/rij;
              Flj[3*i+1] = 4.0*elj*ljtemp*ury/rij;
              Flj[3*i+2] = 4.0*elj*ljtemp*urz/rij;
           }
        }
     }
  }

  /* write forces to external file */
  FILE *ffile = fopen("Force.txt","w");
  fprintf(ffile,"# Fc_x \t Fc_y \t Fc_z \t Flj_x \t Flj_y \t Flj_z \n");
  for(int n=0; n<natoms; ++n)
     fprintf(ffile, " %f %f %f %f %f %f \n",Fc[3*n], Fc[3*n+1], Fc[3*n+2], 
              Flj[3*n], Flj[3*n+1], Flj[3*n+2]);
  fclose(ffile);
  std::cout << " Forces are saved to: Force.txt " << std::endl;

  /* add Coulomb and LJ forces together */
  for(int s=0; s<(3*natoms); ++s)
     forces[s] = Fc[s] + Flj[s];

  free (Fc);
  free (Flj);

  return;
}

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

