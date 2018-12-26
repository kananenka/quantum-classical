#include <cmath>
#include <ctime>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <chrono>
#include "quant.hpp"

void quantum(int tag_mol, int tag_atom, double* sigma, double* eps, 
             double* vij, double* xyz, int atoms_mol, int nmol, 
             int natoms, double* box, double r0, double dr, int Nst,
             double r_cut, double eps_r)
{
/*
    Here we deal with the quantum subsystem. Use DVR basis to
    construct the Hamiltonian
    
    Return:
    -------

    - Two lowest eigenvalues
    - Expectation value <Psi|q|Psi>
    - Hellmann-Feynman forces, for each of the two
      lowest eigenstates

    December 2018

*/  
  std::cout << " Entering QM part " << std::endl;

  /* set indices for tagged atoms */
  int indH = tag_mol*atoms_mol + tag_atom;

  std::cout << " Selected H2O molecule " << tag_mol << std::endl;
  std::cout << " Selected H atom " << indH << std::endl;

  int cindO  = 3*tag_mol*atoms_mol;
  int cindH  = 3*indH;

  /* calculate unit vector in OH direction: dOH{x,y,z} */
  double dOHx = xyz[cindH]   - xyz[cindO];
  double dOHy = xyz[cindH+1] - xyz[cindO+1];
  double dOHz = xyz[cindH+2] - xyz[cindO+2];

  dOHx = minImage(dOHx, box[0]);
  dOHy = minImage(dOHy, box[1]);
  dOHz = minImage(dOHz, box[2]);

  double rOH = sqrt(dOHx*dOHx + dOHy*dOHy + dOHz*dOHz);
 
  dOHx /= rOH;
  dOHy /= rOH;
  dOHz /= rOH;

  /* 
     Build quantum Hamiltonian. Note that DVR used here
     employs atomic units. All energies are first calculated
     in kJ/mol for convenience and then converted to hartrees 
  */
  double *Hqs, *Es, *Esb;
  Hqs = (double *)calloc(Nst*Nst, sizeof(double));
  Es  = (double *)calloc(Nst, sizeof(double));
  Esb = (double *)calloc(Nst, sizeof(double));
 
  for(int n=0; n<(Nst*Nst); ++n)
     Hqs[n] = 0.0;

  /* adding potential energy */
  double ri, lhx, lhy, lhz, rx, ry, rz, rij;
  double Vqs, Vcoul, Vlj;
  double sr, sr2, sr4, sr6, sr12;
  double sr_s, sr2_s, sr4_s, sr6_s, sr12_s, r_s, vt;

  for(int i=0; i<Nst; ++i){
     /* move hydrogen atom ri A away from oxygen atom */ 
     ri = r0 + dr*i;
     /* add quantum subsystem energy */
     Vqs = q_potential(ri, rOH);
     /* 
        adding the coupling terms: V(Coulomb) + V(LJ), but
        since LJ sites are located only at oxygen atoms
        there will be no contribution from LJ interactions
        to the coupling Hamiltonian 
     */
     lhx = xyz[cindO]   + ri*dOHx;
     lhy = xyz[cindO+1] + ri*dOHy;
     lhz = xyz[cindO+2] + ri*dOHz;
     Vcoul = 0.0;
     for(int n=0; n<natoms; ++n){
        vt = vij[indH*natoms+n];
        if(abs(vt) > 1e-6){
           rx = lhx - xyz[3*n];
           ry = lhy - xyz[3*n+1];
           rz = lhz - xyz[3*n+2];
//           std::cout << " distance to: " << n << " " << vij[indH*natoms+n] << std::endl;
           rx = minImage(rx, box[0]);
           ry = minImage(ry, box[1]);
           rz = minImage(rz, box[2]);
//           std::cout << rx << " " << ry << " " << rz << std::endl;
           rij = sqrt(rx*rx + ry*ry + rz*rz);
           Vcoul += vt/(eps_r*rij);
        }
        /* add LJ contribution here if needed */
     }
     Hqs[i*Nst + i] = (Vqs + Vcoul)/kJmH;
     Es[i]  = Vqs/kJmH;
     Esb[i] = Vcoul/kJmH;
  }

  /* add kinetic energy */
  double p;
  double drb = dr*atob;
  double mass = mH*(mD + mO)/(mH + mD + mO); 
  mass *= aum/mass_e;

  for(int i=0; i<Nst; ++i){
     for(int j=0; j<Nst; ++j){
        if(i==j){
           p = M_PI*M_PI/3.0;
        }else{
           p = 2.0/((i-j)*(i-j));
        }
        Hqs[i*Nst + j] += p*pow(-1.0,(i-j))/(2.0*mass*drb*drb);
     }
  }
  
  /* diagonalize the Hamiltonian */
  int info;
  double *evals;
  evals = (double *)calloc(Nst, sizeof(double));

  char jobz = 'V';
  char uplo = 'L';

  info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplo, Nst, Hqs, Nst, evals);

  if(info != 0){
    std::cout << " Error ! Diagonalization has failed. " << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << evals[0] << " " << evals[1] << " " << (evals[0]-evals[1])*219416.0 
            << std::endl;

  /* Calculate <q> expectation value */
  double rexp = 0.0;
  for(int n=0; n<Nst; ++n)
     rexp += (r0 + dr*n)*Hqs[n*Nst+n]*Hqs[n*Nst+n];

  std::cout << " <q> = " << rexp << std::endl;

  /* Write energy components into a file */
  FILE *Efile = fopen("pes.txt", "w");
  for(int n=0; n<Nst; ++n)
     fprintf(Efile, " %f %f %f \n", (r0 + dr*n), Es[n], Esb[n]);

  fclose (Efile);

  /* free arrays */
  free (Hqs);
  free (evals);
  free (Es);
  free (Esb);

  return;
}
