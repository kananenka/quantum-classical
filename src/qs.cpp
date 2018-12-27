#include <cmath>
#include <ctime>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include "qs.hpp"

Subsystem::Subsystem(int ngrid, int tag_mol, int tag_atom,
                     double r0, double dr, int atoms_mol, int natoms):
                     ngrid(ngrid), tag_mol(tag_mol), tag_atom(tag_atom),
                     r0(r0), dr(dr)
{
     /* allocate arrays */
     Es   = (double *)calloc(ngrid, sizeof(double));
     Esb   = (double *)calloc(ngrid, sizeof(double));
     Ek    = (double *)calloc(ngrid*ngrid, sizeof(double));
     Hqs   = (double *)calloc(ngrid*ngrid, sizeof(double));
     evals = (double *)calloc(ngrid, sizeof(double));
     evecs = (double *)calloc(2*ngrid, sizeof(double));
     Fqm   = (double *)calloc(3*natoms, sizeof(double));

     /* set indices for tagged atoms */
     indH = tag_mol*atoms_mol + tag_atom;
     cindO  = 3*tag_mol*atoms_mol;
     cindH  = 3*indH;

     /* set up atomic units for quantum Hamiltonian */
     drb = dr*atob;
     mass = mH*(mD + mO)/(mH + mD + mO);
     mass *= aum/mass_e;

     /* set up kinetic energy matrix */
     kinetic();

     /* Open file for the subsystem energies */
     Efile = fopen("pes.txt", "w");
}

void Subsystem::eval(double *sigma, double *eps, double *vij, 
                       double *xyz, int atoms_mol, int nmols, 
                       int natoms, double *box, double eps_r,
                       int alpha)
{
/*
    Here we deal with the quantum subsystem. Use DVR basis to
    construct the Hamiltonian
    
    Calculate:
    ----------

    - Two lowest eigenvalues and corresponding eigenvectors
    - Expectation value <Psi|q|Psi> for two lowest eigenstates
    - Hellmann-Feynman forces, for a given eigenstate 0 or 1

    December 2018

*/  

  /* calculate unit vector in the H-O direction: dOH{x,y,z} */
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
 
  for(int n=0; n<(ngrid*ngrid); ++n)
     Hqs[n] = 0.0;

  for(int n=0; n<ngrid; ++n){
     evals[n] = 0.0;
     evecs[n] = 0.0;
     evecs[n+ngrid] = 0.0;
  }

  /* adding potential energy */
  double ri, lhx, lhy, lhz, rx, ry, rz, rij, vt;
  double Vqs, Vcoul;

  for(int i=0; i<ngrid; ++i){
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
        /* 
           the tagged atom does not interact with itself
           because vij was set to 0 earlier
        */
        vt = vij[indH*natoms+n];
        if(abs(vt) > 1e-7){
           rx = lhx - xyz[3*n];
           ry = lhy - xyz[3*n+1];
           rz = lhz - xyz[3*n+2];
           rx = minImage(rx, box[0]);
           ry = minImage(ry, box[1]);
           rz = minImage(rz, box[2]);
           rij = sqrt(rx*rx + ry*ry + rz*rz);
           Vcoul += vt/(eps_r*rij);
        }
        /* add LJ contribution here if needed */
     }
     Hqs[i*ngrid + i] = (Vqs + Vcoul)/kJmH;
     Es[i]  = Vqs/kJmH;
     Esb[i] = Vcoul/kJmH;
  }

  /* add kinetic energy */
  for(int n=0; n<ngrid; ++n)
     for(int m=0; m<ngrid; ++m)
        Hqs[n*ngrid+m] += Ek[n*ngrid+m];

  /* diagonalize the Hamiltonian */
  jobz = 'V';
  uplo = 'L';
  info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplo, ngrid, Hqs, ngrid, evals);

  if(info != 0){
    std::cout << " Error ! Diagonalization has failed. " << std::endl;
    exit(EXIT_FAILURE);
  }

  /* Copy eigenvectors corresponding to two lowest eigenvalues */
  for(int n=0; n<ngrid; ++n){
     evecs[n] = Hqs[ngrid*n];
     evecs[ngrid+n] = Hqs[ngrid*n+1];
  }

  /* Calculate <q> expectation value corresponding to two lowest eigenvalues */
  qav0 = 0.0;
  qav1 = 0.0;
  for(int n=0; n<ngrid; ++n){
     qav0 += (r0 + dr*n)*evecs[n]*evecs[n];
     qav1 += (r0 + dr*n)*evecs[ngrid+n]*evecs[ngrid+n];
  }

  /* set up energies converting to kJ/mol and wave numbers where necessary */
  e0  = evals[0]*kJmH;
  e1  = evals[1]*kJmH;
  w01 = (evals[1] - evals[0])*Htocmi;
  w12 = (evals[2] - evals[1])*Htocmi;
  /* calculate anhramonicity in cm-1 */
  anh = w01 - w12;

  /* Write energy components into a file, move this into initialization as well */
  writeEn();

  /* Calculate Helmann--Feynman forces */
  for(int s=0; s<(3*natoms); ++s)
      Fqm[s] = 0.0;

  double tij, urx, ury, urz;

  for(int i=0; i<natoms; ++i){ 
     if(i != indH){
        /* Hellmann--Feynmann force on bath atom i
           due to quantum subsystem */
        tij = vij[indH*natoms+i];
        if(abs(tij) > 1.0e-7){
           for(int n=0; n<ngrid; ++n){
              ri = r0 + dr*n;
              /* location of H atom for a given grid point */
              lhx = xyz[cindO]   + ri*dOHx;
              lhy = xyz[cindO+1] + ri*dOHy;
              lhz = xyz[cindO+2] + ri*dOHz;
              /* distance between bath atom i and tagged atom H */
              rx  = xyz[3*i]   - lhx;
              ry  = xyz[3*i+1] - lhy;
              rz  = xyz[3*i+2] - lhz;
              rx = minImage(rx, box[0]);
              ry = minImage(ry, box[1]);
              rz = minImage(rz, box[2]);
              rij = sqrt(rx*rx + ry*ry + rz*rz);
              /* unit vector in tagged H-i atom direction */
              urx = rx/rij;
              ury = ry/rij;
              urz = rz/rij;
              /* integrate on a grid to get Hellman--Feynman forces 
                 units here are: kJ/mol for the energy and Angstrom for the distance */
              Fqm[3*i]   += dr*evecs[alpha*ngrid+n]*evecs[alpha*ngrid+n]*tij*urx/(eps_r*rij*rij);
              Fqm[3*i+1] += dr*evecs[alpha*ngrid+n]*evecs[alpha*ngrid+n]*tij*ury/(eps_r*rij*rij);
              Fqm[3*i+2] += dr*evecs[alpha*ngrid+n]*evecs[alpha*ngrid+n]*tij*urz/(eps_r*rij*rij);
           }
        }
     }
  }

  return;
}

void Subsystem::kinetic()
{
 /* 
   Build kinetic energy matrix using DVR basis
 */
  
  for(int n=0; n<(ngrid*ngrid); ++n)
     Ek[n] = 0.0;

  double p;

  for(int i=0; i<ngrid; ++i){
     for(int j=0; j<ngrid; ++j){
        if(i==j){
           p = M_PI*M_PI/3.0;
        }else{
           p = 2.0/((i-j)*(i-j));
        }
        Ek[i*ngrid + j] += p*pow(-1.0,(i-j))/(2.0*mass*drb*drb);
     }
  }
  return;
}

void Subsystem::writeEn()
{
/*
  Write energy components for each grid point of a DVR grid
*/
  for(int n=0; n<ngrid; ++n)
     fprintf(Efile, " %f %f %f \n", (r0 + dr*n), Es[n], Esb[n]);

  std::cout << " PES is written to pes.txt file. " << std::endl;
  return;
}

