#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>
#include <stdio.h>
#include "energy.hpp"


void energy(double &ek, double &ep, double &elj, 
            double &ec, double &tempK, double* vel, 
            double* xyz, double* mass, int natoms, 
            int Nconst, double* vij, double* sigma,
            double *eps, double eps_r, double* box,
            double r_cut)
{
/*

  Here we will calculate kinetic energy, potential 
  energy and absolute temperature of the system
  
  December 2018

*/

   /* kinetic energy */
   ek = 0.0;

   for(int s=0; s<natoms; ++s)
      ek += 0.5*mass[s]*(vel[3*s]*vel[3*s] + 
                         vel[3*s+1]*vel[3*s+1] + 
                         vel[3*s+2]*vel[3*s+2]);

   /* units for the kinetic energy above are: a.u.*A^2/ps^2 
      we want to report this energy in kJ/mol */
   ek *= KE_convert;

   std::cout << " Kinetic energy = " << ek << std::endl;

   /* To calculate absolute temperature we use the following
      formula: E_kin = (1/2) * N_{df} * K_B * T,
      where N_{df} is the number of degrees of freedom:
      N_{df} = 3*N - N_{constraints} - N_{com}
      where N is the total number of particles (atoms) in the
      system, N_{constraints} is the number of constraints,
      and N_{com} = 3 additional degrees of freedom must be 
      removed, because the three center-of-mass velocities are 
      constants of the motion, which we set to zero */

   int Ndf = 3*natoms - Nconst - 3;

   tempK = 2.0*ek*1.0e3/(Na*Ndf*Kb);   

   std::cout << " Temperature = " << tempK << std::endl;

   /* potential energy: separate Lennard-Jones and Coulomb terms */
   double rx, ry, rz, rij, r_s;
   double sr, sr2, sr4, sr6, sr12, sr_s, sr2_s, sr4_s, sr6_s, sr12_s;
   double ljtemp;

   elj = 0.0;
   ec  = 0.0;

   for(int i=0; i<natoms; ++i){
      for(int j=0; j<i; ++j){
         /* vij = 0 means the same molecule for which LJ is also zero */
         if(abs(vij[i*natoms+j]) > 0.0){
            rx = xyz[3*i]   - xyz[3*j];
            ry = xyz[3*i+1] - xyz[3*j+1];
            rz = xyz[3*i+2] - xyz[3*j+2];
            rx = minImage(rx, box[0]);
            ry = minImage(ry, box[1]);
            rz = minImage(rz, box[2]);
            rij = sqrt(rx*rx + ry*ry + rz*rz);
            /* Coulomb */
            ec += vij[i*natoms+j]/(eps_r*rij);
            if(rij <= r_cut){
               /* Lennard-Jones truncated and shifted */
               if(abs(eps[i*natoms+j]) > 0.0){
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
                   r_s    = rij/r_cut; 
                   //ljtemp = (sr12 - sr6) - (sr12_s - sr6_s); 
                   ljtemp = sr12 - sr6 
                          + r_s*(6.0*sr12_s - 3.0*sr6_s) 
                          - 7.0*sr12_s + 4.0*sr6_s;
                   elj   += 4.0*eps[i*natoms+j]*ljtemp;
               }
            }
         }
      }
   } 

   std::cout << " Lennard-Jones energy = " << elj << std::endl;   
   std::cout << " Coulomb energy = " << ec << std::endl; 

   return;
}
