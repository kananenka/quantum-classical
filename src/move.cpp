#include <ctime>
#include <cmath>
#include <complex>
#include <string>
#include <chrono>
#include <iostream>
#include "move.hpp"

void move(Subsystem &QS, Settle &Stl, double* xyz, double* vel, double* mass, 
          double* vij, double* sigma, double* eps, double eps_r,
          double* box, double* forces, int natoms, double lj_cut, 
          double dt, int atoms_mol, int nmols, int alpha)
{
/*
   Perform a single velocity-verlet step

   December 2018

*/

  /* save current corrdinates and velocities */
  double *xyzb;

  xyzb = (double *)calloc(natoms*3, sizeof(double));

  for(int s=0; s<(3*natoms); ++s)
     xyzb[s] = xyz[s];

  std::cout << "---- Step ---- " << std::endl;
  std::cout << " Starting xyz and vel " << std::endl;
  std::cout << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
  std::cout << vel[0] << " " << vel[1] << " " << vel[2] << std::endl;

  /* calculate forces */
  QS.eval(sigma, eps, vij, xyz, atoms_mol, nmols, natoms,
          box, eps_r, alpha);

  force_c(QS, forces, vij, sigma, eps, xyz, eps_r, box,
          lj_cut, natoms); 

  /* Step 1: v(t+dt/2) = v(t) + dt*F(t)/(2*m) */
  for(int n=0; n<natoms; ++n){
     vel[3*n  ] += dt*forces[3*n  ]/(2.0*mass[n]);
     vel[3*n+1] += dt*forces[3*n+1]/(2.0*mass[n]);
     vel[3*n+2] += dt*forces[3*n+2]/(2.0*mass[n]);
  }

  /* Step 2: x(t+dt) = x(t) + dt*v(t+dt/2) */
  for(int n=0; n<natoms; ++n){
     xyz[3*n  ] += dt*vel[3*n  ];
     xyz[3*n+1] += dt*vel[3*n+1];
     xyz[3*n+2] += dt*vel[3*n+2];
  }

  /* if atoms moved outside of the box put them back on the other side */
  inbox(xyz, box, natoms);

  /* Step 3: Evaluate forces at a new location */
  QS.eval(sigma, eps, vij, xyz, atoms_mol, nmols, natoms,
          box, eps_r, alpha);

  force_c(QS, forces, vij, sigma, eps, xyz, eps_r, box,
          lj_cut, natoms);

  /* Step 4: v(t+dt) = v(t+dt/2) + dt*F(t+dt)/(2*m) */
  for(int n=0; n<natoms; ++n){
     vel[3*n]   += dt*forces[3*n]/(2.0*mass[n]);
     vel[3*n+1] += dt*forces[3*n+1]/(2.0*mass[n]);
     vel[3*n+2] += dt*forces[3*n+2]/(2.0*mass[n]);
  }

  /* apply contraints here */
  Stl.settle1(xyzb, xyz, vel);

  /* if atoms moved outside of the box put them back on the other side */
  inbox(xyz, box, natoms);

  std::cout << " final xyz and vel " << std::endl;
  std::cout << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
  std::cout << vel[0] << " " << vel[1] << " " << vel[2] << std::endl;
  
  free (xyzb);

  return;
}

void com_v(double* mass, double* vel, int natoms)
{
/*
   This function sets center-of-mass velocity to zero.
   Normally there is no net external force acting on the system and 
   the center-of-mass velocity should remain constant. In practice, 
   however, the update algorithm develops a very slow change in the 
   center-of-mass velocity, and thus in the total kinetic energy of 
   the system, specially when temperature coupling is used.

   December 2018

*/

   double comv_x = 0.0;
   double comv_y = 0.0;
   double comv_z = 0.0;
   double tmass  = 0.0;

   for(int s=0; s<natoms; ++s){
      comv_x += mass[s]*vel[3*s];
      comv_y += mass[s]*vel[3*s+1];
      comv_z += mass[s]*vel[3*s+2];
      tmass  += mass[s];
   }

   comv_x /= tmass;
   comv_y /= tmass;
   comv_z /= tmass;

   for(int s=0; s<natoms; ++s){
      vel[3*s]   -= comv_x;
      vel[3*s+1] -= comv_y;
      vel[3*s+2] -= comv_z;
   }

   std::cout << " Removing center-of-mass velocity: (" 
             << comv_x << ", " << comv_y << ", " << comv_z << ") " << std::endl; 

   return;
}
