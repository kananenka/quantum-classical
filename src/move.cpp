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
  //std::cout << " Starting xyz and vel " << std::endl;
  //std::cout << " XYZ = " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
  //std::cout << " V = " << vel[0] << " " << vel[1] << " " << vel[2] << std::endl;

  /* calculate forces */
  //QS.eval(sigma, eps, vij, xyz, atoms_mol, nmols, natoms,
  //        box, eps_r, alpha);

  force_c(QS, forces, vij, sigma, eps, xyz, eps_r, box,
          lj_cut, natoms); 

  /* Step 1: v(t+dt/2) = v(t) + dt*F(t)/(2*m) */
  for(int n=0; n<natoms; ++n){
     vel[3*n  ] += 0.1*dt*forces[3*n  ]/(2.0*mass[n]);
     vel[3*n+1] += 0.1*dt*forces[3*n+1]/(2.0*mass[n]);
     vel[3*n+2] += 0.1*dt*forces[3*n+2]/(2.0*mass[n]);
  }

  /* Step 2: x(t+dt) = x(t) + dt*v(t+dt/2) */
  for(int n=0; n<natoms; ++n){
     xyz[3*n  ] += dt*vel[3*n  ];
     xyz[3*n+1] += dt*vel[3*n+1];
     xyz[3*n+2] += dt*vel[3*n+2];
  }

  /* if atoms moved outside of the box put them back on the other side */
  //inbox(xyz, box, natoms);

  /* Step 3: Evaluate forces at a new location */
  //QS.eval(sigma, eps, vij, xyz, atoms_mol, nmols, natoms,
  //        box, eps_r, alpha);

  force_c(QS, forces, vij, sigma, eps, xyz, eps_r, box,
          lj_cut, natoms);

  /* Step 4: v(t+dt) = v(t+dt/2) + dt*F(t+dt)/(2*m) */
  for(int n=0; n<natoms; ++n){
     vel[3*n]   += 0.1*dt*forces[3*n]/(2.0*mass[n]);
     vel[3*n+1] += 0.1*dt*forces[3*n+1]/(2.0*mass[n]);
     vel[3*n+2] += 0.1*dt*forces[3*n+2]/(2.0*mass[n]);
  }

  //std::cout << " before settle... " << std::endl;
  //std::cout << " XYZ " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
  //std::cout << " V " << vel[0] << " " << vel[1] << " " << vel[2] << std::endl;

  /* apply contraints here */
  Stl.settle1(xyzb, xyz, vel, box);
 
  /* if atoms moved outside of the box put them back on the other side */
  inbox(xyz, box, natoms);

  //std::cout << " final xyz and vel " << std::endl;
  //std::cout << " XYZ = " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
  //std::cout << " V = " << vel[0] << " " << vel[1] << " " << vel[2] << std::endl;
  
  free (xyzb);

  return;
}

