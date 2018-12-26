#include <ctime>
#include <cmath>
#include <complex>
#include <string>
#include <chrono>
#include <iostream>
#include "move.hpp"

/*
    implement Leapfrog algorithm here....
*/

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
