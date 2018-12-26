#include <cmath>
#include "pes.hpp"

/*
    Use this file to define the potential energy surface
    of and force of a quantum subsystem.

    Units: kJ/mol for energies, kJ/(mol*A) for forces 

    December 2018
*/

/* Define constants below */
double const De  = 555.21;   
double const A   = 2.131;    
double const rOH = 0.98;

double q_potential(double r, double req)
{
/*
   Define potential energy here. 
   Return value of energy in kJ/mol !

*/
   double v = De*pow((1.0 - exp(-A*(r - rOH))),2.0); 
   return v;
}
