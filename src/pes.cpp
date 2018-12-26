#include <cmath>
#include "pes.hpp"

/*
    Use this file to define potential energy surface 
    and forces of a quantum subsystem.
    Please use kJ/mol units for the output energies.

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
   Returning value of energy must be in kJ/mol !

*/
   double v = De*(1.0 - exp(-A*(r - rOH)))*(1.0 - exp(-A*(r - rOH)));
   return v;
}
