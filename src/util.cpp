#include <complex>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "util.hpp"

void inbox(double* xyz, double* box, int natoms)
{
/*
  
  If an atom left the box put it back into the box
  on the other side of it
  
  December 2018

*/
  double ax, ay, az;

  for(int n=0; n<natoms; ++n){
     moveb(xyz[3*n],   box[0]);
     moveb(xyz[3*n+1], box[1]);
     moveb(xyz[3*n+2], box[2]);
  }

  return;
}

void moveb(double &d, double box)
{
  if(d > box)
     d -= box;
  if(d < 0.0)
     d += box;
  return;
}

double minImage(double dist, double boxl)
{
/*
  Minimum image convention
*/
  return dist - boxl*round(dist/boxl);

}


double round(double v){

  return v < 0.0 ? ceil(v-0.5) : floor (v+0.5);

}
