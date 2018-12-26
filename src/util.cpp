#include <complex>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

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
