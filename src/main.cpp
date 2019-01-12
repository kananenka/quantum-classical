#include <cmath>
#include <ctime>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <chrono>
#include "main.hpp"


int main(int argc, char** argv){
/*

  Mixed quantum-classical dynamics of water with
  quantum-mechanical treatment of a tagged OH stretch
  and classical treatment of all other degrees of 
  freedom.  

  OH stretch potential energy surface is modeled as
  Morse oscillator but in the future we will connect this
  code to some electronic structure package

  DVR will be used to obtain transition frequencies of
  OH stretch

  Alexei Kananenka, December 2018 - ...

*/

 /* define files */
 FILE * ET_File;     // temperature and energy file

 /* define arrays */
 double *forces;
 bool *inter;

 /* define other variables */
 double ek, elj, ec, tempK;
 int nconst;

 /* 
    read params here...write function to read everything from the input file 
    temporary define all parameters below
 */
 int nmols = 255; 
 int atoms_mol = 2;
 double eps_r = 1.0;    // relative dielectric constant
 std::string data_file = "data.in";
 double lj_cut = 12.0;   // cut-off LJ interactions
 double r0 = 0.6572;    // OH shortest distance in A
 double dr = 0.02;      // grid step size
 int Ngrid   = 56;
 double dt = 0.001;     // time-step in ps    

 /* 
    end of temporary input 
 */

 int natoms = nmols*atoms_mol;

 System S(natoms, data_file);
 
 S.com_v();

 exit(0);

 /* define quantum subsystem */
 //Subsystem QS(Ngrid, tag_mol, tag_atom, r0, dr, atoms_mol, natoms);

 /* Open files to write some runtime info */
 ET_File = fopen("ET.out","w"); 

 /* create settle */
 //double dOH = 1.0;
 //double dHH = 1.633;
 //Settle Stl(dOH, dHH, nmols, dt); //dOH and dHH for a SPC water model

 /* correct box */
 //inbox(xyz, box, natoms);

 /* remove center of mass velocity */
 //com_v(mass, vel, natoms);

 /* calculate temperature, kinetic and potential energy */
 //energy(ek, elj, ec, tempK, vel, xyz, mass, natoms, nconst, 
 //       vij, sigma, eps, eps_r, box, lj_cut, QS);

 //std::cout << QS.e0 << " " << QS.e1 << " " << QS.w01 << " " << QS.anh 
 //          << " " << QS.qav0 << " " << QS.qav1 << std::endl;

 //std::cout << " Box : " << box[0] << " " << box[1] << " " << box[2] << std::endl;
 //std::cout << " Total energy = " << ek + elj + ec + QS.e0 << std::endl;

 /* main routine */
 //for(int timesteps = 0; timesteps < 100; ++timesteps){
 //   move(QS, Stl, xyz, vel, mass, vij, sigma, eps, eps_r, box, 
 //        forces, natoms, lj_cut, dt, atoms_mol, nmols, alpha);
 //   com_v(mass, vel, natoms);
 //   energy(ek, elj, ec, tempK, vel, xyz, mass, natoms, nconst,
 //          vij, sigma, eps, eps_r, box, lj_cut, QS);
 //   std::cout << " Total energy = " << ek + elj + ec + QS.e0 << std::endl;
 //   fprintf(ET_File," %f %f %f %f %f %f \n ",timesteps*dt,ek,tempK,elj,ec,(ek+elj+ec));
 //}

 /* close files */
 fclose (ET_File);


 return 0;
}

