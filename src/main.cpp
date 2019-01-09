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
 double *xyz, *box, *vel;
 double *charges, *mass, *sigma, *eps, *vij;
 double *forces;
 bool *inter;

 /* define other variables */
 double ek, elj, ec, tempK;
 int nconst;

 /* define constraints */
 const_angles *cang;
 const_bonds  *cbond;

 /* 
    read params here...write function to read everything from the input file 
    temporary define all parameters below
 */
 int nmols = 500; 
 int atoms_mol = 3;
 double eps_r = 1.0;    // relative dielectric constant
 std::string gro_file = "confout.gro";
 std::string ff_file  = "water.ff";
 std::string cst_file = "water.const";
 int tag_mol = 1;
 int tag_atom = 1;      // this can only be 1 or 2 since atom0 is oxygen
 double lj_cut = 9.0;   // cut-off LJ interactions
 double r0 = 0.6572;    // OH shortest distance in A
 double dr = 0.02;      // grid step size
 int Ngrid   = 56;
 int alpha = 0;         // active quantum state
 double dt = 0.001;     // time-step in ps    

 /* 
    end of temporary input 
 */

 int natoms = nmols*atoms_mol;

 /* define quantum subsystem */
 Subsystem QS(Ngrid, tag_mol, tag_atom, r0, dr, atoms_mol, natoms);

 /* allocate arrays */
 cang    = (const_angles *)calloc(1,sizeof(const_angles));
 cbond   = (const_bonds *)calloc(1,sizeof(const_bonds));
 xyz     = (double *)calloc(natoms*3, sizeof(double)); 
 vel     = (double *)calloc(natoms*3, sizeof(double));
 box     = (double *)calloc(3, sizeof(double));
 charges = (double *)calloc(natoms, sizeof(double));
 mass    = (double *)calloc(natoms, sizeof(double));
 sigma   = (double *)calloc(natoms*natoms, sizeof(double));
 eps     = (double *)calloc(natoms*natoms, sizeof(double));
 vij     = (double *)calloc(natoms*natoms, sizeof(double));
 inter   = (bool *)calloc(natoms*natoms, sizeof(double));
 forces  = (double *)calloc(3*natoms, sizeof(double));

 /* Open files to write some runtime info */
 ET_File = fopen("ET.out","w"); 

 /* read coordinates and velocities from *.gro file */
 read_xyz(gro_file, xyz, vel, natoms, box, atoms_mol);

 /* read force field parameters, return sigma_ij and eps_ij mat. */
 read_ff(ff_file, charges, mass, sigma, eps, natoms);

 /* read constraints */
 read_const(cst_file, nconst, natoms, inter, cbond, cang);

 /* create settle */
 double dOH = 1.0;
 double dHH = 1.633;
 Settle Stl(dOH, dHH, nmols, dt); //dOH and dHH for a SPC water model

 /* build coulomb interaction matrix */
 coulomb_m(charges, vij, natoms, inter);

 /* correct box */
 inbox(xyz, box, natoms);

 /* remove center of mass velocity */
 com_v(mass, vel, natoms);

 /* calculate temperature, kinetic and potential energy */
 energy(ek, elj, ec, tempK, vel, xyz, mass, natoms, nconst, 
        vij, sigma, eps, eps_r, box, lj_cut, QS);

 //std::cout << QS.e0 << " " << QS.e1 << " " << QS.w01 << " " << QS.anh 
 //          << " " << QS.qav0 << " " << QS.qav1 << std::endl;

 std::cout << " Box : " << box[0] << " " << box[1] << " " << box[2] << std::endl;
 std::cout << " Total energy = " << ek + elj + ec + QS.e0 << std::endl;

 /* main routine */
 for(int timesteps = 0; timesteps < 10; ++timesteps){
    move(QS, Stl, xyz, vel, mass, vij, sigma, eps, eps_r, box, 
         forces, natoms, lj_cut, dt, atoms_mol, nmols, alpha);
    com_v(mass, vel, natoms);
    energy(ek, elj, ec, tempK, vel, xyz, mass, natoms, nconst,
           vij, sigma, eps, eps_r, box, lj_cut, QS);
    std::cout << " Total energy = " << ek + elj + ec + QS.e0 << std::endl;
 }

 /* close files */
 fclose (ET_File);

 /* free up memory */
 free (cang); free (cbond);
 free (xyz); 
 free (vel);
 free (box);
 free (charges);
 free (mass);
 free (sigma);
 free (eps);
 free (vij);
 free (inter);
 free (forces);

 return 0;
}

