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
 double ek, ep, elj, ec, tempK;
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

 /* 
    end of temporary input 
 */

 int natoms = nmols*atoms_mol;

 /* define quantum subsystem */
 Subsystem QS(Ngrid, tag_mol, tag_atom, r0, dr, atoms_mol);

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

 /* build coulomb interaction matrix */
 coulomb_m(charges, vij, natoms, inter);

 /* correct box */
 //inbox(xyz, natoms, box);

 /* remove center of mass velocity */
 com_v(mass, vel, natoms);

 /* calculate temperature, kinetic and potential energy */
 energy(ek, ep, elj, ec, tempK, vel, xyz, mass, natoms, 
        nconst, vij, sigma, eps, eps_r, box, lj_cut);

 /* calculate properties of a quantum subsystem */
 QS.energy(sigma, eps, vij, xyz, atoms_mol, nmols,
           natoms, box, eps_r);

 /* calculate classical forces */
 force_c(forces, vij, sigma, eps, xyz, eps_r, box, 
         lj_cut, natoms);

 std::cout << QS.e0 << " " << QS.e1 << " " << QS.w01 << " " << QS.anh << " " << QS.qav0 << " " << QS.qav1 << std::endl;

 /* After each integration step the coordinates
of the particles must be examined. If a particle is found to have left
the simulation region, its coordinates must be readjusted to bring it back
inside, which is equivalent to bring in an image particle through the opposite
boundary. Supposing that the simulation region is a rectangular box, this is
done by adding to or subtracting from the affected particle coordinate the
size L of the box along the corresponding direction. May be call it from
verlet/leapfrog propagation routine*/

 /* close files */
 fclose (ET_File);

 /* free up memory */
 free (xyz); 
 free (box);
 free (vel);
 free (charges);
 free (mass);
 free (sigma);
 free (eps);
 free (vij);
 free (inter);
 free (cang); free (cbond);
 free (forces);

 return 0;
}

