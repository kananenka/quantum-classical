#include <string>
#include <cstdlib>
#include <iostream>
#include <cmath> 
#include <fstream>
#include <vector>
#include <sstream>
#include "system.hpp"

System::System(int natoms, std::string data_file):
               natoms(natoms), data_file(data_file)
{
 
 /* allocate arrays */
 xyz   = (double *)calloc(3*natoms, sizeof(double));
 vel   = (double *)calloc(3*natoms, sizeof(double));
 force = (double *)calloc(3*natoms, sizeof(double));
 fc    = (double *)calloc(3*natoms, sizeof(double));
 flj   = (double *)calloc(3*natoms, sizeof(double));
 box   = (double *)calloc(3, sizeof(double));
 mass  = (double *)calloc(natoms, sizeof(double));
 sig   = (double *)calloc(natoms*natoms, sizeof(double));
 eps   = (double *)calloc(natoms*natoms, sizeof(double));
 chg   = (double *)calloc(natoms, sizeof(double));

 cbond = (CONSTR *)calloc(1,sizeof(CONSTR));
 imat  = (INTERM *)calloc(1,sizeof(INTERM));

 /* 
    read everything from input data file: coordinates, velocities,
    force field parameters and constraints 
 */
 read_in();

 /* 
   build LJ and Coulomb interaction matrices
   total N(N-1)/2 interactions possible
   create a structure containing all this information
 */
 build_interaction();

 /* fix box */
 inbox();

 kinetic();

 potential();

}

void System::read_in()
{
/*

   Read coordinates and velocities from the gromacs *.gro file.
   December 21, 2018

*/
  std::ifstream file;
  file.open(data_file.c_str());

  std::string rs;
  std::string buf;

  /* first line is the number of atoms to read */
  {
     getline(file,rs);
     std::stringstream ss(rs);
     std::vector<std::string> strs;
     while(ss >> buf)
        strs.push_back(buf);

     int natoms_file = atoi(strs[0].c_str());
          if(natoms_file != natoms){
       std::cerr << " Number of atoms in input file "
                 << natoms_file << " vs. " << natoms << std::endl;
       exit(EXIT_FAILURE);
     }
  }

  /* 
      start reading coordinates. 
      Format: <name> <name> <id> <element> <X> <Y> <Z> <v_x> <v_y> <v_z> 

  */
  for(int t=0; t<natoms; ++t)
  {
     getline(file,rs);
     std::stringstream ss(rs);
     std::vector<std::string> strs;
     while(ss >> buf)
        strs.push_back(buf);

     for(int n=0; n<3; ++n)
        xyz[3*t+n] = atof(strs[n].c_str());

     for(int n=0; n<3; ++n)
        vel[3*t+n] = atof(strs[n+3].c_str());
  }

  /* read box */
  {
     getline(file,rs);
     std::stringstream ss(rs);
     std::vector<std::string> strs;
     while(ss >> buf)
        strs.push_back(buf);

     box[0] = atof(strs[0].c_str());
     box[1] = atof(strs[1].c_str());
     box[2] = atof(strs[2].c_str());
  }

  /* empty line */
  getline(file,rs);

  /* read number of constraints */
  {
     getline(file,rs);
     std::stringstream ss(rs);
     std::vector<std::string> strs;
     while(ss >> buf)
        strs.push_back(buf);
     ncons = atoi(strs[0].c_str());
  }

  /* check if we can fit all these constraints */
  if(ncons >= nbond_max){
    std::cout << " Too many bond constraints. Change nbond_max in const.hpp to "
              << cbond->N << std::endl;
    exit(EXIT_FAILURE);
  }

  cbond->N = ncons;

  int inda, indb;
  double val;

  cbond->inda = (int *)calloc(ncons, sizeof(int));
  cbond->indb = (int *)calloc(ncons, sizeof(int));
  cbond->val  = (double *)calloc(ncons, sizeof(double));

  for(int t=0; t<cbond->N; ++t)
  {
     getline(file,rs);
     std::stringstream ss(rs);
     std::vector<std::string> strs;
     while(ss >> buf)
        strs.push_back(buf);

     inda = atoi(strs[0].c_str());
     indb = atoi(strs[1].c_str());
     val  = atof(strs[2].c_str());

     cbond->inda[t] = inda;
     cbond->indb[t] = indb;
     cbond->val[t]  = val;
  }

  /* empty line */
  getline(file, rs);

  tmass = 0.0;
  /* number of atoms to read */
  {
     getline(file,rs);
     std::stringstream ss(rs);
     std::vector<std::string> strs;
     while(ss >> buf)
        strs.push_back(buf);

     int natoms_file = atoi(strs[0].c_str());

     if(natoms_file != natoms){
       std::cerr << " Number of interaction parameters in " << data_file 
                 << " does not match the number of atoms in the input file "
                 << natoms_file << " vs. " << natoms << std::endl;
       exit(EXIT_FAILURE);
     }
  }

  /* reading Force field parameters */
  for(int t=0; t<natoms; ++t)
  {
     getline(file,rs);
     std::stringstream ss(rs);
     std::vector<std::string> strs;
     while(ss >> buf)
        strs.push_back(buf);

     mass[t] = atof(strs[0].c_str());
     eps[t]  = atof(strs[1].c_str());
     sig[t]  = atof(strs[2].c_str());
     chg[t]  = atof(strs[3].c_str());

     tmass += mass[t];
  }
  file.close();

  /* calculate and print density */
  vol = box[0]*box[1]*box[2];
  rho = tmass/vol;
  rho *= densf;
  printf(" Density %7.2f [kg/m^3] \n ",rho);
}

void System::build_interaction()
{
/* 
   Build Lennard-Jones and Coulomb interaction matrices
   This code works with diatomic and triatomic molecules
   which can be constrained by setting up bonds. We use
   these constraints to identify molecules and remove
   interactions between atoms of the same molecule

   January 2019
*/
   nmax = natoms*(natoms-1)/2;      

   imat->LJa   = (int *)calloc(nmax, sizeof(int));
   imat->LJb   = (int *)calloc(nmax, sizeof(int));
   imat->eps   = (double *)calloc(nmax, sizeof(double));
   imat->sig   = (double *)calloc(nmax, sizeof(double));

   imat->Ca    = (int *)calloc(nmax, sizeof(int));
   imat->Cb    = (int *)calloc(nmax, sizeof(int));
   imat->vij   = (double *)calloc(nmax, sizeof(double));
   
   /* construct LJ interaction matrix first */
   double enm, snm;
   int ina, inb;
   int counter=0;
   bool connected;
   bool interacting;

   for(int n=0; n<natoms; ++n){
      for(int m=0; m<n; ++m){

         connected=false;

         for(int ib=0; ib<cbond->N; ++ib){
            ina = cbond->inda[ib];
            inb = cbond->indb[ib];     
            
            if(m != ina)
               continue;

            if(inb == n)
               connected = true;
         }

         if(connected)
            continue;

         enm = eps[n]*eps[m];
         snm = sig[n]*sig[m];

         if((abs(enm) < 1.0e-9) && (abs(snm) < 1.0e-9))
            continue;

         imat->LJa[counter] = m;
         imat->LJb[counter] = n;
         imat->eps[counter] = enm;
         imat->sig[counter] = snm;
         counter++;
      }
   }
   imat->nLJ = counter;

   std::cout << " Lennard-Jones interaction pairs: " << imat->nLJ << std::endl;

   double qq;
   counter = 0;

   /* Coulomb list: only charged atoms on different molecules interact */
   for(int n=0; n<natoms; ++n){
      for(int m=0; m<n; ++m){

         connected = false;

         for(int ib=0; ib<cbond->N; ++ib){
            ina = cbond->inda[ib];
            inb = cbond->indb[ib];

            if(m != ina)
              continue;

            if(inb == n)
               connected=true;
         }

         if(connected)
           continue;

         qq = chg[n]*chg[m];

         if(fabs(qq) < 1.0e-8)
           continue;

         imat->Ca[counter] = m;
         imat->Cb[counter] = n;
         imat->vij[counter] = f_electr*qq;
         counter++;
      }
   }
   imat->nC = counter;

   std::cout << " Coulomb interaction pairs: " << imat->nC << std::endl;

}

void System::com_v()
{
/*
   This function sets center-of-mass velocity to zero.
   Normally there is no net external force acting on the system and
   the center-of-mass velocity should remain constant. In practice,
   however, the update algorithm develops a very slow change in the
   center-of-mass velocity, and thus change the total kinetic energy 
   of the system, specially when temperature coupling is used.

   December 2018

*/

   double comv_x = 0.0;
   double comv_y = 0.0;
   double comv_z = 0.0;

   for(int s=0; s<natoms; ++s){
      comv_x += mass[s]*vel[3*s];
      comv_y += mass[s]*vel[3*s+1];
      comv_z += mass[s]*vel[3*s+2];
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
}

void System::inbox()
{
/*

  If an atom left the box put it back into the box
  on the other side of it. Take a look at AB code
  for possible incostistency

  December 2018

*/
  double corr;

  for(int n=0; n<natoms; ++n)
     for(int d=0; d<3; ++d){
        corr = (xyz[3*n+d] > box[d]) ? -box[d] : 0.0;
        xyz[3*n+d]   += corr;
        corr = (xyz[3*n+d] < box[d]) ? box[d] : 0.0;
        xyz[3*n+d] += corr;
     }
}

void System::kinetic()
{
/* 
    Calculate kinetic energy and temperature of the system

    January 2019
*/

  Ek = 0.0;
 
  for(int s=0; s<natoms; ++s)
      Ek += 0.5*mass[s]*(vel[3*s]*vel[3*s] +
                         vel[3*s+1]*vel[3*s+1] +
                         vel[3*s+2]*vel[3*s+2]);

   /* 
      convert units from g*A^2/ps^2 to kJ/mol
    */

   Ek *= KE_convert;

   printf(" Classical Kinetic energy = %9.5f [kJ/mol] \n",Ek);
    
   /* 
      To calculate absolute temperature we use the following
      formula: E_kin = (1/2) * N_{df} * K_B * T,
      where N_{df} is the number of degrees of freedom:
      N_{df} = 3*N - N_{constraints} - N_{com}
      where N is the total number of particles (atoms) in the
      system, N_{constraints} is the number of constraints,
      and N_{com} = 3 additional degrees of freedom must be 
      removed, because the three center-of-mass velocities are 
      constants of the motion, which we set to zero 
   */

   Ndf = 3*natoms - cbond->N - 3;

   Tk = 2.0*Ek/(Ndf*Kb);   

   printf(" Temperature = %5.2f [K] \n",Tk);
}

void System::potential()
{
/*
   Here we will evaluate potential energy and forces

   January 2019
*/

   ELJ = 0.0;
   Ec  = 0.0;

   for(int s=0; s<natoms; ++s){
      force[s] = 0.0;
      fc[s]    = 0.0;
      flj[s]   = 0.0;
   }

   /* Lennard-Jones truncated potential */

}
