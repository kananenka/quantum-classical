#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <vector>
#include <sstream>
#include <stdio.h>
#include "read_write.hpp"
#include "const.hpp"

void read_xyz(std::string xyz_file, double* xyz, double* vel, 
              int natoms, double* box, int natoms_mol)
{
/*
   Read coordinates and velocities from the gromacs *.gro file.
   December 21, 2018
*/
  std::ifstream file;
  file.open(xyz_file.c_str());
  
  std::string rs;
  std::string buf;

  // first line is a comment line
  getline(file,rs);

  // second line is the number of atoms to read
  {
     getline(file,rs);
     std::stringstream ss(rs);
     std::vector<std::string> strs;
     while(ss >> buf)
        strs.push_back(buf);

     int natoms_file = atoi(strs[0].c_str());
   
     if(natoms_file != natoms){
       std::cerr << " Number of atoms in *.gro file does not match the number of atoms in the input file " 
                 << natoms_file << " vs. " << natoms << std::endl;
       exit(EXIT_FAILURE);
     }
  }

  /* start reading coordinates. Format: <name> <name> <id> <element> <X> <Y> <Z> <v_x> <v_y> <v_z> */
  for(int t=0; t<natoms; ++t)
  {
     getline(file,rs);
     std::stringstream ss(rs);
     std::vector<std::string> strs;
     while(ss >> buf)
        strs.push_back(buf); 

     for(int n=0; n<3; ++n)
        xyz[3*t+n] = unitC*atof(strs[n+3].c_str());

     for(int n=0; n<3; ++n)
        vel[3*t+n] = unitC*atof(strs[n+6].c_str());
  } 

  /* read box */
  {
     getline(file,rs);
     std::stringstream ss(rs);
     std::vector<std::string> strs;
     while(ss >> buf)
        strs.push_back(buf);

     box[0] = unitC*atof(strs[0].c_str());
     box[1] = unitC*atof(strs[1].c_str());
     box[2] = unitC*atof(strs[2].c_str());
  }

  file.close();
  return;
}

void read_ff(std::string ff_file, double* charge, double* mass, double* sigma, 
             double* eps, int natoms)
{
/*
   Read force field parameters from file. 

   December 2018

*/
  std::ifstream file;
  file.open(ff_file.c_str());

  std::string rs;
  std::string buf;

  double *e_temp, *s_temp; 

  e_temp = (double *)calloc(natoms, sizeof(double));
  s_temp = (double *)calloc(natoms, sizeof(double));

  /* first line */
  {
     getline(file,rs);
     std::stringstream ss(rs);
     std::vector<std::string> strs;
     while(ss >> buf)
        strs.push_back(buf);

     int natoms_file = atoi(strs[0].c_str());

     if(natoms_file != natoms){
       std::cerr << " Number of atoms in " << ff_file << " does not match the number of atoms in the input file "
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

     mass[t]    = atof(strs[1].c_str());
     charge[t]  = atof(strs[2].c_str());
     s_temp[t]  = atof(strs[3].c_str());
     e_temp[t]  = atof(strs[4].c_str());
  }
  file.close();

  /* 
     build matrix of sigma_i*sigma_j and
     epsilon_i*epsilon_j. This is called OPLS definition
  */

  for(int i=0; i<natoms; ++i){
     for(int j=0; j<natoms; ++j){
        sigma[i*natoms+j] = sqrt(s_temp[i]*s_temp[j]);
        eps[i*natoms+j]   = sqrt(e_temp[i]*e_temp[j]);
     }
  }

  free (e_temp);
  free (s_temp);

  return;
}

void read_const(std::string cst_file, int &nconst,
                int natoms, bool* inter, const_bonds* cbond,
                const_angles* cang)
{
/*

   Read constraints: bonds and angles
   (dihedrals are not implemented)
   December 2018

*/
  std::ifstream file;
  file.open(cst_file.c_str());
  
  std::string rs;
  std::string buf;

  /* first line gives number of bond and angles constraints */
  {
     getline(file,rs);
     std::stringstream ss(rs);
     std::vector<std::string> strs;
     while(ss >> buf)
        strs.push_back(buf);
     cbond->N = atoi(strs[0].c_str()); 
     cang->N  = atoi(strs[1].c_str());
  }

  /* check if we can fit all these constraints */
  if(cbond->N >= nbond_max){
    std::cout << " Too many bond constraints. Change nbond_max in const.hpp to " 
              << cbond->N << std::endl;
    exit(EXIT_FAILURE);
  }

  if(cang->N >= nang_max){
    std::cout << " Too many angle constraints. Change nbond_max in const.hpp to "
              << cang->N << std::endl;
  }

  nconst = cbond->N + cang->N;

  int inda, indb, indc;
  double val;

  /* read bond constraints and build interaction list */
  for(int n=0; n<natoms; ++n)
     for(int m=0; m<natoms; ++m)
        inter[n*natoms+m] = true;

  for(int n=0; n<natoms; ++n)
     inter[n*natoms+n] = false;

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

     inter[indb*natoms + inda] = false;
     inter[inda*natoms + indb] = false;
  } 

  /* read angle constraints */
  for(int t=0; t<cang->N; ++t)
  {
     getline(file,rs);
     std::stringstream ss(rs);
     std::vector<std::string> strs;
     while(ss >> buf)
        strs.push_back(buf); 

     inda = atoi(strs[0].c_str());
     indb = atoi(strs[1].c_str());
     indc = atoi(strs[2].c_str());
     val  = atof(strs[3].c_str());

     cang->inda[t] = inda;
     cang->indb[t] = indb;
     cang->indc[t] = indc;
     cang->val[t]  = val;

     /* atoms inda, indb and indc are connected => remove 
        interactions */
     inter[inda*natoms + indb] = false;
     inter[inda*natoms + indc] = false;
     inter[indb*natoms + indc] = false;
     inter[indb*natoms + inda] = false;
     inter[indc*natoms + inda] = false;
     inter[indc*natoms + indb] = false;

     inter[indb*natoms + inda] = false;
     inter[indc*natoms + inda] = false;
     inter[indc*natoms + indb] = false;
     inter[inda*natoms + indb] = false;
     inter[inda*natoms + indc] = false;
     inter[indb*natoms + indc] = false;
  } 

  return;
}
