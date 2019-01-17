#include <string>
#include <cstdlib>
#include <iostream>
#include <cmath> 
#include <fstream>
#include <vector>
#include <sstream>
#include "system.hpp"

System::System(int natoms, std::string data_file, double rcut, double dt):
               natoms(natoms), data_file(data_file), rcut(rcut), dt(dt)
{
/*
    System's initializer

    December 2018 - ...

*/

 /* set some constants */
 shake_iter_max = 500; 
 shake_tol = 1.0e-7;

 dt2 = dt*dt;
 
 /* allocate arrays */
 xyz   = (double *)calloc(3*natoms, sizeof(double));
 xyz_s = (double *)calloc(3*natoms, sizeof(double));
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

 printe();

 printf(" Step 1 \n");
 vv_vel();
 vv_xyz();
 potential();
 vv_vel();
 kinetic();
 printe();
 save();
 shake();
 kinetic();
 potential();
 printe();

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

  /* allocate arrays for SHAKE */
  dxt   = (double *)calloc(ncons, sizeof(double));
  dyt   = (double *)calloc(ncons, sizeof(double));
  dzt   = (double *)calloc(ncons, sizeof(double));
  dxx   = (double *)calloc(ncons, sizeof(double));
  dyy   = (double *)calloc(ncons, sizeof(double));
  dzz   = (double *)calloc(ncons, sizeof(double));

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

   printf(" Lennard-Jones interaction pairs: %d \n",imat->nLJ);

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

   printf(" Coulomb interaction pairs: %d \n",imat->nC);

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
        corr = (xyz[3*n+d] > box[d]/2.0) ? -box[d] : 0.0;
        xyz[3*n+d]   += corr;
        corr = (xyz[3*n+d] < -box[d]/2.0) ? box[d] : 0.0;
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

}

void System::potential()
{
/*
   Here we will evaluate potential energy and forces

   January 2019
*/


   for(int s=0; s<(3*natoms); ++s){
      force[s] = 0.0;
      fc[s]    = 0.0;
      flj[s]   = 0.0;
   }

   int i, j;
   double ix, iy, iz, jx, jy, jz, dx, dy, dz, r2, r;
   double e, s, er, er3, er6, er12, sr, sr3, sr6;
   double erc, erc3, erc6, erc12, src, src3, src6;
   double ux, uy, uz, gamma;
    
   /* Lennard-Jones truncated potential */
   ELJ = 0.0;
   for(int n=0; n<imat->nLJ; ++n){

      i = imat->LJa[n];
      j = imat->LJb[n];

      ix = xyz[3*i];
      iy = xyz[3*i+1];
      iz = xyz[3*i+2];

      jx = xyz[3*j];
      jy = xyz[3*j+1];
      jz = xyz[3*j+2];

      dx = minImage(ix - jx, box[0]);
      dy = minImage(iy - jy, box[1]);
      dz = minImage(iz - jz, box[2]);

      r2 = dx*dx + dy*dy + dz*dz;
      r  = sqrt(r2);

      if(r < rcut){

         e = imat->eps[n];
         s = imat->sig[n];

         er = e/r;
         er3 = er*er*er;
         er6 = er3*er3;
         er12 = er6*er6;

         sr = s/r;
         sr3 = sr*sr*sr;
         sr6 = sr3*sr3;

         erc = e/rcut;
         erc3 = erc*erc*erc;
         erc6 = erc3*erc3;
         erc12 = erc6*erc6;

         src = s/rcut;
         src3 = src*src*src;
         src6 = src3*src3; 

         ELJ += er12 - sr6 - erc12 + src6;

         ux = dx/r;
         uy = dy/r;
         uz = dz/r;

         gamma = (12.0*er12 - 6.0*sr6)/r;

         flj[3*i]   += gamma*ux;
         flj[3*i+1] += gamma*uy;
         flj[3*i+2] += gamma*uz;

         flj[3*j]   -= gamma*ux;
         flj[3*j+1] -= gamma*uy;
         flj[3*j+2] -= gamma*uz;
      } 
   }

   double qq; 
   double alpha = 1.0/(rcut*rcut);
   double beta  = -1.0/rcut;

   /* Coulomb energy/forces */
   Ec  = 0.0;
   for(int n=0; n<imat->nC; ++n){

      i = imat->Ca[n];
      j = imat->Cb[n];

      ix = xyz[3*i];
      iy = xyz[3*i+1];
      iz = xyz[3*i+2];

      jx = xyz[3*j];
      jy = xyz[3*j+1];
      jz = xyz[3*j+2];

      dx = minImage(ix - jx, box[0]);
      dy = minImage(iy - jy, box[1]);
      dz = minImage(iz - jz, box[2]);

      r2 = dx*dx + dy*dy + dz*dz;
      r  = sqrt(r2);

      if(r < rcut){

         qq = imat->vij[n];

         Ec += qq*(1.0/r + alpha*(r - rcut) + beta);

         ux = dx/r;
         uy = dy/r;
         uz = dz/r;

         gamma = qq*(1.0/r2 - alpha);
 
         fc[3*i]   += gamma*ux;
         fc[3*i+1] += gamma*uy;
         fc[3*i+2] += gamma*uz;

         fc[3*j]   -= gamma*ux;
         fc[3*j+1] -= gamma*uy;
         fc[3*j+2] -= gamma*uz;         
      }
   }

   /* Calculate total energy */
   Ep = ELJ + Ec;
   Et = Ek + Ep;

   /* check net force */
   double fx = 0.0;
   double fy = 0.0;
   double fz = 0.0;

   for(int s=0; s<natoms; ++s){
     fx += flj[3*s];
     fy += flj[3*s+1];
     fz += flj[3*s+2];
   }

   if(abs(fx)>1.0e-8 || abs(fy)>1.0e-8 || abs(fz)>1.0e-8){
      printf(" Net LJ force is not zero: %7.5f %7.5f %7.5f \n",fx,fy,fz);
      exit(EXIT_FAILURE);
   }

   fx = 0.0;
   fy = 0.0;
   fz = 0.0;

   for(int s=0; s<natoms; ++s){
      fx += fc[3*s];
      fy += fc[3*s+1];
      fz += fc[3*s+2];
   }

   if(abs(fx)>1.0e-8 || abs(fy)>1.0e-8 || abs(fz)>1.0e-8){
      printf(" Net Coulomb force is not zero: %7.5f %7.5f %7.5f \n",fx,fy,fz);
      exit(EXIT_FAILURE);
   }

   /* write forces to external file */
   ffile = fopen("Force.txt","w");
   fprintf(ffile,"# Fc_x \t Fc_y \t Fc_z \t Flj_x \t Flj_y \t Flj_z \n");
   for(int n=0; n<natoms; ++n)
      fprintf(ffile, " %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f \n",
             fc[3*n],  fc[3*n+1],  fc[3*n+2], 
             flj[3*n], flj[3*n+1], flj[3*n+2]);
   fclose(ffile);
   printf(" All forces are saved to: Force.txt \n");

   /* add all forces together */
   for(int s=0; s<(3*natoms); ++s)
      force[s] = fc[s] + flj[s];

}

void System::printe()
{
/*
    Print system data: energies, temperature, etc.

    January 2019
*/
   printf("\n------------------- System -------------------\n");
   printf(" Kinetic energy       = %9.2f  [ kJ/mol ] \n",Ek);
   printf(" Lennard-Jones energy = %9.2f  [ kJ/mol ] \n",ELJ);
   printf(" Coulomb energy       = %9.2f  [ kJ/mol ] \n",Ec);
   printf(" Potential energy     = %9.2f  [ kJ/mol ] \n",Ep);
   printf(" Total energy         = %9.2f  [ kJ/mol ] \n",Et);
   printf(" Temperature          = %9.2f  [ K ] \n",Tk);
   printf(" Density              = %9.2f  [ kg/m^3 ] \n ",rho);
   printf("----------------------------------------------\n");
}

void System::vv_vel(){
/*
   Velocity-Verlet algorithm part I.
   v(t+dt/2) = v(t) + dt*F(t)/(2*m)
  
   January 2019
*/
  for(int n=0; n<natoms; ++n){
     vel[3*n  ] += 100.0*dt*force[3*n  ]/(2.0*mass[n]);   
     vel[3*n+1] += 100.0*dt*force[3*n+1]/(2.0*mass[n]);
     vel[3*n+2] += 100.0*dt*force[3*n+2]/(2.0*mass[n]);
  }
}

void System::vv_xyz(){
/*
   Velocity-Verlet algorithm part II.
   x(t+dt) = x(t) + dt*v(t+dt/2)

   January 2019
*/
  for(int n=0; n<natoms; ++n){
     xyz[3*n  ] += dt*vel[3*n  ];
     xyz[3*n+1] += dt*vel[3*n+1];
     xyz[3*n+2] += dt*vel[3*n+2];
  }

  /* correct for box */
  inbox();
}

void System::save(){
/*
   Write coordinates and velocities into a file
*/
  sfile = fopen("xyz.txt","w");
  for(int n=0; n<natoms; ++n)
     fprintf(sfile, " %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f \n",
             xyz[3*n],  xyz[3*n+1],  xyz[3*n+2],
             vel[3*n],  vel[3*n+1],  vel[3*n+2]);
   fclose(sfile);
   printf(" All coordinates and velocities are saved to: xyz.txt \n");
}

void System::shake(){
/*
   Implement SHAKE algorithm for constraints here

   January 2019
*/
   int i, j, k, nsteps;
   double esig, esig1, dx, dy, dz;
   double dis, omega, amti, amtj, gamma, gammi, gammj;
   bool converged = false;

   /* save coordinates */
   for(int s=0; s<(3*natoms); ++s)
      xyz_s[s] = xyz[s];

   /* calculate distances before SHAKE */
   for(int c=0; c<cbond->N; ++c){
      i = cbond->inda[c];
      j = cbond->indb[c];
      dxx[c] = minImage(xyz[3*i] - xyz[3*j],     box[0]);
      dyy[c] = minImage(xyz[3*i+1] - xyz[3*j+1], box[1]);
      dzz[c] = minImage(xyz[3*i+2] - xyz[3*j+2], box[2]);
   }

   for(int iter=0; iter<shake_iter_max; ++iter){   
  
      for(int c=0; c<cbond->N; ++c){
         i = cbond->inda[c];
         j = cbond->indb[c];
         dxt[c] = minImage(xyz[3*i] - xyz[3*j],     box[0]);        
         dyt[c] = minImage(xyz[3*i+1] - xyz[3*j+1], box[1]);
         dzt[c] = minImage(xyz[3*i+2] - xyz[3*j+2], box[2]);
      }

      esig = 0.0;

      for(int c=0; c<cbond->N; ++c){
         dx = dxt[c];
         dy = dyt[c];
         dz = dzt[c];
          
         dis = cbond->val[c];
         esig1 = fabs(dx*dx + dy*dy + dz*dz - dis*dis);

         /* maximum error at this point */
         esig = (esig > esig1) ? esig : esig1 ;
      }

      if (esig < shake_tol)
         converged = true;

      if(converged)
         break;
      else{
 
         for(int c=0; c<cbond->N; ++c){        
            i = cbond->inda[c];
            j = cbond->indb[c];
            dis = cbond->val[c];

            omega = dis*dis;
            amti  = 1.0/mass[i];
            amtj  = 1.0/mass[j];

            dx = minImage(xyz[3*i] - xyz[3*j],     box[0]);
            dy = minImage(xyz[3*i+1] - xyz[3*j+1], box[1]);
            dz = minImage(xyz[3*i+2] - xyz[3*j+2], box[2]);

            gamma = (dx*dx + dy*dy + dz*dz - omega)/
                    (2.0*(amti + amtj)*(dxx[c]*dx + dyy[c]*dy + dzz[c]*dz));

            gammi = gamma*amti;
            xyz[3*i]   = xyz[3*i]   - dxx[c]*gammi;
            xyz[3*i+1] = xyz[3*i+1] - dyy[c]*gammi;
            xyz[3*i+2] = xyz[3*i+2] - dzz[c]*gammi;

            gammj = gamma*amtj;
            xyz[3*i]   = xyz[3*i]   + dxx[c]*gammj;
            xyz[3*i+1] = xyz[3*i+1] + dyy[c]*gammj;
            xyz[3*i+2] = xyz[3*i+2] + dzz[c]*gammj;
         } 
      }

      nsteps = iter + 1;
   }

   if(converged){
      printf(" SHAKE converged in %d steps to accuracy %f \n",nsteps,esig);
   }else{
      printf(" SHAKE did not converged! Nsteps = %d accuracy = %f tolerance = %f \n",
             nsteps,esig,shake_tol);
      exit(EXIT_FAILURE);
   } 

   /* correct for the box */
   inbox();

   /* correct velocities */
   for(int s=0; s<(3*natoms); ++s)
      vel[s] += (xyz[s] - xyz_s[s])/dt;

}
