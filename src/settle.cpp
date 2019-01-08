#include <cmath>
#include <iostream>
#include <fstream>
#include "settle.hpp"

Settle::Settle(double dOH, double dHH, int nmol, double dt):
               nmol(nmol)
{
/* 
   Set parameteres for Settle algorithm:

   dOH - OH distance in A
   dHH - HH distance in A

   January 2019
*/    
   invdt = 1.0/dt;
 
   rc  = dHH/2.0;
   rc2 = dHH;

   mass = mO + 2.0*mH;
   
   ra = 2.0*mH*sqrt(dOH*dOH - rc*rc)/mass;
   rb = sqrt(dOH*dOH - rc*rc) - ra;

   wo = mO/(mO + 2.0*mH);
   wh = mH/(mO + 2.0*mH);

}

void Settle::settle1(double *xyzb, double *xyza, double *vel)
{
/*
   Main routine for SETTLE algorithm
   for a single 3-site water molecule

   Ref.: Miyamoto and Kollman, 
         J. Comput. Chem. 13, 952 (1992)

*/
   int cid;

   double xb0, yb0, zb0, xc0, yc0, zc0;
   double xcom, ycom, zcom;
   double xa1, ya1, za1, xb1, yb1, zb1,xc1,yc1, zc1;
   double xakszd, yakszd, zakszd, xaksxd, yaksxd, zaksxd, xaksyd, yaksyd, zaksyd; 
   double axlng, aylng, azlng;
   double trns11, trns21, trns31, trns12, trns22, trns32, trns13, trns23, trns33;
   double xb0d, yb0d, xc0d, yc0d, xa1d, ya1d, za1d, xb1d, yb1d, zb1d, xc1d, yc1d, zc1d;
   double ya2d, xb2d, t1, t2, yb2d, yc2d;
   double alpha, beta, gamma, al2be2, tmp2, sinthe, costhe;
   double xa3d, ya3d, za3d, xb3d, yb3d, zb3d, xc3d, yc3d, zc3d;
   double xa3, ya3, za3, xb3, yb3, zb3, xc3, yc3, zc3;
   double dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz;

   double sinphi, cosphi, sinpsi, cospsi, temp;
   double tmp;

   bool shake = false;

   //std::cout << ra << " " << rb << " " << rc << " " << rc2 << std::endl;

   ///*-- temporary --*/
   //xyzb[0] = 21.69995;
   //xyzb[1] = 11.80047;
   //xyzb[2] = 3.19961;
   //xyzb[3] = 22.38667;
   //xyzb[4] = 12.47317;
   //xyzb[5] = 2.92414;
   //xyzb[6] = 22.10420;
   //xyzb[7] = 11.15941;
   //xyzb[8] = 3.85201;

   //vel[0] = -1.72800;
   //vel[1] = 0.32400;
   //vel[2] = 0.17400;
   //vel[3] = 10.31800;
   //vel[4] = -11.26300;
   //vel[5] = 2.53900;
   //vel[6] = -16.65600;
   //vel[7] = -2.09400;
   //vel[8] = 6.80200;

   //xyza[0] = 21.69822;
   //xyza[1] = 11.80079;
   //xyza[2] = 3.19979;
   //xyza[3] = 22.39699;
   //xyza[4] = 12.46190;
   //xyza[5] = 2.92668;
   //xyza[6] = 22.08754;
   //xyza[7] = 11.15731;
   //xyza[8] = 3.85881;

   for(int n=0; n<nmol; ++n){
      /* run settle for each water molecule
         correcting positions and velocities

         1. for a configuration before MD step,
         determine OH vectors, here will be
         called a0 and b0 
      */
      //std::cout << " Before array top of SETTLE " << std::endl;
      //std::cout << xyzb[0] << " " << xyzb[1] << " " << xyzb[2] << std::endl;
      //std::cout << xyzb[3] << " " << xyzb[4] << " " << xyzb[5] << std::endl;
      //std::cout << xyzb[6] << " " << xyzb[7] << " " << xyzb[8] << std::endl;

      //std::cout << " After array top of SETTLE " << std::endl; 
      //std::cout << xyza[0] << " " << xyza[1] << " " << xyza[2] << std::endl;
      //std::cout << xyza[3] << " " << xyza[4] << " " << xyza[5] << std::endl;
      //std::cout << xyza[6] << " " << xyza[7] << " " << xyza[8] << std::endl;

      cid = 9*n;

      xb0 = xyzb[cid+3] - xyzb[cid];  // !!!! correct these indices to start from each molecule...
      yb0 = xyzb[cid+4] - xyzb[cid+1];
      zb0 = xyzb[cid+5] - xyzb[cid+2];
      xc0 = xyzb[cid+6] - xyzb[cid]; 
      yc0 = xyzb[cid+7] - xyzb[cid+1];
      zc0 = xyzb[cid+8] - xyzb[cid+2];

      /* for a configuration after MD step,
         determine center-of-mass
      */
      xcom = xyza[0]*wo + (xyza[3] + xyza[6])*wh; 
      ycom = xyza[1]*wo + (xyza[4] + xyza[7])*wh;
      zcom = xyza[2]*wo + (xyza[5] + xyza[8])*wh;

      /*
         subtract COM for configurations after MD step
      */
      xa1 = xyza[0] - xcom;
      ya1 = xyza[1] - ycom;
      za1 = xyza[2] - zcom;
      xb1 = xyza[3] - xcom;
      yb1 = xyza[4] - ycom;
      zb1 = xyza[5] - zcom;
      xc1 = xyza[6] - xcom;
      yc1 = xyza[7] - ycom;
      zc1 = xyza[8] - zcom;

      //std::cout << " xa1 - zc1 array " << std::endl;
      //std::cout << xa1 << " " << ya1 << " " << za1 << std::endl;
      //std::cout << xb1 << " " << yb1 << " " << zb1 << std::endl;
      //std::cout << xc1 << " " << yc1 << " " << zc1 << std::endl;

      /* through series of cross products obtain 
        vectors defining new coordinate systems */
      xakszd = yb0 * zc0 - zb0 * yc0;
      yakszd = zb0 * xc0 - xb0 * zc0;
      zakszd = xb0 * yc0 - yb0 * xc0;
      xaksxd = ya1 * zakszd - za1 * yakszd;
      yaksxd = za1 * xakszd - xa1 * zakszd;
      zaksxd = xa1 * yakszd - ya1 * xakszd;
      xaksyd = yakszd * zaksxd - zakszd * yaksxd;
      yaksyd = zakszd * xaksxd - xakszd * zaksxd;
      zaksyd = xakszd * yaksxd - yakszd * xaksxd;

      /* normalize basis vectors */
      axlng = sqrt(xaksxd*xaksxd + yaksxd*yaksxd + zaksxd*zaksxd);
      aylng = sqrt(xaksyd*xaksyd + yaksyd*yaksyd + zaksyd*zaksyd);
      azlng = sqrt(xakszd*xakszd + yakszd*yakszd + zakszd*zakszd);

      /* calculate unit vectors in each direction, which will serve as
         transformation vectors between two coordinate systems */
      trns11 = xaksxd / axlng;
      trns21 = yaksxd / axlng;
      trns31 = zaksxd / axlng;
      trns12 = xaksyd / aylng;
      trns22 = yaksyd / aylng;
      trns32 = zaksyd / aylng;
      trns13 = xakszd / azlng;
      trns23 = yakszd / azlng;
      trns33 = zakszd / azlng;

      //std::cout << " Transformation matrix " << std::endl;
      //std::cout << trns11 << " " << trns12 << " " << trns13 << std::endl;
      //std::cout << trns21 << " " << trns22 << " " << trns23 << std::endl;
      //std::cout << trns31 << " " << trns32 << " " << trns33 << std::endl;

     /* transform to a new coordinate system (primed X'Y'Z' in the paper) */
      xb0d = trns11 * xb0 + trns21 * yb0 + trns31 * zb0;
      yb0d = trns12 * xb0 + trns22 * yb0 + trns32 * zb0;
      xc0d = trns11 * xc0 + trns21 * yc0 + trns31 * zc0;
      yc0d = trns12 * xc0 + trns22 * yc0 + trns32 * zc0;
      xa1d = trns11 * xa1 + trns21 * ya1 + trns31 * za1;
      ya1d = trns12 * xa1 + trns22 * ya1 + trns32 * za1;
      za1d = trns13 * xa1 + trns23 * ya1 + trns33 * za1;
      xb1d = trns11 * xb1 + trns21 * yb1 + trns31 * zb1;
      yb1d = trns12 * xb1 + trns22 * yb1 + trns32 * zb1;
      zb1d = trns13 * xb1 + trns23 * yb1 + trns33 * zb1;
      xc1d = trns11 * xc1 + trns21 * yc1 + trns31 * zc1;
      yc1d = trns12 * xc1 + trns22 * yc1 + trns32 * zc1;
      zc1d = trns13 * xc1 + trns23 * yc1 + trns33 * zc1;

      //std::cout << " xb0d = " << xb0d << " yb0d = " << yb0d << " xc0d = " << xc0d << " yc0d = " << yc0d << std::endl;
      //std::cout << " xa1d = " << xa1d << " ya1d = " << ya1d << " za1d = " << za1d << " xb1d = " << xb1d << std::endl;
      //std::cout << " yb1d = " << yb1d << " zb1d = " << zb1d << " xc1d = " << xc1d << " yc1d = " << yc1d << " zc1d = " << zc1d << std::endl;

      sinphi = za1d/ra;
      tmp = 1.0 - sinphi*sinphi;
      if(tmp <= 0.0){
         std::cout << " Problem in SETTLE " << std::endl;
         exit(EXIT_FAILURE);
         cosphi = 0.0;
         shake = true;
         // actually use shake just like in gromacs
      }else
        cosphi = sqrt(tmp); 

      sinpsi = (zb1d - zc1d)/(rc2 * cosphi);
      tmp2 = 1.0 - sinpsi*sinpsi;
      if(tmp2 <= 0){
         std::cout << " Problem in SETTLE " << std::endl;
         exit(EXIT_FAILURE);
         cospsi = 0.0;
         shake = true;
         // actually use shake just like in gromacs
      }else
         cospsi = sqrt(tmp2); 

      //std::cout << " sinphi = " << sinphi << " cosphi = " << cosphi << " sinpsi = " << sinpsi << " cospsi = " << cospsi << std::endl;

      if(!shake){
         /* continue normal execution here */
         ya2d = ra * cosphi;
         xb2d = -rc * cospsi;
         t1 = -rb * cosphi;
         t2 = rc * sinpsi * sinphi;
         yb2d = t1 - t2;
         yc2d = t1 + t2; 

         alpha = xb2d * (xb0d - xc0d) + yb0d * yb2d + yc0d * yc2d;
         beta  = xb2d * (yc0d - yb0d) + xb0d * yb2d + xc0d * yc2d;
         gamma = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d;
         al2be2 = alpha * alpha + beta * beta;
         tmp2 = (al2be2 - gamma * gamma);
         sinthe = (alpha * gamma - beta * sqrt(tmp2)) / al2be2;

         //std::cout << " alpha = " << alpha << " beta = " << beta << " gamma = " << gamma << " al2be2 = " << al2be2 << " sinthe = " << sinthe << std::endl;

         tmp2 = 1.0 - sinthe*sinthe;
         costhe = sqrt(tmp2);

         xa3d = -ya2d * sinthe;
         ya3d =  ya2d * costhe;
         za3d =  za1d;
         xb3d =  xb2d * costhe - yb2d * sinthe;
         yb3d =  xb2d * sinthe + yb2d * costhe;
         zb3d =  zb1d;
         xc3d = -xb2d * costhe - yc2d * sinthe;
         yc3d = -xb2d * sinthe + yc2d * costhe;
         zc3d =  zc1d;

         xa3 = trns11 * xa3d + trns12 * ya3d + trns13 * za3d;
         ya3 = trns21 * xa3d + trns22 * ya3d + trns23 * za3d;
         za3 = trns31 * xa3d + trns32 * ya3d + trns33 * za3d;
         xb3 = trns11 * xb3d + trns12 * yb3d + trns13 * zb3d;
         yb3 = trns21 * xb3d + trns22 * yb3d + trns23 * zb3d;
         zb3 = trns31 * xb3d + trns32 * yb3d + trns33 * zb3d;
         xc3 = trns11 * xc3d + trns12 * yc3d + trns13 * zc3d;
         yc3 = trns21 * xc3d + trns22 * yc3d + trns23 * zc3d;
         zc3 = trns31 * xc3d + trns32 * yc3d + trns33 * zc3d;
 
         //std::cout << " xa3 - zc3 vectors " << std::endl;
         //std::cout << xa3 << " " << ya3 << " " << za3 << std::endl;
         //std::cout << xb3 << " " << yb3 << " " << zb3 << std::endl;
         //std::cout << xc3 << " " << yc3 << " " << zc3 << std::endl;

         xyza[0] = xcom + xa3;
         xyza[1] = ycom + ya3;
         xyza[2] = zcom + za3;
         xyza[3] = xcom + xb3;
         xyza[4] = ycom + yb3;
         xyza[5] = zcom + zb3;
         xyza[6] = xcom + xc3;
         xyza[7] = ycom + yc3;
         xyza[8] = zcom + zc3;

         dax = xa3 - xa1;
         day = ya3 - ya1;
         daz = za3 - za1;
         dbx = xb3 - xb1;
         dby = yb3 - yb1;
         dbz = zb3 - zb1;
         dcx = xc3 - xc1;
         dcy = yc3 - yc1;
         dcz = zc3 - zc1;
 
         vel[0] += dax * invdt;
         vel[1] += day * invdt;
         vel[2] += daz * invdt;
         vel[3] += dbx * invdt;
         vel[4] += dby * invdt;
         vel[5] += dbz * invdt;
         vel[6] += dcx * invdt;
         vel[7] += dcy * invdt;
         vel[8] += dcz * invdt;

         //std::cout << " After SETTLE " << std::endl;
         //std::cout << xyza[0] << " " << xyza[1] << " " << xyza[2] << std::endl;
         //std::cout << xyza[3] << " " << xyza[4] << " " << xyza[5] << std::endl;
         //std::cout << xyza[6] << " " << xyza[7] << " " << xyza[8] << std::endl;

         //std::cout << " Velocity " << std::endl;
         //std::cout << vel[0] << " " << vel[1] << " " << vel[2] << std::endl;
         //std::cout << vel[3] << " " << vel[4] << " " << vel[5] << std::endl;
         //std::cout << vel[6] << " " << vel[7] << " " << vel[8] << std::endl;
     }else{
        /* do shake if settle does not work */
        //shake_w();
     } 


   } /* end loop over molecules */
  return; 
}







