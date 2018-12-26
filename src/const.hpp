#ifndef _CONST_
#define _CONST_

/*

 Sizes of arrays of fixed bonds and angles

*/
#define nbond_max  5000
#define nang_max   5000

/* 
    1/(4*pi*eps) factor in Coulomb energy in kJ/(mol*A*e^2)
*/
#define f_electr 1389.35485    

/*
   convert from nm to A
*/
#define unitC    10.0          

/* 
   convert kinetic energy from a.u.*A^2/ps^2 to kJ/mol
   conversion factor = 1.6605e-27 * 1e24 * 1e-20 * 1e-3 / 6.022e23
*/
#define KE_convert 0.0100      

/* 
   Avogadro number
*/
#define Na 6.0221367e23

/*
   Boltzmann constant
*/
#define Kb 1.38064852e-23 

/* 
   Atomic masses
*/
#define mH  1.007825032
#define mO  15.994914620
#define mD  2.014102

/*
   atomic unit of mass
*/
#define aum 1.66054e-27

/*
  Electron mass in kg
*/
#define mass_e 9.10938356e-31

/*
  kJ/mol to Hartree
*/
#define kJmH 2625.5

/*
  Angstrom to bohr
*/
#define atob 1.88973

/*
  Hartree to wavenumbers
*/
#define Htocmi 219416.0


#endif
