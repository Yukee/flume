#include <iostream>
#include <fstream>
#include <math.h>
#include <limits.h>

// Solvers
#include "FD1Solver.h"
#include "timeSolver.h"
#include "EulerSolver.h"
//#include "RK2Solver.h"
#include "RK3Solver.h"

// Containers
#include "Vector.h"
#include "ScalarField.h"
#include "PrescribedField.h"
#include "Flume2DField.h"
#include "Equation.h"

// Functions
#include "Flux.h"
#include "ZeroFlux.h"
#include "Flume2DConvection.h"
#include "Flume2DSource.h"

#include "WriteVectorField.h"

// Flume pb
//#include "Flume3D.h"
#include "Flume2D.h"

using namespace std;
int main(int argc, char *argv[])
{
    // Prescribe the fluxes and source term, store them in an equation
    
   Flux *ptrCF = new Flume2DConvection();
   Flux *ptrDF = new ZeroFlux(2,1);
   Flux *ptrS = new Flume2DSource();
   Equation flume_2d_eq (ptrCF, ptrDF, ptrS);
   
   // Ask the user for the discretisation and timestep infos
   
   Vector<double> dx (2); Vector<double> xi (2); Vector<double> llc (2); double endtime; double timestep; double timebtwfiles; string filename; double sr;
   //cout << "Enter cell width dx:"; cin >> dx[0];
   dx[0] = 0.001;
   
   //cout << "Enter cell heigh dz:"; cin >> dx[1];
   dx[1] = 0.001;
   
   //cout << "Enter domain width:"; cin >> xi[0];
   xi[0] = 4.5;
   xi[1] = 0.25; llc[0] = -xi[0]; llc[1] = 0;
   
   cout << "Enter end time:"; cin >> endtime;
   
   //cout << "Enter time step:"; cin >> timestep;
   timestep = 0.05;
   
   cout << "Enter time between two consecutive file saves:"; cin >> timebtwfiles;
   
   cout << "Enter name for saved file:"; cin >> filename;
   
   cout << "Enter value for the segregation rate:"; cin >> sr;

   FD1Solver *solv_ptr = new FD1Solver (dx, xi, &flume_2d_eq, llc);
   VectorField pos = solv_ptr->get_position();
   Vector<int> xr = solv_ptr->get_nxSteps();

   // Initialise value of the concentration in small particules, and velocity field

   VectorField phi (1, SField (xr));
   SField domain (xr);
   VectorField u0 (2, SField (xr));
   SField dvdy0 (xr);
   
   for(int it=0;it<phi[0].get_size();++it)
   {
    phi[0][it] = phi0(pos[0][it], pos[1][it]);
    domain[it] = boundary(pos[0][it], pos[1][it]);
    u0[0][it] = u(pos[0][it], pos[1][it]);
    u0[1][it] = w(pos[0][it], pos[1][it]);
    dvdy0[it] = dvdy(pos[0][it], pos[1][it]);
   }
   
   phi = domain*phi; u0 = domain*u0; dvdy0 = domain*dvdy0;

   // Set the west and south boundary conditions
   
   ScalarField phiWest (xr.drop(0)); 
   for(int it=0;it<phiWest.get_size();++it) phiWest[it] = phi0( llc[0]-dx[0] , llc[1]+it*dx[1] );
   phi[0].set_bound(0,-1,phiWest);

   ScalarField phiSouth (xr.drop(1)); 
   for(int it=0;it<phiSouth.get_size();++it) phiSouth[it] = phi0( llc[0]+it*dx[0] , llc[1]-dx[1] );
   phi[0].set_bound(1,-1,phiSouth);

   // Sets the parameters of the convection flux: velocity field and segregation rate
   
   SField seg_rate = (sr)*domain; 
   ptrCF->set_parameter(seg_rate); 
   ptrCF->set_parameter(u0); 
   ptrS->set_parameter(dvdy0);

   // Checks

   fstream phiinit;
   phiinit.open("Results/Check/phiinit.tsv",ios::out);
   phi[0].write_in_file(phiinit, dx, llc);
   
   fstream uinit;
   uinit.open("Results/Check/uinit.tsv",ios::out);
   write_VectorField(u0, pos, uinit);

   VectorField df (2); df[0] = ptrCF->get_max_eigenvalue(phi,0); df[1] = ptrCF->get_max_eigenvalue(phi,1);
   fstream duinit;
   duinit.open("Results/Check/duinit.tsv",ios::out);
   write_VectorField(df, pos, duinit);

   // Timestepper
   
   RK3Solver tsolv (timestep, endtime, solv_ptr, phi);
   
   // Compute the time evolution
   
   write_flume_infos(filename, sr);

   //tsolv.get_solution(filename, timebtwfiles);

   return 0;
}
