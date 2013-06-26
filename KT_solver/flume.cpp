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
#include "Equation.h"

// Functions
#include "Flux.h"
#include "ZeroFlux.h"
#include "Flume3DConvection.h"

#include "WriteVectorField.h"

// Flume pb
#include "Flume3D.h"
//#include "Flume2D.h"

using namespace std;
int main(int argc, char *argv[])
{
    // Prescribe the fluxes and source term, store them in an equation
    
   Flux *ptrCF = new Flume3DConvection();
   Flux *ptrDF = new ZeroFlux(3,1);
   Flux *ptrS = new ZeroFlux(3,1);
   Equation flume_2d_eq (ptrCF, ptrDF, ptrS);
   
   // Ask the user for the discretisation and timestep infos
   
   Vector<double> dx (3); Vector<double> xi (3); Vector<double> llc (3); double endtime; double timestep; double timebtwfiles; string filename; double sr;
   //cout << "Enter dx:"; cin >> dx[0];
   dx[0] = 0.05;
   
   //cout << "Enter dy:"; cin >> dx[1];
   dx[1] = 0.05;

   //cout << "Enter dz:"; cin >> dx[2];
   dx[2] = 0.005;
   
   //cout << "Enter domain width:"; cin >> xi[0];
   xi[0] = 4.5;
   xi[1] = 1; 
   xi[2] = 0.25;
   llc[0] = -xi[0]; llc[1] = -0.5; llc[2] = 0;
   
   cout << "Enter end time:"; cin >> endtime;
   
   //cout << "Enter time step:"; cin >> timestep;
   timestep = 0.05;
   
   cout << "Enter time between two consecutive file saves:"; cin >> timebtwfiles;
   
   cout << "Enter name for saved file:"; cin >> filename;
   
   //cout << "Enter value for the segregation rate:"; cin >> sr;
   sr = 0.035;

   FD1Solver *solv_ptr = new FD1Solver (dx, xi, &flume_2d_eq, llc);
   VectorField pos = solv_ptr->get_position();
   Vector<int> xr = solv_ptr->get_nxSteps();

   // Initialise value of the concentration in small particules, and velocity field

   VectorField phi (1, SField (xr));
   SField domain (xr);
   VectorField u0 (3, SField (xr));
   
   for(int it=0;it<phi[0].get_size();++it)
   {
    phi[0][it] = phi0(pos[0][it], pos[1][it], pos[2][it]);
    domain[it] = boundary(pos[0][it], pos[1][it], pos[2][it]);
    u0[0][it] = u(pos[0][it], pos[1][it], pos[2][it]);
    u0[1][it] = w(pos[0][it], pos[1][it], pos[2][it]);
   }
   
   phi = domain*phi; u0 = domain*u0; 

   // Set the west and south boundary conditions
   
   ScalarField phiWest (xr.drop(0)); 
   Vector<int> x (2); Vector<int> y (2); x[0] = 1 ; x[1] = 0; y[0] = 0; y[1] = 1;
   for(int j=0;j<xr[1];j++) for(int k=0;k<xr[2];k++) phiWest(j*x+k*y) = phi0( llc[0]-dx[0], llc[1]+j*dx[1], llc[2]+k*dx[2] );
   phi[0].set_bound(0,-1,phiWest);

   ScalarField phiSouth (xr.drop(2));
   for(int i=0;i<xr[0];i++) for(int j=0;j<xr[1];j++) phiSouth(i*x+j*y) = phi0( llc[0]+i*dx[0], llc[1]+j*dx[1], llc[2]-dx[2] );
   phi[0].set_bound(2,-1,phiSouth);

   // Sets the parameters of the convection flux: velocity field and segregation rate
   
   SField seg_rate = (sr)*domain; 
   ptrCF->set_parameter(seg_rate); 
   ptrCF->set_parameter(u0); 

   // Checks

   fstream phiinit;
   phiinit.open("Results/Check/phiinit.tsv",ios::out);
   phi[0].write_in_file(phiinit, dx, llc);
   
   fstream uinit;
   uinit.open("Results/Check/uinit.tsv",ios::out);
   write_VectorField(u0, pos, uinit);

   VectorField df (3); df[0] = ptrCF->get_max_eigenvalue(phi,0); df[1] = ptrCF->get_max_eigenvalue(phi,1); df[2] = ptrCF->get_max_eigenvalue(phi,2);
   fstream duinit;
   duinit.open("Results/Check/duinit.tsv",ios::out);
   write_VectorField(df, pos, duinit);

   // Timestepper
   
   EulerSolver tsolv (timestep, endtime, solv_ptr, phi);
   
   // Compute the time evolution
   
   write_flume_infos(filename, sr);

   tsolv.get_solution(filename, timebtwfiles);

   return 0;
}
