#include "RK3Solver.h"
#include "FD1Solver.h"
#include <string>
#include <boost/lexical_cast.hpp>
#include "WriteVectorField.h"

using namespace std;

RK3Solver::RK3Solver(double deltaT, double T, FD1Solver *spatialSolver, VectorField initial_conditions) :
  timeSolver(deltaT, T, spatialSolver, initial_conditions)
{
  m_ntSteps = m_T/m_deltaT;
}

void RK3Solver::get_solution(string name, double dt)
{
  timeSolver::get_solution(name, dt);

  Vector<double> deltaX = m_spatialSolver->get_deltaX();
  Vector<double> lowerLeftCorner = m_spatialSolver->get_lowerLeftCorner();
  fstream data;
  string path;


  // used in write_VectorField
  VectorField pos = m_spatialSolver->get_position();

  // Stores the bound to reset it at the end
  ScalarField un_west = m_un[0].get_bound(0,-1);

  VectorField u;
  VectorField df;
  double testDeltaT;
  double newDeltaT = m_deltaT;
  double currenttime = 0;
  int writingCounter = 0;

  SField unity = m_spatialSolver->get_unity();

  for(currenttime=0;currenttime<=m_T;currenttime += newDeltaT)
    {

      newDeltaT = m_deltaT;
      df = m_spatialSolver->get_numerical_flux_gradient(m_un);
      testDeltaT = m_spatialSolver->check_CFL(newDeltaT);
      if(testDeltaT!=newDeltaT) newDeltaT = testDeltaT;

      u = m_un + (-newDeltaT*unity)*df;
      u[0].set_bound(0,-1,un_west);

      df = m_spatialSolver->get_numerical_flux_gradient(u);
      u = (3/4.*unity)*m_un + (1/4.*unity)*(u + (-newDeltaT*unity)*df);
      u[0].set_bound(0,-1,un_west);

      df = m_spatialSolver->get_numerical_flux_gradient(u);
      u = (1/3.*unity)*m_un + (2/3.*unity)*(u + (-newDeltaT*unity)*df);
      u[0].set_bound(0,-1,un_west);

      if( (int)(currenttime/dt) == writingCounter )
        {
	  cout << "writing file number " << writingCounter << endl;
	  path =  "Results/" + name + "_" + boost::lexical_cast<string>(writingCounter) + ".tsv";
	  data.open(path.c_str(), ios::out);
	  u[0].write_in_file(data, deltaX, lowerLeftCorner);
	  data.close();

	  fstream initb;
	  initb.open("Results/Test/phib.tsv",ios::out);
	  initb << u[0].get_bounds();
	  initb.close();

	  writingCounter++;
        }
      m_un = u;
    }
}


