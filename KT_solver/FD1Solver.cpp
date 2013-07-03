#include "FD1Solver.h"
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

FD1Solver::FD1Solver(Vector<double> deltaX, Vector<double> xInterval, Equation *eq, Vector<double> lowerLeftCorner) :
   m_deltaX(deltaX), m_xInterval(xInterval), m_eq(eq), m_lowerLeftCorner(lowerLeftCorner)
{
  // Gets the dimensions of the problem
  m_m = m_eq->get_solved_dimensions(); 
  m_n = m_eq->get_space_dimensions();

  // Gets the base vectors of the physical space
  m_b = ( Vector<int> (m_n) ).get_base_vectors(1,0);

  // Computes the dimensions of each finite-volume cell
  m_nxSteps.resize(m_n);
  for(int i=0;i<m_n;i++) m_nxSteps[i] = m_xInterval[i]/m_deltaX[i];

  // Computes the value of the position field at grid points
  m_position.resize(m_n);
  for(int dir=0;dir<m_n;dir++) m_position[dir].resize_field(m_nxSteps);
  if((int)lowerLeftCorner.size() != m_n) throw invalid_argument("lowerLeftCorner dimension is not equal to the space dimension");
  Vector<int> pIt;
  for(int dir=0;dir<m_n;dir++) 
    {
      for(int it=0;it<m_position[dir].get_size();++it)
	{
	  pIt = m_position[dir].get_pos(it);
	  m_position[dir][it] = m_lowerLeftCorner[dir] + m_deltaX[dir]*pIt[dir];
	}
    }

  // Sets the ranges of the various computation used fields
  upper_right_intermediate_un_values = TensorField (m_n, VectorField (m_m, SField (m_nxSteps ) ) );
  lower_right_intermediate_un_values = TensorField (m_n, VectorField (m_m, SField (m_nxSteps ) ) );
  upper_left_intermediate_un_values = TensorField (m_n, VectorField (m_m, SField (m_nxSteps ) ) );
  lower_left_intermediate_un_values = TensorField (m_n, VectorField (m_m, SField (m_nxSteps ) ) );

  right_convection_flux = TensorField (m_n, VectorField (m_m, SField (m_nxSteps) ) );
  left_convection_flux = TensorField (m_n, VectorField (m_m, SField (m_nxSteps) ) );
  right_diffusion_flux = TensorField (m_n, VectorField (m_m, SField (m_nxSteps) ) );
  left_diffusion_flux = TensorField (m_n, VectorField (m_m, SField (m_nxSteps) ) );

  source_term = VectorField (m_m, SField (m_nxSteps ) );

  un_derivatives = TensorField (m_n, VectorField (m_m, SField (m_nxSteps) ) );

  right_localSpeed = VectorField (m_n, SField (m_nxSteps ) );
  left_localSpeed = VectorField (m_n, SField (m_nxSteps ) );

  // caca
  unity.resize_field (m_nxSteps);
  unity = 1;
}

FD1Solver::~FD1Solver()
{
    if(m_eq) delete m_eq;
    m_eq = NULL;
}

double FD1Solver::check_CFL(double deltaT)
{
    double newDeltaT = deltaT;
    Vector<SField> waveSpeed(m_n);
    Vector<double> maxFrequency(m_n);
    double overallMaxFreq=0;

    for(int i=0;i<m_m;i++)
      {
	for(int dir=0;dir<m_n;dir++)
	  {
	    waveSpeed[dir] = right_localSpeed[dir].max_field(left_localSpeed[dir]);
	    maxFrequency[dir] = waveSpeed[dir].get_max()/m_deltaX[dir];
	    if(maxFrequency[dir]>overallMaxFreq) overallMaxFreq = maxFrequency[dir];
	  }
      }

    if(newDeltaT>1/(8*overallMaxFreq)) newDeltaT = 1/(9*overallMaxFreq);
    return newDeltaT;
}

double FD1Solver::minmod(double a, double b)
{
  double minmod;
  if(a*b<=0) minmod = 0;
  else
    {
      if(a>0) minmod = min(a,b);
      else minmod = max(a,b);
    }
  return minmod;
}

double FD1Solver::three_pts_derivative(int it, int dir, int i)
{
  Vector<int> j = m_un[i].get_pos(it);
  double deriv;
  
  deriv = minmod( (m_un[i][it] - m_un[i](j - m_b[dir]))/m_deltaX[dir], (m_un[i](j + m_b[dir]) - m_un[i][it])/m_deltaX[dir] );
  return deriv;
}// returns the derivative in direction dir for field m_un[i]

void FD1Solver::compute_un_derivatives()
{
  for(int d=0;d<m_n;d++)
    {
      for(int i=0;i<m_m;i++)
	{
	  for(int it=0;it<un_derivatives[d][i].get_size();++it)
	    {
	      un_derivatives[d][i][it] = three_pts_derivative(it,d,i);
	    }
	}
    }
}

void FD1Solver::compute_intermediate_un_values()
{
  Vector<int> j; Vector<int> jp; Vector<int> jm;

  for(int d=0;d<m_n;d++)
    {
      for(int i=0;i<m_m;i++)
	{
	  for(int it=0;it<lower_left_intermediate_un_values[d][i].get_size();++it)
	    {
	      j = lower_left_intermediate_un_values[d][i].get_pos(it);
	      jp = j + m_b[d];
	      jm = j - m_b[d];

	      //un-1/2-
	      lower_left_intermediate_un_values[d][i][it] = m_un[i](jm) + 0.5*m_deltaX[d]*un_derivatives[d][i](jm);
	      //un-1/2+
	      upper_left_intermediate_un_values[d][i][it] = m_un[i][it] - 0.5*m_deltaX[d]*un_derivatives[d][i][it];
	      //un+1/2-
	      lower_right_intermediate_un_values[d][i][it] = m_un[i][it] + 0.5*m_deltaX[d]*un_derivatives[d][i][it];
	      //un+1/2+
	      upper_right_intermediate_un_values[d][i][it] = m_un[i](jp) - 0.5*m_deltaX[d]*un_derivatives[d][i](jp);
	    }
	}
    }
} 

void FD1Solver::compute_localSpeed()
{
  SField lowerSpeed;
  SField upperSpeed;

  for(int dir=0;dir<m_n;dir++)
    {
      upperSpeed = m_eq->get_max_eigenvalue(upper_right_intermediate_un_values[dir], dir);
      lowerSpeed = m_eq->get_max_eigenvalue(lower_right_intermediate_un_values[dir], dir);
      right_localSpeed[dir] = upperSpeed.max_field(lowerSpeed);

      upperSpeed = m_eq->get_max_eigenvalue(upper_left_intermediate_un_values[dir], dir);
      lowerSpeed = m_eq->get_max_eigenvalue(lower_left_intermediate_un_values[dir], dir);
      left_localSpeed[dir] = upperSpeed.max_field(lowerSpeed);
    }
}

void FD1Solver::compute_numerical_convection_flux()
{
  //no particular bc on the flux? -> null flux at West boundary
  VectorField upperFlux;
  VectorField lowerFlux;
	
  for(int d=0;d<m_n;d++)
    {
      upperFlux = m_eq->get_convectionFlux(upper_right_intermediate_un_values[d], d);
      lowerFlux = m_eq->get_convectionFlux(lower_right_intermediate_un_values[d], d);
      right_convection_flux[d] = (0.5*unity)*(upperFlux + lowerFlux) + (-0.5*unity)*right_localSpeed[d]*(upper_right_intermediate_un_values[d] - lower_right_intermediate_un_values[d]);

      upperFlux = m_eq->get_convectionFlux(upper_left_intermediate_un_values[d], d);
      lowerFlux = m_eq->get_convectionFlux(lower_left_intermediate_un_values[d], d);
      left_convection_flux[d] = (0.5*unity)*(upperFlux + lowerFlux) + (-0.5*unity)*left_localSpeed[d]*(upper_left_intermediate_un_values[d] - lower_left_intermediate_un_values[d]);
    }

   // Null flux across south boundary
   for(int i=0;i<m_nxSteps[0];i++) for(int j=0;j<m_nxSteps[1];j++) left_convection_flux[2][0](i*m_b[0]+j*m_b[1]) = 0;
   // Null flux across y=-1, y=1 boundaries
   for(int j=0;j<m_nxSteps[1];j++) for(int k=0;k<m_nxSteps[2];k++) {
	left_convection_flux[1][0](j*m_b[1]+k*m_b[2]) = 0;
	right_convection_flux[1][0](j*m_b[1]+k*m_b[2]+(m_nxSteps[0]-1)*m_b[0]) = 0;
	}
}

void FD1Solver::compute_numerical_diffusion_flux()
{
  //no particular bc on the flux?
  for(int d=0;d<m_n;d++)
    {
      right_diffusion_flux[d] = m_eq->get_diffusionFlux(upper_right_intermediate_un_values[d], d);
      left_diffusion_flux[d] = m_eq->get_diffusionFlux(upper_left_intermediate_un_values[d], d);
    }

}

void FD1Solver::compute_source_term()
{
  source_term = m_eq->get_source_term(m_un);
}

void FD1Solver::compute_south_boundary()
{
  for(int i=0;i<m_nxSteps[0];i++) for(int j=0;j<m_nxSteps[1];j++) m_un[0](i*m_b[0]+j*m_b[1]-m_b[2]) = m_un[0](i*m_b[0]);
} // Copy un value at the south edge into the south boundary; this is specific to the flume pb

//the flux gradient may be infinite, if you take a too large time step.
VectorField FD1Solver::get_numerical_flux_gradient(const VectorField & un)
{
  m_un = un;

  // Copy un value at the south edge into the south boundary; this is specific to the flume pb
  compute_south_boundary();

  compute_un_derivatives();
  compute_intermediate_un_values();
  compute_localSpeed();
  compute_numerical_convection_flux();
  compute_numerical_diffusion_flux();
  compute_source_term();

  // adding everything to form the spatial part of the equation: {convective flux} - {diffusive flux} - {source term} 
  VectorField flux_gradient(m_m, SField (m_nxSteps));
  for(int d=0;d<m_n;d++)
    { 
      flux_gradient = flux_gradient 
	+ ((1./m_deltaX[d])*unity)*( right_convection_flux[d] - left_convection_flux[d]);
      // - right_diffusion_flux[d] + left_diffusion_flux[d] );
    }
  flux_gradient = flux_gradient - source_term;

  return flux_gradient;
}

