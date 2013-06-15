#ifndef FLUMEFLUXNOVEL_H
#define FLUMEFLUXNOVEL_H

#include "Flux.h"

class Flume2DConvectionFluxNoVel: public Flux
{
 private:

  SField m_segregation_rate;
  VectorField m_velocity;

 public:

 Flume2DConvectionFluxNoVel(): Flux(2,1) {m_segregation_rate = 0; m_velocity.resize(2);}

  inline void set_parameter(const SField & sr)
  {
    m_segregation_rate = sr;
  }

  inline void set_parameter(const VectorField & vel)
  {
    m_velocity = vel;
  }

  inline virtual VectorField evaluate(const VectorField & u, const int d)
  {
    switch(d){
    case 0:
      m_evaluated_flux_d[0] = u[0]*m_velocity[0];
      break;

    case 1:
      m_evaluated_flux_d[0] = u[0]*m_velocity[1] - m_segregation_rate*u[0]*(1-u[0]);
      break;

    default:
      throw std::invalid_argument("In FlumeFlux::evaluate dimension must either be 0 or 1");
    }
    return m_evaluated_flux_d;
  }

  inline virtual SField get_max_eigenvalue(const VectorField & u, const int d)
  {
    switch(d){
    case 0:
      m_max_eigenvalue = m_velocity[0];
      break;

    case 1:
      m_max_eigenvalue = m_velocity[1] - m_segregation_rate*(1-2*u[0]);
      break;

    default:
      throw std::invalid_argument("In FlumeFlux::get_max_eigenvalue dimension must either be 0 or 1");
    }
    return m_max_eigenvalue;
  }

};

#endif
