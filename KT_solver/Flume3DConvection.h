#ifndef FLUME3DCONV_H
#define FLUME3DCONV_H

#include "Flux.h"

class Flume3DConvection: public Flux
{
 private:

  SField m_segregation_rate;
  VectorField m_velocity;

 public:

 Flume3DConvection(): Flux(3,1) {m_segregation_rate = 0; m_velocity.resize(3);}

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
      m_evaluated_flux_d[0] = u[0]*m_velocity[1];
      break;

    case 2:
      m_evaluated_flux_d[0] = u[0]*m_velocity[2] - m_segregation_rate*u[0]*(1-u[0]);
      break;

    default:
      throw std::invalid_argument("In FlumeFlux::evaluate dimension must be 0, 1 or 2");
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
      m_max_eigenvalue = m_velocity[1];
      break;

    case 2:
      m_max_eigenvalue = m_velocity[2] - m_segregation_rate*(1-2*u[0]);
      break;

    default:
      throw std::invalid_argument("In FlumeFlux::get_max_eigenvalue dimension must be 0, 1 or 2");
    }
    return m_max_eigenvalue;
  }

};

#endif
