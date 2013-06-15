#ifndef ZEROFLUX_H
#define ZEROFLUX_H

#include "Flux.h"

class ZeroFlux : public Flux
{
 public:

 ZeroFlux(int d, int i) : Flux(d,i){}

  inline VectorField evaluate(const VectorField & u, const int d)
  {
    for(int i=0;i<m_solved_dimensions;i++) m_evaluated_flux_d[i] = 0*u[i];
    return m_evaluated_flux_d;
  }

  inline SField get_max_eigenvalue(const VectorField & u, const int d)
  {
    m_max_eigenvalue = 0*u[0];
    return m_max_eigenvalue;
  }
};

#endif
