#ifndef FLUME2DDIFF_H
#define FLUME2DDIFF_H

#include "Flux.h"

class Flume2DSource: public Flux
{
 private:

  // transverse velocity in the centre plane
  SField m_dvdy;

 public:

 Flume2DSource(): Flux(2,1) {}

  inline void set_parameter(const SField & dvdy)
  {
    m_dvdy = dvdy;
  }

  inline virtual VectorField evaluate(const VectorField & u, const int d)
  {
    switch(d){
    case 0:
      m_evaluated_flux_d[0] = (-1)*u[0]*m_dvdy;
      break;

    default:
      throw std::invalid_argument("In Flume2DSource::evaluate dimension must either be 0");
    }
    return m_evaluated_flux_d;
  }

};

#endif
