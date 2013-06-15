#ifndef EQUATION_H
#define	EQUATION_H

#include "Vector.h"
#include "Flux.h"
#include "ZeroFlux.h"

class Equation
{
protected:
  Flux *m_conv; //convection flux
  Flux *m_diff; //diffusion flux
  Flux *m_source; //source term

  int m_space_dimensions;
  int m_solved_dimensions;

public:
  Equation() {m_conv =  new ZeroFlux(1,1); m_diff = new ZeroFlux(1,1); m_source = new ZeroFlux(1,1);
    m_space_dimensions = m_conv->get_space_dimensions();
    m_solved_dimensions = m_conv->get_solved_dimensions();}

  Equation(Flux *f) : m_conv(f) {
    m_space_dimensions = m_conv->get_space_dimensions();
    m_solved_dimensions = m_conv->get_solved_dimensions();
    m_diff = new ZeroFlux(m_space_dimensions,m_solved_dimensions);
    m_source = new ZeroFlux(m_space_dimensions,m_solved_dimensions);}

  Equation(Flux *f, Flux *diff) : m_conv(f), m_diff(diff){
    m_space_dimensions = m_conv->get_space_dimensions();
    m_solved_dimensions = m_conv->get_solved_dimensions();
    if( m_space_dimensions != m_diff->get_space_dimensions() || m_solved_dimensions != m_diff->get_solved_dimensions())
      throw std::invalid_argument("In Equation(Flux*,Flux*) convective and diffusive fluxes have different dimensions");}

 Equation(Flux *f, Flux *diff, Flux *s) : m_conv(f), m_diff(diff), m_source(s){
    m_space_dimensions = m_conv->get_space_dimensions();
    m_solved_dimensions = m_conv->get_solved_dimensions();
    if( m_space_dimensions != m_diff->get_space_dimensions() || m_solved_dimensions != m_diff->get_solved_dimensions() || m_space_dimensions != m_source->get_space_dimensions() || m_solved_dimensions != m_source->get_solved_dimensions() )
      throw std::invalid_argument("In Equation(Flux*,Flux*,Flux*) convective flux, diffusive flux and source term  have different dimensions");}


  Equation(const Equation &);

  ~Equation();

  VectorField get_convectionFlux(const VectorField &, const int);
  VectorField get_diffusionFlux(const VectorField &, const int);
  SField get_max_eigenvalue(const VectorField &, const int);
  VectorField get_source_term(const VectorField &);

  int get_space_dimensions();
  int get_solved_dimensions();
};

#endif


