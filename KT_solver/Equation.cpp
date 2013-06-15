#include "Equation.h"

Equation::Equation(const Equation & eq)
{
  std::cout << "Copy constructor is called." << std::endl;
  *m_conv = *eq.m_conv;
  *m_diff = *eq.m_diff;
  *m_source = *eq.m_source;
}

Equation::~Equation()
{
  if(m_conv) delete m_conv;
  if(m_diff) delete m_diff;
  if(m_source) delete m_source;
}

VectorField Equation::get_convectionFlux(const VectorField & u, const int d)
{
  return m_conv->evaluate(u,d);
}

VectorField Equation::get_diffusionFlux(const VectorField & u,const int d)
{
  return m_diff->evaluate(u,d);
}

SField Equation::get_max_eigenvalue(const VectorField & u, const int d)
{
  return m_conv->get_max_eigenvalue(u,d);
}

VectorField Equation::get_source_term(const VectorField & u)
{
  return m_source->evaluate(u,0);
}

int Equation::get_space_dimensions()
{
  return m_space_dimensions;
}

int Equation::get_solved_dimensions()
{
  return m_solved_dimensions;
}
