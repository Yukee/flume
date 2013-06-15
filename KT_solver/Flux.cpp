#include "Flux.h"

Flux::Flux() :  m_space_dimensions(1), m_solved_dimensions(1)
{
  m_evaluated_flux_d.resize(m_solved_dimensions);
}

Flux::Flux(int space_dimensions, int solved_dimensions) : m_space_dimensions(space_dimensions), m_solved_dimensions(solved_dimensions)
{
  m_evaluated_flux_d.resize(m_solved_dimensions);
}

int Flux::get_space_dimensions()
{
  return m_space_dimensions;
}

int Flux::get_solved_dimensions()
{
  return m_solved_dimensions;
}

VectorField Flux::evaluate(const VectorField & u, const int d)
{
	for(int i=0;i<m_solved_dimensions;i++) m_evaluated_flux_d[i] = 0*u[i];
	return m_evaluated_flux_d;
}

VectorField Flux::evaluate(const VectorField & u, const VectorField & v, const int d)
{
	return evaluate(u,d);
}

SField Flux::get_max_eigenvalue(const VectorField & u, const int d)
{
	m_max_eigenvalue = 0*u[0];
	return m_max_eigenvalue;
}

void Flux::set_parameter(const SField & param)
{
	m_param = param;
}

void Flux::set_parameter(const VectorField & param)
{
}
