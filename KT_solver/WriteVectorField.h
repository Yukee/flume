#ifndef WRITEVECTOR
#define WRITEVECTOR

#include "Flux.h" // Vector, ScalarField and children
#include <ostream>
#include <stdio.h>

inline void write_VectorField(const VectorField & u, const VectorField & pos, std::ostream & stream)
{

  if(u.size() != pos.size()) throw std::invalid_argument("In write_VectorField, VectorFields position and u must have the same dimension.");
  int n = u.size();

  for(int d=0;d<n;d++) if(!(u[d].get_range() == pos[d].get_range())) throw std::invalid_argument("In write_VectorField, each ScalarField component of VectorFields position and u must have the same range.");

  bool same_sizes = true;
  for(int d=0;d<n-1;d++) same_sizes*=(u[d].get_size() == u[d+1].get_size() && pos[d].get_size() == pos[d+1].get_size());
  if(!same_sizes) throw std::invalid_argument("In write_VectorField, the sizes of each component of the pos and u Fields must be the same.");
  int size = u[0].get_size();

  for(int it=0;it<size;++it)
    {
      for(int d=0;d<n;d++) stream << pos[d][it] << "\t";

      for(int d=0;d<n;d++) 
	{
	  stream << u[d][it];
	  if(d<n-1) stream << "\t";
	}

      stream << std::endl;
    }

}

#endif
