#ifndef FD1SOLVER_H
#define	FD1SOLVER_H

#include "Equation.h"
#include "Vector.h"
#include "ScalarField.h"

enum grid_position {lft, rght};
enum bound {lower, upper};

class FD1Solver
{
 public:
  FD1Solver(Vector<double> deltaX, Vector<double> xInterval, Equation *eq, Vector<double> lowerLeftCorner);

  VectorField get_numerical_flux_gradient(const VectorField & un);

  inline VectorField get_position()
  {
    return m_position;
  }

  double check_CFL(double deltaT);

  ~FD1Solver();

  inline Vector<double> get_deltaX()
  {
    return m_deltaX;
  }
    
  inline Vector<int> get_nxSteps()
  {
    return m_nxSteps;
  }

  inline int get_space_dimensions()
  {
    return m_n;
  }

  inline int get_solved_dimensions()
  {	
    return m_m;
  }

  // returns a 2 vector: first component; left bound, second component: right bound.
  inline Vector< Vector<double> > get_domain_bounds()
  {
    Vector< Vector<double> > dom ( 2, Vector<double> (m_n) );
    dom[0] = m_lowerLeftCorner;
    dom[1] = dom[0] + m_xInterval;
    return dom;
  }

  inline Vector<double> get_lowerLeftCorner()
  {
    return m_lowerLeftCorner;
  }

  inline VectorField get_flux(VectorField u0, int d)
  {
    get_numerical_flux_gradient(u0);
    return m_eq->get_convectionFlux(m_un, d);
  }

  inline VectorField get_intermediate_values(VectorField u0, int d)
  {
    get_numerical_flux_gradient(u0);
    return lower_left_intermediate_un_values[d];
  }

  inline VectorField get_derivatives(VectorField u0, int d)
  {
    get_numerical_flux_gradient(u0);
    return un_derivatives[0];
  }

  inline VectorField get_initial_field(VectorField u0)
  {
    get_numerical_flux_gradient(u0);
    return m_un;
  }

  inline SField get_unity()
  {
    return unity;
  }

  inline Vector< Vector<int> > get_b()
  {
    return m_b;
  }

 private:
  int m_m; //number of solved dimensions
  int m_n; //number of space dimensions

  Vector< Vector<int> > m_b; // base vectors of the physical space

  Vector<double> m_deltaX;
  Vector<double> m_xInterval;
  Equation *m_eq;
  Vector<double> m_lowerLeftCorner;
  Vector<int> m_nxSteps;

  VectorField m_position;

  VectorField m_un;

  TensorField left_convection_flux;
  TensorField right_convection_flux;
  TensorField left_diffusion_flux;
  TensorField right_diffusion_flux;

  TensorField un_derivatives;
  TensorField upper_right_intermediate_un_values;
  TensorField lower_right_intermediate_un_values;
  TensorField upper_left_intermediate_un_values;
  TensorField lower_left_intermediate_un_values;

  VectorField right_localSpeed;
  VectorField left_localSpeed;

  VectorField source_term;

  void compute_south_boundary();

  void compute_un_derivatives();
  void compute_intermediate_un_values();
  void compute_localSpeed();
  void compute_numerical_convection_flux();
  void compute_numerical_diffusion_flux();
  void compute_source_term();

  double minmod(double a, double b);
  double three_pts_derivative(int it, int d, int i); //du/dx at point get_pos(it) in the x or y direction

  // caca
  SField unity;
};


#endif
