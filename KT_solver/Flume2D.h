#ifndef FLUME2D
#define FLUME2D
#include <math.h>
#include <iostream>

/* 2D setup in the centre plane for the flume problem.
Assuming that H = W = 1. If not, multiply by the expression of H and W indicated in the commented text before each function. If "complicated" is indicated, you cannot just multiply by something and you will have to recalculate th expression depending and H and W. If nothing is indicated nothing is to be changed. */

const static double U = 2.3; const static double uF = 2.; const static double alpha=0.5; const static int m = 2; const static int n = 4; const static double shift = 0.0;

double Sech(double x)
{
  return 1/cosh(x);
}

double Csch(double x)
{
  if(x==0) std::cout << "division by zero in Flume2D.h" << std::endl;
  return 1/sinh(x);
}

// complicated
// derivative of the stream function wrt x, ie -1*{ depth averaged v velocity in the travelling frame }
double dpsidx(double x)
{
  return 0;
}

// complicated
// derivative of the stream function wrt y, ie { depth averaged u velocity in the travelling frame }
double dpsidy(double x)
{
  return U*((-1 - 2*n)/((1 + 2*m)*(1 + 2*m + 2*n)) + pow(0,2*m)/pow(-tanh(x),m) - pow(0,2*(m + n))*pow(-tanh(x),-m - n) + (pow(0,2*n)*(1 + 2*n))/((1 + 2*m)*(1 + 2*m + 2*n)*pow(-tanh(x),n)))*tanh(x);
}

// *H/W
double h(double x)
{
  return y0(x);
}

// *W
double y0(double x)
{
  return sqrt(tanh(-x));
}

// *H/W
double dhdx(double x)
{
  return -Sech(x)*Sech(x)/(2*sqrt(tanh(-x)));
}

// *H/W
double dhdy(double x)
{
  return 0;
}

// speed components in the travelling frame

double u(double x, double z)
{
  return -uF + (dpsidy(x)/h(x) + uF)*(alpha + 2*(1 - alpha)*z/h(x));
}

double v(double x, double z)
{
  return (-dpsidx(x)/h(x))*(alpha + 2*(1-alpha)*z/h(x));
}

double w(double x, double z)
{
  return uF*(1-alpha)*(z*z*dhdx(x)/pow(h(x),2)) + (1/pow(h(x),2))*(dhdx(x)*dpsidy(x) - dhdy(x)*dpsidx(x))*(alpha+2*(1-alpha)*z/h(x))*z;
}

// dv/dy in the centre plane

double dvdy(double x, double z)
{
  return ((1 + 2*n)*U*Csch(x)*Sech(x)*(2*z*(-1 + alpha) - alpha*sqrt(-tanh(x))))/((1 + 2*m)*(1 + 2*m + 2*n));
}

// initial concentration of small particules
double phi0(double x, double z)
{
  x+=shift;
  double phi0 = 0;
  if(z < 0.80*h(x) && z >= 0) phi0 = 1;
  return phi0;
}

// returns 1 if within the domain of the flow
double boundary(double x, double z)
{
  x+=shift;
  double boundary = 0;
  if(z < h(x) && z >= 0) boundary = 1;
  return boundary;
}

#endif
