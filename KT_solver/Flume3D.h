#ifndef FLUME3D
#define FLUME3D
#include <math.h>

/* Full 3D setup for the flume problem.
Assuming that H = W = 1. If not, multiply by the expression of H and W indicated in the commented text before each function. If "complicated" is indicated, you cannot just multiply by something and you will have to recalculate th expression depending and H and W. If nothing is indicated nothing is to be changed.
Everything can be further simplified in the centre plane. See Flume2D for the details. */

const static double U = 2.3; const static double uF = 2; const static double alpha=0.5; const static int m = 2; const static int n= 4; const static double shift = 0.0;

double Sech(double x)
{
  return 1/cosh(x);
}

// complicated
// derivative of the stream function wrt x, ie -1*{ depth averaged v velocity in the travelling frame }
double dpsidx(double x, double y)
{
  return (U*y*pow(Sech(x),2)*((1 + 2*m)*(-1 + m + n)*pow(y,2*(m + n)) - (-1 + n)*pow(y,2*n)*pow(-tanh(x),m) - (-1 + m)*(1 + 2*m + 2*n)*pow(y,2*m)*pow(-tanh(x),n) - (1 + 2*n)*pow(-tanh(x),m + n))*pow(-tanh(x),-m - n))/((1 + 2*m)*(1 + 2*m + 2*n));
}

// complicated
// derivative of the stream function wrt y, ie { depth averaged u velocity in the travelling frame }
double dpsidy(double x, double y)
{
  return (U*((1 + 2*m)*(1 + 2*m + 2*n)*pow(y,2*m) - (1 + 2*n)*pow(-tanh(x),m))*
	  (pow(y,2*n) - pow(-tanh(x),n))*pow(-tanh(x),1 - m - n))/((1 + 2*m)*(1 + 2*m + 2*n));
}

// *H/W
double h(double x, double y)
{
  return (pow(y0(x),2*n) - pow(y, 2*n))/pow(y0(x), 2*n - 1);
}

// *W
double y0(double x)
{
  return sqrt(tanh(-x));
}

// *H/W
double dhdx(double x, double y)
{
  return (-Sech(x)*Sech(x)/(2*sqrt(tanh(-x))))*(1+(2*n-1)*pow(y/y0(x), 2*n));
}

// *H/W
double dhdy(double x, double y)
{
  return -2*n*pow(y, 2*n-1);
}

// speed components in the travelling frame

double u(double x, double y, double z)
{
  return -uF + (dpsidy(x,y)/h(x,y) + uF)*(alpha + 2*(1 - alpha)*z/h(x,y));
}

double v(double x, double y, double z)
{
  return (-dpsidx(x,y)/h(x,y))*(alpha + 2*(1-alpha)*z/h(x,y));
}

double w(double x, double y, double z)
{
  return uF*(1-alpha)*(z*z*dhdx(x,y)/pow(h(x,y),2)) + (1/pow(h(x,y),2))*(dhdx(x,y)*dpsidy(x,y) - dhdy(x,y)*dpsidx(x,y))*(alpha+2*(1-alpha)*z/h(x,y))*z;
}

// initial concentration of small particules
double phy0(double x, double y, double z)
{
  x+=shift;
  double phy0 = 0;
  if(z <= 0.95*h(x,y) && z >= 0) phy0 = 1;
  return phy0;
}

// returns 1 if within the domain of the flow
double boundary(double x, double y, double z)
{
  x+=shift;
  double boundary = 0;
  if(z <= h(x,y) && z >= 0) boundary = 1;
  return boundary;
}

#endif
