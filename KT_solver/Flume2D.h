#ifndef FLUME2D
#define FLUME2D
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

/* 2D setup in the centre plane for the flume problem. */

const static double H = 0.25; const static double W = 1;
const static double U = 2.3; const static double uF = 2; const static double alpha=0.5; const static int m = 2; const static int n = 4; const static double shift = 0.;

// write the value of the various quantities in a file
inline void write_flume_infos(std::string path, double sr)
{
	std::fstream infos;
	infos.open( ("Results/" + path + "_flume_infos.tsv").c_str(), std::ios::out);

	infos << "H\t" << "W\t" << "U\t" << "uF\t" << "sr\t" << "alpha\t" << "m\t" << "n\t" << "shift" << std::endl;
	infos << H << "\t" << W << "\t" << U << "\t" << uF << "\t" << sr << "\t" << alpha << "\t" << m << "\t" << n << "\t" << shift << std::endl;

	infos.close();
}	

inline double sech(double x)
{
  return 1/cosh(x);
}

// derivative of the stream function wrt x, ie -1*{ depth averaged v velocity in the travelling frame }
inline double dpsidx(double x)
{
  return 0;
}

// derivative of the stream function wrt y, ie { depth averaged u velocity in the travelling frame }
inline double dpsidy(double x)
{
  return -((H*(1 + 2*n)*U*tanh(x/W))/((1 + 2*m)*(1 + 2*m + 2*n)));
}

inline double y0(double x)
{
  return W*sqrt(tanh(-x/W));
}

inline double dy0dx(double x)
{
  return -0.5*(1 - tanh(x/W)*tanh(x/W))/sqrt(tanh(-x/W));
}

inline double h(double x)
{
  return H*y0(x)/W;
}

inline double dhdx(double x)
{
  return H*dy0dx(x)/W;
}

inline double dhdy(double x)
{
  return 0;
}

// speed components in the travelling frame

inline double u(double x, double z)
{
  return -uF + (dpsidy(x)/h(x) + uF)*(alpha + 2*(1 - alpha)*z/h(x));
}

inline double v(double x, double z)
{
  return (-dpsidx(x)/h(x))*(alpha + 2*(1-alpha)*z/h(x));
}

inline double w(double x, double z)
{
  return uF*(1-alpha)*(z*z*dhdx(x)/pow(h(x),2)) + (1/pow(h(x),2))*(dhdx(x)*dpsidy(x) - dhdy(x)*dpsidx(x))*(alpha+2*(1-alpha)*z/h(x))*z;
}

// dv/dy in the centre plane

inline double dvdy(double x, double z)
{
  return -((U*pow(W,-1 - 2*m - 2*n)*pow(sech(x/W),2)*(2*z*(-1 + alpha) - H*alpha*sqrt(-tanh(x/W)))*(pow(W*sqrt(-tanh(x/W)),2*(m + n)) + 2*n*pow(W,2*(m + n))*pow(-tanh(x/W),m + n))*
       pow(-tanh(x/W),-1 - m - n))/(H*(1 + 2*m)*(1 + 2*m + 2*n)));
}

// initial concentration of small particules
inline double phi0(double x, double z)
{
  x+=shift;
  double phi0 = 0;
  if(z < 0.8*h(x) && z >= 0) phi0 = 1;
  return phi0;
}

// returns 1 if within the domain of the flow
inline double boundary(double x, double z)
{
  x+=shift;
  double boundary = 0;
  if(z < h(x) && z >= 0) boundary = 1;
  return boundary;
}

#endif
