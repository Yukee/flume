#ifndef FLUME3D
#define FLUME3D
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

/* 3D setup in the centre plane for the flume problem. */

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

inline double Sech(double x)
{
  return 1/cosh(x);
}

// derivative of the stream function wrt x, ie -1*{ depth averaged v velocity in the travelling frame }
inline double dpsidx(double x, double y)
{
  return (H*U*y*pow(Sech(x/W),2)*((1 + 2*m)*(-1 + m + n)*pow(y,2*(m + n)) - (-1 + n)*pow(W,2*m)*pow(y,2*n)*pow(-tanh(x/W),m) - 
        (-1 + m)*(1 + 2*m + 2*n)*pow(W,2*n)*pow(y,2*m)*pow(-tanh(x/W),n) - (1 + 2*n)*pow(W,2*(m + n))*pow(-tanh(x/W),m + n)))/
    ((1 + 2*m)*(1 + 2*m + 2*n)*W*pow(W*sqrt(-tanh(x/W)),2*(m + n)));
}

// derivative of the stream function wrt y, ie { depth averaged u velocity in the travelling frame }
inline double dpsidy(double x, double y)
{
  return (H*U*((1 + 2*m)*(1 + 2*m + 2*n)*pow(y,2*m) - (1 + 2*n)*pow(W,2*m)*pow(-tanh(x/W),m))*
      (pow(y,2*n) - pow(W,2*n)*pow(-tanh(x/W),n))*pow(-tanh(x/W),1 - m - n))/((1 + 2*m)*(1 + 2*m + 2*n)*pow(W,2*(m + n)));
}

inline double h(double x, double y)
{
  return (H/W)*(pow(y0(x),2*n)-pow(y,2*n))/pow(y0(x),2*n-1);
}

inline double y0(double x)
{
  return W*sqrt(tanh(-x/W));
}

inline double dy0dx(double x)
{
  return -0.5*(1 - tanh(x/W)*tanh(x/W))/sqrt(tanh(-x/W));
}

inline double dhdx(double x, double y)
{
  return (H/W)*(dy0dx(x) + (2*n-1)*pow(y/y0(x),2*n)*dy0dx(x));
}

inline double dhdy(double x, double y)
{
  return -2*n*(H/W)*pow(y/y0(x), 2*n-1);
}

// speed components in the travelling frame

inline double u(double x, double y, double z)
{
  return -uF + (dpsidy(x,y)/h(x,y) + uF)*(alpha + 2*(1 - alpha)*z/h(x,y));
}

inline double v(double x, double y, double z)
{
  return (-dpsidx(x,y)/h(x,y))*(alpha + 2*(1-alpha)*z/h(x,y));
}

inline double w(double x, double y, double z)
{
  return uF*(1-alpha)*(z*z*dhdx(x,y)/pow(h(x,y),2)) + (1/pow(h(x,y),2))*(dhdx(x,y)*dpsidy(x,y) - dhdy(x,y)*dpsidx(x,y))*(alpha+2*(1-alpha)*z/h(x,y))*z;
}

// initial concentration of small particules
inline double phi0(double x, double y, double z)
{
  x+=shift;
  double phi0 = 0;
  if(z < 0.8*h(x,y) && z >= 0) phi0 = 1;
  return phi0;
}

// returns 1 if within the domain of the flow
inline double boundary(double x, double y, double z)
{
  x+=shift;
  double boundary = 0;
  if(z < h(x,y) && z >= 0) boundary = 1;
  return boundary;
}

#endif
