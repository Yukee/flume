#include "PrescribedField.h"
#include <stdio.h>
#include <math.h>
using namespace std;

#define DEBUG

double & PrescribedField::operator()(Vector<int> component)
{
  for(unsigned int d=0;d<m_r_len;d++) if(component[d] > m_r[d] || component[d] < -1)
					{
					  cout << "received m_r["<<d<<"] = " << m_r[d] << " expected -1 <= m_r[d] <= " << m_r[d] << endl;
					  throw invalid_argument("In PrescribedField::operator()");
					}

  for(unsigned int d=0;d<m_r_len;d++)
    {
      if(component[d] == -1)
	{
	  component = component.drop(d);
	  return m_bounds[2*d](component);
	}
		
      if(component[d] == m_r[d])
	{
	  component = component.drop(d);
	  return m_bounds[2*d+1](component);
	}
    }
  return ScalarField::operator ()(component);
}

void PrescribedField::set_bound(const int d, const int i, const ScalarField & u)
{
	if( i*i != 1 ) throw invalid_argument("In PrescribedField::set_bound(int d, int i, ScalarField u) direction i must be either -1 or 1");
	if( d >= (int)m_r_len || d < 0 ) throw invalid_argument("In PrescribedField::set_bound(int d, int i, ScalarField u) axis number d must be >=0 and <m_r_len");
	if( !(u.get_range() == m_bounds[2*d].get_range()) ) throw invalid_argument("In PrescribedField::set_bound(int d, int i, ScalarField u) invalid u range");
	
	if(i == -1) m_bounds[2*d] = u;
	else m_bounds[2*d+1] = u;
}

ScalarField PrescribedField::get_bound(const int d, const int i)
{
	if( i*i != 1 ) throw invalid_argument("In PrescribedField::get_bound(int d, int i) direction i must be either -1 or 1");
	if( d >= (int)m_r_len || d < 0 ) throw invalid_argument("In PrescribedField::get_bound(int d, int i) axis number d must be >=0 and <m_r_len");

	if(i == -1) return m_bounds[2*d];
	return m_bounds[2*d+1];
}

PrescribedField::~PrescribedField()
{
	if(m_data) delete[] m_data;
	m_data = NULL;
}

void PrescribedField::resize_field(Vector<int> range)
{
    m_r_len = range.size();
    m_r = range;

    m_bounds.resize(2*m_r_len);
    Vector<int> range_surf = m_r;
    //Vector< Vector<int> > b = range_surf.get_base_vectors(1,0); for(unsigned int i=0;i<m_r_len;i++) range_surf = range_surf + 2*b[i];
	 
    for(unsigned int d=0;d<m_r_len;d++)
      {
	m_bounds[2*d].resize_field(range_surf.drop(d)); m_bounds[2*d] = 0;
	m_bounds[2*d+1].resize_field(range_surf.drop(d)); m_bounds[2*d+1] = 0;
      }

    m_data_len = 1;
    for(unsigned int i=0; i<m_r_len; i++) m_data_len*=m_r[i];
    if(m_data) delete[] m_data;
    m_data = new double[m_data_len];
}

Vector<int> PrescribedField::get_pos(int i) const
{
    return ScalarField::get_pos(i);
}

PrescribedField & PrescribedField::operator=(const PrescribedField & u)
{
  m_bounds = u.m_bounds;

    if(!(&u == this))
    {
        m_r_len = u.m_r_len;
        m_r = u.m_r;

        if(m_data_len != u.m_data_len)
        {
            m_data_len = u.m_data_len;
            if(m_data) delete[] m_data;
            m_data = new double[m_data_len];
        }

        for(unsigned int i=0;i<m_data_len;i++) m_data[i] = u.m_data[i];
    }

    return *this;
}

PrescribedField & PrescribedField::operator=(const double & k)
{
  for(unsigned int i=0;i<m_data_len;++i) m_data[i] = k;
  return *this;
} 

PrescribedField operator+(const double & k, const PrescribedField & u)
{
    PrescribedField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = k +u.m_data[i];

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = operator+(k, u.m_bounds[i]);
    return temp;
}

/************************************************/


bool PrescribedField::operator==(const PrescribedField & u)
{
    bool are_equal = false;
    if(m_data_len == u.m_data_len)
    {
        bool same_ranges = true;
        for(unsigned int i=0; i<m_r_len; i++) same_ranges*=(u.m_r[i] == m_r[i]);
        if(same_ranges)
        {
            bool same_data = true;
            for(unsigned int i=0;i<m_data_len;i++) same_data*=(u.m_data[i] == m_data[i]);
            if(same_data) are_equal = true;
        }
    }

    return are_equal;
}


ostream  & operator<<(ostream & output, const PrescribedField & u)
{
  Vector<int> pTemp(u.m_r_len);
    for(unsigned int i=0;i<u.m_data_len;i++)
    {
        pTemp = u.get_pos(i);
        for(unsigned int j=0;j<u.m_r_len;j++) output << pTemp[j] << "\t";
        output << u.m_data[i] << endl;
    }

    return output;
}


PrescribedField operator+(const PrescribedField & u, const PrescribedField & v)
{
#ifdef DEBUG
    if(u.m_r_len != v.m_r_len) throw invalid_argument("In operator+: trying to add two fields of different dimensions");
    bool same_ranges = 1;
    for(unsigned int i=0; i<u.m_r_len; i++) same_ranges*=(u.m_r[i] == v.m_r[i]);
    if(!same_ranges) throw invalid_argument("In operator+: trying to add two fields of different ranges");
#endif

    PrescribedField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = u.m_data[i] + v.m_data[i];

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = operator+(u.m_bounds[i], v.m_bounds[i]);

    return temp;
}

PrescribedField operator-(const double & k, const PrescribedField & u)
{
    PrescribedField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = k - u.m_data[i];

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = operator-(k, u.m_bounds[i]);

    return temp;
}

PrescribedField operator-(const PrescribedField & u, const PrescribedField & v)
{
#ifdef DEBUG
    if(u.m_r_len != v.m_r_len) throw invalid_argument("In operator+: trying to add two fields of different dimensions");
    bool same_ranges = 1;
    for(unsigned int i=0; i<u.m_r_len; i++) same_ranges*=(u.m_r[i] == v.m_r[i]);
    if(!same_ranges) throw invalid_argument("In operator+: trying to add two fields of different ranges");
#endif

    PrescribedField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = u.m_data[i] - v.m_data[i];

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = operator-(u.m_bounds[i], v.m_bounds[i]);

    return temp;
}

PrescribedField operator*(const double & k, const PrescribedField & u)
{
    PrescribedField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = k*u.m_data[i];

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = operator*(k, u.m_bounds[i]);

    return temp;
}

PrescribedField operator*(const PrescribedField & u, const PrescribedField & v)
{
#ifdef DEBUG
    if(u.m_r_len != v.m_r_len) throw invalid_argument("In operator*: trying to multiply two fields of different dimensions");
    bool same_ranges = 1;
    for(unsigned int i=0; i<u.m_r_len; i++) same_ranges*=(u.m_r[i] == v.m_r[i]);
    if(!same_ranges) throw invalid_argument("In operator*: trying to multiply two fields of different ranges");
#endif

    PrescribedField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = u.m_data[i]*v.m_data[i];

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = operator*(u.m_bounds[i], v.m_bounds[i]);

    return temp;
}

double & PrescribedField::operator[](const int & i)
{
  if((unsigned int)i>=m_data_len) {
	  cout << "received i = " << i << " expected 0 <= i <= " << m_data_len - 1 << endl;
	  throw invalid_argument("In PrescribedField::operator[]");
  }
  return m_data[i];
}

PrescribedField PrescribedField::max_field(const PrescribedField u) const
{
#ifdef DEBUG
    if(u.m_r_len != m_r_len) throw invalid_argument("In max_field: trying to compare two fields of different dimensions");
    bool same_ranges = 1;
    for(unsigned int i=0; i<u.m_r_len; i++) same_ranges*=(u.m_r[i] == m_r[i]);
    if(!same_ranges) throw invalid_argument("In max_field: trying to compare two fields of different ranges");
#endif

    PrescribedField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i]=max( fabs(u.m_data[i]), fabs((*this).m_data[i]) );

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = ((*this).m_bounds[i]).ScalarField::max_field(u.m_bounds[i]);

    return temp;
}

double PrescribedField::get_max() const
{
    double maximum=0;
    double temp=0;
    for(unsigned int i=0;i<m_data_len;i++)
    {
        temp = fabs(m_data[i]);
        if(temp>maximum) maximum=temp;
    }
    return maximum;
}

PrescribedField PrescribedField::module() const
{
  for(unsigned int it=0;it<m_data_len;++it) m_data[it] = max( m_data[it], -m_data[it] );

  for(unsigned int i=0;i<2*m_r_len;i++) m_bounds[i] = m_bounds[i].ScalarField::module();

  return *this;
}

void PrescribedField::write_in_file(ostream & output, const Vector<double> deltaX, const Vector<double> lowerLeftCorner)
{
  Vector<int> pTemp(m_r_len);
    for(unsigned int i=0;i<m_data_len;i++)
    {
        pTemp = get_pos(i);
        for(unsigned int j=0;j<m_r_len;j++) output << pTemp[j]*deltaX[j]+lowerLeftCorner[j] << "\t";
        output << m_data[i] << endl;
    }
}

void PrescribedField::write_in_file_matrixform(ostream & output)
{
if(m_r_len != 2) throw invalid_argument("In ScalarField::write_in_file_matrixform: the dimension of the field must be 2");

for(unsigned int it=0;it<m_data_len;++it)
{
output << m_data[it] << "\t";
if( (it + 1) % m_r[0] == 0 ) output << endl;
}
}



