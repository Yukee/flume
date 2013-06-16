#include "Flume2DField.h"
#include <stdio.h>
#include <math.h>
using namespace std;

double & Flume2DField::operator()(Vector<int> component)
{
  for(unsigned int d=0;d<m_r_len;d++) if(component[d] > m_r[d] || component[d] < -1)
					{
					  cout << "received m_r["<<d<<"] = " << m_r[d] << " expected -1 <= m_r[d] <= " << m_r[d] << endl;
					  throw invalid_argument("In Flume2DField::operator()");
					}
		
	if(component[0] == m_r[0])
	{
	  component[0] -= 1;
	}
	
	if(component[1] == -1)
	{
	  component[1] += 1;
	}
				
	if(component[1] == m_r[1])
	{
	  component[1] -= 1;
	}

  return PrescribedField::operator ()(component);
}

// Set west boundary
void Flume2DField::set_bound(const ScalarField & u)
{
	PrescribedField::set_bound(0,-1,u);
}

Flume2DField::~Flume2DField()
{
	if(m_data) delete[] m_data;
	m_data = NULL;
}

void Flume2DField::resize_field(Vector<int> range)
{
    PrescribedField::resize_field(range);
}

Vector<int> Flume2DField::get_pos(int i) const
{
    return ScalarField::get_pos(i);
}

Flume2DField & Flume2DField::operator=(const Flume2DField & u)
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

Flume2DField & Flume2DField::operator=(const double & k)
{
  for(unsigned int i=0;i<m_data_len;++i) m_data[i] = k;
  return *this;
} 

Flume2DField operator+(const double & k, const Flume2DField & u)
{
    Flume2DField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = k +u.m_data[i];

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = operator+(k, u.m_bounds[i]);
    return temp;
}

/************************************************/


bool Flume2DField::operator==(const Flume2DField & u)
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


ostream  & operator<<(ostream & output, const Flume2DField & u)
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


Flume2DField operator+(const Flume2DField & u, const Flume2DField & v)
{
#ifdef DEBUG
    if(u.m_r_len != v.m_r_len) throw invalid_argument("In operator+: trying to add two fields of different dimensions");
    bool same_ranges = 1;
    for(unsigned int i=0; i<u.m_r_len; i++) same_ranges*=(u.m_r[i] == v.m_r[i]);
    if(!same_ranges) throw invalid_argument("In operator+: trying to add two fields of different ranges");
#endif

    Flume2DField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = u.m_data[i] + v.m_data[i];

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = operator+(u.m_bounds[i], v.m_bounds[i]);

    return temp;
}

Flume2DField operator-(const double & k, const Flume2DField & u)
{
    Flume2DField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = k - u.m_data[i];

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = operator-(k, u.m_bounds[i]);

    return temp;
}

Flume2DField operator-(const Flume2DField & u, const Flume2DField & v)
{
#ifdef DEBUG
    if(u.m_r_len != v.m_r_len) throw invalid_argument("In operator+: trying to add two fields of different dimensions");
    bool same_ranges = 1;
    for(unsigned int i=0; i<u.m_r_len; i++) same_ranges*=(u.m_r[i] == v.m_r[i]);
    if(!same_ranges) throw invalid_argument("In operator+: trying to add two fields of different ranges");
#endif

    Flume2DField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = u.m_data[i] - v.m_data[i];

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = operator-(u.m_bounds[i], v.m_bounds[i]);

    return temp;
}

Flume2DField operator*(const double & k, const Flume2DField & u)
{
    Flume2DField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = k*u.m_data[i];

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = operator*(k, u.m_bounds[i]);

    return temp;
}

Flume2DField operator*(const Flume2DField & u, const Flume2DField & v)
{
#ifdef DEBUG
    if(u.m_r_len != v.m_r_len) throw invalid_argument("In operator*: trying to multiply two fields of different dimensions");
    bool same_ranges = 1;
    for(unsigned int i=0; i<u.m_r_len; i++) same_ranges*=(u.m_r[i] == v.m_r[i]);
    if(!same_ranges) throw invalid_argument("In operator*: trying to multiply two fields of different ranges");
#endif

    Flume2DField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = u.m_data[i]*v.m_data[i];

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = operator*(u.m_bounds[i], v.m_bounds[i]);

    return temp;
}

double & Flume2DField::operator[](const int & i)
{
  if((unsigned int)i>=m_data_len) {
	  cout << "received i = " << i << " expected 0 <= i <= " << m_data_len - 1 << endl;
	  throw invalid_argument("In Flume2DField::operator[]");
  }
  return m_data[i];
}

Flume2DField Flume2DField::max_field(const Flume2DField u) const
{
#ifdef DEBUG
    if(u.m_r_len != m_r_len) throw invalid_argument("In max_field: trying to compare two fields of different dimensions");
    bool same_ranges = 1;
    for(unsigned int i=0; i<u.m_r_len; i++) same_ranges*=(u.m_r[i] == m_r[i]);
    if(!same_ranges) throw invalid_argument("In max_field: trying to compare two fields of different ranges");
#endif

    Flume2DField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i]=max( fabs(u.m_data[i]), fabs((*this).m_data[i]) );

    // for(unsigned int i=0;i<2*u.m_r_len;i++) temp.m_bounds[i] = ((*this).m_bounds[i]).ScalarField::max_field(u.m_bounds[i]);

    return temp;
}

double Flume2DField::get_max() const
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

Flume2DField Flume2DField::module() const
{
  for(unsigned int it=0;it<m_data_len;++it) m_data[it] = max( m_data[it], -m_data[it] );

  for(unsigned int i=0;i<2*m_r_len;i++) m_bounds[i] = m_bounds[i].ScalarField::module();

  return *this;
}

void Flume2DField::write_in_file(ostream & output, const Vector<double> deltaX, const Vector<double> lowerLeftCorner)
{
  Vector<int> pTemp(m_r_len);
    for(unsigned int i=0;i<m_data_len;i++)
    {
        pTemp = get_pos(i);
        for(unsigned int j=0;j<m_r_len;j++) output << pTemp[j]*deltaX[j]+lowerLeftCorner[j] << "\t";
        output << m_data[i] << endl;
    }
}

void Flume2DField::write_in_file_matrixform(ostream & output)
{
if(m_r_len != 2) throw invalid_argument("In ScalarField::write_in_file_matrixform: the dimension of the field must be 2");

for(unsigned int it=0;it<m_data_len;++it)
{
output << m_data[it] << "\t";
if( (it + 1) % m_r[0] == 0 ) output << endl;
}
}
