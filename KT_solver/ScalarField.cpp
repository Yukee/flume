#include "ScalarField.h"
#include <stdio.h>
#include <math.h>
using namespace std;

#define DEBUG

ScalarField::ScalarField()
{
  m_r = Vector<int>(2,1);
  m_r_len = m_r.size();
  m_data_len = 1;
  m_data = new double[m_data_len];
  m_data[0] = 0;
}

ScalarField::ScalarField(const ScalarField & u) : m_r(u.m_r), m_r_len(u.m_r_len), m_data_len(u.m_data_len)
{
    m_data = new double[m_data_len];
    for(unsigned int i=0;i<m_data_len;i++) m_data[i] = u.m_data[i];
}

ScalarField::ScalarField(Vector<int> range)
{
    m_r = range;
    m_r_len = m_r.size();
    m_data_len = 1;
    for(unsigned int i=0; i<m_r_len; i++) m_data_len*=m_r[i];
    m_data = new double[m_data_len];
}

ScalarField::~ScalarField()
{
    if(m_data) delete[] m_data;
    m_data = NULL;
}

ScalarField & ScalarField::operator=(const ScalarField & u)
{
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

ScalarField & ScalarField::operator=(const double & k)
{
  for(unsigned int i=0;i<m_data_len;++i) m_data[i] = k;
  return *this;
} 

bool ScalarField::operator==(const ScalarField & u)
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


ostream & operator<<(ostream & output, const ScalarField & u)
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

ScalarField operator+(const double & k, const ScalarField & u)
{
    ScalarField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = k +u.m_data[i];
    return temp;
}

ScalarField operator+(const ScalarField & u, const ScalarField & v)
{
#ifdef DEBUG
    if(u.m_r_len != v.m_r_len) throw invalid_argument("In operator+: trying to add two fields of different dimensions");
    bool same_ranges = 1;
    for(unsigned int i=0; i<u.m_r_len; i++) same_ranges*=(u.m_r[i] == v.m_r[i]);
    if(!same_ranges) throw invalid_argument("In operator+: trying to add two fields of different ranges");
#endif

    ScalarField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = u.m_data[i] + v.m_data[i];
    return temp;
}

ScalarField operator-(const double & k, const ScalarField & u)
{
    ScalarField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = k - u.m_data[i];
    return temp;
}

ScalarField operator-(const ScalarField & u, const ScalarField & v)
{
#ifdef DEBUG
    if(u.m_r_len != v.m_r_len) throw invalid_argument("In operator+: trying to add two fields of different dimensions");
    bool same_ranges = 1;
    for(unsigned int i=0; i<u.m_r_len; i++) same_ranges*=(u.m_r[i] == v.m_r[i]);
    if(!same_ranges) throw invalid_argument("In operator+: trying to add two fields of different ranges");
#endif

    ScalarField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = u.m_data[i] - v.m_data[i];
    return temp;
}

ScalarField operator/(const ScalarField & u, const ScalarField & v)
{
#ifdef DEBUG
    if(u.m_r_len != v.m_r_len) throw invalid_argument("In operator/: trying to divide two fields of different dimensions");
    bool same_ranges = 1;
    for(unsigned int i=0; i<u.m_r_len; i++) same_ranges*=(u.m_r[i] == v.m_r[i]);
    if(!same_ranges) throw invalid_argument("In operator/: trying to divide two fields of different ranges");
#endif

    ScalarField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = u.m_data[i]/v.m_data[i];
    return temp;
}

ScalarField operator*(const double & k, const ScalarField & u)
{
    ScalarField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = k*u.m_data[i];
    return temp;
}

ScalarField operator*(const ScalarField & u, const ScalarField & v)
{
#ifdef DEBUG
    if(u.m_r_len != v.m_r_len) throw invalid_argument("In operator*: trying to multiply two fields of different dimensions");
    bool same_ranges = 1;
    for(unsigned int i=0; i<u.m_r_len; i++) same_ranges*=(u.m_r[i] == v.m_r[i]);
    if(!same_ranges) throw invalid_argument("In operator*: trying to multiply two fields of different ranges");
#endif

    ScalarField temp(u.m_r);
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i] = u.m_data[i]*v.m_data[i];
    return temp;
}

double & ScalarField::operator()(const Vector<int> & component)
{
#ifdef DEBUG
    if(component.size() != m_r_len) throw invalid_argument("In ScalarField::operator(): dimension differs from the field dimension");
    bool within_ranges = 1;
    for(unsigned int i=0; i<m_r_len; i++) within_ranges*=( (component[i] < m_r[i]) && (component[i] >= 0) );
    if(!within_ranges) throw invalid_argument("In ScalarField::operator(): component(s) are outside the ranges of the field");
#endif // DEBUG

    int element = component[m_r_len-1];
    for(unsigned int i=0; i<m_r_len-1; i++)
    {
        element*=m_r[m_r_len-2-i];
        element+=component[m_r_len-2-i];//ex. if m_d=3, element = (component[0] + m_r[0]*(component[1] + m_r[1]*component[2]))
    }
    return m_data[element];
}

double & ScalarField::operator[](const int & i)
{
  if((unsigned int)i>=m_data_len) throw invalid_argument("In ScalarField::operator[]");
  return m_data[i];
}

ScalarField ScalarField::max_field(const ScalarField u) const
{
#ifdef DEBUG
    if(u.m_r_len != m_r_len) throw invalid_argument("In max_field: trying to compare two fields of different dimensions");
    bool same_ranges = 1;
    for(unsigned int i=0; i<u.m_r_len; i++) same_ranges*=(u.m_r[i] == m_r[i]);
    if(!same_ranges) throw invalid_argument("In max_field: trying to compare two fields of different ranges");
#endif

    ScalarField temp;
    temp=u;
    for(unsigned int i=0; i<temp.m_data_len; i++) temp.m_data[i]=max( fabs(u.m_data[i]), fabs((*this).m_data[i]) );
    return temp;
}

double ScalarField::get_max() const
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

ScalarField ScalarField::module() const
{
  for(unsigned int it=0;it<m_data_len;++it) m_data[it] = max( m_data[it], -m_data[it] );
  return *this;
}

void ScalarField::resize_field(Vector<int> range)
{
    m_r_len = range.size();
    m_r = range;

    m_data_len = 1;
    for(unsigned int i=0; i<m_r_len; i++) m_data_len*=m_r[i];
    if(m_data) delete[] m_data;
    m_data = new double[m_data_len];
}

Vector<int> ScalarField::get_pos(int i) const
{
#ifdef DEBUG
  if((unsigned int)i>=m_data_len) throw invalid_argument("In ScalarField::get_pos(int)");
#endif
  Vector<int> p(m_r_len);
    int tempI=i;
    p[0] = i%m_r[0];
    for(unsigned int j=1;j<m_r_len;j++)
    {
        tempI = (tempI - p[j-1])/m_r[j-1];
        p[j] = tempI%m_r[j];
    }
    return p;
}

void ScalarField::write_in_file(ostream & output, const Vector<double> deltaX, const Vector<double> lowerLeftCorner)
{
 Vector<int> pTemp(m_r_len);
    for(unsigned int i=0;i<m_data_len;i++)
    {
        pTemp = get_pos(i);
        for(unsigned int j=0;j<m_r_len;j++) output << pTemp[j]*deltaX[j]+lowerLeftCorner[j] << "\t";
        output << m_data[i] << endl;
    }
}

void ScalarField::write_in_file_matrixform(ostream & output)
{
if(m_r_len != 2) throw invalid_argument("In ScalarField::write_in_file_matrixform: the dimension of the field must be 2");

for(unsigned int it=0;it<m_data_len;++it)
{
output << m_data[it] << "\t";
if( (it + 1) % m_r[0] == 0 ) output << endl;
}
}
