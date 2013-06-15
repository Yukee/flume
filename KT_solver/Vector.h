#ifndef VECTOR_H
#define VECTOR_H

#include <stdexcept>
#include <ostream>

#include <iostream>

template <class T=int>
class Vector
{
public :
  Vector();
  Vector(const unsigned int &);
  Vector(const unsigned int &, const T &);
  Vector(const Vector &);
    ~Vector();
  inline Vector<T> & operator+=(const Vector<T> & right)
  {
    if(right.N != N) throw std::invalid_argument("In Vector::operator+=");
    for(unsigned int i=0;i<N;i++) m_data[i]+=right.m_data[i];
    return *this;
  }

  inline T & operator[](const unsigned int i) const
  {
    if(i>=N) throw std::invalid_argument("In Vector::operator[] const");
    return m_data[i];
  }

  inline T & operator[](const unsigned int i)
  {
    if(i>=N) throw std::invalid_argument("In Vector::operator[]");
    return m_data[i];
  }

  inline friend bool operator==(const Vector<T> & v1, const Vector<T> & v2)
  {
    bool is_eq = (v1.N == v2.N);
    if(is_eq)
      {
	for(unsigned int i=0;i<v1.N;i++) is_eq *= (v1.m_data[i] == v2.m_data[i]);
      }
    return is_eq;
  }

  inline Vector<T> & operator=(const Vector<T> & v)
  {
    if(!(&v == this))
      {
	N = v.N;
	if(m_data) delete[] m_data;
	m_data = 0;
	m_data = new T[N];
	for(unsigned int i=0;i<N;i++) m_data[i] = v.m_data[i];
      }
    return *this;
  }

  inline friend Vector<T> operator+(const Vector<T> & v1, const Vector<T> & v2)
  {
    if(v1.N != v2.N) throw std::invalid_argument("In Vector::operator+");
    Vector<T> temp(v1.N);
    for(unsigned int i=0;i<temp.N;i++) temp.m_data[i] = v1.m_data[i] + v2.m_data[i];
    return temp;
  }

  inline friend Vector<T> operator-(const Vector<T> & left, const Vector<T> & right)
  {
    if(left.N != right.N) throw std::invalid_argument("In Vector::operator-");
    Vector<T> temp(left.N);
    for(unsigned int i=0;i<temp.N;i++) temp.m_data[i] = left.m_data[i] - right.m_data[i];
    return temp;
  }

  inline friend Vector<T> operator*(const T & k, const Vector<T> & v)
  {
    Vector<T> temp(v.N);
    for(unsigned int i=0;i<temp.N;i++) temp.m_data[i] = k*v.m_data[i];
    return temp;
  }

inline friend std::ostream & operator<<(std::ostream & stream, const Vector<T> & right)
{
  for(unsigned int i=0;i<right.N;i++)
    {
      stream << right.m_data[i] << "\t";
      if((i+1)%5==0) stream << std::endl; // writes in lines of 5 components
    }
  stream << std::endl;
  return stream;
}

  unsigned int size() const;

  Vector< Vector<T> > get_base_vectors(const T additiveIdentity, const T multiplicativeIdentity);

  void resize(const unsigned int n);
  
  // drops component j of the Vector, thus returning a Vector of dim N-1
  Vector<T> drop(const int & j);

private :
  T *m_data;
  unsigned int N;
};

template <class T>
Vector<T>::Vector()
{
  N = 2;
  m_data = new T[N];
}

template <class T>
Vector<T>::Vector(const unsigned int & n)
{
  N = n;
  m_data = new T[N];
}

template <class T>
Vector<T>::Vector(const unsigned int & n, const T & value)
{
  N = n;
  m_data = new T[N];
  for(unsigned int i=0;i<N;i++) m_data[i] = value;
}

template <class T>
Vector<T>::Vector(const Vector<T> & v)
{
  N = v.N;
  m_data = new T[N];
  for(unsigned int i=0;i<N;i++) m_data[i] = v.m_data[i];
}

template <class T>
Vector<T>::~Vector()
{
  delete[] m_data;
  m_data = NULL;
}

template <class T>
unsigned int Vector<T>::size() const
{
  return N;
}

template <class T>
Vector< Vector<T> > Vector<T>::get_base_vectors(const T additiveIdentity, const T multiplicativeIdentity)
{
  Vector< Vector<T> > basis(N);
  Vector<T> tempElement(N);
  for(unsigned int d=0;d<N;d++)
    {
      for(unsigned int i=0;i<N;i++)
	{
	  if(i==d) tempElement[i] = additiveIdentity;
	  else tempElement[i] = multiplicativeIdentity;
	}
      basis[d] = tempElement;
    }
  return basis;
}

template <class T>
void Vector<T>::resize(const unsigned int n)
{
    if(n != N)
    {
        N = n;
        if(m_data) delete[] m_data;
        m_data = 0;
        m_data = new T[N];
    }
}

template <class T>
Vector<T> Vector<T>::drop(const int & j)
{
  if(j >= (int)N || j < 0) throw std::invalid_argument("In Vector::drop()");
	
  Vector<T> dropped (N-1);
	
  int it=0;
  for(unsigned int i=0;i<dropped.N;i++)
    {
      if((int)it == j) it++;
      dropped.m_data[i] = m_data[it];
      it++;
    }
	
  return dropped;
}

#endif
