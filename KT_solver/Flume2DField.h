#ifndef F2FIELD_H
#define F2FIELD_H

#include "PrescribedField.h"
#include <stdio.h>

class Flume2DField : public PrescribedField
{
 public:
  virtual double & operator()(Vector<int> );

  Flume2DField() {Flume2DField(Vector<int> (2,1));}
	 
 Flume2DField(Vector<int> range): PrescribedField(range) {
	 // Designed for a 2D problem
	 //if(range.size() != 2) std::throw invalid_argument("In Flume2DField::Flume2DField(range) dimension must be 2");
 }
	 
 Flume2DField(const Flume2DField & u): PrescribedField(u) {}

// set west boundary
void set_bound(const ScalarField & u);

 inline ScalarField get_bound()
 {
	 // West boundary
   return m_bounds[0];
 }

	virtual ~Flume2DField();
     
    Flume2DField & operator=(const Flume2DField &);
    Flume2DField & operator=(const double &);

  /**************************************/

    void resize_field(Vector<int> range);
    inline unsigned int get_space_dimension() const
    {
        return m_r_len;
    }
    inline Vector<int> get_range() const
    {
        return m_r;
    }
    inline int get_size() const
    {
        return m_data_len;
    }
    bool operator==(const Flume2DField &);
    friend std::ostream & operator<<(std::ostream &, const Flume2DField &);
    friend Flume2DField operator+(const double &, const Flume2DField &);
    friend Flume2DField operator+(const Flume2DField &, const Flume2DField &);
    friend Flume2DField operator-(const double &, const Flume2DField &);
    friend Flume2DField operator-(const Flume2DField &, const Flume2DField &);
    friend Flume2DField operator*(const double &, const Flume2DField &);
    friend Flume2DField operator*(const Flume2DField &, const Flume2DField &);
    friend Flume2DField operator/(const Flume2DField &, const Flume2DField &);
    double & operator[](const int &);
    Flume2DField max_field(const Flume2DField) const;
    Flume2DField module() const;
    double get_max() const;
    void write_in_file(std::ostream & output, const Vector<double> deltaX, const Vector<double> lowerLeftCorner);
    void write_in_file_matrixform(std::ostream & output);
    Vector<int> get_pos(int) const;//retrieves the Position of the ith element of m_data  
    
};

#endif
