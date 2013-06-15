#ifndef IMMFIELD_H
#define IMMFIELD_H

#include "ScalarField.h"

class ImmersedField : public ScalarField
{
 public :
  virtual double & operator()(const Vector<int> &);
    
    ImmersedField() : ScalarField() {}
    ImmersedField(Vector<int> range) : ScalarField(range) {}
    ImmersedField(const ImmersedField & u) : ScalarField(u) {}
    virtual ~ImmersedField();
    
    ImmersedField & operator=(const double &);
    
 /*************************/
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
    ImmersedField & operator=(const ImmersedField &);
    bool operator==(const ImmersedField &);
    friend std::ostream & operator<<(std::ostream &, const ImmersedField &);
    friend ImmersedField operator+(const double &, const ImmersedField &);
    friend ImmersedField operator+(const ImmersedField &, const ImmersedField &);
    friend ImmersedField operator-(const double &, const ImmersedField &);
    friend ImmersedField operator-(const ImmersedField &, const ImmersedField &);
    friend ImmersedField operator*(const double &, const ImmersedField &);
    friend ImmersedField operator*(const ImmersedField &, const ImmersedField &);
    friend ImmersedField operator/(const ImmersedField &, const ImmersedField &);
    double & operator[](const int &);
    ImmersedField max_field(const ImmersedField) const;
    ImmersedField module() const;
    double get_max() const;
    void write_in_file(std::ostream & output, const Vector<double> deltaX, const Vector<double> lowerLeftCorner);
    void write_in_file_matrixform(std::ostream & output);
    Vector<int> get_pos(int) const;//retrieves the Position of the ith element of m_data
};

#endif
