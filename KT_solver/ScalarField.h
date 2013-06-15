#ifndef SCALARFIELD_H
#define SCALARFIELD_H

#include "Vector.h"
#include <cstring>
#include <exception>
#include <stdexcept>
#include <iostream>

class ScalarField
{

public:
    // constructs a scalar field on a 2D space, equal to 0 at the origin and not defined elsewhere.
    ScalarField();

    // overloading the copy constructor because the default one let both m_data pointers point the same memory
    ScalarField(const ScalarField &);

    // creates a scalar field on a range.size()D space, defined on a rectangle of lengths range[i]. m_data is NOT initialized.
    ScalarField(Vector<int> range);

    virtual ~ScalarField();
    
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
    ScalarField & operator=(const ScalarField &);
    ScalarField & operator=(const double &);
    bool operator==(const ScalarField &);

    // each line of the stream contains position coordinates followed by field value at that position
    friend std::ostream & operator<<(std::ostream &, const ScalarField &);

    friend ScalarField operator+(const double &, const ScalarField &);
    friend ScalarField operator+(const ScalarField &, const ScalarField &);
    friend ScalarField operator-(const double &, const ScalarField &);
    friend ScalarField operator-(const ScalarField &, const ScalarField &);
    friend ScalarField operator*(const double &, const ScalarField &);
    friend ScalarField operator*(const ScalarField &, const ScalarField &);
    friend ScalarField operator/(const ScalarField &, const ScalarField &);

    // returns field value at position component
    virtual double & operator()(const Vector<int> & component);

    // returns m_data[i]; usefull to do dimension-independant data initialization
    double & operator[](const int & i);

    // costructs a field which is the local maximum of two fields at each position
    ScalarField max_field(const ScalarField) const;

    ScalarField module() const;

    // returns the maximal value of a field
    double get_max() const;

    // each line of the stream contains position coordinates followed by field value at that position
    void write_in_file(std::ostream & output, const Vector<double> deltaX, const Vector<double> lowerLeftCorner);
    // line i column j of the stream is field value at position (i,j)
    void write_in_file_matrixform(std::ostream & output);
    
    // retrieves the Position of the ith element of m_data
    Vector<int> get_pos(int) const;

protected:
    Vector<int> m_r;//has as much components as the number of dimensions of the space, each component is the number of points in a direction of the space
    unsigned int m_r_len;//lenght of the m_r array; corresponds to the dimension of the space where the field lives
    double *m_data;
    unsigned int m_data_len;//lenght of the m_data array; corresponds to the number of points of the field times the dimension of the space
};

#endif // SCALARFIELD_H
