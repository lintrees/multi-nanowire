#ifndef _COORDINATE_H_
#define _COORDINATE_H_

#include <istream>

#include <gsl/gsl_matrix.h>

class Coordinate_2d
{
    public:
        Coordinate_2d(std::istream&);
        ~Coordinate_2d();
        double x(size_t i) const
        {   return gsl_matrix_get(_m, i, 0); }
        double y(size_t i) const
        {   return gsl_matrix_get(_m, i, 1); }
        size_t size() const
        {   return _m->size1; }
        const double* data() const
        {   return _m->data; }
    private:    
        gsl_matrix* _m;
};

class Inner_distance
{
    public:
        Inner_distance(const Coordinate_2d&);
        ~Inner_distance();
        double operator()(size_t i, size_t j) const
        {   return gsl_matrix_get(_m, i, j);}
        size_t size() const
        {   return _m->size1; }
    private:    
        gsl_matrix* _m;
};


#endif //_COORDINATE_H_
