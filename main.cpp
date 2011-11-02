#include <config.h>

#include <iostream>
#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <cstdio>

struct Nanowire
{
    double x, y;  
};

void Readin_coord(std::vector<Nanowire>& nanowires);
gsl_matrix* Coord_matrix(const std::vector<Nanowire>& nanowires);
gsl_matrix* Inner_distance(const gsl_matrix* coord);
gsl_vector* Charge_distri(double phi, double R, const gsl_matrix* inner_dist);

int main()
{
    std::vector<Nanowire> nanowires;
    Readin_coord(nanowires);
    
    gsl_matrix* coord_matrix = Coord_matrix(nanowires);
    gsl_matrix* inner_distance = Inner_distance(coord_matrix);
    
    gsl_vector* Qs = Charge_distri(1, 1e-10, inner_distance);
    
    gsl_vector_fprintf(stdout, Qs, "%g");
}
