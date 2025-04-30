#include "calc.h"

void calc_deviation_vector( const double input_vector0[2], const double input_vector1[2], double output_vector[2])
{
    output_vector[0] = input_vector1[0] - input_vector0[0];
    output_vector[1] = input_vector1[1] - input_vector0[1];
    return;
}

void calc_unit_vector(const double input_vector[2], double output_vector[2])
{
    double n = calc_norm(input_vector);
    
    output_vector[0] = input_vector[0] / n;
    output_vector[1] = input_vector[1] / n;
    
    return;
}

double calc_norm(const double input_vector[2])
{
    double norm = sqrt(input_vector[0] * input_vector[0] + input_vector[1] * input_vector[1]);
    return norm;
}