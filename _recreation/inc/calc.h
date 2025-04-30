#ifndef CALC_H_
#define CALC_H_

#include <math.h>

void calc_deviation_vector( const double input_vector0[2], const double input_vector1[2], double output_vector[2]);
void calc_unit_vector(const double input_vector[2], double output_vector[2]);
double calc_norm(const double input_vector[2]);

#endif