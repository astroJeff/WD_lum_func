#ifndef MATRIX_H
#define MATRIX_H


void invert_3(double XY[3][3], double XY_in[3][3]);
void matrix_zero_3(double XY[3][3]);
void matrix_iden_3(double XY[3][3]);
double determinant_3(double XY[3][3]);
void matrix_mult_3(double A[3], double z[3], double matrix[3][3]);
double interp_3(double x,double A[3]);
void shift_rows_3(double XY[3][3], int start);
void find_coeff_3(double A[3],double x[3],double z[3]);
void print_matrix_3(double matrix_in[3][3]);

#endif
