#include<iostream>
#include"matrix.h"

using namespace std;

void invert_3(double XY[3][3], double XY_in[3][3]){
  int i,j,k;
  int size=3;

  double factor;
  double singular;
  
  singular = determinant_3(XY);
  if(singular == 0 || singular != singular)  cerr << "There is a singular matrix!!" << endl;

  // Diagonalize the matrix
  for(i=0;i<size;i++){
    factor = XY[i][i];

    while(factor==0){  //To deal with singularities
      shift_rows_3(XY,i);
      shift_rows_3(XY_in,i);
      factor = XY[i][i];
    }

    for(k=0;k<size;k++){
      XY[i][k] /= factor;  //Set Diagonal element to 1
      XY_in[i][k] /= factor;  //Do the same thing with the inverse
    }

    for(j=0;j<size;j++){
      if (j==i) continue;  //Only do this for other rows

      factor = XY[j][i]/XY[i][i];

      for(k=0;k<size;k++){
      	XY[j][k] = XY[j][k] - factor*XY[i][k];
      	XY_in[j][k] = (XY_in[j][k]) - (factor)*(XY_in[i][k]);
      }
    }	
  }
}

void matrix_zero_3(double XY[3][3]){
  int i,j;

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      XY[i][j]=0.0;
    }
  }

}

void matrix_iden_3(double XY[3][3]){
  int i,j;

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      if(i==j) XY[i][j]=1.0;
      if(i!=j) XY[i][j] = 0;
    }
  }
}

double determinant_3(double XY[3][3]){
  int i,j,k;
  int count;
  
  double final;
  double det,factor;
  double matrix[3][3];

  final = 1.0;

  for(i=0;i<3;i++) for(j=0;j<3;j++) matrix[i][j] = XY[i][j];
  
  for(i=0;i<3;i++){
    for(j=i+1;j<3;j++){
      if(matrix[j][i] == 0) continue;

      count = 0;
      while(matrix[i][i] == 0){
	shift_rows_3(matrix,i);
	count ++;
	if (count > i) return 0; //Matrix is singular
      }

      factor = matrix[i][i]/matrix[j][i];
      final *= factor;

      for(k=i;k<3;k++){
	matrix[j][k] = matrix[j][k]*factor-matrix[i][k];
      }

    }

  }


  det = 1.0;
  for(i=0;i<3;i++) det *= matrix[i][i];

  det = det/final;

  return det;
}

void matrix_mult_3(double A[3], double z[3], double matrix[3][3]){

  int i,j;
  
  for(i=0;i<3;i++){
    A[i] = 0.0;
    for(j=0;j<3;j++){
      A[i] += z[j]*matrix[i][j];
    }
  }

}


double interp_3(double x,double A[3]){

  double out=0.0;

  out += A[0];
  out += A[1]*x;
  out += A[2]*x*x;

  return out;
}

void shift_rows_3(double XY[3][3], int start){

  int i,j;

  double temp;

  for(i=start;i<2;i++){
    for(j=0;j<3;j++){
      temp = XY[i][j];
      XY[i][j] = XY[i+1][j];
      XY[i+1][j] = temp;
    }
  }

}



void find_coeff_3(double A[3],double x[3],double z[3]){
  int i,j;

  double matrix[3][3];
  double matrix_temp[3][3];

  matrix_iden_3(matrix_temp);
  matrix_zero_3(matrix);

  for(i=0;i<3;i++){
    matrix[i][0] = 1.0;
    matrix[i][1] = x[i];
    matrix[i][2] = x[i]*x[i];
  }

  invert_3(matrix,matrix_temp);
  
  matrix_mult_3(A,z,matrix_temp);
}

void print_matrix_3(double matrix_in[3][3]){
  int i,j;

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      cout << matrix_in[i][j] << "\t";
    }
    cout << endl;
  }
}
