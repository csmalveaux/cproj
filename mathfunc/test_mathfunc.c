#include <stdio.h>
#include <stdint.h>

#include "dsp.h"
#include "complex.h"
#include "matrix.h"
#include "matfileReader.h"
//Header

#define NR (3)
#define NC (3)
#define NR2 (5)
#define NC2 (5)


int main () {

  int n_rows = 4;
  int n_cols = 3;
  //double res[n_cols][n_rows];
  double tol = 1e-5;
  int pass = 1;

  int dim[2] = {n_rows, n_cols};

  // double input[3][2] = {{1 , 2}, {3, 4}, {5, 6}};

  int i, j;
  // double ** input = calloc(n_rows, sizeof(double *));
  // for(i = 0; i < n_rows; i++)
  //   input[i] = calloc(n_cols, sizeof(double));
  // int inc = 1;
  // for(i = 0; i < n_rows; i++)
  //   for(j = 0; j < n_cols; j++){
  //     if(i == j)
  //       input[i][j] = 1.0;
  //     if(i+1 == j)
  //       input[i][j] = -1.0;
  //   }

  // double ** res = calloc(n_cols, sizeof(double *));
  // for(i = 0; i < n_cols; i++)
  //   res[i] = calloc(n_rows, sizeof(double));
  
  // printMatrixDouble("input", dim, input);
  
  // pinvd(n_rows, n_cols, input, res);
  // dim[0] = n_cols;
  // dim[1] = n_rows;
  // printMatrixDouble("res", dim, res);
  // // free(input);
  // // free(res);

  // double ** check = calloc(n_rows, sizeof(double *));
  // for(i = 0; i < n_rows; i++)
  //   check[i] = calloc(n_rows, sizeof(double));

  // multiplyMatrixDouble(n_rows, n_cols, input, n_cols, n_rows, res, check);
  // dim[0] = n_rows;
  // dim[1] = n_rows;
  // printMatrixDouble("check", dim, check);

/// Real Data

  n_cols = 100;
  n_rows = 400;

  int ** BB = calloc(n_rows, sizeof(int *));
  for(i = 0; i < n_rows; i++)
    BB[i] = calloc(n_cols, sizeof(int));

  dim[0] = n_rows;
  dim[1] = n_cols;
  
  matfile2Matrix2Dint("BB.mat", 0, dim, BB);
  //printMatrixInt("BB", dim, BB);

  double ** BBdouble = calloc(n_rows, sizeof(double *));
  for(i = 0; i < n_rows; i++)
    BBdouble[i] = calloc(n_cols, sizeof(double));

  for(i = 0; i < n_rows; i++)
    for(j = 0; j < n_cols; j++)
      BBdouble[i][j] = (double) BB[i][j];

  //printMatrixDouble("BBdouble", dim, BBdouble);
  double ** pinvBB = calloc(n_cols, sizeof(double *));
  for(i = 0; i < n_cols; i++)
    pinvBB[i] = calloc(n_rows, sizeof(double));

  pinvd(n_rows, n_cols, BBdouble, pinvBB);
  dim[0] = n_cols;
  dim[1] = n_rows;
  printMatrixDouble("pinvBB", dim, pinvBB);

    double ** check = calloc(n_rows, sizeof(double *));
  for(i = 0; i < n_rows; i++)
    check[i] = calloc(n_rows, sizeof(double));

  multiplyMatrixDouble(n_rows, n_cols, BBdouble, n_cols, n_rows, pinvBB, check);
  dim[0] = n_rows;
  dim[1] = n_rows;
  printMatrixDouble("check", dim, check);

  return 0;
}

