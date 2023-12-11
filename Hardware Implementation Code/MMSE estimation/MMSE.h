#ifndef MMSE_H_
#define MMSE_H_
#include <iostream>
#include <stdlib.h>
#include <assert.h>
const int N = 24;
const int DIM = 24;

void MatrixMult(double matrix_1[N*N], double matrix_2[N*N], int p, int l, int m, int n, double matrix_product[N*N]);
bool lup(double A[DIM*DIM], double Lo[DIM*DIM], double U[DIM*DIM], double Pe[DIM*DIM]);
void Lower_inv(double Lo[DIM*DIM], double L_inv[DIM*DIM]);
void Upper_inv(double U[DIM*DIM], double U_inv[DIM*DIM]);
void matrix_mult(double U_inv[DIM*DIM], double L_inv[DIM*DIM], double A_inv[DIM*DIM]);
void final_perm(double UL_inv[DIM*DIM], double Pe[DIM*DIM], double A_inv[DIM*DIM]);
bool inverse(double A[DIM*DIM], double A_inv[DIM*DIM]);


void Add(double mat1[DIM*DIM], double mat2[DIM*DIM], int p, int l, double mat12[DIM*DIM]);
void neg(double neg_inv[DIM*DIM], int m, int n, double inv[DIM*DIM]);

//#pragma SDS data mem_attribute(A:PHYSICAL_CONTIGUOUS, Ai:PHYSICAL_CONTIGUOUS,inv: PHYSICAL_CONTIGUOUS,inv_i: PHYSICAL_CONTIGUOUS)
void complex_inverse(double A[N*N], double Ai[N*N], int n, double inv[N*N], double inv_i[N*N]);

void LS(double input[72*24], double input_i[72*24], double output[72*24], double output_i[72*24]);

void Transpose(double matrix[72*24], double matrix_i[72*24], int p, int l, double t_matrix[72*24], double t_matrix_i[72*24]);

// Adding a conatant @diagonal locations:
void Add_Mat_d(double mat[72*24], int p, double Var) ;

void complexMatrixMult(double matrix_1[72*24], double matrix_1_i[72*24], double matrix_2[72*24], double matrix_2_i[72*24], int p, int l, int m, int n, double matrix_product[72*24], double matrix_product_i[72*24]);

void H_MMSE(double SNR_const,double RS_user[72*24],double RS_user_i[72*24],double RS[72*24],double RS_i[72*24], double rec_pilots[72*24], double rec_pilots_i[72*24],double MMSE_Out[72*24],double MMSE_Out_i[72*24]);

#endif
