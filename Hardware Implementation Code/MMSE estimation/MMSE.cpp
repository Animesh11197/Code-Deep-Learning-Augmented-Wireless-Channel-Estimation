#include"MMSE.h"
#include <iostream>
#include <cstdio>
#include<stdio.h>
#include<stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <inttypes.h>
#include <sds_lib.h>
#include "sds_lib.h"
#include "ap_int.h"
#include <ap_fixed.h>
using namespace std;


//LS_hw:
/*
void LS_hw(double input[24*2*2], double output[24*2*2])
{
	data_t local_input [24*2*2];
	data_t local_output[24*2*2];


LOAD:for(int i=0;i<24*2*2;i++)
	{
		#pragma HLS PIPELINE
		local_input[i]=input[i];
	}

#pragma HLS ARRAY_PARTITION variable=local_input complete dim=0
#pragma HLS ARRAY_PARTITION variable=local_output complete dim=0


    for (int i = 0; i < 24*2; i++)
    {
#pragma HLS PIPELINE
        local_output[i] = (local_input[i] + local_input[i + (24 * 2)]);

        local_output[i + (24 * 2)]= ( local_input[i + (24 * 2)] - local_input[i]);

    }

OUT:for(int i=0;i<24*2*2;i++)
    	{
    #pragma HLS PIPELINE
    		output[i]=local_output[i]/2;
    	}


}
*/




//function to multiply Matrices
void MatrixMult(double matrix_1[N*N], double matrix_2[N*N], int p, int l, int m, int n, double matrix_product[N*N]) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {             // not j<M
			matrix_product[i * N + j] = 0;
			for (int k = 0; k < l; k++) {
				matrix_product[i * N + j] += matrix_1[i * N + k] * matrix_2[k * N + j];
			}
		}
	}
	return;

}


bool lup(double Ai[DIM*DIM], double Lo[DIM*DIM], double U[DIM*DIM], double Pe[DIM+1])
{
#pragma HLS inline

	double A[DIM*DIM];
	for(int i=0;i<DIM*DIM;i++)
	{
		A[i]=Ai[i];

	}
PermutMat_Initialize: for (int i = 0; i < DIM; i++)
{
	//#pragma HLS PIPELINE
	Pe[i] = i;
}

int i, j, k;
bool singular = 0;		/*-1 if matrix is singular --> inverse doesn't exist.*/

/*LUP Decomposition is finding L,U,P such that
 * PA=LU, here
 * P-Row Permutation Matrix
 * A-Input Matrix
 * L-Lower Triangular Matrix
 * U-Upper Triangular Matrix
 * Then the inverse for A is-
 * inv(A)=inv(U)*inv(L)*P
 * */
lup_label0: for (i = 0; i < DIM; i++)
{
	//#pragma HLS PIPELINE
						  /* pos-index for maximum value in a column.*/
	int pos = i;

	/*stores the maximum value, initially taken as follows*/
	double max = A[i * DIM + i];

	/*pivoting
	 * Find the maximum value in the column below the diagonal element (for the column)
	 * Swap the row with the one with max value.
	 * */
find_max: for (k = i + 1; k < DIM; k++)
{
	//#pragma HLS PIPELINE
	//#pragma HLS LOOP_TRIPCOUNT min=1 max=3
	double tmp = A[k * DIM + i];
	if (tmp > A[i * DIM + i] & tmp > max)
	{
		pos = k;
		max = tmp;
	}
}
/*check if the matrix is singular
 * */
if (A[pos * DIM + i] == 0.0)
{
	singular = 1; 		//matrix is singular
	return singular;
}
else
{
	if (pos != i)
	{
	swap_row: for (k = 0; k < DIM; k++)
	{

		double tmp = A[pos * DIM + k];
		A[pos * DIM + k] = A[i * DIM + k];
		A[i * DIM + k] = tmp;
	}

	/*update the permutation matrix for the swap*/
	int ind1, ind2;
	for (int i1 = 0; i1 < DIM; i1++)
	{
		if (Pe[i1] == pos)
			ind1 = i1;
		if (Pe[i1] == i)
			ind2 = i1;
	}
	double temp = Pe[ind2];
	Pe[ind2] = Pe[ind1];
	Pe[ind1] = temp;
	}
}

/*extract the L and U Matrices
 * */
lup_label1: for (k = i + 1; k < DIM; k++)
{
	A[k * DIM + i] = A[k * DIM + i] / A[i * DIM + i];
lup_label2: for (j = i + 1; j < DIM; j++)
{
	A[k * DIM + j] = A[k * DIM + j] - A[i * DIM + j] * A[k * DIM + i];
}
}
}
//L matrix: Lower half of A matrix and diagonal elements is 1.
Assign_L0: for (int i = 0; i < DIM; i++)
{
Assign_L1: for (int j = i; j < DIM; j++)
{
	if (i == j)
		Lo[j * DIM + i] = 1;
	else
		Lo[j * DIM + i] = A[j * DIM + i];
}
}

//U matrix: Upper half of A matrix.
Assign_U0: for (int i = 0; i < DIM; i++)
{
Assign_U1: for (int j = i; j < DIM; j++) {
	U[i * DIM + j] = A[i * DIM + j];
}
}

Pe[DIM] = (singular == 1) ? -1 : 0;

return 0;
}

void Lower_inv(double Lo[DIM*DIM], double L_inv[DIM*DIM])
{
#pragma HLS inline
	/* To calculate the inverse of L matrix:
	 * We have,
	 * for i==j, L_inv(i,j)=1/L(i,j)
	 * 						  k=i
	 * for i>j,  L_inv(i,j)=-Summation{L(i,k)*L_inv(k,j)}
	 * 						  k=j
	 * */
	double local_L_inv[DIM*DIM];


linv_label0: for (int i = 0; i < DIM; i++)
	{
		linv_label1: for (int j = 0; j < DIM; j++)
		{

			if (i < j)
				local_L_inv[i * DIM + j] = 0;
			else if (i == j)
				local_L_inv[i * DIM + j] = 1 / Lo[i * DIM + j];
			else
				{
					double sum = 0.0f;
				linv_label2: for (int k = j; k < i; k++)



					//#pragma HLS LOOP_TRIPCOUNT min=1 max=3

					//#pragma HLS PIPELINE
					sum = sum + Lo[i * DIM + k] * local_L_inv[k * DIM + j];

				local_L_inv[i * DIM + j] = -sum;
				}
		}
	}

for(int i=0;i<DIM*DIM;i++)
	{
		L_inv[i]=local_L_inv[i];
	}

}

void Upper_inv(double U[DIM*DIM], double U_inv[DIM*DIM])
{
#pragma HLS inline
	double local_U_inv[DIM*DIM];

	/* Inverse Of Upper Triangular Matrix
	 * We have,
	 * for i==j, U(i,j)=1/U(i,j)
	 * 						  k=i
	 * for i>j   U_inv(j,i)=[-Summation{L(i,k)*L_inv(k,j)} ]/U[i][i]
	 * 						  k=j
	 **/
	uinv_label10: for (int i = 0; i < DIM; i++)
		{
			//#pragma HLS PIPELINE
			uinv_label11: for (int j = 0; j < DIM; j++)
				{
					local_U_inv[i * DIM + j] = 0;
				}
		}

	uinv_label20: for (int i = 0; i < DIM; i++)
		{
			//#pragma HLS PIPELINE
			local_U_inv[i * DIM + i] = (double)1 / U[i * DIM + i];
		}

	univ_label30: for (int i = 0; i < DIM; i++)
		{
			univ_label31: for (int j = 0; j < i; j++)
				{
					//#pragma HLS PIPELINE
					//#pragma HLS LOOP_TRIPCOUNT min=1 max=3
					double sum = 0.0f;
					univ_label32: for (int k = j; k < i; k++)
						{
							//#pragma HLS PIPELINE
							//#pragma HLS LOOP_TRIPCOUNT min=1 max=3
							sum += (U[k * DIM + i] * local_U_inv[j * DIM + k]);
						}
					local_U_inv[j * DIM + i] = -sum / U[i * DIM + i];
				}
		}

	for(int i=0;i<DIM*DIM;i++)
	{
		U_inv[i]=local_U_inv[i];
	}

}

void matrix_mult(double U_inv[DIM*DIM], double L_inv[DIM*DIM], double A_inv[DIM*DIM])
{
#pragma HLS inline
	double sumTemp = 0;
MM_L1: for (int i = 0; i < DIM; i++)
{
MM_L2: for (int j = 0; j < DIM; j++)
{
	//#pragma HLS PIPELINE
	double sumFinal = 0;
MM_L3: for (int k = 0; k < DIM; k++)
{
	sumTemp = U_inv[i * DIM + k] * L_inv[k * DIM + j];
	sumFinal += sumTemp;
}
A_inv[i * DIM + j] = sumFinal;
}
}
}

void final_perm(double UL_inv[DIM*DIM], double Pe[DIM*DIM], double A_inv[DIM*DIM])
{
#pragma HLS inline
	/*Multiplication with P is just the row permuation*/
L1: for (int i = 0; i < DIM; i++)
{
L2: for (int j = 0; j < DIM; j++)
{
	A_inv[j * DIM + i] = UL_inv[j * DIM + (int)Pe[i]];
}
}
}

bool inverse(double A[DIM*DIM], double A_inv[DIM*DIM])
{
#pragma HLS inline
	double Ac [DIM * DIM];
	for (int i = 0; i < DIM * DIM; i++)
	{
		Ac[i] = A[i];
	}
	double Lo [DIM * DIM];
	double U [DIM * DIM];
	double L_inv [DIM * DIM];
	double U_inv [DIM * DIM];
	double UL_inv [DIM * DIM];
	double Pe [DIM + 1];

	// calculate the P,A,L,U such that PA=LU
	bool singular = lup(Ac, Lo, U, Pe);

	if (singular)
	{
		std::cout << "SINGULAR\n";
		for (int i = 0; i < DIM; i++)
		{
			//#pragma HLS PIPELINE
			for (int j = 0; j < DIM; j++)
				A_inv[i * DIM + j] = 0;
		}
		return singular;
	}
	else
	{
		//cout << "///////////////////////////n";
		Lower_inv(Lo, L_inv);
		Upper_inv(U, U_inv);
		matrix_mult(U_inv, L_inv, UL_inv);
		final_perm(UL_inv, Pe, A_inv);


		/*printf("\n\nPrinting L_inverse Matrix\n");
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				printf("%f ", L_inv[i*DIM+j]);
			}
			printf("\n");
		}
		printf("\n\nPrinting U_inverse Matrix\n");
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				printf("%f ", U_inv[i*DIM+j]);
			}
			printf("\n");
		}
		printf("\n\nPrinting LU_inverse Matrix\n");
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				printf("%f ", UL_inv[i*DIM+j]);
			}
			printf("\n");
		}*/
		/*printf("\n\nPrinting A_inverse Matrix\n");
		for (int i = 0; i < DIM; i++)
		{
			for (int j = 0; j < DIM; j++)
			{
				printf("%f ", A_inv[i*DIM + j]);
			}
			printf("\n");
		}*/

	}

	return 0;
}

// Adds matrix
void Add(double mat1[DIM*DIM], double mat2[DIM*DIM], int p, int l, double mat12[DIM*DIM])
{

	for (int i = 0; i < DIM * DIM; i++)//for (int i = 0; i < p * l; i++)
	{
		mat12[i] = mat2[i] + mat1[i];
		//cout<<mat12[i]<<"  ";
	}
}

// Neg matrix function
void neg(double neg_inv[DIM*DIM], int m, int n, double inv[DIM*DIM])
{
	for (int i = 0; i < N * N; i++)//for (int i = 0; i < m * n; i++)
	{
		inv[i] = -1 * neg_inv[i];
	}
}

//finds the complex inverse of matrix
void complex_inverse(double port_A[72*24], double port_Ai[72*24], int n, double inv[72*24], double inv_i[72*24])
{
	double A[N*N],Ai[N*N];
	for(int i=0;i<DIM*DIM;i++)
	{
		A[i]=port_A[i];
		Ai[i]=port_Ai[i];
	}
	double local_inv[N*N],local_inv_i[N*N];

	double neg_inv_i[N * N];
	double A_i[N * N];
	double BA_i_B[N * N];
	double A_BA_i_B [N * N];
	double BA_i [N * N];

	inverse(A, A_i);
	//display(A_i, n, n);
//display(A,n,n);
	MatrixMult(A_i, Ai, N, N, N, N, BA_i);

	MatrixMult(Ai, BA_i, N, N, N, N, BA_i_B);
	Add(A, BA_i_B, N, N, A_BA_i_B);
	inverse(A_BA_i_B, local_inv);
	//display(A_BA_i_B,n,n);
	//display(inv, n, n);
	MatrixMult(BA_i, local_inv, N, N, N, N, neg_inv_i);
	neg(neg_inv_i, N, N, local_inv_i);
	//display(inv_i, n, n);
	for(int i=0;i<DIM*DIM;i++)
	{
		inv[i]=local_inv[i];
		inv_i[i]=local_inv_i[i];
	}
}



void LS(double input[72*24], double input_i[72*24], double output[72*24], double output_i[72*24])
{
	for (int i = 0; i < 24 * 1; i++)
	{
		output[i] = 0.5 * (input[i] + input_i[i]);

		output_i[i] = 0.5 * (input_i[i] - input[i]);

	}
}

//function to transpose Matrix
void Transpose(double matrix[72*24], double matrix_i[72*24], int p, int l, double t_matrix[72*24], double t_matrix_i[72*24]) {
	for (int i = 0; i < p; i++) {
		for (int j = 0; j < l; j++)
		{
			t_matrix[j * p + i] = matrix[i * l + j];
			t_matrix_i[j * p + i] = -1 * matrix_i[i * l + j];
		}
		//cout << endl;
	}
	return;
}


// Adding a conatant @diagonal locations:
void Add_Mat_d(double mat[72*24], int p, double Var) // for sqaure Matrix
{
	for (int i = 0; i < p; i++)
	{
		mat[i * p + i] = mat[i * p + i]+ Var;
		//mat_i[i * p + i] += Var_i;
	}

}

//complex matrix multiplication function
void complexMatrixMult(double matrix_1[72*24], double matrix_1_i[72*24], double matrix_2[72*24], double matrix_2_i[72*24], int p, int l, int m, int n, double matrix_product[72*24], double matrix_product_i[72*24]) {
	for (int i = 0; i < p; i++) {
		for (int j = 0; j < n; j++) {             // not j<M
			matrix_product[i * n + j] = 0;
			matrix_product_i[i * n + j] = 0;
			for (int k = 0; k < l; k++) {
				matrix_product[i * n + j] += (matrix_1[i * l + k] * matrix_2[k * n + j]) - (matrix_1_i[i * l + k] * matrix_2_i[k * n + j]);
				matrix_product_i[i * n + j] += (matrix_1[i * l + k] * matrix_2_i[k * n + j]) + (matrix_1_i[i * l + k] * matrix_2[k * n + j]);;
			}
		}
	}
	return;

}



void H_MMSE(double SNR_const,double RS_user[72*24],double RS_user_i[72*24],double RS[72*24],double RS_i[72*24], double rec_pilots[72*24], double rec_pilots_i[72*24],double MMSE_Out[72*24],double MMSE_Out_i[72*24])
{
	//double* H_ref, * H_ref_i;
	//H_ref = (double*)malloc((24 * 1) * sizeof(double));//H_ref = Rx symbol without noise
	//H_ref_i = (double*)malloc((24 * 1) * sizeof(double));
	double H_ref[72*24],H_ref_i[72*24];


	//H_Pilot = H_Ref;	//H_Ref=Recieved_pilot./Pliot_value_user

/*	double* Rhh, * Rhh_i, * Rs, * R_i2;	// Rhh = H_Pilot * H_Pilot'
	Rhh = (double*)malloc((24 * 24) * sizeof(double));
	Rhh_i = (double*)malloc((24 * 24) * sizeof(double));*/
	double Rhh[72*24],Rhh_i[72*24];


/*	double H_LS[24 * 1], H_LS_i[24 * 1];*/
	double H_LS[72*24],H_LS_i[72*24];


	LS(RS_user, RS_user_i, H_ref, H_ref_i);

////	double t_H_ref[1 * 24], t_H_ref_i[1 * 24];
	double t_H_ref[72 * 24], t_H_ref_i[72 * 24];
	Transpose(H_ref, H_ref_i, 24, 1, t_H_ref, t_H_ref_i);
	//Transpose(H_ref_i, 24, 1, t_H_ref_i);

	double Mult_Out[72 * 24], Mult_Out_i[72 * 24];

	complexMatrixMult(RS, RS_i, t_H_ref, t_H_ref_i, 72, 1, 1, 24, Mult_Out, Mult_Out_i);// Find the First Term (autocorelation): CIR * sym_wo_noise

	complexMatrixMult(H_ref, H_ref_i, t_H_ref, t_H_ref_i, 24, 1, 1, 24, Rhh, Rhh_i);		//Finds Rhh
	Add_Mat_d(Rhh, 24, SNR_const);

/*	double Rhh_inv[24 * 24], Rhh_inv_i[24 * 24];*/
	double Rhh_inv[72 * 24], Rhh_inv_i[72 * 24];

	complex_inverse(Rhh, Rhh_i, 24, Rhh_inv, Rhh_inv_i);// Find the middle Term: INV[Rhh+sigma^2.I]

	double comb_out[72 * 24], comb_out_i[72 * 24];
	complexMatrixMult(Mult_Out, Mult_Out_i, Rhh_inv, Rhh_inv_i, 72, 24, 24, 24, comb_out, comb_out_i);

	LS(rec_pilots, rec_pilots_i, H_LS, H_LS_i); // The 3rd Term : H_LS


	complexMatrixMult(comb_out, comb_out_i, H_LS, H_LS_i, 72, 24, 24, 1, MMSE_Out, MMSE_Out_i);



}

