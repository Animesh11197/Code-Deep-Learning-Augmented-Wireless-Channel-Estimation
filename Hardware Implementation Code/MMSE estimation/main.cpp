// including required Header Files and macros
//#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <sds_lib.h>
#include "sds_lib.h"
#include <math.h>
#include <time.h>
#include "mmse.h"
#include <cmath>

int num_of_frames_each_SNR = 10;
//SNR_const[i]=(1 / 10^(SNR / 10))
double SNR_const[7] = { 3.1623, 1, 0.3162, 0.1, 0.0316, 0.0100, 0.0032 };

using namespace std;


void interpolate(double* original_img, double* output)
{
    //printf("\nTHIS FUNCTION WILL PERFORM BILINEAR INTERPOLATION.!\n");
    //input dimension is 24*2
    //ouput dimension is 72*14
    double old_h = 72;
    double old_w = 2;
    double new_h = 72;
    double new_w = 14;
    double h_scale_factor = old_h / new_h;//old_h/new_h
    double w_scale_factor = old_w / new_w; //old_w / new_w
    for (int channel = 0; channel < 2; channel++)
    {

        for (int i = 0; i < new_h; i++)
        {
            for (int j = 0; j < new_w-3; j++)
            {
                //printf("\n*********************************************\n");
                //cout<<"i="<< i<< "\tj="<< j<<endl;
                //map the coordinates back to the original image
                double x = i * h_scale_factor;
                double y = j * w_scale_factor;
                //cout<<"x="<<x <<"\ty="<< y<<endl;
                //calculate the cordinates of 4 surrounding pixels
                int x_floor = floor(x);
                int x_ceil = min(old_h - 1, ceil(x));
                int y_floor = floor(y);
                int y_ceil = min(old_w - 1, ceil(y));
                double q;
                //cout<<"x_floor= "<<x_floor<<"\tx_ceil="<<x_ceil<<"\ty_floor="<<y_floor<<"\ty_ceil="<<y_ceil ;

                if ((x_ceil == x_floor) && (y_ceil == y_floor))
                {
                    q = original_img[(channel * 24 * 2) + int(x) * 2 + int(y)];
                    //printf("\nq =%f\n", q);
                }
                else if (x_ceil == x_floor)
                {
                    double q1 = original_img[(channel * 24 * 2) + int(x) * 2 + int(y_floor)];
                    double q2 = original_img[(channel * 24 * 2) + int(x) * 2 + int(y_ceil)];
                    q = q1 * (y_ceil - y) + q2 * (y - y_floor);
                    //cout<<"\nq1="<<q1<<"\tq2="<<q2<<"\tq="<<q;
                }
                else if (y_ceil == y_floor)
                {
                    double q1 = original_img[(channel * 24 * 2) + int(x_floor) * 2 + int(y)];
                    double q2 = original_img[(channel * 24 * 2) + int(x_ceil) * 2 + int(y)];
                    q = (q1 * (x_ceil - x)) + (q2 * (x - x_floor));
                    //cout<<"\nq1="<<q1<<"\tq2="<<q2<<"\tq="<< q;
                }
                else
                {
                    double v1 = original_img[(channel * 24 * 2) + x_floor * 2 + y_floor];
                    double v2 = original_img[(channel * 24 * 2) + x_ceil * 2 + y_floor];
                    double v3 = original_img[(channel * 24 * 2) + x_floor * 2 + y_ceil];
                    double v4 = original_img[(channel * 24 * 2) + x_ceil * 2 + y_ceil];
                    //cout<<"\nv1="<<v1<<"\tv2="<<v2<<"\tv3="<<v3<<"\tv4="<<v4;

                    double q1 = v1 * (x_ceil - x) + v2 * (x - x_floor);
                    double q2 = v3 * (x_ceil - x) + v4 * (x - x_floor);
                    q = q1 * (y_ceil - y) + q2 * (y - y_floor);
                    //cout<<"\nq1="<<q1<<"\tq2="<<q2<<"\tq="<<q;

                }

                output[((i) * 14) + j + 3] = q;

            }
        }


        for (int i = new_h ; i >= 0; i--)
        {
            for (int j = 3; j > 0; j--)
            {
                output[(i * 14) + j - 1] = output[i * 14 + j];
                //output[(i * 14) + j - 1] = 0;

            }
        }


    }
    //return output;
}



int main()
{
	//initialising required Matrices
	double* RS, * RS_i, * RS_user, * RS_user_i;	// RS = CIR: Channel impulse Response!
	RS = (double*)malloc((72 * 24) * sizeof(double));
	RS_i = (double*)malloc((72 * 24) * sizeof(double));
	RS_user = (double*)malloc((24 * 72) * sizeof(double));//RS_User = Rx symbol signal without noise
	RS_user_i = (double*)malloc((24 * 72) * sizeof(double));

	double* H_ref, * H_ref_i;
	H_ref = (double*)malloc((24 * 1) * sizeof(double));//H_ref = Rx symbol without noise
	H_ref_i = (double*)malloc((24 * 1) * sizeof(double));


	//H_Pilot = H_Ref;	//H_Ref=Recieved_pilot./Pliot_value_user
	double* Rhh, * Rhh_i, * Rs, * R_i2;	// Rhh = H_Pilot * H_Pilot'
	Rhh = (double*)malloc((24 * 24) * sizeof(double));
	Rhh_i = (double*)malloc((24 * 24) * sizeof(double));

	double* rec_pilots, * rec_pilots_i;// Unrecovered Signal
	rec_pilots = (double*)malloc((24 * 72) * sizeof(double));		//Urecover_sig = Unrecovered Signal
	rec_pilots_i = (double*)malloc((24 * 72) * sizeof(double));
	double H_LS[24*1],H_LS_i[24*1];


	//int* golden_out;
	//golden_out = (int*)malloc(14 * sizeof(int));

	////int *golden_out = new int[14];
	//int iter, data_points;
	//int* out, * out2;
	//out = (int*)malloc(N * sizeof(int));
	//out2 = (int*)malloc(N * sizeof(int));


	//int out[1 * N];int out2[1 * N];
	double count = 0; double count2 = 0;

	//File Management:
	FILE* fin, * fin_i, * fin2, * fin2_i, * fin3, * fin4, * fin5;
	double inputsample, inputsample_i, inputsample2, inputsample2_i;
	double inputsample3, inputsample4, inputsample5;

	fin = fopen("RS_72x2_real.dat", "r");
	if (!fin) {
		printf("could not open file \a\n");
		exit(101);
	}
	fin_i = fopen("RS_72x2_imag.dat", "r");
	if (!fin_i) {
		printf("could not open file \a\n");
		exit(101);
	}
	fin2 = fopen("RS_user_24x2_real.dat", "r");
	 if (!fin2) {
	printf("could not open file \a\n");
	exit(101);
	}
	 fin2_i = fopen("RS_user_24x2_imag.dat", "r");
	 if (!fin2_i) {
		 printf("could not open file \a\n");
		 exit(101);
	 }
	 fin3 = fopen("Received_pilot_LS_24x2_real.dat", "r");
	 if (!fin3) {
		 printf("could not open file \a\n");
		 exit(101);
	 }
	 fin4 = fopen("Received_pilot_LS_24x2_imag.dat", "r");
	 if (!fin4) {
		 printf("could not open file \a\n");
		 exit(101);
	 }
		fin5 = fopen("test_golden_output.dat", "r");
		 if (!fin5) {
			 printf("could not open file \a\n");
			 exit(101);
		 }

	 double golden_out[72 * 14 * 2];

	 double MMSE_Out[72 * 2]; double MMSE_Out_i[72 * 2];

	 //Counter for calculation of execution time
	     uint64_t cnt_start = 0;
	     uint64_t cnt_stop = 0;
	     //uint64_t frequency = sds_clock_frequency();


	     //Separate counter for SW and HW functions
	     uint64_t total_count_sw = 0;
	     uint64_t total_count_hw = 0;

	     double sum_sw,sum_hw;




	 // for (int snr = -5; snr <= -5; snr += 5) {
	 for (int snr = 0; snr <= 6; snr += 1) {
		//double sum_sw = 0;
		//double  sum_hw = 0;
	 	double sum = 0;
	 	for (int frame = 0; frame < num_of_frames_each_SNR; frame++)
	 	{
	 		 //Reading Golden Out
	 		 for (int j = 0; j < 72; j++)
	 		 {
	 			 for (int k = 0; k < 2; k++) {
	 				 for (int i = 0; i < 14; i++) {
	 					 fscanf(fin5, "%f", &inputsample5);
	 					 golden_out[i + j * 14 + k * (14 * 72)] = inputsample5;
	 					 //weight_buffer[l * (3 * 2 * 8)] = inputsample2;
	 					 //printf("%f\t", inputsample1);
	 				 	 }
	 				 //std::cout << "\n";
	 			 	 }
	 			 //std::cout << "\n";
	 			 //std::cout << "\n------\n";
	 		 }

			 for (int pilot_sym = 0; pilot_sym < 2; pilot_sym++)
				 {

					 //  fin2_i = fopen("data_R_i.dat", "r"); fin3 = fopen("Truth.dat", "r"); fin4 = fopen("K.dat", "r");
					 //fin5 = fopen("data_points.dat", "r");

					 //Reading RS
					 for (int i = 0; i < 72 * 1; i++)
					 {
						 fscanf(fin, "%f", &inputsample); fscanf(fin_i, "%f", &inputsample_i);
						 RS[i] = inputsample;		RS_i[i] = inputsample_i;

						 //cout << inputsample << "\t" << inputsample_i<<endl;
						 //A2[i] = inputsample;		A_i2[i] = inputsample_i;

					 }
					 //Reading RS_user
					 for (int i = 0; i < 24 * 1; i++)
					 {
						 fscanf(fin2, "%f", &inputsample2); fscanf(fin2_i, "%f", &inputsample2_i);
						 RS_user[i] = inputsample2;		RS_user_i[i] = inputsample2_i;

						 fscanf(fin3, "%f", &inputsample3); fscanf(fin4, "%f", &inputsample4);
						 rec_pilots[i] = inputsample3;		rec_pilots_i[i] = inputsample4;

						 //cout << rec_pilots[i] << "\t" << rec_pilots_i[i] <<endl;
						 //A2[i] = inputsample;		A_i2[i] = inputsample_i;

					 }

					 //LS(RS_user, RS_user_i, H_ref, H_ref_i);

					 //double t_H_ref[1 * 24], t_H_ref_i[1 * 24];
					 //Transpose(H_ref, H_ref_i, 24, 1, t_H_ref, t_H_ref_i);
					 ////Transpose(H_ref_i, 24, 1, t_H_ref_i);

					 //double Mult_Out[72 * 24], Mult_Out_i[72 * 24];

					 //complexMatrixMult(RS, RS_i, t_H_ref, t_H_ref_i, 72, 1, 1, 24, Mult_Out, Mult_Out_i);// Find the First Term (autocorelation): CIR * sym_wo_noise

					 //complexMatrixMult(H_ref, H_ref_i, t_H_ref, t_H_ref_i, 24, 1, 1, 24, Rhh, Rhh_i);		//Finds Rhh
					 //Add_Mat_d(Rhh, 24, 3.1623);
					 //
					 //double Rhh_inv[24 * 24], Rhh_inv_i[24 * 24];
					 //complex_inverse(Rhh, Rhh_i, 24,Rhh_inv,Rhh_inv_i);// Find the middle Term: INV[Rhh+sigma^2.I]

					 //double comb_out[72*24], comb_out_i[72*24];
					 //complexMatrixMult(Mult_Out, Mult_Out_i, Rhh_inv, Rhh_inv_i, 72, 24, 24, 24, comb_out, comb_out_i);

					 //LS(rec_pilots, rec_pilots_i, H_LS, H_LS_i); // The 3rd Term : H_LS
					 //
					 //double MMSE_Out[72 * 1], MMSE_Out_i[72 * 1];
					 //complexMatrixMult(comb_out, comb_out_i, H_LS,H_LS_i,72, 24, 24, 1, MMSE_Out, MMSE_Out_i);

			/*		 double H_MMSE_Out[72 * 1], H_MMSE_Out_i[72 * 1];*/
					 double H_MMSE_Out[72 * 24], H_MMSE_Out_i[72 * 24];

				cnt_start = sds_clock_counter(); //storing the time stamp of starting the sw function  call
					 H_MMSE(SNR_const[snr],RS_user, RS_user_i, RS, RS_i, rec_pilots, rec_pilots_i, H_MMSE_Out, H_MMSE_Out_i);
					 for (int i = 0; i < 72; i++)
					 {
						 MMSE_Out[i*2+pilot_sym] = H_MMSE_Out[i];
						 MMSE_Out_i[i * 2+ pilot_sym] = H_MMSE_Out_i[i];
					 }
				cnt_stop = sds_clock_counter(); //storing the time stamp of stoping the sw function call
				total_count_sw = total_count_sw + (cnt_stop - cnt_start);

				 }

				 double Out[72 * 14]; double Out_i[72 * 14];
				cnt_start = sds_clock_counter(); //storing the time stamp of starting the sw function  call
					interpolate(MMSE_Out,Out);
					interpolate(MMSE_Out_i, Out_i);
				cnt_stop = sds_clock_counter(); //storing the time stamp of stoping the sw function call
				total_count_sw = total_count_sw + (cnt_stop - cnt_start);


					 //CALCULATING MSE :
					 double squared_sum = 0;
					 double diff = 0; double diff2 = 0;
					 double mse = 0;
					 for (int i = 0; i < 72 * 14; i++)
					 {
						 diff = (golden_out[i] - Out[i]); diff2 = (golden_out[72 * 14 + i] - Out_i[i]);
						 diff = diff * diff; diff2 = diff2 * diff2;
						 squared_sum += diff; squared_sum += diff2;
					 }
					 mse = squared_sum / (72 * 14);

					 sum = sum + mse;//Acculumating MSE for All frames



				//DISPLAY H_ref
				for (int i = 0; i < 72; i++)
				{
					for (int j = 0; j < 14; j++)
					{
						//cout << Out[i * 14 + j] << "+i" << Out_i[i * 14 + j] << "\t";
					}
					//cout << endl<<endl;
				}

				for (int i = 0; i < 72; i++)
				{
					for (int j = 0; j < 24; j++)
					{
						//cout << comb_out[i * 24 + j] << "+i" << comb_out_i[i * 24 + j] << "\t";
					}
					//cout << endl << endl;
				}
	 	}
		 printf("\n\ MMSE  MSE_over_SNR for SNR %d: %f", (-5 + (snr*5)), sum / num_of_frames_each_SNR);
	 }


	    printf("\nExecution time in seconds: %f\n", (double)total_count_sw / (double)sds_clock_frequency());
	    //printf("Execution time of HW in seconds: %f\n", (double)total_count_hw / (double)sds_clock_frequency());
	//cout << (RS[1]*H_ref[1])- (RS_i[1] * H_ref_i[1]) << endl;


	//double A[16] = {
	//0.8444 ,   0.4219 ,   0.0684   , 0.9236,
	//0.1951,    0.0459  ,  0.4534  ,  0.3664,
	//0.0257 ,   0.5861   , 0.1941  ,  0.0198,
	//0.0209  ,  0.9590 ,   0.5401,    0.9453
	//};
	//double A_i[16] = {
	//0.0358,    0.6851,    0.7305,    0.7307,
	//0.8869,    0.9694,    0.2914,    0.1592,
	//0.9246,    0.1560,    0.0077,    0.3031,
	//0.1469,    0.8181,    0.9617,    0.4906
	//};
	//double A_inv[16]; double A_inv_i[16];
	//complex_inverse(A, A_i, 4, A_inv, A_inv_i);
	//for (int i = 0; i < 4; i++)
	//{
	//	for (int j = 0; j < 4; j++)
	//	{
	//		cout << A_inv[i * 4 + j] << "+i" << A_inv_i[i * 4 + j] << "\t";
	//	}
	//	cout << endl;
	//}

	//Transpose(H_ref, 24, 1, t_H_ref);


}
