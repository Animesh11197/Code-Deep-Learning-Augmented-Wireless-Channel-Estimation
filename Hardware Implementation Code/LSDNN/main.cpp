//#define _CRT_SECURE_NO_WARNINGS
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
#include "DNN.h"

#define num_of_frames_each_SNR 10

//#pragma warning(disable : 4996)
using namespace std;

void LS(float* input, float* output)
{
    for (int i = 0; i < 24*2; i++)
    {
        output[i] = 0.5 * (input[i] + input[i + (24 * 2)]);

        output[i + (24 * 2)]= 0.5 * ( input[i + (24 * 2)] - input[i]);

    }

}


void interpolate(float* original_img, float* output)
{
    //printf("\nTHIS FUNCTION WILL PERFORM BILINEAR INTERPOLATION.!\n");
    //input dimension is 24*2
    //ouput dimension is 72*14
    float old_h = 24;
    float old_w = 2;
    float new_h = 72;
    float new_w = 14;
    float h_scale_factor = old_h / new_h;//old_h/new_h
    float w_scale_factor = old_w / new_w; //old_w / new_w
    for (int channel = 0; channel < 8; channel++)
    {

        for (int i = 0; i < new_h - 1; i++)
        {
            for (int j = 0; j < new_w - 3; j++)
            {
                //printf("\n*********************************************\n");
                //cout<<"i="<< i<< "\tj="<< j<<endl;
                //map the coordinates back to the original image
                float x = i * h_scale_factor;
                float y = j * w_scale_factor;
                //cout<<"x="<<x <<"\ty="<< y<<endl;
                //calculate the cordinates of 4 surrounding pixels
                int x_floor = floor(x);
                int x_ceil = min(old_h - 1, ceil(x));
                int y_floor = floor(y);
                int y_ceil = min(old_w - 1, ceil(y));
                float q;
                //cout<<"x_floor= "<<x_floor<<"\tx_ceil="<<x_ceil<<"\ty_floor="<<y_floor<<"\ty_ceil="<<y_ceil ;

                if ((x_ceil == x_floor) && (y_ceil == y_floor))
                {
                    q = original_img[(channel * 24 * 2) + int(x) * 2 + int(y)];
                    //printf("\nq =%f\n", q);
                }
                else if (x_ceil == x_floor)
                {
                    float q1 = original_img[(channel * 24 * 2) + int(x) * 2 + int(y_floor)];
                    float q2 = original_img[(channel * 24 * 2) + int(x) * 2 + int(y_ceil)];
                    q = q1 * (y_ceil - y) + q2 * (y - y_floor);
                    //cout<<"\nq1="<<q1<<"\tq2="<<q2<<"\tq="<<q;
                }
                else if (y_ceil == y_floor)
                {
                    float q1 = original_img[(channel * 24 * 2) + int(x_floor) * 2 + int(y)];
                    float q2 = original_img[(channel * 24 * 2) + int(x_ceil) * 2 + int(y)];
                    q = (q1 * (x_ceil - x)) + (q2 * (x - x_floor));
                    //cout<<"\nq1="<<q1<<"\tq2="<<q2<<"\tq="<< q;
                }
                else
                {
                    float v1 = original_img[(channel * 24 * 2) + x_floor * 2 + y_floor];
                    float v2 = original_img[(channel * 24 * 2) + x_ceil * 2 + y_floor];
                    float v3 = original_img[(channel * 24 * 2) + x_floor * 2 + y_ceil];
                    float v4 = original_img[(channel * 24 * 2) + x_ceil * 2 + y_ceil];
                    //cout<<"\nv1="<<v1<<"\tv2="<<v2<<"\tv3="<<v3<<"\tv4="<<v4;

                    float q1 = v1 * (x_ceil - x) + v2 * (x - x_floor);
                    float q2 = v3 * (x_ceil - x) + v4 * (x - x_floor);
                    q = q1 * (y_ceil - y) + q2 * (y - y_floor);
                    //cout<<"\nq1="<<q1<<"\tq2="<<q2<<"\tq="<<q;

                }

                output[(channel * 72 * 14) + ((i + 1) * 14) + j + 3] = q;

            }
        }


        for (int i = new_h - 1; i > 0; i--)
        {
            for (int j = 3; j > 0; j--)
            {
                output[(channel * 72 * 14) + (i * 14) + j - 1] = output[(channel * 72 * 14) + (i * 14 + j)];
                //output[(i * 14) + j - 1] = 0;

            }
        }

        for (int j = new_w - 1; j > -1; j--)
        {
            output[(channel * 72 * 14) + (0 * 14) + j] = output[(channel * 72 * 14) + (1 * 14 + j)];
            //output[(0 * 14) + j ] = 0;
        }

    }


    //test[0] = 10;
    //return output;
}


int main()
{
    //float input[24 * 2 * 2];
    float* input = (float*)sds_alloc(sizeof(float*) * 96);

    float input2[24 * 2 * 2];
    //float DNN_input[96];
    float* DNN_input = (float*)sds_alloc(sizeof(float*) * 96);
    float interpolate_out[72 * 14 * 2];
    //float golden_out[72 * 14 * 2];
    float* golden_out = (float*)sds_alloc(sizeof(float*) * 2016);
    //float recieved_pilots[24 * 2 * 2];
    float* recieved_pilots = (float*)sds_alloc(sizeof(float*) * 96);

    float golden_out_DNN[72 * 14 * 2];
    //float* golden_out = (float*)sds_alloc(sizeof(float*) * 2016);


    FILE* fina, * finb, * fin5;

    float inputsample1, inputsample3, inputsample4, inputsample5;
/*    fina = fopen("Biliner_input.dat", "r");
    if (!fina) {
        printf("could not open file \a\n");
        exit(101);
    }*/
    finb = fopen("recieved_pilots.dat", "r");
    if (!finb) {
        printf("could not open file \a\n");
        exit(101);
    }


    fin5 = fopen("test_golden_output.dat", "r");
    if (!fin5) {
        printf("could not open file \a\n");
        exit(101);
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    float* input_file_buffer = (float*)sds_alloc(sizeof(float*) * 96);
    float* w_dense1 = (float*)sds_alloc(sizeof(float*) * 96 * 48);
    float* b_dense1 = (float*)sds_alloc(sizeof(float*) * 48);
    float* w_dense2 = (float*)sds_alloc(sizeof(float*) * 48 * (72 * 14 * 2));
    float* b_dense2 = (float*)sds_alloc(sizeof(float*) * 72 * 14 * 2);
    float* output_buffer = (float*)sds_alloc(sizeof(float*) * 2016);


    FILE* fin, * fin2, * fin3, * fin4;
    float inputsample, inputsample_i, inputsample2, inputsample2_i;
    //int inputsample3, inputsample4,inputsample5;

    fin = fopen("w_dense1.dat", "r");
    fin2 = fopen("w_dense2.dat", "r");
    fin3 = fopen("b_dense1.dat", "r");
    fin4 = fopen("b_dense2.dat", "r");



    for (int i = 0; i < 96; i++)
    {
        //fscanf(fin, "%f", &inputsample);// fscanf(fin_i, "%f", &inputsample_i);
        input_file_buffer[i] = 1;		//A_i[i] = inputsample_i;
        //A2[i] = inputsample;		A_i2[i] = inputsample_i;

    }
    //w_DNN1
    for (int i = 0; i < 96 * 48; i++)
    {
        fscanf(fin, "%f", &inputsample2);// fscanf(fin_i, "%f", &inputsample_i);
        w_dense1[i] = inputsample2;		//A_i[i] = inputsample_i;
        //A2[i] = inputsample;		A_i2[i] = inputsample_i;

    }
    //b_DNN1
    for (int i = 0; i < 48; i++)
    {
        fscanf(fin3, "%f", &inputsample2);// fscanf(fin_i, "%f", &inputsample_i);
        b_dense1[i] = inputsample2;		//A_i[i] = inputsample_i;
        //A2[i] = inputsample;		A_i2[i] = inputsample_i;

    }
    //w_DNN2
    for (int i = 0; i < 48 * 2016; i++)
    {
        fscanf(fin2, "%f", &inputsample2);
        w_dense2[i] = inputsample2;// inputsample2;
    }
    //b_DNN2
    for (int i = 0; i < 2016; i++)
    {
        fscanf(fin4, "%f", &inputsample2);
        b_dense2[i] = inputsample2;// inputsample2;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Counter for calculation of execution time
    uint64_t cnt_start = 0;
    uint64_t cnt_stop = 0;
    //uint64_t frequency = sds_clock_frequency();


    //Separate counter for SW and HW functions
    uint64_t total_count_sw = 0;
    uint64_t total_count_hw = 0;

    float sum_sw,sum_hw;



    for (int snr = -5; snr <= 25; snr += 5) {

        float sum_LS = 0;
        float sum_LS_DNN = 0;
        for (int frame = 0; frame < num_of_frames_each_SNR; frame++)
        {

            //READING INPUT FROM .DAT file

                    //READING recieved_pilots INPUT FROM .DAT file

            for (int j = 0; j < 24; j++)
            {
                for (int k = 0; k < 2; k++) {
                    for (int i = 0; i < 2; i++) {
                        fscanf(finb, "%f", &inputsample3);
                        recieved_pilots[i + j * 2 + k * (2 * 24)] = inputsample3;
                        //weight_buffer[l * (3 * 2 * 8)] = inputsample2;
                        //printf("%f\t", inputsample1);
                    }
                    //cout << "\n";
                }
                //cout << "\n";
            //cout << "\n------\n";
            }



            //for (int j = 0; j < 24; j++)
            //{
            //    for (int k = 0; k < 2; k++) {
            //        for (int i = 0; i < 2; i++) {
            //            fscanf(fin1, "%f", &inputsample1);
            //            input[i + j * 2 + k * (2 * 24)] = inputsample1;
            //            //weight_buffer[l * (3 * 2 * 8)] = inputsample2;
            //            //printf("%d)\t%f\t", (i + j * 2 + k * (2 * 24)), inputsample1);
            //        }
            //        cout << "\n";
            //    }
            //    //cout << "\n";
            //cout << "\n------\n";
            //}

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


            //reading Golden out for DNN output:
            for (int i = 0,l=0; i < 14; i++)
            {
                for (int k = 0; k < 2; k++) {
                    for (int j = 0; j < 72; j++) {
                        //fscanf(fin5, "%f", &inputsample5);
                        golden_out_DNN[l] = golden_out[i + j * 14 + k * (14 * 72)];
                        l++;
                        //weight_buffer[l * (3 * 2 * 8)] = inputsample2;
                        //printf("%f\t", inputsample1);
                    }
                    //std::cout << "\n";
                }
                //std::cout << "\n";
            //std::cout << "\n------\n";
            }

   	cnt_start = sds_clock_counter(); //storing the time stamp of starting the sw function  call
    LS(recieved_pilots, input);
	cnt_stop = sds_clock_counter(); //storing the time stamp of stoping the sw function call
	total_count_sw = total_count_sw + (cnt_stop - cnt_start);


    //interpolate(input, interpolate_out);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////LS_DNN/////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Transform recieved pilots for DNN_input.!
    for (int i = 0, k = 0; k < 96; i = i + 2, k++)
    {
        DNN_input[k] = input[i];
        if (k == 23)
        {
            i = 46;
        }
        if (k == 47)
        {
            i = -1;
        }
        if (k == 71)
        {
            i = 47;
        }
        //cout << recieved_pilots[i]<<endl;
    }

    //    for (int i = 0; i < 96; i++)
    //    {
    //        std::cout << i<<"]"<<DNN_input[i] << endl;
    //}
    //float hid_out[48];
    float* hid_out = (float*)sds_alloc(sizeof(float*) * 48);

    //for (int i = 0; i < 96; i++)
    //{
    //    DNN_input[i] = 1;
    //}

   	cnt_start = sds_clock_counter(); //storing the time stamp of starting the sw function  call
   	perform_dense_L2(DNN_input,hid_out, w_dense1, b_dense1, 96, 48, 1);
   	perform_dense_L2(hid_out, output_buffer, w_dense2, b_dense2, 48, 2016, 0);
	cnt_stop = sds_clock_counter(); //storing the time stamp of stoping the sw function call
	total_count_sw = total_count_sw + (cnt_stop - cnt_start);


    //Display output of LS_DNN
    //for (int i = 0; i < 2016; i++)
    //{

    //    cout << i+1 << "]" << golden_out_DNN[i] << "  " << output_buffer[i] << endl;
    //}

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    //for (int j = 0; j < 24; j++)
    //{
    //    for (int k = 0; k < 2; k++) {
    //        for (int i = 0; i < 2; i++) {
    //            fscanf(fin1, "%f", &inputsample1);
    //            input2[i + j * 2 + k * (2 * 24)] = inputsample1;
    //            //weight_buffer[l * (3 * 2 * 8)] = inputsample2;
    //            printf("%d)\t%f\t%f\t", (i + j * 2 + k * (2 * 24)), input[i + j * 2 + k * (2 * 24)], input2[i + j * 2 + k * (2 * 24)]);
    //        }
    //        cout << "\n";
    //    }
    //    //cout << "\n";
    //    cout << "\n------\n";
    //}

    //cout << "OUTPUT: \n";
    //for (int i = 0; i < 72 * 2; i++) {
    //    if (i % 72 == 0)
    //    {
    //        printf("-------\n");
    //    }
    //    for (int j = 0; j < 14; j++)
    //    {
    //        cout << interpolate_out[i * 14 + j] << "\t";
    //        if (j == 13) {
    //            j == 0;
    //        }
    //    }
    //    cout << "\n";
    //}


    //CALCULATING MSE :
/*
    float squared_sum = 0;
    float diff = 0;
    float mse = 0;
    for (int i = 0; i < 72 * 14 * 2; i++)
    {
        diff = (golden_out[i] - interpolate_out[i]);
        diff = diff * diff;
        squared_sum += diff;
    }
    mse = squared_sum / (72 * 14 );

    sum_LS = sum_LS + mse;//Acculumating MSE for All frames
*/


    //CALCULATING MSE :

    float squared_sum = 0;
    float diff = 0;
    float mse = 0;
    for (int i = 0; i < 72 * 14 * 2; i++)
    {
        diff = (golden_out_DNN[i] - output_buffer[i]);
        diff = diff * diff;
        squared_sum += diff;
    }
    mse = squared_sum / (72 * 14);

    sum_LS_DNN = sum_LS_DNN + mse;//Acculumating MSE for All frames


    //std::cout << "\n\n" << "\tGOLDEN_MSE :\t\t" << 0.00026377;
    //std::cout << "\n" << "Interpolated_Resnet_MSE :\t" << mse << endl << "\n";



        }


        //printf("\n\n\ LS  MSE_over_SNR for SNR %d: %f", snr, sum_LS / num_of_frames_each_SNR);

        printf("\n\n\ LS-DNN  MSE_over_SNR for SNR %d: %f", snr, sum_LS_DNN / num_of_frames_each_SNR);


    }
    cout << endl;
       printf("Execution time of SW in seconds: %f\n", (double)total_count_sw / (double)sds_clock_frequency());
    return 0;

}
