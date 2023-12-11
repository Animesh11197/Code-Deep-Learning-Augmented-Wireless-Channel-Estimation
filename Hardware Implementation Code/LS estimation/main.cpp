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
#include "LS_hw.h"
#define num_of_frames_each_SNR 10

//#pragma warning(disable : 4996)
using namespace std;

void LS(float input[24*2*2], float output[24*2*2])
{
    for (int i = 0; i < 24*2; i++)
    {
        output[i] = 0.5 * (input[i] + input[i + (24 * 2)]);

        output[i + (24 * 2)]= 0.5 * ( input[i + (24 * 2)] - input[i]);

    }

}



void interpolate(float original_img[24*2*2], float out[72*14*8])
{
    //printf("\nTHIS FUNCTION WILL PERFORM BILINEAR INTERPOLATION.!\n");
    //input dimension is 24*2
    //ouput dimension is 72*14
	float output[72*8*14];
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

    for (int i=0;i<72*14*8;i++)
    	out[i]=output[i];

    //test[0] = 10;
    //return output;
}


int main()
{
    float* recieved_pilots=(float*)sds_alloc(sizeof(float*)* 24 * 2 * 2);
    float* input=(float*)sds_alloc(sizeof(float*)* 24 * 2 * 2);

    //float input[24 * 2 * 2];
    float input2[24 * 2 * 2];
    float interpolate_out[72 * 14 * 8];
    float golden_out[72 * 14 * 2];
    //float recieved_pilots[24 * 2 * 2];

    FILE* fin1, * fin2, * fin3;
    FILE* fin4_w, * fin4_b;
    FILE* fin5;
    float inputsample1, inputsample2, inputsample3, inputsample4, inputsample5;
    fin1 = fopen("Biliner_input.dat", "r");
    if (!fin1) {
        printf("could not open file \a\n");
        exit(101);
    }
    fin3 = fopen("recieved_pilots.dat", "r");
    if (!fin3) {
        printf("could not open file \a\n");
        exit(101);
    }


    fin5 = fopen("test_golden_output.dat", "r");
    if (!fin5) {
        printf("could not open file \a\n");
        exit(101);
    }


    //Counter for calculation of execution time
    uint64_t cnt_start = 0;
    uint64_t cnt_stop = 0;
    //uint64_t frequency = sds_clock_frequency();


    //Separate counter for SW and HW functions
    uint64_t total_count_sw = 0;
    uint64_t total_count_hw = 0;

    float sum_sw,sum_hw;


    for (int snr = -5; snr <= 25; snr += 5) {

        float sum = 0;
        for (int frame = 0; frame < num_of_frames_each_SNR; frame++)
        {

    //READING INPUT FROM .DAT file

            //READING recieved_pilots INPUT FROM .DAT file

            for (int j = 0; j < 24; j++)
            {
                for (int k = 0; k < 2; k++) {
                    for (int i = 0; i < 2; i++) {
                        fscanf(fin3, "%f", &inputsample3);
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

	cnt_start = sds_clock_counter(); //storing the time stamp of starting the sw function  call
		LS(recieved_pilots, input);
		interpolate(input, interpolate_out);
	cnt_stop = sds_clock_counter(); //storing the time stamp of stoping the sw function call
	total_count_sw = total_count_sw + (cnt_stop - cnt_start);


	cnt_start = sds_clock_counter(); //storing the time stamp of starting the sw function  call
		LS_hw(recieved_pilots, input);
		interpolate(input, interpolate_out);
	cnt_stop = sds_clock_counter(); //storing the time stamp of stoping the sw function call
	total_count_hw = total_count_hw + (cnt_stop - cnt_start);



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
    float squared_sum = 0;
    float diff = 0;
    float mse = 0;
    for (int i = 0; i < 72 * 14 * 2; i++)
    {
        diff = (golden_out[i] - interpolate_out[i]);
        diff = diff * diff;
        squared_sum += diff;
    }
    mse = squared_sum / (72 * 14);

    sum = sum + mse;//Acculumating MSE for All frames


    //std::cout << "\n\n" << "\tGOLDEN_MSE :\t\t" << 0.00026377;
    //std::cout << "\n" << "Interpolated_Resnet_MSE :\t" << mse << endl << "\n";



        }


        printf("\n\ LS  MSE_over_SNR for SNR %d: %f", snr, sum / num_of_frames_each_SNR);

    }
    cout << endl;
       printf("Execution time of SW in seconds: %f\n", (double)total_count_sw / (double)sds_clock_frequency());
       printf("Execution time of HW in seconds: %f\n", (double)total_count_hw / (double)sds_clock_frequency());
    return 0;

}
