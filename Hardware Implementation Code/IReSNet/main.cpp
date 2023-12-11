//#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstdio>
#include<stdio.h>
#include<stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include<cmath>
#include <inttypes.h>
#include <sds_lib.h>
#include "sds_lib.h"
#include "hw_functions.h"
//#include "hw_functions.cpp"
#define num_of_frames_each_SNR 10

//#pragma warning(disable : 4996)
using namespace std;
/*
M   INPUT channels
N   OUTPUT channels
I   input size assuming I*I
O   output size assuming O*O
K   Filter size K*K
*/

//void perform_conv(float* input, float* output, float* weight, float* bias, int M, int n, int i, int o, int k)
//stride set at 1
//void perform_conv(float* input, float* output, const float* weight, const float* bias, int M, int N, int In_h, int In_w, int K, int pad)
//{
//    int I_w = In_w + 2 * pad; //For including padding into the input buffer
//    int I_h = In_h + 2 * pad;
//    int O_w = I_w - K + 1; //Here, we have taken stride=1 if stride is NOT 1, we need to write it as O=((In-K+2*pad)/stride)+1
//    int O_h = I_h - K + 1;
//
//    float local_input_buffer[8][26][4];
//    float local_output_buffer[8][26][4];
//    float local_weight_buffer[8][8][3][3];
//
//    for (int p = 0; p < M; p++) //for all input channels
//    {
//        for (int iter = 0, iter2 = 0, i = 0, j = 0; iter < I_h * I_w; iter++, j++) //For complete image I*I
//        {
//            if (j == I_w)
//            {
//                j = 0;
//                i++;
//            }
//            if (i < pad || i >= In_h + pad || j < pad || j >= In_w + pad)
//            {
//                local_input_buffer[p][i][j] = 0;
//            }
//            else
//            {
//                local_input_buffer[p][i][j] = input[iter2 + p * In_h * In_w]; //store the image on chip.!!
//                iter2++;
//            }
//        }
//        cout << "INPUT stored\n";
//    }
//
//    /*/ PRINTS INPUT BUFFER
//    cout<<"INPUT BUFFER : \n";
//    for(int i=0;i<I_h;i++){
//        for(int j=0;j<I_w;j++)
//        {
//            cout<<local_input_buffer[0][i][j]<<"\t";
//        }
//        cout<<"\n";
//    }*/
//
//    for (int p = 0; p < N; p++) //for all output channels
//    {
//        for (int iter = 0, i = 0, j = 0; iter < O_h * O_w; iter++, j++) //For complete image I*I
//        {
//            if (j == O_w)
//            {
//                j = 0;
//                i++;
//            }
//            local_output_buffer[p][i][j] = 0.0; //store the image on chip.!!
//
//        }
//        cout << "OUTPUT stored\n";
//    }
//
//    for (int o = 0; o < N; o++)    //For all Output channels
//    {
//        for (int p = 0; p < M; p++)    //For all input channels
//        {
//            for (int iter = 0, i = 0, j = 0; iter < K * K; iter++, j++)  //For filter size K*K
//            {
//                if (j == K)
//                {
//                    j = 0;
//                    i++;
//                }
//                local_weight_buffer[o][p][i][j] = weight[iter + p * K * K + o * M * K * K];   //Read the weights and store in onchip memory
//
//            }
//        }
//        cout << "WEIGHTS stored\n";
//    }
//
//    for (int trr = 0; trr < O_h; trr++) //For Output Image Width
//    {
//        for (int tcc = 0; tcc < O_w; tcc++) //For O/p image Height
//        {
//            for (int too = 0; too < N; too++)  // For all output channels
//            {
//                for (int tii = 0; tii < M; tii++)  // For all input channels
//                {
//                    for (int i = 0; i < K; i++)  // For filter Width
//                    {
//                        for (int j = 0; j < K; j++)    // For filter height
//                        {
//                            //Acculumate the product of weight and input value for the output pixel [too][trr][tcc]
//                            //cout<<"in loop\n";
//                            local_output_buffer[too][trr][tcc] += (local_weight_buffer[too][tii][i][j] * local_input_buffer[tii][trr + i][tcc + j]);
//                            //cout<<"YESSSS\nIts working fine.!";
//                        }
//                    }
//                }
//            }
//            //cout<<"-->"<<local_output_buffer[0][trr][tcc];
//        }
//        //cout<<"\n";
//    }
//
//
//    for (int n = 0; n < N; n++)
//    {
//        for (int x = 0; x < O_h; x++)
//        {
//            for (int y = 0; y < O_w; y++)
//            {
//                float biased = local_output_buffer[n][x][y] + bias[n];
//                //local_output_buffer[n][x][y] = (biased > float(0)) ? biased : float(0);//Performing bias + ReLu action//removing ReLu activation
//                local_output_buffer[n][x][y] = biased;
//            }
//        }
//    }
//
//    /*cout<<"LOCAL OUTPUT BUFFER: dimensions are "<<O_w<<"*"<<O_h<<" \n";
//    for(int i=0;i<O_h;i++){
//        for(int j=0;j<O_w;j++)
//        {
//            cout<<local_output_buffer[0][i][j]<<"\t";
//        }
//        cout<<"\n";
//    }*/
//
//
//
//    for (int p = 0; p < N; p++)
//    {
//        for (int iter = 0, i = 0, j = 0; iter < O_w * O_h; iter++, j++)
//        {
//            if (j == O_w)
//            {
//                j = 0;
//                i++;
//            }
//            output[iter + p * O_w * O_h] = local_output_buffer[p][i][j];
//            //cout<<"output assigned at "<<iter+ p*O_h*O_w<<"\t"<<local_output_buffer[p][i][j]<<"\t i="<<i<<" j="<<j<<"\n";
//        }
//    }
//
//}
void perform_conv(float input[72*14*8], float output[72*14*8], const float weight[36*7*8*8], const float bias[1*1*8], int M, int N, int In_h, int In_w, int K_h, int K_w, int pad_left, int pad_right, int pad_above, int pad_below)
{
    int I_w = In_w + pad_right + pad_left; //For including padding into the input buffer
    int I_h = In_h + pad_above + pad_below;
    int O_w = I_w - K_w + 1; //Here, we have taken stride=1 if stride is NOT 1, we need to write it as O=((In-K+2*pad)/stride)+1
    int O_h = I_h - K_h + 1;
/*    cout << "I_w)" << I_w << endl;
    cout << "I_h)" << I_h << endl;
    cout << "K_w)" << K_w << endl;
    cout << "K_h)" << K_h << endl;
    cout << "O_w)" << O_w << endl;
    cout << "O_h)" << O_h << endl;*/
    float local_input_buffer[8][72 + 17 + 18][14 + 3 + 3];
    float local_output_buffer[8][72][14];
    float local_weight_buffer[8][8][36][7];

    //#pragma HLS ARRAY_PARTITION variable = local_output_buffer complete
    //#pragma HLS ARRAY_PARTITION variable=local_input_buffer dim=0 complete
    //#pragma HLS ARRAY_PARTITION variable=local_weight_buffer dim=0 complete
    //#pragma HLS ARRAY_PARTITION variable=bias dim=0 complete

    for (int p = 0; p < M; p++) //for all input channels
    {
        for (int iter = 0, iter2 = 0, i = 0, j = 0; iter < I_h * I_w; iter++, j++) //For complete image I*I
        {
            if (j == I_w)
            {
                j = 0;
                i++;
            }
            if (i < pad_above || i >= In_h + pad_above || j < pad_left || j >= In_w + pad_left)
            {
                local_input_buffer[p][i][j] = 0;
            }
            else
            {
                local_input_buffer[p][i][j] = input[iter2 + p * In_h * In_w]; //store the image on chip.!!
                iter2++;
            }
        }
        //cout << "INPUT stored\n";
    }

    // PRINTS INPUT BUFFER
/*    cout<<"INPUT BUFFER : \n";
    for(int i=0;i<I_h;i++){
        for(int j=0;j<I_w;j++)
        {
            cout<<local_input_buffer[0][i][j]<<"\t";
        }
        cout<<"\n";
    }*/

    // 0 out output buffer
    for (int p = 0; p < N; p++) //for all output channels
    {
        for (int iter = 0, i = 0, j = 0; iter < O_h * O_w; iter++, j++) //For complete image I*I
        {
            if (j == O_w)
            {
                j = 0;
                i++;
            }
            local_output_buffer[p][i][j] = 0.0; //store the image on chip.!!

        }
        //cout << "OUTPUT stored\n";
    }


    //load weights to weight buffer
    for (int o = 0; o < N; o++)    //For all Output channels
    {
        for (int p = 0; p < M; p++)    //For all input channels
        {
            for (int iter = 0, i = 0, j = 0; iter < K_w * K_h; iter++, j++)  //For filter size K*K
            {
                if (j == K_w)
                {
                    j = 0;
                    i++;
                }
                local_weight_buffer[o][p][i][j] = weight[iter + p * K_w * K_h + o * M * K_w * K_h];   //Read the weights and store in onchip memory

            }
        }
        //cout << "WEIGHTS stored\n";
    }

    for (int trr = 0; trr < O_h; trr++) //For Output Image Width
    {
        for (int tcc = 0; tcc < O_w; tcc++) //For O/p image Height
        {
            for (int too = 0; too < N; too++)  // For all output channels
            {
                for (int tii = 0; tii < M; tii++)  // For all input channels
                {
                    for (int i = 0; i < K_h; i++)  // For filter Width
                    {
                        for (int j = 0; j < K_w; j++)    // For filter height
                        {
                            //Acculumate the product of weight and input value for the output pixel [too][trr][tcc]
                            //cout<<"in loop\n";
                            local_output_buffer[too][trr][tcc] += (local_weight_buffer[too][tii][i][j] * local_input_buffer[tii][trr + i][tcc + j]);
                            //cout<<"YESSSS\nIts working fine.!";
                        }
                    }
                }
            }
            //cout<<"-->"<<local_output_buffer[0][trr][tcc];
        }
        //cout<<"\n";
    }


    for (int n = 0; n < N; n++)
    {
        for (int x = 0; x < O_h; x++)
        {
            for (int y = 0; y < O_w; y++)
            {
                float biased = local_output_buffer[n][x][y] + bias[n];
                //local_output_buffer[n][x][y] = (biased > float(0)) ? biased : float(0);//Performing bias + ReLu action//removing ReLu activation
                local_output_buffer[n][x][y] = biased;
            }
        }
    }

    //cout<<"LOCAL OUTPUT BUFFER: dimensions are "<<O_w<<"*"<<O_h<<" \n";
    //for(int i=0;i<O_h;i++){
    //    for(int j=0;j<O_w;j++)
    //    {
    //        cout<<local_output_buffer[0][i][j]<<"\t";
    //    }
    //    cout<<"\n";
    //}



    for (int p = 0; p < N; p++)
    {
        for (int iter = 0, i = 0, j = 0; iter < O_w * O_h; iter++, j++)
        {
            if (j == O_w)
            {
                j = 0;
                i++;
            }
            output[iter + p * O_w * O_h] = local_output_buffer[p][i][j];
            //cout<<"output assigned at "<<iter+ p*O_h*O_w<<"\t"<<local_output_buffer[p][i][j]<<"\t i="<<i<<" j="<<j<<"\n";
        }
    }

}


void ReLU(float* input, float* ReLU_out) {
    for (int i = 0; i < 24 * 2 * 8;i++)
    {
        ReLU_out[i]= (input[i] > float(0)) ? input[i] : float(0);
    }
}

void Add_layer(float input1[72*14*8], float input2[72*14*8], float add_out[72*14*8])
{
    for (int i = 0; i < 24 * 2 * 8; i++)
    {
#pragma HLS PIPELINE II=1
        add_out[i] = input1[i]+input2[i];
    }
}

void Add_layer5(float input1[72*14*8], float input2[72*14*8], float input3[72*14*8], float input4[72*14*8], float input5[72*14*8], float input6[72*14*8], float add_out[72*14*8])
{
    for (int i = 0; i < 24 * 2 * 8; i++)
    {
#pragma HLS PIPELINE II=1
        add_out[i] = input1[i] + input2[i] + input3[i] + input4[i] + input5[i] + input6[i];
    }
}

//void interpolate(float* original_img, float* output)
//{
//    printf("\nTHIS FUNCTION WILL PERFORM BILINEAR INTERPOLATION.!\n");
//    //input dimension is 24*2
//    //ouput dimension is 72*14
//    float old_h = 24;
//    float old_w = 2;
//    float new_h = 72;
//    float new_w = 14;
//    float h_scale_factor = old_h / new_h;//old_h/new_h
//    float w_scale_factor = old_w / new_w; //old_w / new_w
//    for (int i = 0; i < new_h - 1; i++)
//    {
//        for (int j = 0; j < new_w - 3; j++)
//        {
//            //printf("\n*********************************************\n");
//            //cout<<"i="<< i<< "\tj="<< j<<endl;
//            //map the coordinates back to the original image
//            float x = i * h_scale_factor;
//            float y = j * w_scale_factor;
//            //cout<<"x="<<x <<"\ty="<< y<<endl;
//            //calculate the cordinates of 4 surrounding pixels
//            int x_floor = floor(x);
//            int x_ceil = min(old_h - 1, ceil(x));
//            int y_floor = floor(y);
//            int y_ceil = min(old_w - 1, ceil(y));
//            float q;
//            //cout<<"x_floor= "<<x_floor<<"\tx_ceil="<<x_ceil<<"\ty_floor="<<y_floor<<"\ty_ceil="<<y_ceil ;
//
//            if ((x_ceil == x_floor) && (y_ceil == y_floor))
//            {
//                q = original_img[int(x) * 2 + int(y)];
//                //printf("\nq =%f\n", q);
//            }
//            else if (x_ceil == x_floor)
//            {
//                float q1 = original_img[int(x) * 2 + int(y_floor)];
//                float q2 = original_img[int(x) * 2 + int(y_ceil)];
//                q = q1 * (y_ceil - y) + q2 * (y - y_floor);
//                //cout<<"\nq1="<<q1<<"\tq2="<<q2<<"\tq="<<q;
//            }
//            else if (y_ceil == y_floor)
//            {
//                float q1 = original_img[int(x_floor) * 2 + int(y)];
//                float q2 = original_img[int(x_ceil) * 2 + int(y)];
//                q = (q1 * (x_ceil - x)) + (q2 * (x - x_floor));
//                //cout<<"\nq1="<<q1<<"\tq2="<<q2<<"\tq="<< q;
//            }
//            else
//            {
//                float v1 = original_img[x_floor * 2 + y_floor];
//                float v2 = original_img[x_ceil * 2 + y_floor];
//                float v3 = original_img[x_floor * 2 + y_ceil];
//                float v4 = original_img[x_ceil * 2 + y_ceil];
//                //cout<<"\nv1="<<v1<<"\tv2="<<v2<<"\tv3="<<v3<<"\tv4="<<v4;
//
//                float q1 = v1 * (x_ceil - x) + v2 * (x - x_floor);
//                float q2 = v3 * (x_ceil - x) + v4 * (x - x_floor);
//                q = q1 * (y_ceil - y) + q2 * (y - y_floor);
//                //cout<<"\nq1="<<q1<<"\tq2="<<q2<<"\tq="<<q;
//
//            }
//
//            output[((i + 1) * 14) + j + 3] = q;
//
//        }
//    }
//
//
//    for (int i = new_h - 1; i > 0; i--)
//    {
//        for (int j = 3; j > 0; j--)
//        {
//            output[(i * 14) + j - 1] = output[i * 14 + j];
//            //output[(i * 14) + j - 1] = 0;
//
//        }
//    }
//
//    for (int j = new_w - 1; j > -1; j--)
//    {
//        output[(0 * 14) + j] = output[1 * 14 + j];
//        //output[(0 * 14) + j ] = 0;
//    }
//    //test[0] = 10;
//    //return output;
//}

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
    cout << "Initializing Model" << endl;

    //float* input_buffer = (float*)malloc(sizeof(float*) * 27 * 2 * 2);
    //float* padded_input_buffer = (float*)malloc(sizeof(float*) * (612 + 10) * (14 + 10) * 64);//for Padding of max 5 for 612*14 input
    //float* output_buffer = (float*)malloc(sizeof(float*) * 612 * 14 * 64);

    float* input_file_buffer = (float*)sds_alloc(sizeof(float*) * 1);
    //float* weight_buffer = (float*)malloc(sizeof(float*) * 3 * 3 * 2 * 8);
    float* input=(float*)sds_alloc(sizeof(float*) * 72 * 14 * 8);//[24 * 2 * 2];// = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
    float* out_conv_1=(float*)sds_alloc(sizeof(float*)* 72 * 14 * 8);
    float* out_conv_5=(float*)sds_alloc(sizeof(float*)* 72 * 14 * 8);//[24 * 2 * 8];
    float* out_conv_7=(float*)sds_alloc(sizeof(float*)* 72 * 14 * 8);//[24 * 2 * 8];
    float* out_conv_9=(float*)sds_alloc(sizeof(float*)* 72 * 14 * 8);//[24 * 2 * 8];
    float* add_out_1=(float*)sds_alloc(sizeof(float*)* 72 * 14 * 8);//[24 * 2 * 8];
    float* add_out_2=(float*)sds_alloc(sizeof(float*)* 72 * 14 * 8);//[24 * 2 * 8];
    float* add_out_3=(float*)sds_alloc(sizeof(float*)* 72 * 14 * 8);//[24 * 2 * 8];
    float* add_out_4=(float*)sds_alloc(sizeof(float*)* 72 * 14 * 8);//[24 * 2 * 8];
    float* add_out_5=(float*)sds_alloc(sizeof(float*)* 72 * 14 * 8);//[24 * 2 * 8];
    float* interpolate_out=(float*)sds_alloc(sizeof(float*)* 72 * 14 * 8);

    float golden_out[72 * 14 * 2];

    float* output=(float*)sds_alloc(sizeof(float*)* 72 * 14 * 8);//[24*2*8];
    float* relu_out=(float*)sds_alloc(sizeof(float*)* 72 * 14 * 8);//[24 * 2 * 8];
    float* weight_buffer = (float*)sds_alloc(sizeof(float*) * 36 * 7 * 8 * 8);//3 * 3 * 2 * 8
    float* bias_buffer = (float*)sds_alloc(sizeof(float*) * 8);
    float* weight_2 = (float*)sds_alloc(sizeof(float*) * 36 * 7 * 8 * 8);
    float* weight_3 = (float*)sds_alloc(sizeof(float*) * 36 * 7 * 8 * 8);
    float* weight_4 = (float*)sds_alloc(sizeof(float*) * 36 * 7 * 8 * 8);
    float* weight_5 = (float*)sds_alloc(sizeof(float*) * 36 * 7 * 8 * 8);
    float* weight_6 = (float*)sds_alloc(sizeof(float*) * 36 * 7 * 8 * 8);
    float* weight_7 = (float*)sds_alloc(sizeof(float*) * 36 * 7 * 8 * 8);
    float* weight_8 = (float*)sds_alloc(sizeof(float*) * 36 * 7 * 8 * 8);
    float* weight_9 = (float*)sds_alloc(sizeof(float*) * 36 * 7 * 8 * 8);
    float* weight_10 = (float*)sds_alloc(sizeof(float*) * 36 * 7 * 8 * 8);
    float* bias_2 = (float*)sds_alloc(sizeof(float*) * 8);
    float* bias_3 = (float*)sds_alloc(sizeof(float*) * 8);
    float* bias_4 = (float*)sds_alloc(sizeof(float*) * 8);
    float* bias_5 = (float*)sds_alloc(sizeof(float*) * 8);
    float* bias_6 = (float*)sds_alloc(sizeof(float*) * 8);
    float* bias_7 = (float*)sds_alloc(sizeof(float*) * 8);
    float* bias_8 = (float*)sds_alloc(sizeof(float*) * 8);
    float* bias_9 = (float*)sds_alloc(sizeof(float*) * 8);
    float* bias_10 = (float*)sds_alloc(sizeof(float*) * 8);
    //////layer 22:
    float* Weight_22 = (float*)sds_alloc(sizeof(float*) * 36 * 7 * 8 * 8);
    float* Bias_22 = (float*)sds_alloc(sizeof(float*) *  8);
    float* Out_22 = (float*)sds_alloc(sizeof(float*) * 72*14*8);
    //float* bias_ = (float*)malloc(sizeof(float*) * 8);
    //float weight[9] = { 1,1,1,1,1,1,1,1,1 };
    //float bias[9] = { 0,1,0,0,0,0,0,0,0 };

    FILE * fin1,* fin2,* fin3;
    FILE* fin4_w, * fin4_b;
    FILE* fin5;
    float inputsample1,inputsample2,inputsample3,inputsample4, inputsample5;
    fin1 = fopen("test_input.dat", "r");
    if (!fin1) {
        printf("could not open file \a\n");
        exit(101);
    }
    fin2 = fopen("test_conv_1.dat", "r");
    if (!fin2) {
        printf("could not open file \a\n");
        exit(101);
    }
    fin3 = fopen("test_bias_1.dat", "r");
    if (!fin3) {
        printf("could not open file \a\n");
        exit(101);
    }
    fin4_w = fopen("test_layer22_Weight.dat", "r");
    if (!fin2) {
        printf("could not open file \a\n");
        exit(101);
    }
    fin4_b = fopen("test_layer22_Bias.dat", "r");
    if (!fin3) {
        printf("could not open file \a\n");
        exit(101);
    }
    fin5 = fopen("test_golden_output.dat", "r");
    if (!fin5) {
        printf("could not open file \a\n");
        exit(101);
    }
    //for (int l = 0; l < 8; l++)
    //{   for (int k = 0; k < 2; k++) {
    //        for (int j = 0; j < 3; j++) {
    //            for (int i = 0; i < 3; i++) {
    //                fscanf(fin2, "%f", &inputsample2);
    //                weight_buffer[i * 3 * 3 * 2 + j];
    //                printf("%f\t", inputsample2);
    //            }
    //            cout << "\n";
    //        }
    //        cout << "\n";
    //    }
    //cout << "\n------\n";
    //}



    //int out_channels = 2;
        //READING WEIGHTS FROM .DAT file
    for (int j = 0; j < 3; j++) {
        for (int l = 0; l < 8; l++) {
            for (int k = 0; k < 2; k++) {
                for (int i = 0; i < 3; i++) {
                    fscanf(fin2, "%f", &inputsample2);
                    weight_buffer[i + j * 3 + k * (3 * 3) + l * (3 * 3 * 2)] = inputsample2;
                    //weight_buffer[l * (3 * 2 * 8)] = inputsample2;
                    //printf("%f\t", inputsample2);
                }
                //cout << "\n";
            }
            //cout << "\n";
        }
        //cout << "\n------\n";
    }


    for (int layer = 0; layer < 10; layer++) {
        for (int j = 0; j < 3; j++){
            for (int l = 0; l < 8; l++) {
                for (int k = 0; k < 8; k++) {
                    for (int i = 0; i < 3; i++) {
                        fscanf(fin2, "%f", &inputsample2);
                        switch (layer) {
                        case 0:
                            weight_2[i + j * 3 + k * (3 * 3) + l * (3 * 3 * 8)] = inputsample2;
                            break;
                        case 1:
                            weight_3[i + j * 3 + k * (3 * 3) + l * (3 * 3 * 8)] = inputsample2;
                            break;
                        case 2:
                            weight_4[i + j * 3 + k * (3 * 3) + l * (3 * 3 * 8)] = inputsample2;
                            break;
                        case 3:
                            weight_5[i + j * 3 + k * (3 * 3) + l * (3 * 3 * 8)] = inputsample2;
                            break;
                        case 4:
                            weight_6[i + j * 3 + k * (3 * 3) + l * (3 * 3 * 8)] = inputsample2;
                            break;
                        case 5:
                            weight_7[i + j * 3 + k * (3 * 3) + l * (3 * 3 * 8)] = inputsample2;
                            break;
                        case 6:
                            weight_8[i + j * 3 + k * (3 * 3) + l * (3 * 3 * 8)] = inputsample2;
                            break;
                        case 7:
                            weight_9[i + j * 3 + k * (3 * 3) + l * (3 * 3 * 8)] = inputsample2;
                            break;
                        case 8:
                            weight_10[i + j * 3 + k * (3 * 3) + l * (3 * 3 * 8)] = inputsample2;
                            break;

                        }

                        //weight_buffer[l * (3 * 2 * 8)] = inputsample2;
                        //printf("%f\t", inputsample2);
                    }
                    //cout << "\n";
                }
                //cout << "\n";
            }
            //cout << "\n------\n";
        }

     }


    //PRINTS WEIGHTS:
    /*cout << "WEIGHTS:\n";
    for (int l = 0; l < 8; l++)
    {   for (int k = 0; k < 8; k++) {
            for (int j = 0; j < 3; j++) {
                for (int i = 0; i < 3; i++) {
                    //fscanf(fin2, "%f", &inputsample2);
                    cout<<weight_3[i + j * 3 + k * (3 * 3) + l * (3 * 3 * 8)]<<"\t";
                    //printf("%f\t", inputsample2);
                }
                cout << "\n";
            }
            cout << "\n";
        }
    cout << "\n------\n";
    }*/
    ///Layer 22 Weights Reading:
    for (int j = 0; j < 36; j++) {
        for (int l = 0; l < 2; l++) {
            for (int k = 0; k < 8; k++) {
                for (int i = 0; i < 7; i++) {
                    fscanf(fin4_w, "%f", &inputsample4);
                    Weight_22[i + j * 7 + k * (36 * 7) + l * (36 * 7 * 8)] = inputsample4;
                    //weight_buffer[l * (3 * 2 * 8)] = inputsample2;
                    //printf("%f\t", inputsample2);
                }
                //cout << "\n";
            }
            //cout << "\n";
        }
        //cout << "\n------\n";
    }


        //READING BIAS:
    //for (int layer = 0; layer < 2; layer++) {
        for (int i = 0; i < 8; i++) {
            fscanf(fin3, "%f", &inputsample3);
            //if(layer==0)
                bias_buffer[i] = inputsample3;
            //if(layer==1)
                //bias_2[i] = inputsample3;
        }
     //}

       // cout << "\nBias:\n";
    for (int layer = 0; layer < 10; layer++) {
        for (int i = 0; i < 8; i++) {
            fscanf(fin3, "%f", &inputsample3);
            switch (layer) {
                case 0:
                    bias_2[i] = inputsample3;
                    break;
                case 1:
                    bias_3[i] = inputsample3;
                    break;
                case 2:
                    bias_4[i] = inputsample3;
                    break;
                case 3:
                    bias_5[i] = inputsample3;
                    break;
                case 4:
                    bias_6[i] = inputsample3;
                    break;
                case 5:
                    bias_7[i] = inputsample3;
                    break;
                case 6:
                    bias_8[i] = inputsample3;
                    break;
                case 7:
                    bias_9[i] = inputsample3;
                    break;
                case 8:
                    bias_10[i] = inputsample3;
                    break;

            }

        }
    }
    /*for (int i = 0; i < 8; i++)
    {
        std::cout << bias_3[i] << "\n";
    }*/

    ///layer 22 Bias read:
    for (int i = 0; i < 2; i++) {
        fscanf(fin4_b, "%f", &inputsample4);
        //if(layer==0)
        Bias_22[i] = inputsample4;
        //if(layer==1)
            //bias_2[i] = inputsample3;
    }



    //Counter for calculation of execution time
    uint64_t cnt_start = 0;
    uint64_t cnt_stop = 0;
    //uint64_t frequency = sds_clock_frequency();


    //Separate counter for SW and HW functions
    uint64_t total_count_sw = 0;
    uint64_t total_count_hw = 0;

    float sum_sw,sum_hw;

for (int snr=-5;snr<=25;snr+=5){

sum_sw = 0;
sum_hw = 0;
	for (int frame = 0; frame < num_of_frames_each_SNR; frame++)
	{


		//READING INPUT FROM .DAT file

		for (int j = 0; j < 24; j++)
		{
			for (int k = 0; k < 2; k++) {
				for (int i = 0; i < 2; i++) {
					fscanf(fin1, "%f", &inputsample1);
					input[i + j * 2 + k * (2 * 24)] = inputsample1;
					//weight_buffer[l * (3 * 2 * 8)] = inputsample2;
					//printf("%f\t", inputsample1);
				}
				//cout << "\n";
			}
			//cout << "\n";
		//cout << "\n------\n";
		}


			////PRINTING INPUTS: //ONLY FOR DEBUGGING//
			//cout << "INPUTS: \n";
			//for (int k = 0; k < 2; k++) {
			//    for (int j = 0; j < 24; j++)
			//    {
			//        for (int i = 0; i < 2; i++)
			//            cout<<input[i+j*2+k*(2*24)]<<"\t";
			//        cout << "\n";
			//    }
			//    cout << "\n----\n";
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
			perform_conv(input, out_conv_1, weight_buffer, bias_buffer, 2, 8, 24, 2, 3,3,1,1,1,1);
				perform_conv(out_conv_1, output, weight_2, bias_2, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
				ReLU(output, relu_out);
				perform_conv(relu_out, output, weight_3, bias_3, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
					Add_layer(out_conv_1, output,add_out_1);
				perform_conv(add_out_1, output, weight_4, bias_4, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
				ReLU(output, relu_out);
				perform_conv(relu_out, out_conv_5, weight_5, bias_5, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
					Add_layer(out_conv_5, add_out_1, add_out_2);
				perform_conv(add_out_2, output, weight_6, bias_6, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
				ReLU(output, relu_out);
				perform_conv(relu_out, out_conv_7, weight_7, bias_7, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
					Add_layer(out_conv_7, add_out_2,add_out_3);
				perform_conv(add_out_3, output, weight_8, bias_8, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
				ReLU(output, relu_out);
				perform_conv(relu_out, out_conv_9, weight_9, bias_9, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
					Add_layer(out_conv_9, add_out_3,add_out_4);
				perform_conv(add_out_4, output, weight_10, bias_10, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
			Add_layer5(out_conv_1, add_out_1, add_out_2, add_out_3, add_out_4, output,add_out_5);
			interpolate(add_out_5,interpolate_out);
			perform_conv(interpolate_out, Out_22, Weight_22, Bias_22, 8, 2, 72, 14, 36, 7, 3, 3, 17, 18);
			cnt_stop = sds_clock_counter(); //storing the time stamp of stoping the sw function call
			total_count_sw = total_count_sw + (cnt_stop - cnt_start);


			//CALCULATING MSE :
			   float squared_sum = 0;
			   float diff = 0;
			   float mse_sw = 0;
			   //function for calculating MSE:
			   for (int i = 0; i < 72 * 14 * 2; i++)
			   {
				   diff = (golden_out[i] - Out_22[i]);
				   diff = diff * diff;
				   squared_sum += diff;
			   }
			   mse_sw = squared_sum / (72 * 14);
			   sum_sw += mse_sw;
		//////////////////////////////////////////////////////////////////////////////////////////////////////////



			/////////////////////////////////////////////////////////////////////////////////////////////////////
		cnt_start = sds_clock_counter(); //storing the time stamp of starting the sw function  call
			perform_conv_hw(input, out_conv_1, weight_buffer, bias_buffer, 2, 8, 24, 2, 3,3,1,1,1,1);
				perform_conv_hw(out_conv_1, output, weight_2, bias_2, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
				ReLU(output, relu_out);
				perform_conv_hw(relu_out, output, weight_3, bias_3, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
					Add_layer(out_conv_1, output,add_out_1);
				perform_conv_hw(add_out_1, output, weight_4, bias_4, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
				ReLU(output, relu_out);
				perform_conv_hw(relu_out, out_conv_5, weight_5, bias_5, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
					Add_layer(out_conv_5, add_out_1, add_out_2);
				perform_conv_hw(add_out_2, output, weight_6, bias_6, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
				ReLU(output, relu_out);
				perform_conv_hw(relu_out, out_conv_7, weight_7, bias_7, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
					Add_layer(out_conv_7, add_out_2,add_out_3);
				perform_conv_hw(add_out_3, output, weight_8, bias_8, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
				ReLU(output, relu_out);
				perform_conv_hw(relu_out, out_conv_9, weight_9, bias_9, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
					Add_layer(out_conv_9, add_out_3,add_out_4);
				perform_conv_hw(add_out_4, output, weight_10, bias_10, 8, 8, 24, 2, 3, 3, 1, 1, 1, 1);
			Add_layer5(out_conv_1, add_out_1, add_out_2, add_out_3, add_out_4, output,add_out_5);
			interpolate(add_out_5,interpolate_out);
			perform_conv_hw(interpolate_out, Out_22, Weight_22, Bias_22, 8, 2, 72, 14, 36, 7, 3, 3, 17, 18);
			cnt_stop = sds_clock_counter(); //storing the time stamp of stoping the sw function call
			total_count_hw = total_count_hw + (cnt_stop - cnt_start);

			//FOR DEBUGGING FILE INPUTS:
			//for (int l = 0; l < 8; l++)
			//{
			//    for (int k = 0; k < 2; k++) {
			//        for (int j = 0; j < 3; j++) {
			//            for (int i = 0; i < 3; i++) {
			//                //fscanf(fin2, "%f", &inputsample2);
			//                cout<<weight_buffer[i + (j * 3) + (k * (3 * 3) )+ (l * (3 * 3 *2))]<<"\t";
			//                //printf("%f\t", );
			//            }
			//            cout << "\n";
			//        }
			//        cout << "\n";
			//    }
			//    cout << "\n------\n";
			//}



			//cout << "OUTPUT: \n";
			//for (int i = 0; i < 24*8; i++) {
			//    if (i % 24 == 0)
			//    {
			//        printf("-------\n");
			//    }
			//    for (int j = 0; j < 2; j++)
			//    {
			//        cout << add_out_5[i * 2 + j] << "\t";
			//        if (j == 1) {
			//            j == 0;
			//        }
			//    }
			//    cout << "\n";
			//}

			//cout << "OUTPUT: \n";
			//for (int i = 0; i < 72 * 8; i++) {
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

		//// PRINTING OUTPUT
		//    std::cout << "OUTPUT: \n\n";
		//    for (int channel = 0; channel < 2; channel++)
		//    {
		//        for (int i = 0; i < 72; i++)
		//        {
		//            for (int j = 0; j < 14; j++)
		//            {
		//                //fscanf(fin2, "%f", &inputsample2);
		//                //in[i * N + j] = inputsample2;
		//                //input[i][j] = inputsample2;
		//                //printf("%f ", out[i*14+j]);
		//                std::cout << (channel * 72 * 14) + (i * 14) + j + 1 << ")" << golden_out[(channel * 72 * 14) + (i * 14 + j)] << " ";
		//                std::cout << Out_22[(channel * 72 * 14) + (i * 14 + j)] << "\t";
		//                //fprintf(fin2_out, "%f\n", out[(channel * 72 * 14) + (i * 14 + j)]);
		//            }
		//            printf("\n\n");
		//        }
		//        std::cout << "-----" << channel << "----\n";
		//    }


			//CALCULATING MSE :
			squared_sum = 0;
			diff = 0;
			float mse_hw = 0;
			//function for calculating MSE:
			for (int i = 0; i < 72 * 14 * 2; i++)
			{
				diff = (golden_out[i] - Out_22[i]);
				diff = diff * diff;
				squared_sum += diff;
			}
			mse_hw = squared_sum / (72 * 14);
			sum_hw+=mse_hw;

	}

	printf("\n\nInterpolated_ResNet MSE_over_SNR for SNR %d on SW: %f", snr, sum_sw / num_of_frames_each_SNR);
	printf("\n\nInterpolated_ResNet MSE_over_SNR for SNR %d on HW: %f", snr, sum_hw / num_of_frames_each_SNR);

}
    //std::cout << "\n\n" << "\tGOLDEN_MSE :\t\t" << 0.00026377;
    //std::cout << "\n" << "Interpolated_Resnet_MSE_SW :\t" << mse_sw<<endl<<"\n";
    //std::cout << "\n" << "Interpolated_Resnet_MSE_HW :\t" << mse_hw<<endl<<"\n";

    printf("Execution time of SW in seconds: %f\n", (double)total_count_sw / (double)sds_clock_frequency());
    printf("Execution time of HW in seconds: %f\n", (double)total_count_hw / (double)sds_clock_frequency());

    return 0;

}


/*
// PRINTS INPUT
cout<<"INPUT: \n";
for(int i=0;i<5;i++){
    for(int j=0;j<5;j++)
    {
        cout<<input[i*5+j]<<"\t";
    }
    cout<<"\n";
}
*/

//PRINTS WEIGHTS
//cout << "WEIGHTS: \n";
//for (int i = 0; i < 3; i++) {
//    for (int j = 0; j < 3; j++)
//    {
//        cout << weight[i * 3 + j] << "\t";
//    }
//    cout << "\n";
//}
//
//perform_conv(input, output, weight, bias, 1, 1, 4, 5, 3, 2);
//
//////prints output matrix
//cout << "OUTPUT: \n";
//for (int i = 0; i < 6; i++) {
//    for (int j = 0; j < 7; j++)
//    {
//        cout << output[i * 7 + j] << "\t";
//    }
//    cout << "\n";
//}
//
