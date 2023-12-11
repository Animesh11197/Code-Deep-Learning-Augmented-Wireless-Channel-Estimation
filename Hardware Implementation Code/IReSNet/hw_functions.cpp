#include"hw_functions.h"
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
#include "hls_half.h"
#include <math.h>
#include "ap_int.h"
#include <ap_fixed.h>
typedef ap_fixed<16,3,AP_RND> data_t;
//typedef float data_t;
//typedef half data_t;

/*
void perform_conv_hw(float input[72*14*8], float output[72*14*8], const float weight[36*7*8*8], const float bias[1*1*8], int M, int N, int In_h, int In_w, int K_h, int K_w, int pad_left, int pad_right, int pad_above, int pad_below)
{
    int I_w = In_w + pad_right + pad_left; //For including padding into the input buffer
    int I_h = In_h + pad_above + pad_below;
    int O_w = I_w - K_w + 1; //Here, we have taken stride=1 if stride is NOT 1, we need to write it as O=((In-K+2*pad)/stride)+1
    int O_h = I_h - K_h + 1;
    std::cout << "I_w)" << I_w << std::endl;
    std::cout << "I_h)" << I_h << std::endl;
    std::cout << "K_w)" << K_w << std::endl;
    std::cout << "K_h)" << K_h << std::endl;
    std::cout << "O_w)" << O_w << std::endl;
    std::cout << "O_h)" << O_h << std::endl;

    float local_input_buffer[8][72 + 17 + 18][14 + 3 + 3];
    float local_output_buffer[8][72][14];
    float local_weight_buffer[8][8][36][7];



//#pragma HLS ARRAY_PARTITION variable=local_output_buffer complete
//#pragma HLS ARRAY_PARTITION variable=local_input_buffer complete
//#pragma HLS ARRAY_PARTITION variable=local_weight_buffer complete
//#pragma HLS ARRAY_PARTITION variable=bias complete


    LOAD_INPUT:
	for (int p = 0; p < M; p++) //for all input channels
    {
//#pragma HLS loop_tripcount min=1 max=8
//#pragma HLS pipeline
        for (int iter = 0, iter2 = 0, i = 0, j = 0; iter < I_h * I_w; iter++, j++) //For complete image I*I
        {
//#pragma HLS loop_tripcount min=1 max=72*14
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
        //std::cout << "INPUT stored\n";
    }

    // PRINTS INPUT BUFFER
    std::cout<<"INPUT BUFFER : \n";
    for(int i=0;i<I_h;i++){
        for(int j=0;j<I_w;j++)
        {
            std::cout<<local_input_buffer[0][i][j]<<"\t";
        }
        std::cout<<"\n";
    }

    // 0 out output buffer
    LOAD_OUTPUT:
	for (int p = 0; p < N; p++) //for all output channels
    {
//#pragma HLS loop_tripcount min=1 max=8
//#pragma HLS pipeline
        for (int iter = 0, i = 0, j = 0; iter < O_h * O_w; iter++, j++) //For complete image I*I
        {
//#pragma HLS loop_tripcount min=1 max=72*14
            if (j == O_w)
            {
                j = 0;
                i++;
            }
            local_output_buffer[p][i][j] = 0.0; //store the image on chip.!!

        }
        //std::cout << "OUTPUT stored\n";
    }


    //load weights to weight buffer
    LOAD_WEIGHTS:
	for (int o = 0; o < N; o++)    //For all Output channels
    {
//#pragma HLS loop_tripcount min=1 max=8
//#pragma HLS pipeline
        for (int p = 0; p < M; p++)    //For all input channels
        {
//#pragma HLS loop_tripcount min=1 max=8
//#pragma HLS pipeline
            for (int iter = 0, i = 0, j = 0; iter < K_w * K_h; iter++, j++)  //For filter size K*K
            {
//#pragma HLS loop_tripcount min=1 max=36*7
                if (j == K_w)
                {
                    j = 0;
                    i++;
                }
                local_weight_buffer[o][p][i][j] = weight[iter + p * K_w * K_h + o * M * K_w * K_h];   //Read the weights and store in onchip memory

            }
        }
        //std::cout << "WEIGHTS stored\n";
    }

    CONVOLVE_LOOP:
	for (int trr = 0; trr < O_h; trr++) //For Output Image Width
    {
//#pragma HLS loop_tripcount min=1 max=14
//#pragma HLS UNROLL
        for (int tcc = 0; tcc < O_w; tcc++) //For O/p image Height
        {
//#pragma HLS loop_tripcount min=1 max=72
            for (int too = 0; too < N; too++)  // For all output channels
            {
//#pragma HLS loop_tripcount min=1 max=8
//#pragma HLS UNROLL
            	for (int tii = 0; tii < M; tii++)  // For all input channels
                {
//#pragma HLS loop_tripcount min=1 max=8
                    for (int i = 0; i < K_h; i++)  // For filter Width
                    {
//#pragma HLS loop_tripcount min=1 max=14
//#pragma HLS UNROLL
                        for (int j = 0; j < K_w; j++)    // For filter height
                        {
//#pragma HLS loop_tripcount min=1 max=72
//#pragma HLS PIPELINE
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


	 BIAS_LOOP:for (int n = 0; n < N; n++)
    {
#pragma HLS UNROLL
#pragma HLS loop_tripcount min=1 max=8
        for (int x = 0; x < O_h; x++)
        {
#pragma HLS loop_tripcount min=1 max=72
            for (int y = 0; y < O_w; y++)
            {
#pragma HLS loop_tripcount min=1 max=14
#pragma HLS PIPELINE II=1
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


OUT_LOOP:
	for (int p = 0; p < N; p++)
    {
#pragma HLS loop_tripcount min=1 max=8
//#pragma HLS Pipeline
        for (int iter = 0, i = 0, j = 0; iter < O_w * O_h; iter++, j++)
        {
#pragma HLS loop_tripcount min=1 max=72*14
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
*/


void perform_conv_hw(float input[72*14*8], float output[72*14*8], const float weight[36*7*8*8],  float bias[1*1*8], int M, int N, int In_h, int In_w, int K_h, int K_w, int pad_left, int pad_right, int pad_above, int pad_below)
{
    int I_w = In_w + pad_right + pad_left; //For including padding into the input buffer
    int I_h = In_h + pad_above + pad_below;
    int O_w = I_w - K_w + 1; //Here, we have taken stride=1 if stride is NOT 1, we need to write it as O=((In-K+2*pad)/stride)+1
    int O_h = I_h - K_h + 1;


    //I_w = In_w + pad_right + pad_left;
    //I_h = In_h + pad_above + pad_below;
    //O_w=In_w + pad_right + pad_left - K_w + 1;
    //O_h=In_h + pad_above + pad_below - K_h + 1;
    /*std::cout << "I_w)" << I_w << endl;
    std::cout << "I_h)" << I_h << endl;
    std::cout << "K_w)" << K_w << endl;
    std::cout << "K_h)" << K_h << endl;
    std::cout << "O_w)" << O_w << endl;
    std::cout << "O_h)" << O_h << endl;*/
    data_t  local_input_buffer[8][72 + 17 + 18][14 + 3 + 3];
    data_t local_output_buffer[8][72][14];
    data_t local_weight_buffer[8][8][36][7];
    data_t local_bias[8];

//#pragma HLS ARRAY_PARTITION variable=local_output_buffer complete dim=0
//#pragma HLS ARRAY_PARTITION variable=local_input_buffer complete dim=0
//#pragma HLS ARRAY_PARTITION variable=local_weight_buffer complete dim=0
//#pragma HLS ARRAY_PARTITION variable=bias complete dim=0


   LOAD_BIAS:for(int i=0;i<8;i++)
    {
#pragma HLS PIPELINE
    	local_bias[i]=bias[i];
    }

    LOAD_INPUT:for (int p = 0; p < M; p++) //for all input channels
    {
#pragma HLS loop_tripcount min=1 max=8
#pragma HLS pipeline
        for (int iter = 0, iter2 = 0, i = 0, j = 0; iter < (In_h + pad_above + pad_below) * (In_w + pad_right + pad_left); iter++, j++) //For complete image I*I
        {
#pragma HLS loop_tripcount min=1 max=72*14
            if (j == (In_w + pad_right + pad_left))
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
        //std::cout << "INPUT stored\n";
    }

    /*/ PRINTS INPUT BUFFER
    std::cout<<"INPUT BUFFER : \n";
    for(int i=0;i<I_h;i++){
        for(int j=0;j<I_w;j++)
        {
            std::cout<<local_input_buffer[0][i][j]<<"\t";
        }
        std::cout<<"\n";
    }*/

    // 0 out output buffer
    LOAD_OUTPUT:for (int p = 0; p < N; p++) //for all output channels
    {
#pragma HLS loop_tripcount min=1 max=8
#pragma HLS pipeline
        for (int iter = 0, i = 0, j = 0; iter < (In_h + pad_above + pad_below - K_h + 1) * (In_w + pad_right + pad_left - K_w + 1); iter++, j++) //For complete image I*I
        {
#pragma HLS loop_tripcount min=1 max=72*14
            if (j == (In_w + pad_right + pad_left - K_w + 1))
            {
                j = 0;
                i++;
            }
            local_output_buffer[p][i][j] = 0.0; //store the image on chip.!!

        }
        //std::cout << "OUTPUT stored\n";
    }


    //load weights to weight buffer
    LOAD_WEIGHTS:for (int o = 0; o < N; o++)    //For all Output channels
    {
#pragma HLS loop_tripcount min=1 max=8
#pragma HLS pipeline
        for (int p = 0; p < M; p++)    //For all input channels
        {
#pragma HLS loop_tripcount min=1 max=8
#pragma HLS pipeline
            for (int iter = 0, i = 0, j = 0; iter < K_w * K_h; iter++, j++)  //For filter size K*K
            {
#pragma HLS loop_tripcount min=1 max=36*7
                if (j == K_w)
                {
                    j = 0;
                    i++;
                }
                local_weight_buffer[o][p][i][j] = weight[iter + p * K_w * K_h + o * M * K_w * K_h];   //Read the weights and store in onchip memory

            }
        }
        //std::cout << "WEIGHTS stored\n";
    }

    CONV_LOOP:for (int trr = 0; trr < O_h; trr++) //For Output Image Width
    {
#pragma HLS loop_tripcount min=1 max=14
//#pragma HLS PIPELINE
        for (int tcc = 0; tcc < (In_w + pad_right + pad_left - K_w + 1); tcc++) //For O/p image Height
        {
#pragma HLS loop_tripcount min=1 max=72
            for (int too = 0; too < N; too++)  // For all output channels
            {
#pragma HLS loop_tripcount min=1 max=8
                for (int tii = 0; tii < M; tii++)  // For all input channels
                {
#pragma HLS loop_tripcount min=1 max=8
                    for (int i = 0; i < K_h; i++)  // For filter Width
                    {
#pragma HLS loop_tripcount min=1 max=7
                        for (int j = 0; j < K_w; j++)    // For filter height
                        {
#pragma HLS loop_tripcount min=1 max=36

                            //Acculumate the product of weight and input value for the output pixel [too][trr][tcc]
                            //std::cout<<"in loop\n";
                            local_output_buffer[too][trr][tcc] += (local_weight_buffer[too][tii][i][j] * local_input_buffer[tii][trr + i][tcc + j]);
                            //std::cout<<"YESSSS\nIts working fine.!";
                        }
                    }
                }
            }
            //std::cout<<"-->"<<local_output_buffer[0][trr][tcc];
        }
        //std::cout<<"\n";
    }


    BIAS_ADD_LOOP:for (int n = 0; n < N; n++)
    {
#pragma HLS loop_tripcount min=1 max=8
#pragma HLS PIPELINE
        for (int x = 0; x < O_h; x++)
        {
#pragma HLS loop_tripcount min=1 max=72
            for (int y = 0; y < (In_w + pad_right + pad_left - K_w + 1); y++)
            {
#pragma HLS loop_tripcount min=1 max=14
#pragma HLS PIPELINE
            	float biased = local_output_buffer[n][x][y] + local_bias[n];
                //local_output_buffer[n][x][y] = (biased > float(0)) ? biased : float(0);//Performing bias + ReLu action//removing ReLu activation
                local_output_buffer[n][x][y] = biased;
            }
        }
    }

    //std::cout<<"LOCAL OUTPUT BUFFER: dimensions are "<<O_w<<"*"<<O_h<<" \n";
    //for(int i=0;i<O_h;i++){
    //    for(int j=0;j<O_w;j++)
    //    {
    //        std::cout<<local_output_buffer[0][i][j]<<"\t";
    //    }
    //    std::cout<<"\n";
    //}



    OUT_LOOP:for (int p = 0; p < N; p++)
    {
#pragma HLS loop_tripcount min=1 max=8
#pragma HLS PIPELINE
        for (int iter = 0, i = 0, j = 0; iter < (In_w + pad_right + pad_left - K_w + 1) * O_h; iter++, j++)
        {
#pragma HLS loop_tripcount min=1 max=72*14
            if (j == (In_w + pad_right + pad_left - K_w + 1))
            {
                j = 0;
                i++;
            }
            output[iter + p * (In_w + pad_right + pad_left - K_w + 1) * O_h] = local_output_buffer[p][i][j];
            //std::cout<<"output assigned at "<<iter+ p*O_h*O_w<<"\t"<<local_output_buffer[p][i][j]<<"\t i="<<i<<" j="<<j<<"\n";
        }
    }

}
