#include "DNN.h"
#include <math.h>
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
typedef ap_fixed<32,12,AP_RND> dtype;
//typedef float dtype;
// M --> # of input channles/Neurons
// N --> # of output channels/neurons
// O --> Output Image Size
// I --> Input Image Size
// K --> Filter Size

// FC layer y = WX + b
void perform_dense_sw(float* input, float* output, const float* weight, const float* bias, int M, int N, int Relu)
{

    for (int n = 0; n < N; n++) { // for all outputs
        float temp = bias[n];
        for (int m = 0; m < M; m++) { // for all inputs
            int w_index = m + n * M; // Since weight is a 1-D Pointer, calculating the pointer index.
            temp += input[m] * weight[w_index]; //WX
        }
        //temp = temp + bias[n]; // + bias
        if (Relu == 1)
            temp = (temp > float(0)) ? temp : float(0); // ReLU
        //else if (Relu == 0)
            //temp = 1 / float(1 + exp(-temp)); // Sigmoid
        output[n] = temp;


    }
}

void perform_dense_L2(float input[48], float output[2016], const float weight[48*2016], const float bias[2016], int M, int N, int Relu)
{


	dtype local_input[48];
	dtype local_weight[48][2016];
#pragma HLS array_partition variable=local_input cyclic factor=4
#pragma HLS array_partition variable=local_weight dim=2 cyclic factor=4

	for(int i=0;i<48;i++)
	{
		#pragma HLS PIPELINE
		local_input[i]=input[i];
	}

	for(int i=0,k=0,j=0;i<48*2016;i++,j++)
	{
		#pragma HLS PIPELINE
		local_weight[k][j]=weight[i];
		if(j==47)
		{
			j=-1;
			k++;
		}
	}
/*	for(int k=0;k<2016;k++)
	{
#pragma HLS PIPELINE
		for (int j=0;j<48;j++)
		{

		local_weight[k][j]=weight[j+k*48];
		}
	}*/


    for (int n = 0; n < 2016; n++) { // for all outputs
//#pragma HLS LOOP_TRIPCOUNT max=2016
#pragma HLS PIPELINE
        dtype temp = bias[n];
        for (int m = 0; m < 48; m++) { // for all inputs
//#pragma HLS PIPELINE
            //int w_index = m + n * 48; // Since weight is a 1-D Pointer, calculating the pointer index.
            temp += local_input[m] * local_weight[n][m]; //WX
        }
        //temp = temp + bias[n]; // + bias
        if (Relu ==1)
            temp = (temp > dtype(0)) ? temp : dtype(0); // ReLU
        //else if (Relu == 0)
            //temp = 1 / float(1 + exp(-temp)); // Sigmoid
        output[n] = temp;


    }
}
