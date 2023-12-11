#include"LS_hw.h"
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
#include "hls_half.h"
#include <math.h>
#include "ap_int.h"
#include <ap_fixed.h>
typedef float data_t;
//typedef half data_t;
//typedef ap_fixed<6,4,AP_RND> data_t;


void LS_hw(float input[24*2*2], float output[24*2*2])
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
