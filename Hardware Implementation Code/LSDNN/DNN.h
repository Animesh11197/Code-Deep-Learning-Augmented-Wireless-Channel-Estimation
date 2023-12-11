#ifndef DNN_H_
#define DNN_H_
#include <iostream>
#include <stdlib.h>
#include<assert.h>
//#define dtype float



void perform_dense_sw(float* input, float* output, const float* weight, const float* bias, int M, int N, int Relu);
#pragma SDS data zero_copy(input[0:48],output[0:2016],weight[0:48*2016], bias[0:2016])
#pragma SDS data mem_attribute(weight:PHYSICAL_CONTIGUOUS, bias:PHYSICAL_CONTIGUOUS,input: PHYSICAL_CONTIGUOUS,output: PHYSICAL_CONTIGUOUS)
//void perform_dense(float* input, float* output, const float* weight, const float* bias, int Relu);
void perform_dense_L2(float input[48], float output[2016], const float weight[48*2016], const float bias[2016], int M, int N, int Relu);


#endif
