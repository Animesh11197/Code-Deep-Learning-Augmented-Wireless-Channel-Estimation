#ifndef LS_HW_H_
#define LS_HW_H_
#include <iostream>
#include <stdlib.h>
#include<assert.h>
#include<cmath>


//#pragma SDS data mem_attribute(weight:PHYSICAL_CONTIGUOUS, bias:PHYSICAL_CONTIGUOUS,input: PHYSICAL_CONTIGUOUS,output: PHYSICAL_CONTIGUOUS)
//void perform_conv_hw(float input[72*14*8], float output[72*14*8], const float weight[36*7*8*8],  float bias[1*1*8], int M, int N, int In_h, int In_w, int K_h, int K_w, int pad_left, int pad_right, int pad_above, int pad_below);

#pragma SDS data mem_attribute(input: PHYSICAL_CONTIGUOUS,output: PHYSICAL_CONTIGUOUS)
void LS_hw(float input[24*2*2], float output[24*2*2]);



#endif
