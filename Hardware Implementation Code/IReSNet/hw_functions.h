#ifndef HW_FUNCTIONS_H_
#define HW_FUNCTIONS_H_
#include <iostream>
#include <stdlib.h>
#include<assert.h>


#pragma SDS data mem_attribute(weight:PHYSICAL_CONTIGUOUS, bias:PHYSICAL_CONTIGUOUS,input: PHYSICAL_CONTIGUOUS,output: PHYSICAL_CONTIGUOUS)
void perform_conv_hw(float input[72*14*8], float output[72*14*8], const float weight[36*7*8*8],  float bias[1*1*8], int M, int N, int In_h, int In_w, int K_h, int K_w, int pad_left, int pad_right, int pad_above, int pad_below);





#endif
