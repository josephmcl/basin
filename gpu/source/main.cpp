/*
#include "timing.h"

int main(int argc, char **argv) {
    return 0;
}
*/

#include <hip/hip_runtime.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>

#define SUCCESS 0
#define FAILURE 1

using namespace std;

__global__ void helloworld(char* in, char* out)
{
	int num = hipThreadIdx_x + hipBlockDim_x * hipBlockIdx_x;
	out[num] = in[num] + 1;
}

int main(int argc, char **argv) {

    hipDeviceProp_t devProp;
    hipGetDeviceProperties(&devProp, 0);
    std::cout 
        << "Device information. " << std::endl << "Minor ["
        << devProp.minor << "] Major [" << devProp.major << "] Name [" 
        << devProp.name  << "]" << std::endl;

	/* Initial input,output for the host and create memory objects for the kernel*/
	const char* input = "GdkknVnqkc";
	size_t strlength = strlen(input);
	cout << "input string:" << endl;
	cout << input << endl;
	char *output = (char*) malloc(strlength + 1);

	char* inputBuffer;
	char* outputBuffer;
	hipMalloc((void**)&inputBuffer, (strlength + 1) * sizeof(char));
    hipMalloc((void**)&outputBuffer, (strlength + 1) * sizeof(char));

    hipMemcpy(inputBuffer, input, (strlength + 1) * sizeof(char), hipMemcpyHostToDevice);

	hipLaunchKernelGGL(helloworld,
                  dim3(1),
                  dim3(strlength),
                  0, 0,
                  inputBuffer ,outputBuffer );

	hipMemcpy(output, outputBuffer,(strlength + 1) * sizeof(char), hipMemcpyDeviceToHost);

    hipFree(inputBuffer);
    hipFree(outputBuffer);

	output[strlength] = '\0';	//Add the terminal character to the end of output.
	cout << "\noutput string:" << endl;
	cout << output << endl;

	free(output);

	std::cout<<"Passed!\n";
	return SUCCESS;
}