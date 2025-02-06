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
#include <vector>
#include <numbers>

#include "rocalution.hpp"
#include "rocsolver.h"


#include "definitions.h"

#include "data.h"
#include "connect.h"

using real_t = type::real_t;


int main() {

  std::size_t vl_n = 4;
  std::size_t el_n = 3;

  const std::size_t l_blocks = el_n;
  const std::size_t n_blocks = l_blocks * l_blocks;

  auto g = data<real_t> {};

  vv<std::size_t> interfaces;
    std::size_t n_interfaces = make_connectivity(interfaces, l_blocks); 

  auto span = 1. / static_cast<double>(l_blocks);
  auto n_points_x = vl_n;
  auto n = vl_n;

  /*
  // Generate ranges for the x and y of each block, given the number of 
  // blocks in each dimension.
  auto block_grid = range_t(0, 1., l_blocks + 1);
  std::vector<range_t> grids; 
  for (std::size_t i = 0; i != block_grid.size() - 1; ++i)  {
    grids.push_back(range_t(*block_grid[i], *block_grid[i + 1], n));
  }

  rocblas_int M;
  rocblas_int N;
  rocblas_int lda;

  // here is where you would initialize M, N and lda with desired values

  rocblas_handle handle;
  rocblas_create_handle(&handle);

  size_t size_A = size_t(lda) * N;          // the size of the array for the matrix
  size_t size_piv = size_t(std::min(M, N)); // the size of array for the Householder scalars

  std::vector<double> hA(size_A);      // creates array for matrix in CPU
  std::vector<double> hIpiv(size_piv); // creates array for householder scalars in CPU

  double *dA, *dIpiv;
  hipMalloc(&dA, sizeof(double)*size_A);      // allocates memory for matrix in GPU
  hipMalloc(&dIpiv, sizeof(double)*size_piv); // allocates memory for scalars in GPU

  // here is where you would initialize matrix A (array hA) with input data
  // note: matrices must be stored in column major format,
  //       i.e. entry (i,j) should be accessed by hA[i + j*lda]

  // copy data to GPU
  hipMemcpy(dA, hA.data(), sizeof(double)*size_A, hipMemcpyHostToDevice);
  // compute the QR factorization on the GPU
  rocsolver_dgeqrf(handle, M, N, dA, lda, dIpiv);


    rocsolver_dcsrrf_solve()

  // copy the results back to CPU
  hipMemcpy(hA.data(), dA, sizeof(double)*size_A, hipMemcpyDeviceToHost);
  hipMemcpy(hIpiv.data(), dIpiv, sizeof(double)*size_piv, hipMemcpyDeviceToHost);

  // the results are now in hA and hIpiv, so you can use them here

  hipFree(dA);                        // de-allocate GPU memory
  hipFree(dIpiv);
  rocblas_destroy_handle(handle);     // destroy handle
  */
}


/*
int main()
{
    rocblas_int n = 10240;
    float alpha = 10.0;

    vector<float> hx(n);
    vector<float> hz(n);
    float* dx;

    rocblas_handle handle;
    rocblas_create_handle(&handle);

    // allocate memory on device
    hipMalloc(&dx, n * sizeof(float));

    // Initial Data on CPU,
    srand(1);
    for( int i = 0; i < n; ++i )
    {
        hx[i] = rand() % 10 + 1;  //generate a integer number between [1, 10]
    }

    // copy array from host memory to device memory
    hipMemcpy(dx, hx.data(), sizeof(float) * n, hipMemcpyHostToDevice);

    // call rocBLAS function
    rocblas_status status = rocblas_sscal(handle, n, &alpha, dx, 1);

    // check status for errors
    if(status == rocblas_status_success)
    {
        cout << "status == rocblas_status_success" << endl;
    }
    else
    {
        cout << "rocblas failure: status = " << status << endl;
    }

    // copy output from device memory to host memory
    hipMemcpy(hx.data(), dx, sizeof(float) * n, hipMemcpyDeviceToHost);

    hipFree(dx);
    rocblas_destroy_handle(handle);
    return 0;
}


__global__ void helloworld(char* in, char* out)
{
	int num = hipThreadIdx_x + hipBlockDim_x * hipBlockIdx_x;
	out[num] = in[num] + 1;
}

int main(int argc, char **argv) {

    hipDeviceProp_t devProp;
    hipGetDeviceProperties(&devProp, 0);
    std::cout 
        << "Device information [Minor [" << devProp.minor 
		<< "] Major [" << devProp.major << "] Name [" 
        << devProp.name  << "]]" << std::endl;

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
*/