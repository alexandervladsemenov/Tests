#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cublas_v2.h"
#include "cusolverDn.h"
#include <iostream>
#include <iomanip>
inline cublasStatus_t addingCublas(cublasHandle_t handle, float al, float* dev_a, float* dev_b, int size)
{
    return cublasSaxpy(handle, size, &al, dev_b, 1, dev_a, 1);
}
inline cublasStatus_t addingCublas(cublasHandle_t handle, double al, double* dev_a, double* dev_b, int size)
{
    return cublasDaxpy(handle, size, &al, dev_b, 1, dev_a, 1);
}
inline cublasStatus_t addingCublas(cublasHandle_t handle, cuComplex al, cuComplex* dev_a, cuComplex* dev_b, int size)
{
    return cublasCaxpy(handle, size, &al, dev_b, 1, dev_a, 1);
}
inline cublasStatus_t addingCublas(cublasHandle_t handle, cuDoubleComplex al, cuDoubleComplex* dev_a, cuDoubleComplex* dev_b, int size)
{
    return cublasZaxpy(handle, size, &al, dev_b, 1, dev_a, 1);
}
template<typename T>
cudaError_t addwithBlas(T* c, T* a, T* b, T al, int size)
{
    cudaError_t cudaStatus;
    cublasHandle_t handle; // CUBLAS context
    T* dev_a = 0;
    T* dev_b = 0;
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        std::cout << stderr << "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?";
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(T));
    if (cudaStatus != cudaSuccess) {
        std::cout << stderr << "cudaMalloc failed for array a!";
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(T));
    if (cudaStatus != cudaSuccess) {
        std::cout << stderr << "cudaMalloc failed for array b!";
        goto Error;
    }
    cublasCreate_v2(&handle);
    auto stat = cublasSetVector(size, sizeof(*a), a, 1, dev_a, 1);
    if (stat != CUBLAS_STATUS_SUCCESS)
        std::cout << "something went with copying to device operation for vector a";
    stat = cublasSetVector(size, sizeof(*b), b, 1, dev_b, 1);
    if (stat != CUBLAS_STATUS_SUCCESS)
        std::cout << "something went with copying to device operation for vector b";
    stat = addingCublas(handle, al, dev_a, dev_b, size);
    if (stat != CUBLAS_STATUS_SUCCESS)
        std::cout << "something went with addition operation";
    stat = cublasGetVector(size, sizeof(T), dev_a, 1, c, 1);
    if (stat != CUBLAS_STATUS_SUCCESS)
        std::cout << "something went with copying to host operation";
Error:
    cudaFree(dev_a);
    cudaFree(dev_b);
    cublasDestroy(handle);
    return cudaStatus;
}
void printArray(float* inpArray, int size)
{
    for (int i = 0; i< size; i++)
    {
        std::cout << std::setprecision(5) << inpArray[i] << ' ';
    }
}
void printArray(double* inpArray, int size)
{
    for (int i = 0; i < size; i++)
        std::cout << std::setprecision(10) << inpArray[i] << ' ';
}
void printArray(int* inpArray, int size)
{
    std::cout << std::setw(10);
    for (int i = 0; i < size; i++)
        std::cout << inpArray[i] << ' ';
}
int main()
{
    const int arraySize = 5;
    float aa[arraySize] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    float bb[arraySize] = { 10.0, 20.0, 30.0, 40.0, 50.0 };
    float cc[arraySize] = {0.0};
    float multiplier1 = 1/7.0;
    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    auto cudaStatus = addwithBlas<float>(cc, aa, bb, multiplier1, arraySize);
//    cudaStatus = addwithBlas<int>(c, a, b, multiplier, arraySize);
    std::cout << "COMPUTING with CUBLAS\n";
    printArray(cc, arraySize);
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        std::cout << stderr<< "cudaDeviceReset failed!";
        return 1;
    }
    return 0;
}
