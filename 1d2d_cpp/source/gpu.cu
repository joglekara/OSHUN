#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <assert.h>

#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include "gpu.h"

#define BLOCKSIZE_x 32
#define BLOCKSIZE_y 4

/********************/
/* CUDA ERROR CHECK */
/********************/
// --- Credit to http://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %d, %s %s %d\n", code, cudaGetErrorString(code), file, line);
      if (abort) { exit(code); }
   }
}

extern "C" void gpuErrchk(cudaError_t ans) { gpuAssert((ans), __FILE__, __LINE__); }

/***************************/
/* CUSPARSE ERROR CHECKING */
/***************************/
static const char *_cusparseGetErrorEnum(cusparseStatus_t error)
{
    switch (error)
    {

        case CUSPARSE_STATUS_SUCCESS:
            return "CUSPARSE_STATUS_SUCCESS";

        case CUSPARSE_STATUS_NOT_INITIALIZED:
            return "CUSPARSE_STATUS_NOT_INITIALIZED";

        case CUSPARSE_STATUS_ALLOC_FAILED:
            return "CUSPARSE_STATUS_ALLOC_FAILED";

        case CUSPARSE_STATUS_INVALID_VALUE:
            return "CUSPARSE_STATUS_INVALID_VALUE";

        case CUSPARSE_STATUS_ARCH_MISMATCH:
            return "CUSPARSE_STATUS_ARCH_MISMATCH";

        case CUSPARSE_STATUS_MAPPING_ERROR:
            return "CUSPARSE_STATUS_MAPPING_ERROR";

        case CUSPARSE_STATUS_EXECUTION_FAILED:
            return "CUSPARSE_STATUS_EXECUTION_FAILED";

        case CUSPARSE_STATUS_INTERNAL_ERROR:
            return "CUSPARSE_STATUS_INTERNAL_ERROR";

        case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
            return "CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED";

        case CUSPARSE_STATUS_ZERO_PIVOT:
            return "CUSPARSE_STATUS_ZERO_PIVOT";
    }

    return "<unknown>";
}

inline void __cusparseSafeCall(cusparseStatus_t err, const char *file, const int line)
{
    if(CUSPARSE_STATUS_SUCCESS != err) {
        fprintf(stderr, "CUSPARSE error in file '%s', line %Ndims \n objs %s\nerror %Ndims: %s\nterminating! \n objs",__FILE__, __LINE__,err, \
                                _cusparseGetErrorEnum(err)); \
        cudaDeviceReset(); assert(0); \
    }
}

extern "C" void cusparseSafeCall(cusparseStatus_t err) { __cusparseSafeCall(err, __FILE__, __LINE__); }

/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/

void GPU_interface_routines::AllocateMatrixSystemOnHost(int totalsize, 
                double& ld, double &dd, double &ud, double &fin)
{
    double a;

    cudaHostAlloc((void**)ld, totalsize*sizeof(a),cudaHostAllocDefault);
    cudaHostAlloc((void**)dd, totalsize*sizeof(a),cudaHostAllocDefault);
    cudaHostAlloc((void**)ud, totalsize*sizeof(a),cudaHostAllocDefault);
    cudaHostAlloc((void**)fin, totalsize*sizeof(a),cudaHostAllocDefault);

}
void GPU_interface_routines::TDsolve( int calculations_per_loop, int n_systems,
                            double *ld, double *dd, double *ud, double *fin, 
                            int device)
{
    cudaSetDevice(device);
    // --- Initialize cuSPARSE
    cusparseHandle_t handle;    cusparseSafeCall(cusparseCreate(&handle));

    const int N     =  n_systems*calculations_per_loop;        // --- Size of the linear system

    // std::cout << "\n 10 \n";
    // // --- Lower diagonal, diagonal and upper diagonal of the system matrix
    double *d_ld;   gpuErrchk(cudaMalloc(&d_ld, N * sizeof(double)));
    double *d_d;    gpuErrchk(cudaMalloc(&d_d,  N * sizeof(double)));
    double *d_ud;   gpuErrchk(cudaMalloc(&d_ud, N * sizeof(double)));

    // double *d_ld;   double *d_d;  double *d_ud; double *d_x;
    // cudaGetSymbolAddress((void **)&d_ld, lowerdiagonal);
    // cudaGetSymbolAddress((void **)&d_d, diagonal);
    // cudaGetSymbolAddress((void **)&d_ud, upperdiagonal);
    // cudaGetSymbolAddress((void **)&d_x, solution);

    gpuErrchk(cudaMemcpy(d_ld, ld, N * sizeof(double), cudaMemcpyHostToDevice));//std::cout << "\n 11 \n";
    gpuErrchk(cudaMemcpy(d_d,  dd, N * sizeof(double), cudaMemcpyHostToDevice));//std::cout << "\n 12 \n";
    gpuErrchk(cudaMemcpy(d_ud, ud, N * sizeof(double), cudaMemcpyHostToDevice));//std::cout << "\n 13 \n";

    // --- Allocating and defining dense device data vectors
    double *d_x;        gpuErrchk(cudaMalloc(&d_x, N * sizeof(double)));   
    gpuErrchk(cudaMemcpy(d_x, fin, N * sizeof(double), cudaMemcpyHostToDevice));//std::cout << "\n 14 \n";

    // --- Solve for solution
    cusparseSafeCall(cusparseDgtsvStridedBatch(handle, calculations_per_loop, d_ld, d_d, d_ud, d_x, n_systems, calculations_per_loop));

    // --- Copy back into host
    cudaMemcpy(fin, d_x, N * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_ld);cudaFree(d_ud);cudaFree(d_d);cudaFree(d_x);
    cusparseSafeCall(cusparseDestroy(handle));
}
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
FokkerPlanckOnGPU::FokkerPlanckOnGPU()
{
    calc_per_loop = 0;
    n_sys = 0;

}
void FokkerPlanckOnGPU::initialize(int calculations_per_loop, int n_systems, int device)
        // : calc_per_loop(calculations_per_loop), n_sys(n_systems)
{
    calc_per_loop = calculations_per_loop;
    n_sys = n_systems;

    cudaSetDevice(device);
    gpuErrchk(cudaMalloc(&d_ld, calc_per_loop * n_sys * sizeof(double)));
    gpuErrchk(cudaMalloc(&d_d,  calc_per_loop * n_sys * sizeof(double)));
    gpuErrchk(cudaMalloc(&d_ud, calc_per_loop * n_sys * sizeof(double)));
    gpuErrchk(cudaMalloc(&d_x, calc_per_loop * n_sys  * sizeof(double)));   

}
/********//********//********//********//********//********//********//********//********//********/
void FokkerPlanckOnGPU::destroy(int device)
{
    cudaSetDevice(device);
    cudaFree(d_ld);cudaFree(d_ud);cudaFree(d_d);cudaFree(d_x);
}
/********//********//********//********//********//********//********//********//********//********/
void FokkerPlanckOnGPU::SolveTridiagonal(double *ld, double *dd, double *ud, double *fin, int device)
{
    cudaSetDevice(device);
    // --- Initialize cuSPARSE
    cusparseHandle_t handle;    cusparseSafeCall(cusparseCreate(&handle));

    const int N     =  n_sys*calc_per_loop;        // --- Size of the linear system

    gpuErrchk(cudaMemcpy(d_ld, ld, N * sizeof(double), cudaMemcpyHostToDevice));//std::cout << "\n 11 \n";
    gpuErrchk(cudaMemcpy(d_d,  dd, N * sizeof(double), cudaMemcpyHostToDevice));//std::cout << "\n 12 \n";
    gpuErrchk(cudaMemcpy(d_ud, ud, N * sizeof(double), cudaMemcpyHostToDevice));//std::cout << "\n 13 \n";

    // --- Allocating and defining dense device data vectors
    gpuErrchk(cudaMemcpy(d_x, fin, N * sizeof(double), cudaMemcpyHostToDevice));//std::cout << "\n 14 \n";

    // --- Solve for solution
    cusparseSafeCall(cusparseDgtsvStridedBatch(handle, calc_per_loop, d_ld, d_d, d_ud, d_x, n_sys, calc_per_loop));

    // --- Copy back into host
    cudaMemcpy(fin, d_x, N * sizeof(double), cudaMemcpyDeviceToHost);

    cusparseSafeCall(cusparseDestroy(handle));
}

/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
