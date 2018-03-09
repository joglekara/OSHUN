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
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
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
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
// -------------
// column (nump) derivative
// -------------
/********//********//********//********//********//********//********//********//********//********/
__global__ void e_times_derivative_p(double *f, double *df, int nump, int offset)
{  
    __shared__ double s_f[4][BLOCKSIZE_x + 4]; // 2-wide halo

    int ip   = threadIdx.x;
    int ix   = blockIdx.x*blockDim.y + threadIdx.y;

    int si = ip + 2;       // local i for shared memory access + halo offset
    int sj = threadIdx.y; // local j for shared memory access

    int globalIdx = offset + ix * nump + ip;

    s_f[sj][si] = f[globalIdx];

    __syncthreads();

    // fill in periodic images in shared memory array 
    if (ix < 2) 
    {
        s_f[sj][si-2]  = s_f[sj][si+BLOCKSIZE_x-3];
        s_f[sj][si+BLOCKSIZE_x] = s_f[sj][si+1];   
    }

    __syncthreads();

    df[globalIdx] = 
                        ( (4./3.) * ( s_f[sj][si+1] - s_f[sj][si-1] )
                        - (1./6.) * ( s_f[sj][si+2] - s_f[sj][si-2] ) );
}
/********//********//********//********//********//********//********//********//********//********/
// -------------
// column (numx) derivative
// -------------
/********//********//********//********//********//********//********//********//********//********/
/*
__global__ void derivative_x(double *f, double *df)
{
  __shared__ double s_f[BLOCKSIZE_y+4][4];

  int ip  = blockIdx.x*blockDim.x + threadIdx.x;
  int ix  = threadIdx.y;
  
  int si = threadIdx.x;
  int sj = ix + 2;

  int globalIdx =  ix * BLOCKSIZE_x + ip;

  s_f[sj][si] = f[globalIdx];

  __syncthreads();

  if (j < 2) {
    s_f[sj-2][si]  = s_f[sj+BLOCKSIZE_y-3][si];
    s_f[sj+BLOCKSIZE_y][si] = s_f[sj+1][si];
  }

  __syncthreads();

  df[globalIdx] = 
    ( 4./3. * ( s_f[sj+1][si] - s_f[sj-1][si] )
    - 1./6. * ( s_f[sj+2][si] - s_f[sj-2][si] ) );
} */
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
__global__ void v_times_derivative_x(double *f, double *df, int numx, int offset)
{
    __shared__ double s_f[BLOCKSIZE_y+4][32];

    int ip  = blockIdx.x*blockDim.x + threadIdx.x;
    int si = threadIdx.x;

    for (int ix = threadIdx.y; ix < BLOCKSIZE_y; ix += blockDim.y) {
    int globalIdx = offset + ix * BLOCKSIZE_x + ip;
    int sj = ix + 2;
    s_f[sj][si] = f[globalIdx];
    }

    __syncthreads();

    int sj = threadIdx.y + 2;
    if (sj < 2) 
    {
        s_f[sj-2][si]  = s_f[sj+BLOCKSIZE_y-3][si];
        s_f[sj+BLOCKSIZE_y][si] = s_f[sj+1][si];   
    }

    __syncthreads();

    for (int ix = threadIdx.y; ix < BLOCKSIZE_y; ix += blockDim.y) 
    {
        int globalIdx = offset + ix * BLOCKSIZE_x + ip;
        int sj = ix + 2;
        df[globalIdx] = ( (4./3.) * ( s_f[sj+1][si] - s_f[sj-1][si] )
                            - (1./6.) * ( s_f[sj+2][si] - s_f[sj-2][si] ) );
    }
}
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
void GPU_interface_routines::setupTDsolve(double *d_ld, double *d_d, double *d_ud, double *d_x, int device)
{
    cudaSetDevice(device);

    gpuErrchk(cudaMalloc(&d_ld, N * sizeof(double)));
    gpuErrchk(cudaMalloc(&d_d,  N * sizeof(double)));
    gpuErrchk(cudaMalloc(&d_ud, N * sizeof(double)));
}


void GPU_interface_routines::TDsolve( int calculations_per_loop, int n_systems,
                            double *ld, 
                            double *dd, 
                                  double *ud,      
                                  double *fin,// int device)
                                  double *d_ld, double *d_d, double *d_ud, double *d_x)
{
    cudaSetDevice(device);
    // --- Initialize cuSPARSE
    cusparseHandle_t handle;    cusparseSafeCall(cusparseCreate(&handle));

    const int N     =  n_systems*calculations_per_loop;        // --- Size of the linear system

    // // --- Lower diagonal, diagonal and upper diagonal of the system matrix
    // double *d_ld;   gpuErrchk(cudaMalloc(&d_ld, N * sizeof(double)));
    // double *d_d;    gpuErrchk(cudaMalloc(&d_d,  N * sizeof(double)));
    // double *d_ud;   gpuErrchk(cudaMalloc(&d_ud, N * sizeof(double)));

    gpuErrchk(cudaMemcpy(d_ld, ld, N * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_d,  dd, N * sizeof(double), cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(d_ud, ud, N * sizeof(double), cudaMemcpyHostToDevice));

    // --- Allocating and defining dense device data vectors
    // double *d_x;        gpuErrchk(cudaMalloc(&d_x, N * sizeof(double)));   
    gpuErrchk(cudaMemcpy(d_x, fin, N * sizeof(double), cudaMemcpyHostToDevice));

    // --- Solve for solution
    cusparseSafeCall(cusparseDgtsvStridedBatch(handle, calculations_per_loop, d_ld, d_d, d_ud, d_x, n_systems, calculations_per_loop));

    // --- Copy back into host
    cudaMemcpy(fin, d_x, N * sizeof(double), cudaMemcpyDeviceToHost);

    // cudaFree(d_ld);cudaFree(d_ud);cudaFree(d_d);cudaFree(d_x);
    cusparseSafeCall(cusparseDestroy(handle));
}
void GPU_interface_routines::setupTDsolve(double *d_ld, double *d_d, double *d_ud, double *d_x)
{
    cudaFree(d_ld);cudaFree(d_ud);cudaFree(d_d);cudaFree(d_x);
}

/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
void GPU_interface_routines::calc_fieldxdf(int numx, int nump, int numdist, double *fin, 
                                            // double *dxf, double *dvf, 
                                            // double *ex, double *vtemp, 
                                            int device)
{
    cudaSetDevice(device);

    /// ------------------------------
    /// Row and Column Sizing
    /// ------------------------------
    // const int nump = _nump;
    // const int numx = _numx;
    // const int numdist = _numdist;
    const int totalsize = (2*numdist)*numx*nump;
    
    /// ------------------------------
    /// Allocate and initialize Arrays
    /// ------------------------------
    double *d_dxf;     gpuErrchk(cudaMalloc(&d_dxf,    totalsize  * sizeof(double)));
    gpuErrchk(cudaMemset(d_dxf, 0.,             totalsize  * sizeof(double)));

    double *d_dvf;     gpuErrchk(cudaMalloc(&d_dvf,    totalsize  * sizeof(double)));
    gpuErrchk(cudaMemset(d_dvf, 0.,             totalsize  * sizeof(double)));

    // double *d_vtemp;   gpuErrchk(cudaMalloc(&d_vtemp,  nump    * sizeof(double)));
    // gpuErrchk(cudaMemcpy(d_vtemp, vtemp,    nump  * sizeof(double), cudaMemcpyHostToDevice));

    // double *d_ex;      gpuErrchk(cudaMalloc(&d_ex,     numx    * sizeof(double)));
    // gpuErrchk(cudaMemcpy(d_ex,  ex,          numx  * sizeof(double), cudaMemcpyHostToDevice));

    double *d_fin;        gpuErrchk(cudaMalloc(&d_fin, totalsize * sizeof(double)));   
    gpuErrchk(cudaMemcpy(d_fin, fin, totalsize * sizeof(double), cudaMemcpyHostToDevice));
    
    /// ------------------------------
    /// Parallelization Grid on GPU
    /// ------------------------------
    dim3 threadsPerBlock(BLOCKSIZE_x, BLOCKSIZE_y);
    dim3 numBlocks(nump / threadsPerBlock.x, numx / threadsPerBlock.y);

    /// ------------------------------
    /// vgradf or Edfdv
    /// ------------------------------
    // e_times_derivative_p<<<numBlocks, threadsPerBlock>>>(d_fin, d_ex,       d_dvf, nump);
    int baseidx;
    for (int id(0); id < 2*numdist; ++id)
    {
        baseidx = id*numx*nump;
        
        // gpuErrchk(cudaMemcpy(d_dxftemp, d_fin+baseidx, nump*numx * sizeof(double), cudaMemcpyDeviceToDevice));
        v_times_derivative_x<<<numBlocks, threadsPerBlock>>>(d_fin, d_dxf, nump, baseidx);
        // gpuErrchk(cudaMemcpy(d_fin+baseidx, d_dxf, nump*numx * sizeof(double), cudaMemcpyDeviceToDevice));

    }
    
    /// ------------------------------
    /// Back to CPU
    /// ------------------------------
    // thrust::copy(d_dvf,d_dvf+totalsize,dvf.begin());
    gpuErrchk(cudaMemcpy(fin, d_dxf, totalsize * sizeof(double), cudaMemcpyDeviceToHost));


    /// ------------------------------
    /// Free resources
    /// ------------------------------
    cudaFree(d_fin);cudaFree(d_dxf);cudaFree(d_dvf);
    // cudaFree(d_ex);cudaFree(d_vtemp);

}

/********//********//********//********//********//********//********//********//********//********/
/********//********//********//********//********//********//********//********//********//********/
