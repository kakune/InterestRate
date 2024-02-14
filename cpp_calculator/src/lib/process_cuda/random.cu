#include <curand_kernel.h>

#include <iostream>

#include "process_cuda/random.hpp"

__global__ void genBrownAntitheticKernel( double* inResult, int inNPath,
                                          int inNTerm, double* inTerms,
                                          unsigned long long seed = 5556 )
{
    std::size_t lIdx    = threadIdx.x + blockDim.x * blockIdx.x;
    std::size_t lStride = blockDim.x * gridDim.x;

    curandState state;
    for ( std::size_t iTerm = lIdx + 1; iTerm < inNTerm; iTerm += lStride )
    {
        curand_init( seed, iTerm, 0, &state );
        std::size_t lIndexStart = iTerm * inNPath;
        for ( std::size_t iPath = 0; iPath < inNPath; iPath += 2 )
        {
            double rnd = curand_normal_double( &state ) *
                         ( inTerms[iTerm] - inTerms[iTerm - 1] );
            inResult[lIndexStart + iPath]     = rnd;
            inResult[lIndexStart + iPath + 1] = -rnd;
        }
    }
}

__global__ void genBrownPlainKernel( double* inResult, int inNPath, int inNTerm,
                                     double* inTerms,
                                     unsigned long long seed = 5556 )
{
    std::size_t lIdx    = threadIdx.x + blockDim.x * blockIdx.x;
    std::size_t lStride = blockDim.x * gridDim.x;

    curandState state;
    for ( std::size_t iTerm = lIdx + 1; iTerm < inNTerm; iTerm += lStride )
    {
        curand_init( seed, iTerm, 0, &state );
        std::size_t lIndexStart = iTerm * inNPath;
        for ( std::size_t iPath = 0; iPath < inNPath; ++iPath )
        {
            double rnd = curand_normal_double( &state ) *
                         ( inTerms[iTerm] - inTerms[iTerm - 1] );
            inResult[lIndexStart + iPath] = rnd;
        }
    }
}

void ProcessCUDA::Random::PathBrownPlain::makeRandomVals(
    unsigned long long inSeed )
{
    cudaMalloc( (void**)&mpcRandomValues, mNPath * mNTerm * sizeof( double ) );
    int threadsPerBlock = 256;
    int blocksPerGrid =
        ( mNPath * mNTerm + threadsPerBlock - 1 ) / threadsPerBlock;
    genBrownPlainKernel<<<blocksPerGrid, threadsPerBlock>>>(
        mpcRandomValues, mNPath, mNTerm, mpcTerms, inSeed );

    // std::vector<double> h_result( mNPath * mNTerm );
    // cudaMemcpy( h_result.data(), mpcRandomValues,
    //             mNPath * mNTerm * sizeof( double ), cudaMemcpyDeviceToHost );

    // for ( int i = 0; i < mNPath; ++i )
    // {
    //     for ( int j = 0; j < mNTerm; ++j )
    //     {
    //         std::cout << h_result[i * mNTerm + j] << " ";
    //     }
    //     std::cout << "\n";
    // }
}

void ProcessCUDA::Random::PathBrownAntithetic::makeRandomVals(
    unsigned long long inSeed )
{
    cudaMalloc( (void**)&mpcRandomValues, mNPath * mNTerm * sizeof( double ) );
    int threadsPerBlock = 256;
    int blocksPerGrid =
        ( mNPath * mNTerm + threadsPerBlock - 1 ) / threadsPerBlock;
    genBrownAntitheticKernel<<<blocksPerGrid, threadsPerBlock>>>(
        mpcRandomValues, mNPath, mNTerm, mpcTerms, inSeed );

    // std::vector<double> h_result( mNPath * mNTerm );
    // cudaMemcpy( h_result.data(), mpcRandomValues,
    //             mNPath * mNTerm * sizeof( double ), cudaMemcpyDeviceToHost );

    // for ( int i = 0; i < mNPath; ++i )
    // {
    //     for ( int j = 0; j < mNTerm; ++j )
    //     {
    //         std::cout << h_result[j * mNPath + i] << " ";
    //     }
    //     std::cout << "\n";
    // }
}
