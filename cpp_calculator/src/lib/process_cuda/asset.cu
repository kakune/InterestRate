/**
 * @file asset.cu
 * @brief This implements calculator using CUDA of assets.
 * @author kakune
 * @date 2/14/2024
 */

#define THREADS_PER_BLOCK 256

#include <curand_kernel.h>

#include "process_cuda/asset.hpp"

namespace ProcessCUDA
{
namespace Asset
{
__global__ void priceCallOptionKernel( double* inResult, std::size_t inNPath,
                                       double inStrike, double* inPrice )
{
    extern __shared__ double lData[];
    std::size_t lTid      = threadIdx.x;
    std::size_t lInd      = blockIdx.x * ( blockDim.x * 2 ) + threadIdx.x;
    std::size_t lSizeGrid = blockDim.x * 2 * gridDim.x;
    lData[lTid]           = 0;

    while ( lInd < inNPath )
    {
        lData[lTid] += ( lInd < inNPath && inPrice[lInd] - inStrike > 0.0 )
                           ? inPrice[lInd] - inStrike
                           : 0.0;
        lData[lTid] += ( lInd + blockDim.x < inNPath &&
                         inPrice[lInd + blockDim.x] - inStrike > 0.0 )
                           ? inPrice[lInd + blockDim.x] - inStrike
                           : 0.0;
        lInd += lSizeGrid;
    }
    __syncthreads();
    for ( std::size_t i = blockDim.x / 2; i > 0; i >>= 1 )
    {
        if ( lTid < i ) { lData[lTid] += lData[lTid + i]; }
        __syncthreads();
    }
    if ( lTid == 0 ) { inResult[blockIdx.x] = lData[0]; }
}

__global__ void pricePutOptionKernel( double* inResult, std::size_t inNPath,
                                      double inStrike, double* inPrice )
{
    extern __shared__ double lData[];
    std::size_t lTid      = threadIdx.x;
    std::size_t lInd      = blockIdx.x * ( blockDim.x * 2 ) + threadIdx.x;
    std::size_t lSizeGrid = blockDim.x * 2 * gridDim.x;
    lData[lTid]           = 0;

    while ( lInd < inNPath )
    {
        lData[lTid] += ( lInd < inNPath && inStrike - inPrice[lInd] > 0.0 )
                           ? inStrike - inPrice[lInd]
                           : 0.0;
        lData[lTid] += ( lInd + blockDim.x < inNPath &&
                         inStrike - inPrice[lInd + blockDim.x] > 0.0 )
                           ? inStrike - inPrice[lInd + blockDim.x]
                           : 0.0;
        lInd += lSizeGrid;
    }
    __syncthreads();
    for ( std::size_t i = blockDim.x / 2; i > 0; i >>= 1 )
    {
        if ( lTid < i ) { lData[lTid] += lData[lTid + i]; }
        __syncthreads();
    }
    if ( lTid == 0 ) { inResult[blockIdx.x] = lData[0]; }
}

__global__ void calcBlackEachForwardKernel( double* inResult,
                                            std::size_t inNPath,
                                            std::size_t inNTerm, double inVol,
                                            double inInitPrice, double* inSqDts,
                                            unsigned long long inSeed )
{
    std::size_t lIdx    = threadIdx.x + blockDim.x * blockIdx.x;
    std::size_t lStride = blockDim.x * gridDim.x;

    curandState state;
    std::size_t lHalfPath = inNPath / 2;
    for ( std::size_t iPath = lIdx; iPath < lHalfPath; iPath += lStride )
    {
        std::size_t lInd1 = iPath * 2;
        std::size_t lInd2 = lInd1 + 1;
        inResult[lInd1]   = inInitPrice;
        inResult[lInd2]   = inInitPrice;
        curand_init( inSeed, iPath, 0, &state );
        for ( std::size_t iTerm = 1; iTerm < inNTerm; ++iTerm )
        {
            double lRnd =
                inVol * inSqDts[iTerm] * curand_normal_double( &state );
            inResult[lInd1] += inResult[lInd1] * lRnd;
            inResult[lInd2] -= inResult[lInd2] * lRnd;
        }
    }
}

__global__ void calcBlackEachForwardWithLogKernel(
    double* inResult, std::size_t inNPath, std::size_t inNTerm, double inVol,
    double inInitPrice, double* inSqDts, unsigned long long inSeed )
{
    std::size_t lIdx    = threadIdx.x + blockDim.x * blockIdx.x;
    std::size_t lStride = blockDim.x * gridDim.x;

    curandState state;
    std::size_t lHalfPath = inNPath / 2;
    double lLogInitPrice  = log( inInitPrice );
    double lVol2          = -0.5 * inVol * inVol;
    for ( std::size_t iPath = lIdx; iPath < lHalfPath; iPath += lStride )
    {
        std::size_t lInd1 = iPath * 2;
        std::size_t lInd2 = lInd1 + 1;
        inResult[lInd1]   = lLogInitPrice;
        inResult[lInd2]   = lLogInitPrice;
        curand_init( inSeed, iPath, 0, &state );
        for ( std::size_t iTerm = 1; iTerm < inNTerm; ++iTerm )
        {
            double lRnd =
                inVol * inSqDts[iTerm] * curand_normal_double( &state );
            double lDrift = lVol2 * inSqDts[iTerm] * inSqDts[iTerm];
            inResult[lInd1] += lRnd + lDrift;
            inResult[lInd2] += -lRnd + lDrift;
        }
        inResult[lInd1] = exp( inResult[lInd1] );
        inResult[lInd2] = exp( inResult[lInd2] );
    }
}

__global__ void calcSABREachForwardKernel(
    double* inResult, std::size_t inNPath, std::size_t inNTerm,
    double inInitPrice, double inInitVol, double inCorr, double inExponent,
    double inVolvol, double* inSqDts, unsigned long long inSeed )
{
    std::size_t lIdx    = threadIdx.x + blockDim.x * blockIdx.x;
    std::size_t lStride = blockDim.x * gridDim.x;

    curandState state;
    double lAuxCorr       = sqrt( 1.0 - inCorr * inCorr );
    std::size_t lHalfPath = inNPath / 2;
    for ( std::size_t iPath = lIdx; iPath < lHalfPath; iPath += lStride )
    {
        std::size_t lInd1 = iPath * 2;
        std::size_t lInd2 = lInd1 + 1;
        inResult[lInd1]   = inInitPrice;
        inResult[lInd2]   = inInitPrice;
        double lVol1      = inInitVol;
        double lVol2      = inInitVol;
        curand_init( inSeed, iPath, 0, &state );
        for ( std::size_t iTerm = 1; iTerm < inNTerm; ++iTerm )
        {
            double lBrown1 = inSqDts[iTerm] * curand_normal_double( &state );
            double lBrown2 = inSqDts[iTerm] * curand_normal_double( &state );
            lBrown2        = inCorr * lBrown1 + lAuxCorr * lBrown2;
            inResult[lInd1] +=
                lVol1 * pow( inResult[lInd1], inExponent ) * lBrown1;
            inResult[lInd2] -=
                lVol2 * pow( inResult[lInd2], inExponent ) * lBrown1;
            lVol1 += inVolvol * lVol1 * lBrown2;
            lVol2 -= inVolvol * lVol2 * lBrown2;
        }
    }
}
__global__ void calcSABREachForwardWithLogKernel(
    double* inResult, std::size_t inNPath, std::size_t inNTerm,
    double inInitPrice, double inInitVol, double inCorr, double inExponent,
    double inVolvol, double* inSqDts, unsigned long long inSeed )
{
    std::size_t lIdx    = threadIdx.x + blockDim.x * blockIdx.x;
    std::size_t lStride = blockDim.x * gridDim.x;

    curandState state;
    double lAuxCorr       = sqrt( 1.0 - inCorr * inCorr );
    double lLogInitPrice  = log( inInitPrice );
    double lLogInitVol    = log( inInitVol );
    double lVolvol2       = inVolvol * inVolvol;
    double lExponentM     = inExponent - 1.0;
    std::size_t lHalfPath = inNPath / 2;
    for ( std::size_t iPath = lIdx; iPath < lHalfPath; iPath += lStride )
    {
        std::size_t lInd1 = iPath * 2;
        std::size_t lInd2 = lInd1 + 1;
        inResult[lInd1]   = lLogInitPrice;
        inResult[lInd2]   = lLogInitPrice;
        double lLogVol1   = lLogInitVol;
        double lLogVol2   = lLogInitVol;
        curand_init( inSeed, iPath, 0, &state );
        for ( std::size_t iTerm = 1; iTerm < inNTerm; ++iTerm )
        {
            double lBrown1 = inSqDts[iTerm] * curand_normal_double( &state );
            double lBrown2 = inSqDts[iTerm] * curand_normal_double( &state );
            lBrown2        = inCorr * lBrown1 + lAuxCorr * lBrown2;

            double lFac1 =
                exp( lLogVol1 ) * pow( exp( inResult[lInd1] ), lExponentM );
            double lFac2 =
                exp( lLogVol2 ) * pow( exp( inResult[lInd2] ), lExponentM );

            double lHalfDt = -0.5 * inSqDts[iTerm] * inSqDts[iTerm];
            inResult[lInd1] += lFac1 * ( lBrown1 + lFac1 * lHalfDt );
            inResult[lInd2] += lFac2 * ( -lBrown1 + lFac2 * lHalfDt );

            lLogVol1 += inVolvol * lBrown2 + lVolvol2 * lHalfDt;
            lLogVol2 += -inVolvol * lBrown2 + lVolvol2 * lHalfDt;
        }
        inResult[lInd1] = exp( inResult[lInd1] );
        inResult[lInd2] = exp( inResult[lInd2] );
    }
}

double ModelForwardAbstract::priceCallOption( double inStrike )
{
    std::size_t threadsPerBlock = THREADS_PER_BLOCK;
    std::size_t blocksPerGrid =
        ( mNPath + threadsPerBlock - 1 ) / threadsPerBlock;
    std::size_t smemSize = threadsPerBlock * sizeof( double );
    double* lpcResult;
    cudaMalloc( (void**)&lpcResult, blocksPerGrid * sizeof( double ) );
    double* lpResult = new double[blocksPerGrid];
    priceCallOptionKernel<<<blocksPerGrid, threadsPerBlock, smemSize>>>(
        lpcResult, mNPath, inStrike, mpcForwardPrice );
    cudaMemcpy( lpResult, lpcResult, blocksPerGrid * sizeof( double ),
                cudaMemcpyDeviceToHost );
    double lResult = 0;
    for ( std::size_t i = 0; i < blocksPerGrid; ++i )
    {
        lResult += lpResult[i];
    }
    cudaFree( lpcResult );
    delete[] lpResult;

    return lResult / mNPath;
}

double ModelForwardAbstract::pricePutOption( double inStrike )
{
    std::size_t threadsPerBlock = THREADS_PER_BLOCK;
    std::size_t blocksPerGrid =
        ( mNPath + threadsPerBlock - 1 ) / threadsPerBlock;
    std::size_t smemSize = threadsPerBlock * sizeof( double );
    double* lpcResult;
    cudaMalloc( (void**)&lpcResult, blocksPerGrid * sizeof( double ) );
    double* lpResult = new double[blocksPerGrid];
    pricePutOptionKernel<<<blocksPerGrid, threadsPerBlock, smemSize>>>(
        lpcResult, mNPath, inStrike, mpcForwardPrice );
    cudaMemcpy( lpResult, lpcResult, blocksPerGrid * sizeof( double ),
                cudaMemcpyDeviceToHost );
    double lResult = 0;
    for ( std::size_t i = 0; i < blocksPerGrid; ++i )
    {
        lResult += lpResult[i];
    }
    cudaFree( lpcResult );
    delete[] lpResult;

    return lResult / mNPath;
}

void BlackSholesForward::calcEachForwardPrice()
{
    cudaMalloc( (void**)&mpcForwardPrice, mNPath * sizeof( double ) );
    std::size_t threadsPerBlock = THREADS_PER_BLOCK;
    std::size_t blocksPerGrid =
        ( mNPath + threadsPerBlock - 1 ) / threadsPerBlock;
    calcBlackEachForwardWithLogKernel<<<blocksPerGrid, threadsPerBlock>>>(
        mpcForwardPrice, mNPath, mNTerm, mVol, mInitPrice, mpcSqDts );
}

void SABRForward::calcEachForwardPrice()
{
    cudaMalloc( (void**)&mpcForwardPrice, mNPath * sizeof( double ) );
    std::size_t threadsPerBlock = THREADS_PER_BLOCK;
    std::size_t blocksPerGrid =
        ( mNPath + threadsPerBlock - 1 ) / threadsPerBlock;
    calcSABREachForwardWithLogKernel<<<blocksPerGrid, threadsPerBlock>>>(
        mpcForwardPrice, mNPath, mNTerm, mInitPrice, mInitVol, mCorr, mExponent,
        mVolvol, mpcSqDts );
}

}  // namespace Asset
}  // namespace ProcessCUDA