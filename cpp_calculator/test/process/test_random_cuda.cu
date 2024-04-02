#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include <iostream>
#include <vector>

#include "process/random.hpp"

__global__ void myKernel( Process::Random::StdBrownPlain* obj,
                          double* inpResult )
{
    obj->dGenerateRandomVal( inpResult );
}

int main( int argc, char* argv[] )
{
    std::size_t lNTerms = 5;
    std::size_t lNPath  = 30;
    double dt           = 0.1;
    std::vector<double> lTerms( lNTerms, 0 );
    for ( std::size_t iTerm = 1; iTerm < lNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + dt;
    }
    auto lsTerms = std::make_shared<std::vector<double>>( lTerms );
    Process::Random::StdBrownPlain lPathObj( lNPath, lsTerms );
    Process::Random::StdBrownPlain* lpdPathObj;
    cudaMalloc( &lpdPathObj, sizeof( Process::Random::StdBrownPlain ) );
    cudaMemcpy( lpdPathObj, &lPathObj, sizeof( Process::Random::StdBrownPlain ),
                cudaMemcpyHostToDevice );

    double* ldResult;
    cudaMalloc( &ldResult, sizeof( double ) );

    myKernel<<<1, 1>>>( lpdPathObj, ldResult );
    cudaDeviceSynchronize();
    double lResult;
    cudaMemcpy( &lResult, ldResult, sizeof( double ), cudaMemcpyDeviceToHost );
    std::cout << lResult << std::endl;

    // // デバイスからホストへ更新されたオブジェクトをコピー
    // cudaMemcpy( &h_obj, d_obj, sizeof( MyDeviceClass ),
    //             cudaMemcpyDeviceToHost );

    // // リソースを解放
    // cudaFree( d_obj );

    return 0;
}