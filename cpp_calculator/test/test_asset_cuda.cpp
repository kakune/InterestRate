#include <curand_kernel.h>

#include <iostream>
#include <vector>

#include "process_cuda/asset.hpp"

int main( int argc, char* argv[] )
{
    std::size_t lNTerm = 50;
    std::size_t lNPath = 10000000;
    double dt          = 1.0 / double( lNTerm - 1 );

    std::vector<double> lSqDtsVec( lNTerm, std::sqrt( dt ) );

    double* lSqDts;
    cudaMalloc( (void**)&lSqDts, lNTerm * sizeof( double ) );
    cudaMemcpy( lSqDts, lSqDtsVec.data(), lNTerm * sizeof( double ),
                cudaMemcpyHostToDevice );

    double lVol       = 0.1;
    double lInitPrice = 100.0;
    ProcessCUDA::Asset::BlackSholesForwardBuilder lBSBuilder;
    lBSBuilder.setNPath( lNPath );
    lBSBuilder.setSqDts( lNTerm, lSqDts, dt * ( lNTerm - 1 ) );
    lBSBuilder.setVol( lVol );
    lBSBuilder.setInitPrice( lInitPrice );

    auto lBSObj = lBSBuilder.build();
    std::cout << lBSObj.priceCallOption( 100.0 ) << std::endl;
    std::cout << lBSObj.impliedVolatility( 100.0 ) << std::endl;

    double lVolvol   = 0.3;
    double lInitVol  = 0.2;
    double lCorr     = -0.2;
    double lExponent = 1.0;
    ProcessCUDA::Asset::SABRForwardBuilder lSABRBuilder;
    lSABRBuilder.setNPath( lNPath );
    lSABRBuilder.setSqDts( lNTerm, lSqDts, dt * ( lNTerm - 1 ) );
    lSABRBuilder.setInitPrice( lInitPrice );
    lSABRBuilder.setInitVol( lInitVol );
    lSABRBuilder.setExponent( lExponent );
    lSABRBuilder.setVolvol( lVolvol );
    lSABRBuilder.setCorr( lCorr );

    auto lSABRObj = lSABRBuilder.build();
    std::cout << lSABRObj.priceCallOption( 100.0 ) << std::endl;
    std::cout << lSABRObj.impliedVolatility( 100.0 ) << std::endl;

    cudaFree( lSqDts );

    // std::cout << std::endl;
    // for ( double lStrike = 40.0; lStrike < 160.0; lStrike += 1.0 )
    // {
    //     std::cout << lStrike << " " << std::setprecision( 12 )
    //               << lBSObj.impliedVolatility( lStrike, lNTerm - 1 )
    //               << std::endl;
    // }
    // std::cout << std::endl;

    return 0;
}