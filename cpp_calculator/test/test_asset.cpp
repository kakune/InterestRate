#include <iomanip>
#include <iostream>
#include <vector>

#include "analytical/Black_Scholes.hpp"
#include "process/asset.hpp"
#include "process/random.hpp"

int main( int argc, char* argv[] )
{
    std::size_t lNTerms = 100;
    std::size_t lNPath  = 100000;
    double lDt          = 0.1;
    double lVol         = 0.5;
    double lInitPrice   = 100.0;
    std::vector< double > lTerms( lNTerms, 0 );
    for ( std::size_t iTerm = 1; iTerm < lNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    auto lsTerms = std::make_shared< std::vector< double > >( lTerms );

    auto luRandom = std::make_unique< Process::Random::PathBrownAntithetic >(
        lNPath, lsTerms );

    Process::Asset::BlackSholesForwardBuilder lBuilder;
    lBuilder.setNPath( lNPath );
    lBuilder.setTerms( lsTerms );
    lBuilder.setVol( lVol );
    lBuilder.setRandomPath( std::move( luRandom ) );
    lBuilder.setInitPrice( lInitPrice );

    auto lObj = lBuilder.build();

    std::cout << std::endl;
    for ( double lStrike = 80.0; lStrike < 120.0; lStrike += 1.0 )
    {
        std::cout << lStrike << " " << std::setprecision( 12 )
                  << lObj.impliedVolatility( lStrike, lNTerms - 1 )
                  << std::endl;
    }
    std::cout << std::endl;
    return 0;
}