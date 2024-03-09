#include <iostream>
#include <memory>
#include <vector>

#include "process/market.hpp"
#include "process/random.hpp"
#include "process/short_rate.hpp"

int main( int argc, char* argv[] )
{
    std::size_t lNTerms = 100;
    std::size_t lNPath  = 100000;
    double lDt          = 0.01;
    std::vector<double> lTerms( lNTerms, 0 );
    std::vector<double> lZCB( lNTerms, 1.0 );
    for ( std::size_t iTerm = 1; iTerm < lNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
        lZCB[iTerm]   = lZCB[iTerm - 1] - 0.004;
    }
    auto lsTerms  = std::make_shared<std::vector<double> >( lTerms );
    auto luRandom = std::make_unique<Process::Random::PathBrownAntithetic>(
        lNPath, lsTerms );

    Process::Market::Data lMarketData( lsTerms );
    lMarketData.setZCB( lZCB );

    Process::ShortRate::HoLeeBuilder lBuilder;
    lBuilder.setTerms( lsTerms );
    lBuilder.setMarketData(
        std::make_shared<Process::Market::Data>( lMarketData ) );
    lBuilder.setNPath( lNPath );
    lBuilder.setVol( 0.2 );
    lBuilder.setRandom( std::move( luRandom ) );

    Process::ShortRate::HoLee lHoLee = lBuilder.build();

    lHoLee.build();
    for ( std::size_t iTerm = 0; iTerm < lNTerms; ++iTerm )
    {
        std::cout << lZCB.at( iTerm ) << ":"
                  << lHoLee.priceZCB( 0.0, lsTerms->at( iTerm ) ) << std::endl;
    }

    // Process::ShortRate::ConstantRate lObj( lsTerms, nullptr, lRate );
    // lObj.build();

    // for ( std::size_t iTerm = 0; iTerm < lNTerms; ++iTerm )
    // {
    //     std::cout << lTerms[iTerm] << " ";
    // }
    // std::cout << std::endl;
    // double x = 0.0;
    // for ( std::size_t iTerm = 0; iTerm < lNTerms; ++iTerm )
    // {
    //     std::cout << lObj.priceZCB( 0.0, lTerms[iTerm] ) << " ";
    // }
    // std::cout << std::endl;
    return 0;
}