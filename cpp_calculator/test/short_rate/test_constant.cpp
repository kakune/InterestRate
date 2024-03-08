#include <iostream>
#include <vector>

#include "process/short_rate.hpp"

int main( int argc, char* argv[] )
{
    std::size_t lNTerms = 10;
    std::size_t lNPath  = 10;
    double lDt          = 0.1;
    double lRate        = 0.05;
    std::vector<double> lTerms( lNTerms, 0 );
    for ( std::size_t iTerm = 1; iTerm < lNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    auto lsTerms = std::make_shared<std::vector<double> >( lTerms );

    Process::ShortRate::ConstantRate lObj( lsTerms, nullptr, lRate );
    lObj.build();

    for ( std::size_t iTerm = 0; iTerm < lNTerms; ++iTerm )
    {
        std::cout << lTerms[iTerm] << " ";
    }
    std::cout << std::endl;
    double x = 0.0;
    for ( std::size_t iTerm = 0; iTerm < lNTerms; ++iTerm )
    {
        std::cout << lObj.priceZCB( 0.0, lTerms[iTerm] ) << " ";
    }
    std::cout << std::endl;
    return 0;
}