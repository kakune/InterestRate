#include <iostream>
#include <vector>

#include "process/market.hpp"

int main( int argc, char* argv[] )
{
    std::size_t lNTerms = 10;
    std::size_t lNPath  = 10;
    double lDt          = 0.1;
    double lRate        = 0.05;
    std::vector<double> lTerms( lNTerms, 0 );
    std::vector<double> lZCB( lNTerms, 1.0 );
    for ( std::size_t iTerm = 1; iTerm < lNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
        lZCB[iTerm]   = lZCB[iTerm - 1] - 0.04;
    }
    auto lsTerms = std::make_shared<std::vector<double> >( lTerms );

    Process::Market::Data lObj( lsTerms );
    lObj.setZCB( lZCB );
    std::vector<double> lFR( lNTerms );
    for ( std::size_t iTerm = 0; iTerm < lNTerms; ++iTerm )
    {
        lFR.at( iTerm ) = lObj.mInterpInstantaneousForwardRate( lTerms[iTerm] );
        std::cout << lFR.at( iTerm ) << std::endl;
    }

    lObj.setForwardRate( lFR );
    for ( std::size_t iTerm = 0; iTerm < lNTerms; ++iTerm )
    {
        std::cout << lObj.mInterpZCB( lTerms[iTerm] ) << std::endl;
    }
    return 0;
}