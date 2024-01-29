#include <iostream>
#include <vector>

#include "process/random.hpp"

int main( int argc, char* argv[] )
{
    std::size_t lNTerms = 10;
    std::size_t lNPath  = 10;
    double dt           = 0.1;
    std::vector< double > lTerms( lNTerms, 0 );
    for ( std::size_t iTerm = 1; iTerm < lNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + dt;
    }
    auto lsTerms = std::make_shared< std::vector< double > >( lTerms );
    Process::Random::PathBrownAntithetic lPathObj( lNPath, lsTerms );
    lPathObj.makePath();
    for ( std::size_t iTerm = 0; iTerm < lNTerms; ++iTerm )
    {
        std::cout << lTerms[iTerm] << " ";
    }
    std::cout << std::endl;
    double x = 0.0;
    for ( std::size_t iPath = 0; iPath < lNPath; ++iPath )
    {
        for ( std::size_t iTerm = 0; iTerm < lNTerms; ++iTerm )
        {
            std::cout << lPathObj[iPath][iTerm] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << lPathObj( 0, dt * 0.1 ) << std::endl;
    return 0;
}