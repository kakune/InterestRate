#include <iostream>
#include <vector>

#include "process/random.hpp"
#include "process_cuda/random.hpp"

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
    auto lsTerms = std::make_shared<std::vector<double> >( lTerms );
    ProcessCUDA::Random::PathBrownAntithetic lPathObj( lNPath, lsTerms );
    lPathObj.makeRandomVals();

    return 0;
}