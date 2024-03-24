#include <chrono>
#include <memory>
#include <vector>

#include "process/random.hpp"
int main()
{
    std::vector<double> lTerms{ 0.0, 1.0, 2.0 };
    auto lsTerms = std::make_shared<std::vector<double> >( lTerms );
    Process::RandomVec::PathBrownPlain lRandom( 100, lsTerms, 100000000 );
    lRandom.setIndexTime( 1 );
    auto lVec1  = lRandom.generateRandomVal();
    auto lVec2  = lRandom.generateRandomVal();
    auto start  = std::chrono::high_resolution_clock::now();
    double lRes = dot( lVec1, lVec2 );
    auto end    = std::chrono::high_resolution_clock::now();
    std::cout << lRes << std::endl;
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;
}
