#include <iostream>

#include "analytical/Black_Scholes.hpp"

int main( int argc, char* argv[] )
{
    double lInitS  = 1.0;
    double lStrike = 1.0;
    double lRate   = 0.0;
    double lVol    = 0.01;
    double lTMat   = 1.0;
    std::cout << Analytical::BlackScholes::europeanCallOptionPrice(
                     lInitS, lStrike, lRate, lVol, lTMat )
              << std::endl;
    std::cout << Analytical::BlackScholes::europeanPutOptionPrice(
                     lInitS, lStrike, lRate, lVol, lTMat )
              << std::endl;
    return 0;
}