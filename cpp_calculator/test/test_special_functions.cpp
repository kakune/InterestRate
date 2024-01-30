#include <iomanip>
#include <iostream>
#include <vector>

#include "math/special_functions.hpp"

int main( int argc, char* argv[] )
{
    std::cout << Math::SpecialFunctions::normalCDF( -0.0 ) << std::endl;
    return 0;
}