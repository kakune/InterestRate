#include <iostream>
#include <memory>

#include "math/interpolate_1d.hpp"

int main( int argc, char* argv[] )
{
    std::vector<double> lX{ 0.0, 2.0, 4.0, 6.0, 8.0 };
    std::vector<double> lY{ 0.3, 0.7, 0.1, 0.6, 0.4 };
    Math::Interpolate1d::NewtonSpline lObj( 2 );
    lObj.build( std::make_shared<std::vector<double> >( lX ),
                std::make_shared<std::vector<double> >( lY ) );
    for ( std::size_t i = 0; i <= 10; ++i )
    {
        double lTmpX = 8.0 / 10.0 * i;
        std::cout << lTmpX << " : " << lObj( lTmpX ) << ", "
                  << lObj.deriv( lTmpX, 1 ) << std::endl;
    }
    return 0;
}