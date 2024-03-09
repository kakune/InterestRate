#include <iostream>
#include <memory>

#include "math/interpolate_1d.hpp"

int main( int argc, char* argv[] )
{
    std::vector<double> lX{ 0.0, 2.0, 4.0, 6.0, 8.0 };
    std::vector<double> lY{ 0.0, 4.0, 16.0, 36.0, 64.0 };
    Math::Interpolate1d::NewtonSpline lObj( 3 );
    lObj.build( std::make_shared<std::vector<double> >( lX ),
                std::make_shared<std::vector<double> >( lY ) );
    for ( std::size_t i = 0; i <= 10; ++i )
    {
        double lTmpX = 8.0 / 10.0 * i;
        std::cout << lTmpX << " : " << lObj( lTmpX ) << ", "
                  << lObj.deriv( lTmpX, 1 ) << std::endl;
    }
    lObj.buildIntegral();
    std::cout << lObj.integral( 7.99 ) + lObj.integral( 7.99, 7.999 )
              << std::endl;

    return 0;
}