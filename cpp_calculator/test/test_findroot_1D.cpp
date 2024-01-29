#include <iostream>

#include "math/findroot_1d.hpp"

double testFunction( double x ) { return x * x - 2.0; }
class TestClass
{
private:
    double a;

public:
    TestClass( double inA ) : a( inA ) {}
    double testFunc( double x ) { return x * x - a; }
};
int main( int argc, char* argv[] )
{
    std::cout << Math::FindRoot1D::Brent( testFunction, 0.0, 2.0 ) << std::endl;
    TestClass Obj( 9.0 );
    auto boundFunc = [&Obj]( double x ) -> double { return Obj.testFunc( x ); };
    std::cout << Math::FindRoot1D::Brent( boundFunc, 0.0, 10.0 ) << std::endl;
    return 0;
}