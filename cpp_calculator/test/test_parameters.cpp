#include <iostream>

#include "utils/parameters.hpp"

int main( int argc, char* argv[] )
{
    std::string lFilePath = argv[1];
    Utils::Parameters lParams;
    lParams.readParameters( lFilePath );
    lParams.setCommonSectionName( "COMMON" );
    lParams.setCurrentSectionName( "PARAM1" );
    std::cout << lParams( "a" ) << std::endl;
    return 0;
}