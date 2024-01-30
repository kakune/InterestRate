/**
 * @file parameters.cpp
 * @brief This is implementation of functions in Utils::Parameters.
 * @author kakune
 * @date 1/28/2024
 */
#include "parameters.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

namespace Utils
{

void Parameters::readParameters( const std::string& inFilePath )
{
    std::ifstream lParamFile( inFilePath );
    if ( !lParamFile )
    {
        std::cerr << "Could not open parameter file: " << inFilePath
                  << std::endl;
        return;
    }

    std::string lLine;
    std::string lSection = "DEFAULT";
    while ( std::getline( lParamFile, lLine ) )
    {
        std::istringstream iss( lLine );
        std::string lName, lDelim;
        double lValue;
        if ( lLine[0] == '[' && lLine[lLine.size() - 1] == ']' )
        {
            lSection = lLine.substr( 1, lLine.size() - 2 );
        }
        if ( !( iss >> lName >> lDelim >> lValue ) ) { continue; }

        mData[lSection][lName] = lValue;
    }

    lParamFile.close();
}

}  // namespace Utils