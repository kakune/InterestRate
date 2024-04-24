#include "utils/csv.hpp"

#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace Utils::CSV
{

std::map<std::string, std::vector<double>> readFile(
    const std::string& inPathCSV )
{
    std::map<std::string, std::vector<double>> lMapResult;

    std::ifstream file( inPathCSV );
    if ( !file.is_open() )
    {
        std::cerr << "Error: Utils::CSV::readFile()" << std::endl
                  << "Failed to open file: " << inPathCSV << std::endl;
        return lMapResult;
    }

    std::string line;
    std::getline( file, line );

    std::stringstream lStreamHeader( line );
    std::string lCellHeader;
    std::vector<std::string> lHeaderRow;
    while ( std::getline( lStreamHeader, lCellHeader, ',' ) )
    {
        lHeaderRow.push_back( lCellHeader );
    }

    while ( std::getline( file, line ) )
    {
        std::stringstream lStream( line );
        std::string lCell;
        std::vector<std::string> lRow;

        while ( std::getline( lStream, lCell, ',' ) )
        {
            lRow.push_back( lCell );
        }

        for ( size_t i = 0; i < lHeaderRow.size(); ++i )
        {
            try
            {
                double lValue = std::stod( lRow[i] );
                lMapResult[lHeaderRow[i]].push_back( lValue );
            }
            catch ( const std::invalid_argument& e )
            {
                std::cerr << "Error: Utils::CSV::readFile()" << std::endl
                          << "Invalid argument: " << lRow[i] << std::endl;
                lMapResult[lHeaderRow[i]].push_back(
                    std::numeric_limits<double>::quiet_NaN() );
            }
            catch ( const std::out_of_range& e )
            {
                std::cerr << "Error: Utils::CSV::readFile()" << std::endl
                          << "Out of range: " << lRow[i] << std::endl;
                lMapResult[lHeaderRow[i]].push_back(
                    std::numeric_limits<double>::quiet_NaN() );
            }
        }
    }

    file.close();
    return lMapResult;
}

void prepareZCBColumn( std::map<std::string, std::vector<double>>& inMapMarket,
                       const Parameters& inParams )
{
    if ( inMapMarket.find( "ZCB" ) != inMapMarket.end() ) { return; }
    if ( inMapMarket.find( "ZCBRate" ) != inMapMarket.end() )
    {
        double lFactor   = inParams( "ZCBRateFactor" );
        double lCompound = inParams( "Compound" );
        if ( lFactor <= 0.0 ) { lFactor = 1.0; }
        if ( lCompound <= 0.0 ) { lCompound = 1.0; }
        for ( std::size_t i = 0; i < inMapMarket["Maturity"].size(); ++i )
        {
            inMapMarket["ZCB"].push_back(
                std::pow( 1.0 + lFactor * inMapMarket["ZCBRate"][i],
                          -inMapMarket["Maturity"][i] / lCompound ) );
        }
        return;
    }
    std::cerr << "Error: Utils::CSV::prepareZCBColumn()" << std::endl
              << "Cannot find corresponding column." << std::endl;
    return;
}

}  // namespace Utils::CSV