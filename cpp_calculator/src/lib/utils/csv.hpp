#ifndef UTILS_CSV_HPP
#define UTILS_CSV_HPP

#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Utils
{
namespace CSV
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
    std::getline( file, line );  // 最初の行を読み飛ばす（ヘッダー行）

    // ヘッダー行から列名を取得
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

        // 各列の値をmapのvectorに追加
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
}  // namespace CSV
}  // namespace Utils

#endif