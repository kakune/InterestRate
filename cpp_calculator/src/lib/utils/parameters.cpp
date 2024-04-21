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

template <> std::string castData( const DataType& inVal )
{
    return std::visit(
        []( auto&& inArg ) -> std::string
        {
            using Type = std::decay_t<decltype( inArg )>;
            if constexpr ( std::is_same_v<Type, std::string> ) { return inArg; }
            if constexpr ( std::is_same_v<Type, int> ||
                           std::is_same_v<Type, double> )
            {
                return std::to_string( inArg );
            }
            throw std::invalid_argument( "castData() mismatch." );
        },
        inVal );
}
template <> int castData( const DataType& inVal )
{
    return std::visit(
        []( auto&& inArg ) -> int
        {
            using Type = std::decay_t<decltype( inArg )>;
            if constexpr ( std::is_same_v<Type, int> ) { return inArg; }
            throw std::invalid_argument( "castData() mismatch." );
        },
        inVal );
}
template <> double castData( const DataType& inVal )
{
    return std::visit(
        []( auto&& inArg ) -> double
        {
            using Type = std::decay_t<decltype( inArg )>;
            if constexpr ( std::is_same_v<Type, double> ) { return inArg; }
            if constexpr ( std::is_same_v<Type, int> )
            {
                return static_cast<double>( inArg );
            }
            throw std::invalid_argument( "castData() mismatch." );
        },
        inVal );
}
template <> std::vector<int> castData( const DataType& inVal )
{
    return std::visit(
        []( auto&& inArg ) -> std::vector<int>
        {
            using Type = std::decay_t<decltype( inArg )>;
            if constexpr ( std::is_same_v<Type, std::vector<int>> )
            {
                return inArg;
            }
            if constexpr ( std::is_same_v<Type, int> )
            {
                return std::vector<int>( 1, inArg );
            }
            throw std::invalid_argument( "castData() mismatch." );
        },
        inVal );
}
template <> std::vector<double> castData( const DataType& inVal )
{
    return std::visit(
        []( auto&& inArg ) -> std::vector<double>
        {
            using Type = std::decay_t<decltype( inArg )>;
            if constexpr ( std::is_same_v<Type, std::vector<double>> )
            {
                return inArg;
            }
            if constexpr ( std::is_same_v<Type, std::vector<int>> )
            {
                return std::vector<double>( inArg.begin(), inArg.end() );
            }
            if constexpr ( std::is_same_v<Type, int> ||
                           std::is_same_v<Type, double> )
            {
                return std::vector<double>( 1, inArg );
            }
            throw std::invalid_argument( "castData() mismatch." );
        },
        inVal );
}
template <> std::vector<std::vector<int>> castData( const DataType& inVal )
{
    return std::visit(
        []( auto&& inArg ) -> std::vector<std::vector<int>>
        {
            using Type = std::decay_t<decltype( inArg )>;
            if constexpr ( std::is_same_v<Type, std::vector<std::vector<int>>> )
            {
                return inArg;
            }
            if constexpr ( std::is_same_v<Type, std::vector<int>> )
            {
                return std::vector<std::vector<int>>( 1, inArg );
            }
            if constexpr ( std::is_same_v<Type, int> )
            {
                return std::vector<std::vector<int>>(
                    1, std::vector<int>( inArg ) );
            }
            throw std::invalid_argument( "castData() mismatch." );
        },
        inVal );
}
template <> std::vector<std::vector<double>> castData( const DataType& inVal )
{
    return std::visit(
        []( auto&& inArg ) -> std::vector<std::vector<double>>
        {
            using Type = std::decay_t<decltype( inArg )>;
            if constexpr ( std::is_same_v<Type,
                                          std::vector<std::vector<double>>> )
            {
                return inArg;
            }
            if constexpr ( std::is_same_v<Type, std::vector<std::vector<int>>> )
            {
                std::vector<std::vector<double>> lResult;
                for ( const auto& subVec : inArg )
                {
                    std::vector<double> tempVec( subVec.begin(), subVec.end() );
                    lResult.push_back( tempVec );
                }
                return lResult;
            }
            if constexpr ( std::is_same_v<Type, std::vector<double>> )
            {
                return std::vector<std::vector<double>>( 1, inArg );
            }
            if constexpr ( std::is_same_v<Type, std::vector<int>> )
            {
                return std::vector<std::vector<double>>(
                    1, std::vector<double>( inArg.begin(), inArg.end() ) );
            }
            if constexpr ( std::is_same_v<Type, int> ||
                           std::is_same_v<Type, double> )
            {
                return std::vector<std::vector<double>>(
                    1, std::vector<double>( inArg ) );
            }
            throw std::invalid_argument( "castData() mismatch." );
        },
        inVal );
}

DataType parseNumber( const std::string& text )
{
    if ( text.find( '.' ) != std::string::npos ) { return std::stod( text ); }
    return std::stoi( text );
}

std::vector<int> parseIntVector( const std::string& text )
{
    std::vector<int> result;
    std::regex number_regex( "[-+]?[0-9]*\\.?[0-9]+" );
    auto numbers_begin =
        std::sregex_iterator( text.begin(), text.end(), number_regex );
    auto numbers_end = std::sregex_iterator();
    for ( std::sregex_iterator i = numbers_begin; i != numbers_end; ++i )
    {
        std::smatch match = *i;
        result.push_back( std::stoi( match.str() ) );
    }
    return result;
}

std::vector<double> parseDoubleVector( const std::string& text )
{
    std::vector<double> result;
    std::regex number_regex( "[-+]?[0-9]*\\.?[0-9]+" );
    auto numbers_begin =
        std::sregex_iterator( text.begin(), text.end(), number_regex );
    auto numbers_end = std::sregex_iterator();
    for ( std::sregex_iterator i = numbers_begin; i != numbers_end; ++i )
    {
        std::smatch match = *i;
        result.push_back( std::stod( match.str() ) );
    }
    return result;
}
std::vector<std::vector<int>> parseIntMatrix( const std::string& text )
{
    std::vector<std::vector<int>> result;
    std::regex vector_regex( "\\{([^}]*)\\}" );
    auto vectors_begin =
        std::sregex_iterator( text.begin(), text.end(), vector_regex );
    auto vectors_end = std::sregex_iterator();
    for ( std::sregex_iterator i = vectors_begin; i != vectors_end; ++i )
    {
        std::smatch match = *i;
        result.push_back( parseIntVector( match[1].str() ) );
    }
    return result;
}
std::vector<std::vector<double>> parseDoubleMatrix( const std::string& text )
{
    std::vector<std::vector<double>> result;
    std::regex vector_regex( "\\{([^}]*)\\}" );
    auto vectors_begin =
        std::sregex_iterator( text.begin(), text.end(), vector_regex );
    auto vectors_end = std::sregex_iterator();
    for ( std::sregex_iterator i = vectors_begin; i != vectors_end; ++i )
    {
        std::smatch match = *i;
        result.push_back( parseDoubleVector( match[1].str() ) );
    }
    return result;
}

bool isNumber( const std::string& s )
{
    std::regex number_pattern( "^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$" );
    return std::regex_match( s, number_pattern );
}

void Parameters::readParameters( const std::string& inFilePath )
{
    std::ifstream lParamFile( inFilePath );
    if ( !lParamFile )
    {
        std::cerr << "Could not open parameter file: " << inFilePath
                  << std::endl;
        return;
    }

    std::string lLine, lCurrentData;
    std::string lSection = "DEFAULT";
    std::regex lPattern( "^([^=]+)=\\s*(.*)$" );
    std::smatch lMatches;
    while ( std::getline( lParamFile, lLine ) )
    {
        std::istringstream lIss( lLine );
        lLine = std::regex_replace( lLine, std::regex( "\\s+" ), "" );
        if ( lLine.empty() ) { continue; }
        if ( lLine[0] == '[' && lLine[lLine.size() - 1] == ']' )
        {
            lSection     = lLine.substr( 1, lLine.size() - 2 );
            lCurrentData = "";
            continue;
        }
        lCurrentData += lLine;
        int lNOpen =
            std::count( lCurrentData.begin(), lCurrentData.end(), '{' );
        int lNClose =
            std::count( lCurrentData.begin(), lCurrentData.end(), '}' );
        if ( lNOpen != lNClose ) { continue; }

        if ( std::regex_match( lCurrentData, lMatches, lPattern ) )
        {
            std::string lKey = lMatches[1].str();
            std::string lVal = lMatches[2].str();

            DataType lData;
            if ( lVal.front() == '{' && lVal.back() == '}' )
            {
                if ( lVal.find( "{{" ) == 0 )
                {
                    if ( lVal.find( '.' ) != std::string::npos )
                    {
                        lData = parseDoubleMatrix( lVal );
                    }
                    else { lData = parseIntMatrix( lVal ); }
                }
                else
                {
                    if ( lVal.find( '.' ) != std::string::npos )
                    {
                        lData = parseDoubleVector( lVal );
                    }
                    else { lData = parseIntVector( lVal ); }
                }
            }
            else if ( isNumber( lVal ) ) { lData = parseNumber( lVal ); }
            else { lData = lVal; }

            mData[lSection][lKey] = lData;
            std::cout << "Read:" << lSection << ":" << lKey << "." << std::endl;
            lCurrentData = "";
        }
    }

    lParamFile.close();
}

}  // namespace Utils