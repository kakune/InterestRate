/**
 * @file parameters.hpp
 * @brief define struct "Utils::Parameters" that organize numerical parameters
 * @author kakune
 * @date 1/28/2024
 */
#ifndef UTILS_PARAMETERS_HPP
#define UTILS_PARAMETERS_HPP

#include <map>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace Utils
{
using DataType =
    std::variant<std::string, int, double, std::vector<int>,
                 std::vector<double>, std::vector<std::vector<int>>,
                 std::vector<std::vector<double>>>;
template <typename T> T castData( const DataType& inVal )
{
    throw std::invalid_argument( "castData() mismatch." );
}
template <> std::string castData( const DataType& inVal );
template <> int castData( const DataType& inVal );
template <> double castData( const DataType& inVal );
template <> std::vector<double> castData( const DataType& inVal );
template <> std::vector<int> castData( const DataType& inVal );
template <> std::vector<std::vector<int>> castData( const DataType& inVal );
template <> std::vector<std::vector<double>> castData( const DataType& inVal );
/**
 * @brief This stores parameters in map.
 * @details Parameters[S][T] is a double variable T in section S.
 * S and T are std::string.
 * Parameters that are not in any section are considered to be in section
 * "DEFAULT".
 */
class Parameters
{
private:
    //! data read from ini file
    std::map<std::string, std::map<std::string, DataType>> mData;
    //! name of section that used operator()
    std::string mNameCurrentSection;
    //! name of section that used if current section has no corresponding
    //! parameters
    std::string mNameCommonSection = "DEFAULT";

public:
    /**
     * @brief This reads ini file.
     * @param inFilePath The path to ini file to read
     */
    void readParameters( const std::string& inFilePath );

    /**
     * @brief This accesses the map of each section.
     * @param inSectionName Name of section
     * @return const std::map< std::string, double >&
     */
    const std::map<std::string, DataType>& operator[](
        const std::string& inSectionName ) const
    {
        return mData.at( inSectionName );
    }

    /**
     * @brief This accesses the map of each section.
     * @param inSectionName Name of section
     * @return std::map< std::string, double >&
     */
    std::map<std::string, DataType>& operator[](
        const std::string& inSectionName )
    {
        return mData[inSectionName];
    }

    /**
     * @brief This accesses the parameter using mNameCurrentSection and
     * mNameCommonSection.
     * @details Return parameters in mNameCurrentSection.
     * If there is no parameter corresponding to that section,
     * return the parameter in mNameCommonSection.
     * @param inParameterName
     * @return double
     */
    template <typename T = double>
    T operator()( const std::string& inParameterName )
    {
        if ( mData[mNameCurrentSection].find( inParameterName ) ==
             mData[mNameCurrentSection].end() )
        {
            return castData<T>( mData[mNameCommonSection][inParameterName] );
        }
        return castData<T>( mData[mNameCurrentSection][inParameterName] );
    }

    /**
     * @brief This accesses the parameter using mNameCurrentSection and
     * mNameCommonSection.
     * @details Return parameters in mNameCurrentSection.
     * If there is no parameter corresponding to that section,
     * return the parameter in mNameCommonSection.
     * @param inParameterName
     * @return const double
     */
    template <typename T = double>
    const T operator()( const std::string& inParameterName ) const
    {
        if ( mData.at( mNameCurrentSection ).find( inParameterName ) ==
             mData.at( mNameCurrentSection ).end() )
        {
            return castData<T>(
                mData.at( mNameCommonSection ).at( inParameterName ) );
        }
        return castData<T>(
            mData.at( mNameCurrentSection ).at( inParameterName ) );
    }

    /**
     * @brief This sets the mNameCommonSection
     * @param inNameCommonSection
     */
    void setNameCommonSection( const std::string& inNameCommonSection )
    {
        if ( mData.find( inNameCommonSection ) == mData.end() )
        {
            throw std::invalid_argument(
                std::string( "Error: setNameCommonSecion()\n" ) +
                std::string( "There is no section named " ) +
                inNameCommonSection + std::string( "." ) );
        }
        mNameCommonSection = inNameCommonSection;
    }
    /**
     * @brief This sets the mNameCurrentSection
     * @param inNameCurrentSection
     */
    void setNameCurrentSection( const std::string& inNameCurrentSection )
    {
        if ( mData.find( inNameCurrentSection ) == mData.end() )
        {
            throw std::invalid_argument(
                std::string( "Error: setNameCurrentSecion()\n" ) +
                std::string( "There is no section named " ) +
                inNameCurrentSection + std::string( "." ) );
        }
        mNameCurrentSection = inNameCurrentSection;
    }
};

}  // namespace Utils

#endif