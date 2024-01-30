/**
 * @file parameters.hpp
 * @brief define struct "Utils::Parameters" that organize numerical parameters
 * @author kakune
 * @date 1/28/2024
 */
#ifndef UTILS_PARAMETERS_HPP
#define UTILS_PARAMETERS_HPP

#include <map>
#include <string>

namespace Utils
{

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
    std::map< std::string, std::map< std::string, double > > mData;
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
    const std::map< std::string, double >& operator[](
        const std::string& inSectionName ) const
    {
        return mData.at( inSectionName );
    }

    /**
     * @brief This accesses the map of each section.
     * @param inSectionName Name of section
     * @return std::map< std::string, double >&
     */
    std::map< std::string, double >& operator[](
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
    double operator()( const std::string& inParameterName )
    {
        if ( mData[mNameCurrentSection].find( inParameterName ) ==
             mData[mNameCurrentSection].end() )
        {
            return mData[mNameCommonSection][inParameterName];
        }
        return mData[mNameCurrentSection][inParameterName];
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
    const double operator()( const std::string& inParameterName ) const
    {
        if ( mData.at( mNameCurrentSection ).find( inParameterName ) ==
             mData.at( mNameCurrentSection ).end() )
        {
            return mData.at( mNameCommonSection ).at( inParameterName );
        }
        return mData.at( mNameCurrentSection ).at( inParameterName );
    }

    /**
     * @brief This sets the mNameCommonSection
     * @param inNameCommonSection
     */
    void setNameCommonSection( const std::string& inNameCommonSection )
    {
        mNameCommonSection = inNameCommonSection;
    }
    /**
     * @brief This sets the mNameCurrentSection
     * @param inNameCurrentSection
     */
    void setNameCurrentSection( const std::string& inNameCurrentSection )
    {
        mNameCurrentSection = inNameCurrentSection;
    }
};

}  // namespace Utils

#endif