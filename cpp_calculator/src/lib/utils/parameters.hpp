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
    std::string mCurrentSectionName;
    //! name of section that used if current section has no corresponding
    //! parameters
    std::string mCommonSectionName = "DEFAULT";

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
     * @brief This accesses the parameter using mCurrentSectionName and
     * mCommonSectionName.
     * @details Return parameters in mCurrentSectionName.
     * If there is no parameter corresponding to that section,
     * return the parameter in mCommonSectionName.
     * @param inParameterName
     * @return double
     */
    double operator()( const std::string& inParameterName )
    {
        if ( mData[mCurrentSectionName].find( inParameterName ) ==
             mData[mCurrentSectionName].end() )
        {
            return mData[mCommonSectionName][inParameterName];
        }
        return mData[mCurrentSectionName][inParameterName];
    }

    /**
     * @brief This accesses the parameter using mCurrentSectionName and
     * mCommonSectionName.
     * @details Return parameters in mCurrentSectionName.
     * If there is no parameter corresponding to that section,
     * return the parameter in mCommonSectionName.
     * @param inParameterName
     * @return const double
     */
    const double operator()( const std::string& inParameterName ) const
    {
        if ( mData.at( mCurrentSectionName ).find( inParameterName ) ==
             mData.at( mCurrentSectionName ).end() )
        {
            return mData.at( mCommonSectionName ).at( inParameterName );
        }
        return mData.at( mCurrentSectionName ).at( inParameterName );
    }

    /**
     * @brief This sets the mCommonSectionName
     * @param inCommonSectionName
     */
    void setCommonSectionName( const std::string& inCommonSectionName )
    {
        mCommonSectionName = inCommonSectionName;
    }
    /**
     * @brief This sets the mCurrentSectionName
     * @param inCommonSectionName
     */
    void setCurrentSectionName( const std::string& inCurrentSectionName )
    {
        mCurrentSectionName = inCurrentSectionName;
    }
};

}  // namespace Utils

#endif