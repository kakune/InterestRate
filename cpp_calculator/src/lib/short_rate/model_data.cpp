/**
 * @file model_data.cpp
 * @brief This implements classes for model data.
 * @author kakune
 * @date 4/21/2024
 */

#include "short_rate/model_data.hpp"

namespace ShortRate
{

std::vector<std::vector<double>> calcDFFromSpotRates(
    const Process::MarketData::Terms& inTerms,
    std::shared_ptr<const std::vector<std::vector<double>>> insDataSpotRate )
{
    std::size_t lNPath = insDataSpotRate->size();
    std::vector<std::vector<double>> lResult(
        lNPath, std::vector<double>( inTerms.size(), 1.0 ) );
    for ( std::size_t iTerm = 1; iTerm < inTerms.size(); ++iTerm )
    {
        double lHalfTmpDt = 0.5 * ( inTerms[iTerm] - inTerms[iTerm - 1] );
        for ( std::size_t iPath = 0; iPath < insDataSpotRate->size(); ++iPath )
        {
            lResult[iPath][iTerm] =
                lResult[iPath][iTerm - 1] *
                std::exp( -( ( *insDataSpotRate )[iPath][iTerm - 1] +
                             ( *insDataSpotRate )[iPath][iTerm] ) *
                          lHalfTmpDt );
        }
    }
    return lResult;
}

Process::MarketData::ZCB createZCBFromSpotRates(
    const Process::MarketData::Terms& inTerms,
    const std::vector<std::vector<double>>& inDFs, std::size_t inDeg )
{
    std::vector<double> lZCB( inTerms.size(), 1.0 );
    for ( std::size_t iTerm = 1; iTerm < inTerms.size(); ++iTerm )
    {
        lZCB[iTerm] = 0.0;
        for ( std::size_t iPath = 0; iPath < inDFs.size(); ++iPath )
        {
            lZCB[iTerm] += inDFs[iPath][iTerm];
        }
        lZCB[iTerm] /= inDFs.size();
    }
    return Process::MarketData::ZCB( inTerms, lZCB, inDeg );
}

SpotRates::SpotRates(
    const Process::MarketData::Terms& inTerms,
    std::shared_ptr<const std::vector<std::vector<double>>> insDataSpotRate,
    std::size_t inDegZCB ) :
    mTerms( inTerms ),
    mNPath( insDataSpotRate->size() ),
    msDataSpotRate( insDataSpotRate ),
    msDataDF( std::make_shared<const std::vector<std::vector<double>>>(
        calcDFFromSpotRates( inTerms, insDataSpotRate ) ) ),
    msZCB( std::make_shared<const Process::MarketData::ZCB>(
        createZCBFromSpotRates( inTerms, *msDataDF, inDegZCB ) ) )
{
}
SpotRates::SpotRates( const Process::MarketData::Terms& inTerms,
                      std::vector<std::vector<double>> inDataSpotRate,
                      std::size_t inDegZCB ) :
    SpotRates( inTerms,
               std::make_shared<const std::vector<std::vector<double>>>(
                   inDataSpotRate ),
               inDegZCB )
{
}
const std::vector<double>& SpotRates::operator[]( std::size_t inIndex ) const
{
    return msDataSpotRate->operator[]( inIndex );
}
double SpotRates::term( std::size_t inIndex ) const { return mTerms[inIndex]; }
const Process::MarketData::Terms& SpotRates::getTerms() const { return mTerms; }
const std::vector<std::vector<double>>& SpotRates::getDF() const
{
    return *msDataDF;
}
std::size_t SpotRates::sizeTerms() const { return mTerms.size(); }
std::size_t SpotRates::sizePath() const { return mNPath; }
const Process::MarketData::ZCB& SpotRates::getZCB() const { return ( *msZCB ); }

double SpotRates::calcMeanMMA( std::size_t inIndStart,
                               std::size_t inIndEnd ) const
{
    double lResult = 0.0;
    for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
    {
        double lTmp = 1.0;
        for ( std::size_t iTerm = inIndStart + 1; iTerm <= inIndEnd; ++iTerm )
        {
            lTmp *= ( *msDataDF )[iPath][iTerm];
        }
        lResult += 1.0 / lTmp;
    }
    return lResult / mNPath;
}

double SpotRates::calcCaplet( double inStrike, std::size_t inIndStart,
                              std::size_t inIndEnd ) const
{
    return std::max( 0.0, calcMeanMMA( inIndStart, inIndEnd ) -
                              std::exp( inStrike * ( mTerms[inIndEnd] -
                                                     mTerms[inIndStart] ) ) );
}
double SpotRates::calcFloorlet( double inStrike, std::size_t inIndStart,
                                std::size_t inIndEnd ) const
{
    return std::max(
        0.0, std::exp( inStrike * ( mTerms[inIndEnd] - mTerms[inIndStart] ) -
                       calcMeanMMA( inIndStart, inIndEnd ) ) );
}

double SpotRates::calcCap( double inStrike,
                           std::vector<std::size_t> inTenor ) const
{
    double lResult = 0;
    for ( std::size_t iIndTenor = 0; iIndTenor < inTenor.size() - 1;
          ++iIndTenor )
    {
        lResult +=
            calcCaplet( inStrike, inTenor[iIndTenor], inTenor[iIndTenor + 1] );
    }
    return lResult;
}
double SpotRates::calcFloor( double inStrike,
                             std::vector<std::size_t> inTenor ) const
{
    double lResult = 0;
    for ( std::size_t iIndTenor = 0; iIndTenor < inTenor.size() - 1;
          ++iIndTenor )
    {
        lResult += calcFloorlet( inStrike, inTenor[iIndTenor],
                                 inTenor[iIndTenor + 1] );
    }
    return lResult;
}

}  // namespace ShortRate