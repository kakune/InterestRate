/**
 * @file model_data.cpp
 * @brief This implements classes for model data.
 * @author kakune
 * @date 3/30/2024
 */

#include "process/model_data.hpp"

#include <cmath>
#include <memory>
#include <vector>

namespace Process
{
namespace ModelData
{

std::vector<std::vector<double>> SpotRates::calcDFFromSpotRate(
    const MarketData::Terms& inTerms,
    std::shared_ptr<const std::vector<std::vector<double>>> insDataSpotRate )
{
    std::size_t lNPath = insDataSpotRate->size();
    std::vector<std::vector<double>> lResult(
        lNPath, std::vector<double>( inTerms.size(), 1.0 ) );
    for ( std::size_t iTerm = 1; iTerm < inTerms.size(); ++iTerm )
    {
        double lHalfTmpDt = 0.5 * ( inTerms[iTerm] - inTerms[iTerm - 1] );
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lResult[iPath][iTerm] =
                lResult[iPath][iTerm - 1] *
                std::exp( -( insDataSpotRate->operator[]( iPath )[iTerm - 1] +
                             insDataSpotRate->operator[]( iPath )[iTerm] ) *
                          lHalfTmpDt );
        }
    }
    return lResult;
}

SpotRates::SpotRates(
    const MarketData::Terms& inTerms,
    std::shared_ptr<const std::vector<std::vector<double>>> insDataSpotRate ) :
    mTerms( inTerms ),
    mNPath( insDataSpotRate->size() ),
    msDataSpotRate( insDataSpotRate ),
    msDataDF( std::make_shared<const std::vector<std::vector<double>>>(
        calcDFFromSpotRate( inTerms, insDataSpotRate ) ) )
{
}
SpotRates::SpotRates( const MarketData::Terms& inTerms,
                      std::vector<std::vector<double>> inDataSpotRate ) :
    SpotRates( inTerms,
               std::make_shared<const std::vector<std::vector<double>>>(
                   inDataSpotRate ) )
{
}
const std::vector<double>& SpotRates::operator[]( std::size_t inIndex ) const
{
    return msDataSpotRate->operator[]( inIndex );
}
double SpotRates::term( std::size_t inIndex ) const { return mTerms[inIndex]; }
const MarketData::Terms& SpotRates::getTerms() const { return mTerms; }
const std::vector<std::vector<double>>& SpotRates::getDF() const
{
    return *msDataDF;
}
std::size_t SpotRates::sizeTerms() const { return mTerms.size(); }
std::size_t SpotRates::sizePath() const { return mNPath; }

MarketData::ZCB SpotRates::createZCB( std::size_t inDeg ) const
{
    std::vector<double> lZCB( mTerms.size(), 1.0 );
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
    {
        lZCB[iTerm] = 0.0;
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lZCB[iTerm] += msDataDF->operator[]( iPath )[iTerm];
        }
        lZCB[iTerm] /= mNPath;
    }
    return MarketData::ZCB( mTerms, lZCB, inDeg );
}

ForwardRates::ForwardRates(
    const MarketData::Terms& inTerms,
    const std::vector<std::size_t>& inIndTenor,
    std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
        insDataForwardRate ) :
    mTerms( inTerms ),
    mNPath( insDataForwardRate->size() ),
    mIndTenor( inIndTenor ),
    msDataForwardRate( insDataForwardRate )
{
}
ForwardRates::ForwardRates(
    const MarketData::Terms& inTerms,
    const std::vector<std::size_t>& inIndTenor,
    std::vector<std::vector<Math::Vec>> inDataForwardRate ) :
    ForwardRates( inTerms, inIndTenor,
                  std::make_shared<const std::vector<std::vector<Math::Vec>>>(
                      inDataForwardRate ) )
{
}
const std::vector<Math::Vec>& ForwardRates::operator[](
    std::size_t inIndex ) const
{
    return msDataForwardRate->operator[]( inIndex );
}
double ForwardRates::term( std::size_t inIndex ) const
{
    return mTerms[inIndex];
}
const MarketData::Terms& ForwardRates::getTerms() const { return mTerms; }
std::size_t ForwardRates::sizeTerms() const { return mTerms.size(); }
std::size_t ForwardRates::sizePath() const { return mNPath; }

}  // namespace ModelData
}  // namespace Process
