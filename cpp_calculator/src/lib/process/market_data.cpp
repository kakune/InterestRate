/**
 * @file interpolate_1d.cpp
 * @brief This implements classes for market data.
 * @author kakune
 * @date 3/21/2024
 */

#include "process/market_data.hpp"

#include <cmath>
#include <memory>
#include <vector>

namespace Process
{
namespace MarketData
{

Terms::Terms( std::shared_ptr<const std::vector<double>> insData ) :
    msData( insData )
{
}
Terms::Terms( std::vector<double> inData ) :
    Terms( std::make_shared<const std::vector<double>>( inData ) )
{
}
double Terms::operator[]( std::size_t inIndex ) const
{
    return msData->operator[]( inIndex );
}
double Terms::at( std::size_t inIndex ) const { return msData->at( inIndex ); }
std::size_t Terms::size() const { return msData->size(); }
const std::shared_ptr<const std::vector<double>> Terms::ptr() const
{
    return msData;
}
double Terms::front() const { return msData->front(); }
double Terms::back() const { return msData->back(); }

std::vector<std::vector<double>> SpotRates::calcDFFromSpotRate(
    const Terms& inTerms,
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
    const Terms& inTerms,
    std::shared_ptr<const std::vector<std::vector<double>>> insDataSpotRate ) :
    mTerms( inTerms ),
    mNPath( insDataSpotRate->size() ),
    msDataSpotRate( insDataSpotRate ),
    msDataDF( std::make_shared<const std::vector<std::vector<double>>>(
        calcDFFromSpotRate( inTerms, insDataSpotRate ) ) )
{
}
SpotRates::SpotRates( const Terms& inTerms,
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
const Terms& SpotRates::getTerms() const { return mTerms; }
const std::vector<std::vector<double>>& SpotRates::getDF() const
{
    return *msDataDF;
}
std::size_t SpotRates::sizeTerms() const { return mTerms.size(); }
std::size_t SpotRates::sizePath() const { return mNPath; }

std::vector<double> ZCB::calcZCBFromSpotRates(
    const SpotRates& inSpotRates ) const
{
    std::vector<double> lResults( inSpotRates.sizeTerms(), 1.0 );
    for ( std::size_t iTerm = 1; iTerm < inSpotRates.sizeTerms(); ++iTerm )
    {
        lResults[iTerm] = 0.0;
        for ( std::size_t iPath = 0; iPath < inSpotRates.sizePath(); ++iPath )
        {
            lResults[iTerm] += inSpotRates.getDF()[iPath][iTerm];
        }
        lResults[iTerm] /= double( inSpotRates.sizePath() );
    }
    return lResults;
}

ZCB::ZCB( const Terms& inTerms,
          std::shared_ptr<const std::vector<double>> insData,
          std::size_t inDeg ) :
    mTerms( inTerms ),
    msData( insData ),
    msSpline( std::make_shared<Math::Interpolate1D::NewtonSpline>(
        inTerms.ptr(), insData, inDeg ) )
{
}
ZCB::ZCB( const Terms& inTerms, std::vector<double> inData,
          std::size_t inDeg ) :
    ZCB( inTerms, std::make_shared<const std::vector<double>>( inData ), inDeg )
{
}
ZCB::ZCB( const SpotRates& inSpotRate, std::size_t inDeg ) :
    ZCB( inSpotRate.getTerms(), calcZCBFromSpotRates( inSpotRate ), inDeg )
{
}
double ZCB::operator[]( std::size_t inIndex ) const
{
    return msData->operator[]( inIndex );
}
double ZCB::term( std::size_t inIndex ) const { return mTerms[inIndex]; }
const Terms& ZCB::getTerms() const { return mTerms; }
std::size_t ZCB::sizeTerms() const { return mTerms.size(); }

double ZCB::operator()( double inTime ) const
{
    return msSpline->operator()( inTime );
}

double ZCB::operator()( double inStartTime, double inMaturityTime ) const
{
    return operator()( inMaturityTime ) / operator()( inStartTime );
}

double ZCB::forwardRate( double inStartTime, double inTerminalTime ) const
{
    return ( std::log( operator()( inStartTime ) ) -
             std::log( operator()( inTerminalTime ) ) ) /
           ( inTerminalTime - inStartTime );
}

double ZCB::instantaneousForwardRate( double inTime ) const
{
    return -msSpline->deriv( inTime, 1 ) / operator()( inTime );
}
double ZCB::derivInstantaneousForwardRate( double inTime ) const
{
    double lInvP = 1.0 / operator()( inTime );
    double lP1   = msSpline->deriv( inTime, 1 );
    double lP2   = msSpline->deriv( inTime, 2 );
    return ( -lP2 + lP1 * lP1 * lInvP ) * lInvP;
}
double ZCB::initialSpotRate() const
{
    return instantaneousForwardRate( 0.95 * mTerms[0] + 0.05 * mTerms[1] );
}

}  // namespace MarketData
}  // namespace Process
