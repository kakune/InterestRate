/**
 * @file market_data.cpp
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
std::vector<double> Terms::calcDifTime(
    std::shared_ptr<const std::vector<double>> insTime )
{
    std::vector<double> lResult( insTime->size(), 0.0 );
    for ( std::size_t i = 1; i < insTime->size(); ++i )
    {
        lResult[i] = insTime->operator[]( i ) - insTime->operator[]( i - 1 );
    }
    return lResult;
}
Terms::Terms( std::shared_ptr<const std::vector<double>> insData ) :
    msData( insData ), mDifTime( calcDifTime( insData ) )
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
double Terms::difTime( std::size_t inIndex ) const { return mDifTime[inIndex]; }
std::size_t Terms::size() const { return msData->size(); }
const std::shared_ptr<const std::vector<double>> Terms::ptr() const
{
    return msData;
}
double Terms::front() const { return msData->front(); }
double Terms::back() const { return msData->back(); }

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
    return ( operator()( inStartTime ) / operator()( inTerminalTime ) - 1.0 ) /
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
