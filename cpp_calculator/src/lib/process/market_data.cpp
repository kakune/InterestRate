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

#include "math/findroot_1d.hpp"

namespace Process
{
namespace MarketData
{

static std::vector<double> calcDifTime(
    std::shared_ptr<const std::vector<double>> insTime )
{
    std::vector<double> lResult( insTime->size(), 0.0 );
    for ( std::size_t i = 1; i < insTime->size(); ++i )
    {
        lResult[i] = ( *insTime )[i] - ( *insTime )[i - 1];
    }
    return lResult;
}
static std::vector<double> calcSqrtDifTime(
    std::shared_ptr<const std::vector<double>> insTime )
{
    std::vector<double> lResult( insTime->size(), 0.0 );
    for ( std::size_t i = 1; i < insTime->size(); ++i )
    {
        lResult[i] = std::sqrt( ( *insTime )[i] - ( *insTime )[i - 1] );
    }
    return lResult;
}
Terms::Terms( std::shared_ptr<const std::vector<double>> insData ) :
    msData( insData ),
    msDifTime(
        std::make_shared<const std::vector<double>>( calcDifTime( insData ) ) ),
    msSqrtDifTime( std::make_shared<const std::vector<double>>(
        calcSqrtDifTime( insData ) ) )
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
double Terms::difTime( std::size_t inIndex ) const
{
    return ( *msDifTime )[inIndex];
}
double Terms::sqrtDifTime( std::size_t inIndex ) const
{
    return ( *msSqrtDifTime )[inIndex];
}
std::size_t Terms::size() const { return msData->size(); }
const std::shared_ptr<const std::vector<double>> Terms::ptr() const
{
    return msData;
}
double Terms::front() const { return msData->front(); }
double Terms::back() const { return msData->back(); }

static Math::Vec calcTau(
    const Terms& inTerms,
    std::shared_ptr<const std::vector<std::size_t>> insData )
{
    Math::Vec lResult( insData->size() - 1 );
    for ( std::size_t i = 0; i < insData->size() - 1; ++i )
    {
        lResult( i ) = inTerms[( *insData )[i + 1]] - inTerms[( *insData )[i]];
    }
    return lResult;
}
static std::vector<std::size_t> calcMinIndAtEachTime(
    const Terms& inTerms,
    std::shared_ptr<const std::vector<std::size_t>> insData )
{
    std::vector<std::size_t> lResult( inTerms.size() );
    std::size_t iTmp = 0;
    for ( std::size_t iTerm = 0; iTerm < inTerms.size(); ++iTerm )
    {
        while ( iTmp < insData->size() && ( *insData )[iTmp] < iTerm )
        {
            ++iTmp;
        }
        lResult[iTerm] = iTmp;
    }
    return lResult;
}

Tenor::Tenor( const Terms& inTerms,
              const std::shared_ptr<const std::vector<std::size_t>> insData ) :
    mNSize( insData->size() - 1 ),
    mTerms( inTerms ),
    msData( insData ),
    msTau( std::make_shared<const Math::Vec>( calcTau( inTerms, insData ) ) ),
    msMinIndAtEachTime( std::make_shared<const std::vector<std::size_t>>(
        calcMinIndAtEachTime( inTerms, insData ) ) )
{
}

Tenor::Tenor( const Terms& inTerms, const std::vector<std::size_t>& inData ) :
    Tenor( inTerms, std::make_shared<const std::vector<std::size_t>>( inData ) )
{
}

std::size_t Tenor::operator[]( std::size_t inIndex ) const
{
    return ( *msData )[inIndex];
}
double Tenor::term( std::size_t inIndex ) const
{
    return mTerms[( *msData )[inIndex]];
}
double Tenor::tau( std::size_t inIndex ) const { return ( *msTau )( inIndex ); }
const Math::Vec& Tenor::getTauVec() const { return *msTau; }
std::size_t Tenor::minIndex( std::size_t inIndex ) const
{
    return ( *msMinIndAtEachTime )[inIndex];
}
std::size_t Tenor::size() const { return mNSize; }

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

static std::vector<double> calcZCBVec( const Terms& inTerms, const ZCB& inZCB )
{
    std::vector<double> lResults( inTerms.size() - 1 );
    for ( std::size_t i = 0; i < lResults.size(); ++i )
    {
        lResults[i] = inZCB( inTerms[i + 1] );
    }
    return lResults;
}
Caplets::Caplets(
    const Terms& inTerms, std::shared_ptr<const std::vector<double>> insStrikes,
    std::shared_ptr<const std::vector<std::vector<double>>> insCaplets,
    const ZCB& inZCB ) :
    mTerms( inTerms ),
    msStrikes( insStrikes ),
    msCaplets( insCaplets ),
    msBlack( std::make_shared<Analytical::Black76::Model>( *inTerms.ptr() ) )
{
    ( *msBlack ).setInitZCB( calcZCBVec( inTerms, inZCB ) );
}

Caplets::Caplets( const Terms& inTerms, const std::vector<double>& inStrikes,
                  const std::vector<std::vector<double>>& inCaplets,
                  const ZCB& inZCB ) :
    Caplets(
        inTerms, std::make_shared<const std::vector<double>>( inStrikes ),
        std::make_shared<const std::vector<std::vector<double>>>( inCaplets ),
        inZCB )
{
}
double Caplets::strike( std::size_t inIndex ) const
{
    return ( *msStrikes )[inIndex];
}
const std::vector<double>& Caplets::operator[]( std::size_t inIndex ) const
{
    return ( *msCaplets )[inIndex];
}

double Caplets::impliedBlackVol( std::size_t inIndexStrike,
                                 std::size_t inIndexTerm ) const
{
    double lStrike      = ( *msStrikes )[inIndexStrike];
    double lPriceCaplet = ( *msCaplets )[inIndexStrike][inIndexTerm];
    auto lFuncDif       = [this, inIndexTerm, lStrike,
                     lPriceCaplet]( double inVol ) -> double
    {
        ( *msBlack ).setVol( inVol );
        return lPriceCaplet - ( msBlack->priceCaplet )( lStrike, inIndexTerm );
    };
    return Math::FindRoot1D::Brent( lFuncDif, 1e-10, 1e2 );
}

}  // namespace MarketData
}  // namespace Process
