/**
 * @file payoff.cpp
 * @brief This implements payoff function objects.
 * @author kakune
 * @date 4/21/2024
 */
#include "LIBOR/forward/payoff.hpp"

#include <iostream>

namespace LIBOR::Forward::Payoff
{

CapletFloorlet::CapletFloorlet(
    const Process::MarketData::Tenor& inTenor,
    std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
        insDataForwardRates,
    std::size_t inIndTenorPay, std::size_t inIndTenorFR, double inStrike,
    bool inIsCaplet ) :
    msDataForwardRates( insDataForwardRates ),
    mIndTenorPay( inIndTenorPay ),
    mIndTenorFR( inIndTenorFR ),
    mIndTermsFR( inTenor[inIndTenorFR] ),
    mTau( inTenor.tau( inIndTenorFR ) ),
    mStrike( inStrike ),
    mIsCaplet( inIsCaplet )
{
}
std::size_t CapletFloorlet::getIndexTenorPay() const { return mIndTenorPay; }
double CapletFloorlet::operator()( std::size_t inIndPath ) const
{
    if ( mIsCaplet )
    {
        return std::max(
            0.0,
            mTau *
                ( ( *msDataForwardRates )[inIndPath][mIndTermsFR][mIndTenorFR] -
                  mStrike ) );
    }
    return std::max(
        0.0,
        mTau *
            ( mStrike -
              ( *msDataForwardRates )[inIndPath][mIndTermsFR][mIndTenorFR] ) );
}

static std::vector<std::size_t> transformIndsFromTenorToTerms(
    const Process::MarketData::Tenor& inTenor,
    const std::vector<std::size_t>& inIndsTenorPay )
{
    std::vector<std::size_t> lResult( inIndsTenorPay.size() );
    for ( std::size_t i = 0; i < lResult.size(); ++i )
    {
        lResult[i] = inTenor[inIndsTenorPay[i]];
    }
    return lResult;
}

static Math::Vec calcTaus( const Process::MarketData::Tenor& inTenor,
                           std::size_t inIndTenorStart,
                           std::size_t inIndTenorLast )
{
    Math::Vec lResult = Math::makeVec( inIndTenorLast - inIndTenorStart );
    for ( std::size_t i = 0; i < lResult.size(); ++i )
    {
        lResult[i] = inTenor.tau( inIndTenorStart + i );
    }
    return lResult;
}

Swaption::Swaption( const Process::MarketData::Tenor& inTenor,
                    std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
                        insDataForwardRates,
                    std::size_t inIndTenorStart, std::size_t inIndTenorLast,
                    double inStrike, bool inIsPayer ) :
    mTenor( inTenor ),
    msDataForwardRates( insDataForwardRates ),
    mIndTenorStart( inIndTenorStart ),
    mIndTermsStart( inTenor[inIndTenorStart] ),
    mIndTenorLast( inIndTenorLast ),
    mIndTermsLast( inTenor[inIndTenorLast] ),
    mTaus( calcTaus( inTenor, inIndTenorStart, inIndTenorLast ) ),
    mStrike( inStrike ),
    mIsPayer( inIsPayer )
{
}
std::size_t Swaption::getIndexTenorFirstPay() const
{
    return mIndTenorStart + 1;
}
std::size_t Swaption::getIndexTenorLastPay() const { return mIndTenorLast; }

static Math::Vec calcZCBFromForwardRate(
    const Math::Vec& inForwardRate, const Process::MarketData::Tenor& inTenor,
    std::size_t inIndTenorStart, std::size_t inIndTenorLast )
{
    Math::Vec lZCB = Math::makeVec( inIndTenorLast - inIndTenorStart );
    for ( std::size_t i = 0; i < lZCB.size(); ++i )
    {
        lZCB[i] = ( i == 0 ) ? 1.0 : lZCB[i - 1];
        lZCB[i] /= 1.0 + inTenor.tau( inIndTenorStart + i ) *
                             inForwardRate[inIndTenorStart + i];
    }
    return lZCB;
}

static Math::Vec calcRealizedTauForwardRate(
    const std::vector<Math::Vec>& inForwardRates,
    const Process::MarketData::Tenor& inTenor, std::size_t inIndTenorStart,
    std::size_t inIndTenorLast )
{
    Math::Vec lForwardRate = Math::makeVec( inIndTenorLast - inIndTenorStart );
    for ( std::size_t i = 0; i < lForwardRate.size(); ++i )
    {
        lForwardRate[i] = 1.0 + inTenor.tau( inIndTenorStart + i ) *
                                    inForwardRates[inTenor[inIndTenorStart + i]]
                                                  [inIndTenorStart + i];
    }
    return lForwardRate - 1.0;
}

Math::Vec Swaption::operator()( std::size_t inIndPath ) const
{
    Math::Vec lZCB = calcZCBFromForwardRate(
        ( *msDataForwardRates )[inIndPath][mIndTermsStart], mTenor,
        mIndTenorStart, mIndTenorLast );
    double lSwapRate =
        ( 1.0 - lZCB[lZCB.size() - 1] ) / Math::dot( lZCB, mTaus );
    Math::Vec lResult = Math::makeVec( lZCB.size(), 0.0 );
    if ( mIsPayer )
    {
        if ( lSwapRate <= mStrike ) { return lResult; }
        lResult =
            calcRealizedTauForwardRate( ( *msDataForwardRates )[inIndPath],
                                        mTenor, mIndTenorStart, mIndTenorLast );
        return lResult - mTaus * mStrike;
    }
    if ( lSwapRate >= mStrike ) { return lResult; }
    lResult =
        calcRealizedTauForwardRate( ( *msDataForwardRates )[inIndPath], mTenor,
                                    mIndTenorStart, mIndTenorLast );
    return mTaus * mStrike - lResult;
}

}  // namespace LIBOR::Forward::Payoff
