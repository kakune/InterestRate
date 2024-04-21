/**
 * @file payoff.cpp
 * @brief This implements payoff function objects.
 * @author kakune
 * @date 4/21/2024
 */

#include "LIBOR/forward/payoff.hpp"

namespace LIBOR::Forward
{

CapletFloorletPayoff::CapletFloorletPayoff(
    const Process::MarketData::Tenor& inTenor,
    std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
        insDataForwardRates,
    double inStrike, std::size_t inIndTenorFR, std::size_t inIndTenorPay,
    bool inIsCaplet ) :
    msDataForwardRates( insDataForwardRates ),
    mStrike( inStrike ),
    mIndTermsFR( inTenor[inIndTenorFR] ),
    mIndTenorFR( inIndTenorFR ),
    mIndTenorPay( inIndTenorPay ),
    mTau( inTenor.tau( inIndTenorFR ) ),
    mIsCaplet( inIsCaplet )
{
}
std::size_t CapletFloorletPayoff::getIndexTenorPay() const
{
    return mIndTenorPay;
}
double CapletFloorletPayoff::operator()( std::size_t inIndPath ) const
{
    if ( mIsCaplet )
    {
        return std::max(
            0.0, mTau * ( ( *msDataForwardRates )[inIndPath][mIndTermsFR](
                              mIndTenorFR ) -
                          mStrike ) );
    }
    return std::max(
        0.0, mTau * ( mStrike - ( *msDataForwardRates )[inIndPath][mIndTermsFR](
                                    mIndTenorFR ) ) );
}

}  // namespace LIBOR::Forward
