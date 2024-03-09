/**
 * @file market.hpp
 * @brief This implements the construction of ZCB and Forward rate in the market
 * data
 * @author kakune
 * @date 3/8/2024
 */

#include "process/market.hpp"

#include <cmath>

#include "math/integral_1d.hpp"

namespace Process
{
namespace Market
{

void Data::setZCB( std::vector<double> inZCB )
{
    mInterpZCB.build( msTerms, std::make_shared<std::vector<double>>( inZCB ) );
    std::vector<double> lInstantaneousForwardRate( mFineTerms.size() );
    for ( std::size_t iTerm = 0; iTerm < mFineTerms.size(); ++iTerm )
    {
        lInstantaneousForwardRate.at( iTerm ) =
            -mInterpZCB.deriv( mFineTerms.at( iTerm ), 1 ) /
            mInterpZCB( mFineTerms.at( iTerm ) );
    }
    mInterpInstantaneousForwardRate.build(
        std::make_shared<std::vector<double>>( mFineTerms ),
        std::make_shared<std::vector<double>>( lInstantaneousForwardRate ) );
}

void Data::setForwardRate( std::vector<double> inForwardRate )
{
    mInterpInstantaneousForwardRate.build(
        msTerms, std::make_shared<std::vector<double>>( inForwardRate ) );

    std::vector<double> lFineRates( mFineTerms.size() );
    for ( std::size_t iTerm = 0; iTerm < mFineTerms.size(); ++iTerm )
    {
        lFineRates.at( iTerm ) =
            mInterpInstantaneousForwardRate( mFineTerms.at( iTerm ) );
    }
    std::vector<double> lIntegralRates =
        Math::Integral::eachTrapezoidal( mFineTerms, lFineRates );
    for ( std::size_t iTerm = 0; iTerm < mFineTerms.size(); ++iTerm )
    {
        lIntegralRates.at( iTerm ) = std::exp( -lIntegralRates.at( iTerm ) );
    }
    mInterpZCB.build( std::make_shared<std::vector<double>>( mFineTerms ),
                      std::make_shared<std::vector<double>>( lIntegralRates ) );
}

}  // namespace Market
}  // namespace Process