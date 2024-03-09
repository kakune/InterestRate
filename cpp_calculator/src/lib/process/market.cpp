/**
 * @file market.hpp
 * @brief This implements the construction of ZCB and Forward rate in the market
 * data
 * @author kakune
 * @date 3/8/2024
 */

#include "process/market.hpp"

#include <cmath>

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
    Math::Interpolate1d::NewtonSpline lSplineRates( mNDimSpline );
    lSplineRates.build( std::make_shared<std::vector<double>>( mFineTerms ),
                        std::make_shared<std::vector<double>>( lFineRates ) );
    lSplineRates.buildIntegral();

    std::vector<double> lZCB( mFineTerms.size() );
    for ( std::size_t iTerm = 0; iTerm < mFineTerms.size(); ++iTerm )
    {
        lZCB.at( iTerm ) =
            std::exp( -lSplineRates.integral( mFineTerms.at( iTerm ) ) );
    }
    mInterpZCB.build( std::make_shared<std::vector<double>>( mFineTerms ),
                      std::make_shared<std::vector<double>>( lZCB ) );
}

}  // namespace Market
}  // namespace Process