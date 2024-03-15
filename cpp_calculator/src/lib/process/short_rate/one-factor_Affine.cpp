/**
 * @file one-factor_Affine.cpp
 * @brief This implements one-factor Affine short-rate models. (see 3.2.4 of
 * Brigo)
 * @author kakune
 * @date 3/14/2024
 */

#include "process/short_rate/one-factor_Affine.hpp"

#include <iostream>
#include <memory>
#include <vector>

#include "math/ODE.hpp"
#include "process/market.hpp"

namespace Process
{
namespace ShortRate
{

double ConstantAffine::driftCoeff( std::size_t inIndPath,
                                   std::size_t inIndTerm ) const
{
    return mLambda * mSpotRates.at( inIndPath ).at( inIndTerm - 1 ) + mEta;
}
double ConstantAffine::volCoeff( std::size_t inIndPath,
                                 std::size_t inIndTerm ) const
{
    return std::sqrt( mGamma * mSpotRates.at( inIndPath ).at( inIndTerm - 1 ) +
                      mDelta );
}

void ConstantAffine::buildAB( double inTimeMesh )
{
    std::vector<double> lInitAB = { 1.0, 0.0 };
    auto ODEForAB               = [this]( double inTime,
                            std::vector<double> inAB ) -> std::vector<double>
    {
        std::vector<double> lResult( 2 );
        double B2  = inAB[1] * inAB[1];
        lResult[1] = -this->mLambda * inAB[1] + 0.5 * this->mGamma * B2 - 1.0;
        lResult[0] =
            inAB[0] * ( this->mEta * inAB[1] - 0.5 * this->mDelta * B2 );
        return lResult;
    };
    double lMaxDif = ( msTerms->back() - msTerms->front() ) / inTimeMesh;
    std::vector<double> lStartTimes = { 0.0 };
    std::vector<double> lEndTimes   = { 0.0 };
    std::vector<double> lAValues    = { 1.0 };
    std::vector<double> lBValues    = { 0.0 };
    for ( std::size_t iTerm = 1; iTerm < msTerms->size(); ++iTerm )
    {
        auto lResults = Math::ODE::solveSIMLRungeKutta45(
            ODEForAB, msTerms->operator[]( iTerm ), lInitAB,
            msTerms->front() - lMaxDif, 1e-6, lMaxDif );

        lEndTimes.resize( lEndTimes.size() + lResults[0].size(),
                          msTerms->operator[]( iTerm ) );
        lStartTimes.insert( lStartTimes.end(),
                            std::make_move_iterator( lResults[0].begin() ),
                            std::make_move_iterator( lResults[0].end() ) );
        lAValues.insert( lAValues.end(),
                         std::make_move_iterator( lResults[1].begin() ),
                         std::make_move_iterator( lResults[1].end() ) );
        lBValues.insert( lBValues.end(),
                         std::make_move_iterator( lResults[2].begin() ),
                         std::make_move_iterator( lResults[2].end() ) );
    }
    std::vector<std::shared_ptr<const std::vector<double> > > lsTimes = {
        std::make_shared<const std::vector<double> >( lStartTimes ),
        std::make_shared<const std::vector<double> >( lEndTimes ) };
    mInterpA.build( lsTimes,
                    std::make_shared<std::vector<double> >( lAValues ) );
    mInterpB.build( lsTimes,
                    std::make_shared<std::vector<double> >( lBValues ) );
}

double ConstantAffine::priceZCBByAB( double inMaturityTime ) const
{
    return mInterpA( { msTerms->front(), inMaturityTime } ) *
           std::exp( -mInitSpotRate *
                     mInterpB( { msTerms->front(), inMaturityTime } ) );
}

}  // namespace ShortRate
}  // namespace Process
