/**
 * @file one-factor_Affine.cpp
 * @brief This implements one-factor Affine short-rate models. (see 3.2.4 of
 * Brigo)
 * @author kakune
 * @date 3/14/2024
 */

#include "process/short_rate/one-factor_Affine.hpp"

#include <iostream>
#include <limits>
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

void AffineAbstract::buildAB( double inTimeMesh )
{
    std::vector<double> lInitAB = { 1.0, 0.0 };
    auto lFunc                  = [this]( double inTime,
                         std::vector<double> inAB ) -> std::vector<double>
    { return this->ODEForAB( inTime, inAB ); };
    double lMaxDif = ( msTerms->back() - msTerms->front() ) / inTimeMesh;
    std::vector<double> lStartTimes = { 0.0 };
    std::vector<double> lEndTimes   = { 0.0 };
    std::vector<double> lAValues    = { 1.0 };
    std::vector<double> lBValues    = { 0.0 };
    for ( std::size_t iTerm = 1; iTerm < msTerms->size(); ++iTerm )
    {
        auto lResults = Math::ODE::solveSIMLRungeKutta45(
            lFunc, msTerms->operator[]( iTerm ), lInitAB,
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

double AffineAbstract::priceZCBByAB( double inMaturityTime ) const
{
    return mInterpA( { msTerms->front(), inMaturityTime } ) *
           std::exp( -mInitSpotRate *
                     mInterpB( { msTerms->front(), inMaturityTime } ) );
}

double AffineAbstract::priceZCBByAB( double inStartTime,
                                     double inMaturityTime ) const
{
    return priceZCBByAB( inMaturityTime ) / priceZCBByAB( inStartTime );
}
double AffineAbstract::priceZCBByAB( std::size_t inIndStartTime,
                                     std::size_t inIndMaturityTime ) const
{
    if ( inIndStartTime >= msTerms->size() ||
         inIndMaturityTime >= msTerms->size() )
    {
        std::cerr << "Error: Process::ShortRate::ConstantAffine::priceZCBByAB()"
                  << std::endl
                  << "Argument is out of range." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    return priceZCBByAB( msTerms->operator[]( inIndStartTime ),
                         msTerms->operator[]( inIndMaturityTime ) );
}

double AffineAbstract::instantaneousForwardRateByAB( double inTime ) const
{
    return instantaneousForwardRateByAB( inTime, msTerms->front(),
                                         mInitSpotRate );
}

double AffineAbstract::instantaneousForwardRateByAB(
    std::size_t inIndTime ) const
{
    if ( inIndTime >= msTerms->size() )
    {
        std::cerr << "Error: Process::ShortRate::ConstantAffine::priceZCBByAB()"
                  << std::endl
                  << "Argument is out of range." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    return instantaneousForwardRateByAB( msTerms->operator[]( inIndTime ) );
}

double AffineAbstract::instantaneousForwardRateByAB( double inTime,
                                                     double inObservingTime,
                                                     double inInitRate ) const
{
    if ( inObservingTime > inTime )
    {
        std::cerr << "Error: "
                     "Process::ShortRate::ConstantAffine::"
                     "instantaneousForwardRateByAB()"
                  << std::endl
                  << "ObservingTime must be smaller than inTime." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    return -mInterpA.deriv( { inObservingTime, inTime }, 1, 1 ) /
               mInterpA( { inObservingTime, inTime } ) +
           inInitRate * mInterpB.deriv( { inObservingTime, inTime }, 1, 1 );
}

double AffineAbstract::instantaneousForwardRateByAB(
    std::size_t inIndTime, std::size_t inIndObservingTime,
    double inInitRate ) const
{
    if ( inIndTime >= msTerms->size() || inIndObservingTime >= msTerms->size() )
    {
        std::cerr << "Error: Process::ShortRate::ConstantAffine::priceZCBByAB()"
                  << std::endl
                  << "Argument is out of range." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    return instantaneousForwardRateByAB(
        msTerms->operator[]( inIndTime ),
        msTerms->operator[]( inIndObservingTime ), inInitRate );
}

std::vector<double> ConstantAffine::ODEForAB( double inTime,
                                              std::vector<double> inAB ) const
{
    std::vector<double> lResult( 2 );
    double B2  = inAB[1] * inAB[1];
    lResult[1] = -this->mLambda * inAB[1] + 0.5 * this->mGamma * B2 - 1.0;
    lResult[0] = inAB[0] * ( this->mEta * inAB[1] - 0.5 * this->mDelta * B2 );
    return lResult;
}

}  // namespace ShortRate
}  // namespace Process
