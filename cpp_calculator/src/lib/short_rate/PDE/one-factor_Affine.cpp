/**
 * @file one-factor_Affine.hpp
 * @brief This implements classes for one-factor Affine model calculation with
 * PDE.
 * @author kakune
 * @date 3/21/2024
 */

#include "short_rate/PDE/one-factor_Affine.hpp"

#include <memory>
#include <vector>

#include "math/ODE.hpp"

namespace ShortRate
{
namespace PDE
{

void AffineAbstract::buildAB()
{
    std::vector<double> lInitAB = { 1.0, 0.0 };
    auto lFunc                  = [this]( double inTime,
                         std::vector<double> inAB ) -> std::vector<double>
    { return this->ODEForAB( inTime, inAB ); };
    std::vector<double> lStartTimes = { mStartTime };
    std::vector<double> lEndTimes   = { mStartTime };
    std::vector<double> lAValues    = { 1.0 };
    std::vector<double> lBValues    = { 0.0 };
    for ( double lTmpTime = mStartTime; lTmpTime <= mEndTime + mStepTime;
          lTmpTime += mStepTime )
    {
        auto lResults = Math::ODE::solveSIMLRungeKutta45(
            lFunc, lTmpTime, lInitAB, mStartTime - mStepTime, 1e-6, mStepTime );

        lEndTimes.resize( lEndTimes.size() + lResults[0].size(), lTmpTime );
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
    std::vector<std::shared_ptr<const std::vector<double>>> lsTimes = {
        std::make_shared<const std::vector<double>>( lStartTimes ),
        std::make_shared<const std::vector<double>>( lEndTimes ) };
    mInterpA = std::make_unique<Math::InterpolateMulti::RBFGaussian>(
        lsTimes, std::make_shared<std::vector<double>>( lAValues ), 0.0,
        1.0e-6 * ( mEndTime - mStartTime ) );
    mInterpB = std::make_unique<Math::InterpolateMulti::RBFGaussian>(
        lsTimes, std::make_shared<std::vector<double>>( lBValues ), 0.0,
        1.0e-6 * ( mEndTime - mStartTime ) );
}

double AffineAbstract::priceZCB( double inMaturityTime,
                                 double inInitSpotRate ) const
{
    return priceZCBObservedAtArbitraryTime( mStartTime, inMaturityTime,
                                            inInitSpotRate );
}

double AffineAbstract::priceZCB( double inStartTime, double inMaturityTime,
                                 double inInitSpotRate ) const
{
    return priceZCB( inMaturityTime, inInitSpotRate ) /
           priceZCB( inStartTime, inInitSpotRate );
}
double AffineAbstract::priceZCBObservedAtArbitraryTime(
    double inObservingTime, double inMaturityTime, double inInitSpotRate ) const
{
    return mInterpA->operator()( { inObservingTime, inMaturityTime } ) *
           std::exp(
               -inInitSpotRate *
               mInterpB->operator()( { inObservingTime, inMaturityTime } ) );
}

double AffineAbstract::priceZCBObservedAtArbitraryTime(
    double inObservingTime, double inStartTime, double inMaturityTime,
    double inInitSpotRate ) const
{
    return priceZCBObservedAtArbitraryTime( inObservingTime, inMaturityTime,
                                            inInitSpotRate ) /
           priceZCBObservedAtArbitraryTime( inObservingTime, inStartTime,
                                            inInitSpotRate );
}

double AffineAbstract::instantaneousForwardRate( double inTime,
                                                 double inInitSpotRate ) const
{
    return instantaneousForwardRateAtArbitraryTime( mStartTime, inTime,
                                                    inInitSpotRate );
}

double AffineAbstract::instantaneousForwardRateAtArbitraryTime(
    double inObservingTime, double inTime, double inInitRate ) const
{
    return -mInterpA->deriv( { inObservingTime, inTime }, 1, 1 ) /
               mInterpA->operator()( { inObservingTime, inTime } ) +
           inInitRate * mInterpB->deriv( { inObservingTime, inTime }, 1, 1 );
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

}  // namespace PDE
}  // namespace ShortRate
