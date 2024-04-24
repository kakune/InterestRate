/**
 * @file one-factor_Affine.hpp
 * @brief This defines classes for one-factor Affine model calculation with PDE.
 * @author kakune
 * @date 3/21/2024
 */

#ifndef SHORT_RATE_PDE_ONE_FACTOR_AFFINE_HPP
#define SHORT_RATE_PDE_ONE_FACTOR_AFFINE_HPP

#include <memory>
#include <vector>

#include "math/interpolate_multi.hpp"
#include "short_rate/PDE/core.hpp"

namespace ShortRate
{
namespace PDE
{

class AffineAbstract : public ModelAbstract
{
protected:
    std::unique_ptr<const Math::InterpolateMulti::RBFGaussian> mInterpA,
        mInterpB;  //! Affine A,B defined as $P = Ae^{-rB}$

    virtual std::vector<double> ODEForAB( double inTime,
                                          std::vector<double> inAB ) const = 0;
    /**
     * @brief This calculates mInterpA and mInterpB by ODE.
     */
    void buildAB();

public:
    AffineAbstract( double inStartTime, double inEndTime,
                    double inMinStepTime ) :
        ModelAbstract( inStartTime, inEndTime, inMinStepTime ),
        mInterpA( nullptr ),
        mInterpB( nullptr )
    {
    }

    /**
     * @brief This calculates ZCB starting at
     * mStartTime using mInterpA and mInterpB.
     * @param inMaturityTime maturity time of ZCB
     * @param inInitSpotRate spot rate at mStartTime
     * @return double P(mStartTime, inMaturityTime)
     */
    double priceZCB( double inMaturityTime, double inInitSpotRate ) const;
    /**
     * @brief This calculates ZCB price within arbitrary interval observed at
     * mStartTime using mInterpA and mInterpB.
     * @param inStartTime start time of ZCB
     * @param inMaturityTime maturity time of ZCB
     * @param inInitSpotRate spot rate at mStartTime
     * @return double P(inStartTime, inMaturityTime)
     */
    double priceZCB( double inStartTime, double inMaturityTime,
                     double inInitSpotRate ) const;
    /**
     * @brief This calculates ZCB starting at
     * inObservingTime using mInterpA and mInterpB.
     * @param inObservingTime time observing spot rate
     * @param inMaturityTime maturity time of ZCB
     * @param inInitSpotRate spot rate at inObservingTime
     * @return double P(inObserveTime, inMaturityTime)
     */
    double priceZCBObservedAtArbitraryTime( double inObservingTime,
                                            double inMaturityTime,
                                            double inInitSpotRate ) const;
    /**
     * @brief This calculates ZCB starting at
     * inObservingTime using mInterpA and mInterpB.
     * @param inObservingTime time observing spot rate
     * @param inStartTime start time of ZCB
     * @param inMaturityTime maturity time of ZCB
     * @param inInitSpotRate spot rate at inObservingTime
     * @return double P(inObserveTime, inMaturityTime)
     */
    double priceZCBObservedAtArbitraryTime( double inObservingTime,
                                            double inStartTime,
                                            double inMaturityTime,
                                            double inInitSpotRate ) const;
    /**
     * @brief This calculates instantaneous forward rate at arbitary time
     * observed at Terms[0] using mInterpA and mInterpB.
     * @param inTime time of forward rate
     * @param inInitSpotRate spot rate at mStartTime
     * @return double f(inTime)
     */
    double instantaneousForwardRate( double inTime,
                                     double inInitSpotRate ) const;

    /**
     * @brief This calculates instantaneous forward rate at arbitary time
     * observed at Terms[0] using mInterpA and mInterpB.
     * @param inObservingTime time observing forward rate
     * @param inTime time of forward rate
     * @param inInitRate short rate at inObservingTime
     * @return double f(inTime)
     */
    double instantaneousForwardRateAtArbitraryTime(
        double inObservingTime, double inTime, double inInitSpotRate ) const;
};

/**
 * @brief This is affine short-rate model with constant coefficients.
 */
class ConstantAffine : public AffineAbstract
{
private:
    double mLambda, mEta;   //! drift coefficients ( lambda x + eta )dt
    double mGamma, mDelta;  //! vol coefficients ( gamma x + delta )dW
    virtual std::vector<double> ODEForAB(
        double inTime, std::vector<double> inAB ) const override;

public:
    ConstantAffine( double inStartTime, double inEndTime, double inStepTime,
                    double inLambda, double inEta, double inGamma,
                    double inDelta ) :
        AffineAbstract( inStartTime, inEndTime, inStepTime ),
        mLambda( inLambda ),
        mEta( inEta ),
        mGamma( inGamma ),
        mDelta( inDelta )
    {
        buildAB();
    }
};

class ConstantAffineBuilder : public ModelAbstractBuilder
{
private:
    double mLambda, mEta;   //! drift coefficients ( lambda x + eta )dt
    double mGamma, mDelta;  //! vol coefficients ( gamma x + delta )dW

public:
    ConstantAffineBuilder& setDrift( double inLambda, double inEta )
    {
        mLambda = inLambda;
        mEta    = inEta;
        return *this;
    }
    ConstantAffineBuilder& setVol( double inGamma, double inDelta )
    {
        mGamma = inGamma;
        mDelta = inDelta;
        return *this;
    }
    ConstantAffine build()
    {
        return ConstantAffine( mStartTime, mEndTime, mStepTime, mLambda, mEta,
                               mGamma, mDelta );
    }
};

}  // namespace PDE
}  // namespace ShortRate

#endif