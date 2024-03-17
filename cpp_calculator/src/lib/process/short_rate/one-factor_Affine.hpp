/**
 * @file one-factor_Affine.hpp
 * @brief This defines one-factor Affine short-rate models. (see 3.2.4 of
 * Brigo)
 * @author kakune
 * @date 3/14/2024
 */

#ifndef PROCESS_SHORT_RATE_ONE_FACTOR_AFFINE_HPP
#define PROCESS_SHORT_RATE_ONE_FACTOR_AFFINE_HPP

#include <iostream>
#include <memory>
#include <vector>

#include "math/interpolate_multi.hpp"
#include "process/market.hpp"
#include "process/short_rate/core.hpp"

namespace Process
{
namespace ShortRate
{

class AffineAbstract : public OneFactorAbstract
{
protected:
    Math::InterpolateMulti::RBFGaussian mInterpA,
        mInterpB;  //! Affine A,B defined as $P = Ae^{-rB}$
    virtual std::vector<double> ODEForAB( double inTime,
                                          std::vector<double> inAB ) const = 0;

public:
    AffineAbstract(
        std::size_t inNPath,
        std::shared_ptr<const std::vector<double>> insTerms,
        std::shared_ptr<const Market::Data> insMarketData,
        double inInitSpotRate,
        std::unique_ptr<Process::Random::PathAbstract> inuRandomPath ) :
        OneFactorAbstract( inNPath, insTerms, insMarketData, inInitSpotRate,
                           std::move( inuRandomPath ) ),
        mInterpA( 0.0, 1.0e-6 * ( insTerms->back() - insTerms->front() ) ),
        mInterpB( 0.0, 1.0e-6 * ( insTerms->back() - insTerms->front() ) )
    {
    }
    /**
     * @brief This calculates mInterpA and mInterpB by ODE.
     * @param inTimeMesh step of ODE is bounded by $T / inTimeMesh$.
     */
    void buildAB( double inTimeMesh = 10.0 );
    /**
     * @brief This calculates ZCB starting at
     * Terms[0] using mInterpA and mInterpB.
     * @param inMaturityTime maturity time of ZCB
     * @return double P(Terms[0], inMaturityTime)
     */
    double priceZCBByAB( double inMaturityTime ) const;
    /**
     * @brief This calculates ZCB price within arbitrary interval observed at
     * Terms[0] using mInterpA and mInterpB.
     * @param inStartTime start time of ZCB
     * @param inMaturityTime maturity time of ZCB
     * @return double P(inStartTime, inMaturityTime)
     */
    double priceZCBByAB( double inStartTime, double inMaturityTime ) const;
    double priceZCBByAB( std::size_t inIndStartTime,
                         std::size_t inIndMaturityTime ) const;
    /**
     * @brief This calculates instantaneous forward rate at arbitary time
     * observed at Terms[0] using mInterpA and mInterpB.
     * @param inTime time of forward rate
     * @return double f(inTime)
     */
    double instantaneousForwardRateByAB( double inTime ) const;
    double instantaneousForwardRateByAB( std::size_t inIndTime ) const;

    /**
     * @brief This calculates instantaneous forward rate at arbitary time
     * observed at Terms[0] using mInterpA and mInterpB.
     * @param inTime time of forward rate
     * @param inObservingTime time observing forward rate
     * @param inInitRate short rate at inObservingTime
     * @return double f(inTime)
     */
    double instantaneousForwardRateByAB( double inTime, double inObservingTime,
                                         double inInitRate ) const;
    double instantaneousForwardRateByAB( std::size_t inIndTime,
                                         std::size_t inIndObservingTime,
                                         double inInitRate ) const;
};

/**
 * @brief This is affine short-rate model with constant coefficients.
 */
class ConstantAffine : public AffineAbstract
{
private:
    double mLambda, mEta;   //! drift coefficients ( lambda x + eta )dt
    double mGamma, mDelta;  //! vol coefficients ( gamma x + delta )dW
    double driftCoeff( std::size_t inIndPath,
                       std::size_t inIndTerm ) const override;
    double volCoeff( std::size_t inIndPath,
                     std::size_t inIndTerm ) const override;
    virtual std::vector<double> ODEForAB(
        double inTime, std::vector<double> inAB ) const override;

public:
    ConstantAffine(
        std::size_t inNPath,
        std::shared_ptr<const std::vector<double>> insTerms,
        std::shared_ptr<const Market::Data> insMarketData,
        double inInitSpotRate,
        std::unique_ptr<Process::Random::PathAbstract> inuRandomPath,
        double inLambda, double inEta, double inGamma, double inDelta ) :
        AffineAbstract( inNPath, insTerms, insMarketData, inInitSpotRate,
                        std::move( inuRandomPath ) ),
        mLambda( inLambda ),
        mEta( inEta ),
        mGamma( inGamma ),
        mDelta( inDelta )
    {
    }
};

class ConstantAffineBuilder : public OneFactorAbstractBuilder
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
        return ConstantAffine( mNPath, msTerms, msMarketData, mInitSpotRate,
                               std::move( muRandomPath ), mLambda, mEta, mGamma,
                               mDelta );
    }
};

}  // namespace ShortRate
}  // namespace Process

#endif