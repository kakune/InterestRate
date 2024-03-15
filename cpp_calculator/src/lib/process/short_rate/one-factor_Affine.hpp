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

/**
 * @brief This is affine short-rate model with constant coefficients.
 */
class ConstantAffine : public OneFactorAbstract
{
private:
    double mLambda, mEta;   //! drift coefficients ( lambda x + eta )dt
    double mGamma, mDelta;  //! vol coefficients ( gamma x + delta )dW
    double driftCoeff( std::size_t inIndPath,
                       std::size_t inIndTerm ) const override;
    double volCoeff( std::size_t inIndPath,
                     std::size_t inIndTerm ) const override;
    Math::InterpolateMulti::RBFGaussian mInterpA, mInterpB;

public:
    ConstantAffine(
        std::size_t inNPath,
        std::shared_ptr<const std::vector<double>> insTerms,
        std::shared_ptr<const Market::Data> insMarketData,
        double inInitSpotRate,
        std::unique_ptr<Process::Random::PathAbstract> inuRandomPath,
        double inLambda, double inEta, double inGamma, double inDelta ) :
        OneFactorAbstract( inNPath, insTerms, insMarketData, inInitSpotRate,
                           std::move( inuRandomPath ) ),
        mLambda( inLambda ),
        mEta( inEta ),
        mGamma( inGamma ),
        mDelta( inDelta ),
        mInterpA( 0.0, 1.0e-6 * ( insTerms->back() - insTerms->front() ) ),
        mInterpB( 0.0, 1.0e-6 * ( insTerms->back() - insTerms->front() ) )
    {
    }
    void buildAB( double inTimeMesh = 10.0 );
    double priceZCBByAB( double inMaturityTime ) const;
    // double analyticalPriceZCB( std::size_t inIndStartTime,
    //                            std::size_t inIndMaturityTime ) const;
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