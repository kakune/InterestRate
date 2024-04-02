/**
 * @file Affine.hpp
 * @brief This defines one-factor Affine short-rate models. (see 3.2.4 of
 * Brigo)
 * @author kakune
 * @date 3/14/2024
 */

#ifndef SHORT_RATE_ONE_FACTOR_AFFINE_HPP
#define SHORT_RATE_ONE_FACTOR_AFFINE_HPP

#include <memory>
#include <vector>

#include "short_rate/one-factor/core.hpp"

namespace ShortRate
{
namespace OneFactor
{

/**
 * @brief This is affine short-rate model with constant coefficients.
 */
class ConstantAffine : public OneFactorAbstract
{
private:
    const double mLambda, mEta;   //! drift coefficients ( lambda x + eta )dt
    const double mGamma, mDelta;  //! vol coefficients ( gamma x + delta )dW

    double driftCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots ) const override;
    double volCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots ) const override;

public:
    ConstantAffine(
        std::size_t inNPath, const Process::MarketData::Terms& inTerms,
        double inInitSpotRate,
        std::unique_ptr<Process::Random::StdBrownAbstract> inuRandomPath,
        double inLambda, double inEta, double inGamma, double inDelta ) :
        OneFactorAbstract( inNPath, inTerms, inInitSpotRate,
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
        return ConstantAffine( mNPath, *muTerms, mInitSpotRate,
                               std::move( muStdBrown ), mLambda, mEta, mGamma,
                               mDelta );
    }
};

}  // namespace OneFactor
}  // namespace ShortRate

#endif