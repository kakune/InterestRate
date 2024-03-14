/**
 * @file one-factor_Gauss.hpp
 * @brief This defines one-factor Affine short-rate models. (see 3.2.4 of
 * Brigo)
 * @author kakune
 * @date 3/12/2024
 */

#ifndef PROCESS_SHORT_RATE_ONE_FACTOR_GAUSS_HPP
#define PROCESS_SHORT_RATE_ONE_FACTOR_GAUSS_HPP

#include <iostream>
#include <memory>
#include <vector>

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

public:
    ConstantAffine( std::size_t inNPath,
           std::shared_ptr<const std::vector<double>> insTerms,
           std::shared_ptr<const Market::Data> insMarketData,
           double inInitSpotRate,
           std::unique_ptr<Process::Random::PathAbstract> inuRandomPath,
           double inLambda, double inEta, double inGamma, double inDelta ) :
        OneFactorAbstract( inNPath, insTerms, insMarketData, inInitSpotRate,
                           std::move( inuRandomPath ) ),
        mVol( inVol ),
        mVol2( inVol * inVol )
    {
    }
    double analyticalPriceZCB( double inStartTime,
                               double inMaturityTime ) const;
    double analyticalPriceZCB( std::size_t inIndStartTime,
                               std::size_t inIndMaturityTime ) const;
};

class HoLeeBuilder : public OneFactorAbstractBuilder
{
private:
    double mVol;

public:
    HoLeeBuilder& setVol( double inVol )
    {
        mVol = inVol;
        return *this;
    }
    HoLee build()
    {
        return HoLee( mNPath, msTerms, msMarketData, mInitSpotRate,
                      std::move( muRandomPath ), mVol );
    }
};


}  // namespace ShortRate
}  // namespace Process

#endif