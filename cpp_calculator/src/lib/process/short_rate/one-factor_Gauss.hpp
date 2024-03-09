/**
 * @file one-factor_Gauss.hpp
 * @brief This defines one-factor Gaussian short-rate models. (see Chap.10 of
 * Andersen & Piterbarg)
 * @author kakune
 * @date 3/8/2024
 */

#ifndef PROCESS_SHORT_RATE_ONE_FACTOR_GAUSS_HPP
#define PROCESS_SHORT_RATE_ONE_FACTOR_GAUSS_HPP

#include <memory>

#include "process/market.hpp"
#include "process/short_rate/core.hpp"

namespace Process
{
namespace ShortRate
{

class HoLee : public OneFactorAbstract
{
private:
    double mVol, mVol2;
    double driftCoeff( std::size_t inIndPath,
                       std::size_t inIndTerm ) const override;
    double volCoeff( std::size_t inIndPath,
                     std::size_t inIndTerm ) const override;

public:
    HoLee( std::size_t inNPath,
           std::shared_ptr<const std::vector<double> > insTerms,
           std::shared_ptr<const Market::Data> insMarketData,
           double inInitSpotRate,
           std::unique_ptr<Process::Random::PathAbstract> inuRandomPath,
           double inVol ) :
        OneFactorAbstract( inNPath, insTerms, insMarketData, inInitSpotRate,
                           std::move( inuRandomPath ) ),
        mVol( inVol ),
        mVol2( inVol * inVol )
    {
    }
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

class Vasicek : public OneFactorAbstract
{
private:
    double mKappa, mMean, mVol;
    double mMeanFactor1, mMeanFactor2;
    double driftCoeff( std::size_t inIndPath,
                       std::size_t inIndTerm ) const override;
    double volCoeff( std::size_t inIndPath,
                     std::size_t inIndTerm ) const override;

public:
    Vasicek( std::size_t inNPath,
             std::shared_ptr<const std::vector<double> > insTerms,
             std::shared_ptr<const Market::Data> insMarketData,
             double inInitSpotRate,
             std::unique_ptr<Process::Random::PathAbstract> inuRandomPath,
             double inVol, double inKappa, double inMean ) :
        OneFactorAbstract( inNPath, insTerms, insMarketData,
                           insMarketData == nullptr
                               ? inInitSpotRate
                               : insMarketData->mInterpInstantaneousForwardRate(
                                     insTerms->at( 0 ) ),
                           std::move( inuRandomPath ) ),
        mVol( inVol ),
        mMean( inMean ),
        mKappa( inKappa ),
        mMeanFactor1( 1.0 / inKappa ),
        mMeanFactor2( inVol * inVol / ( 2.0 * inKappa * inKappa ) )
    {
    }
};

class VasicekBuilder : public OneFactorAbstractBuilder
{
private:
    double mVol, mMean, mKappa;

public:
    VasicekBuilder& setVol( double inVol )
    {
        mVol = inVol;
        return *this;
    }
    VasicekBuilder& setMean( double inMean )
    {
        mMean = inMean;
        return *this;
    }
    VasicekBuilder& setKappa( double inKappa )
    {
        mKappa = inKappa;
        return *this;
    }
    Vasicek build()
    {
        return Vasicek( mNPath, msTerms, msMarketData, mInitSpotRate,
                        std::move( muRandomPath ), mVol, mKappa, mMean );
    }
};

}  // namespace ShortRate
}  // namespace Process

#endif