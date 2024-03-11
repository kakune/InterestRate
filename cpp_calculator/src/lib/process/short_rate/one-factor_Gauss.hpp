/**
 * @file one-factor_Gauss.hpp
 * @brief This defines one-factor Gaussian short-rate models. (see Chap.10 of
 * Andersen & Piterbarg)
 * @author kakune
 * @date 3/8/2024
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
 * @brief This is HoLee short-rate model.
 */
class HoLee : public OneFactorAbstract
{
private:
    double mVol;   //! volatility of rate
    double mVol2;  //! square of mVol
    double driftCoeff( std::size_t inIndPath,
                       std::size_t inIndTerm ) const override;
    double volCoeff( std::size_t inIndPath,
                     std::size_t inIndTerm ) const override;

public:
    HoLee( std::size_t inNPath,
           std::shared_ptr<const std::vector<double>> insTerms,
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

/**
 * @brief This is Vasicek short-rate model.
 */
class Vasicek : public OneFactorAbstract
{
private:
    double mVol, mKappa, mMean;
    double mMeanFactor1, mMeanFactor2;
    double driftCoeff( std::size_t inIndPath,
                       std::size_t inIndTerm ) const override;
    double volCoeff( std::size_t inIndPath,
                     std::size_t inIndTerm ) const override;

public:
    Vasicek( std::size_t inNPath,
             std::shared_ptr<const std::vector<double>> insTerms,
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
    double analyticalPriceZCB( double inMaturityTime ) const;
    double analyticalPriceZCB( double inStartTime,
                               double inMaturityTime ) const;
    double analyticalPriceZCB( std::size_t inIndStartTime,
                               std::size_t inIndMaturityTime ) const;
};

class VasicekBuilder : public OneFactorAbstractBuilder
{
private:
    double mVol, mKappa, mMean;

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

/**
 * @brief This is General One-Factor Gaussian Short Rate Model.
 */
class GSR : public OneFactorAbstract
{
private:
    Math::Interpolate1d::NewtonSpline mInterpVol, mInterpKappa, mInterpMean;
    std::vector<std::vector<double>> mFactors;
    double driftCoeff( std::size_t inIndPath,
                       std::size_t inIndTerm ) const override;
    double volCoeff( std::size_t inIndPath,
                     std::size_t inIndTerm ) const override;
    double factorCoeff( std::size_t inIndPath, std::size_t inIndTerm ) const;
    void buildSpot() override;

public:
    GSR( std::size_t inNPath,
         std::shared_ptr<const std::vector<double>> insTerms,
         std::shared_ptr<const Market::Data> insMarketData,
         double inInitSpotRate,
         std::unique_ptr<Process::Random::PathAbstract> inuRandomPath,
         Math::Interpolate1d::NewtonSpline inInterpVol,
         Math::Interpolate1d::NewtonSpline inInterpKappa,
         Math::Interpolate1d::NewtonSpline inInterpMean ) :
        OneFactorAbstract( inNPath, insTerms, insMarketData,
                           insMarketData == nullptr
                               ? inInitSpotRate
                               : insMarketData->mInterpInstantaneousForwardRate(
                                     insTerms->at( 0 ) ),
                           std::move( inuRandomPath ) ),
        mInterpVol( inInterpVol ),
        mInterpKappa( inInterpKappa ),
        mInterpMean( inInterpMean )
    {
        if ( insMarketData != nullptr )
        {
            mFactors.resize( mNPath,
                             std::vector<double>( msTerms->size(), 0.0 ) );
        }
    }
};

class GSRBuilder : public OneFactorAbstractBuilder
{
private:
    Math::Interpolate1d::NewtonSpline mInterpVol, mInterpKappa, mInterpMean;

public:
    GSRBuilder( std::size_t inNDegVol = 3, std::size_t inNDegKappa = 3,
                std::size_t inNDegMean = 3 ) :
        mInterpVol( inNDegVol ),
        mInterpKappa( inNDegKappa ),
        mInterpMean( inNDegMean )
    {
    }
    GSRBuilder& setInterpVol( Math::Interpolate1d::NewtonSpline inInterpVol )
    {
        mInterpVol = inInterpVol;
        return *this;
    }
    GSRBuilder& setInterpVol( std::vector<double> inVols )
    {
        mInterpVol.build(
            msTerms, std::make_shared<const std::vector<double>>( inVols ) );
        return *this;
    }
    GSRBuilder& setInterpVol( std::vector<double> inTerms,
                              std::vector<double> inVols )
    {
        mInterpVol.build(
            std::make_shared<const std::vector<double>>( inTerms ),
            std::make_shared<const std::vector<double>>( inVols ) );
        return *this;
    }
    GSRBuilder& setInterpKappa(
        Math::Interpolate1d::NewtonSpline inInterpKappa )
    {
        mInterpKappa = inInterpKappa;
        return *this;
    }
    GSRBuilder& setInterpKappa( std::vector<double> inKappas )
    {
        mInterpKappa.build(
            msTerms, std::make_shared<const std::vector<double>>( inKappas ) );
        return *this;
    }
    GSRBuilder& setInterpKappa( std::vector<double> inTerms,
                                std::vector<double> inKappas )
    {
        mInterpKappa.build(
            std::make_shared<const std::vector<double>>( inTerms ),
            std::make_shared<const std::vector<double>>( inKappas ) );
        return *this;
    }
    GSRBuilder& setInterpMean( Math::Interpolate1d::NewtonSpline inInterpMean )
    {
        mInterpMean = inInterpMean;
        return *this;
    }
    GSRBuilder& setInterpMean( std::vector<double> inMeans )
    {
        mInterpMean.build(
            msTerms, std::make_shared<const std::vector<double>>( inMeans ) );
        return *this;
    }
    GSRBuilder& setInterpMean( std::vector<double> inTerms,
                               std::vector<double> inMeans )
    {
        mInterpMean.build(
            std::make_shared<const std::vector<double>>( inTerms ),
            std::make_shared<const std::vector<double>>( inMeans ) );
        return *this;
    }
    GSR build()
    {
        return GSR( mNPath, msTerms, msMarketData, mInitSpotRate,
                    std::move( muRandomPath ), mInterpVol, mInterpKappa,
                    mInterpMean );
    }
};

}  // namespace ShortRate
}  // namespace Process

#endif