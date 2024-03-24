/**
 * @file Gauss.hpp
 * @brief This defines one-factor Gaussian short-rate models. (see Chap.10 of
 * Andersen & Piterbarg)
 * @author kakune
 * @date 3/8/2024
 */

#ifndef PROCESS_SHORT_RATE_MC_ONE_GAUSS_HPP
#define PROCESS_SHORT_RATE_MC_ONE_GAUSS_HPP

#include <iostream>
#include <memory>
#include <vector>

#include "process/market_data.hpp"
#include "process/short_rate_MC_one/core.hpp"

namespace Process
{
namespace ShortRateMCOne
{

/**
 * @brief This is HoLee short-rate model.
 */
class HoLee : public OneFactorAbstract
{
private:
    double mVol;  //! volatility of rate
    double driftCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots ) const override;
    double volCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots ) const override;

public:
    HoLee( std::size_t inNPath, const MarketData::Terms inTerms,
           double inInitSpotRate,
           std::unique_ptr<Process::Random::PathAbstract> inuRandomPath,
           double inVol ) :
        OneFactorAbstract( inNPath, inTerms, inInitSpotRate,
                           std::move( inuRandomPath ) ),
        mVol( inVol )
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
        return HoLee( mNPath, *muTerms, mInitSpotRate,
                      std::move( muRandomPath ), mVol );
    }
};

/**
 * @brief This is HoLee short-rate model quoting market data.
 */
class HoLeeWithMarket : public OneFactorAbstract
{
private:
    double mVol;   //! volatility of rate
    double mVol2;  //! square of mVol
    double driftCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots ) const override;
    double volCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots ) const override;
    MarketData::ZCB mMarketZCB;

public:
    HoLeeWithMarket(
        std::size_t inNPath, const MarketData::Terms inTerms,
        std::unique_ptr<Process::Random::PathAbstract> inuRandomPath,
        double inVol, const MarketData::ZCB& inMarketZCB ) :
        OneFactorAbstract( inNPath, inTerms, inMarketZCB.initialSpotRate(),
                           std::move( inuRandomPath ) ),
        mVol( inVol ),
        mVol2( inVol * inVol ),
        mMarketZCB( inMarketZCB )
    {
    }
};

class HoLeeWithMarketBuilder : public OneFactorAbstractBuilder
{
private:
    double mVol;
    std::unique_ptr<MarketData::ZCB> muMarketZCB;

public:
    HoLeeWithMarketBuilder& setVol( double inVol )
    {
        mVol = inVol;
        return *this;
    }
    HoLeeWithMarketBuilder& setMarketZCB( const MarketData::ZCB& inMarketZCB )
    {
        muMarketZCB = std::make_unique<MarketData::ZCB>( inMarketZCB );
        return *this;
    }
    HoLeeWithMarket build()
    {
        return HoLeeWithMarket( mNPath, *muTerms, std::move( muRandomPath ),
                                mVol, *muMarketZCB );
    }
};

/**
 * @brief This is Vasicek short-rate model quoting market data.
 */
class Vasicek : public OneFactorAbstract
{
private:
    double mVol, mKappa, mMean;
    double driftCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots ) const override;
    double volCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots ) const override;

public:
    Vasicek( std::size_t inNPath, const MarketData::Terms inTerms,
             double inInitSpotRate,
             std::unique_ptr<Process::Random::PathAbstract> inuRandomPath,
             double inVol, double inKappa, double inMean ) :
        OneFactorAbstract( inNPath, inTerms, inInitSpotRate,
                           std::move( inuRandomPath ) ),
        mVol( inVol ),
        mMean( inMean ),
        mKappa( inKappa )
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
        return Vasicek( mNPath, *muTerms, mInitSpotRate,
                        std::move( muRandomPath ), mVol, mKappa, mMean );
    }
};

/**
 * @brief This is Vasicek short-rate model.
 */
class VasicekWithMarket : public OneFactorAbstract
{
private:
    double mVol, mKappa;
    double mMeanFactor1, mMeanFactor2;
    double driftCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots ) const override;
    double volCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots ) const override;
    MarketData::ZCB mMarketZCB;

public:
    VasicekWithMarket(
        std::size_t inNPath, const MarketData::Terms inTerms,
        std::unique_ptr<Process::Random::PathAbstract> inuRandomPath,
        double inVol, double inKappa, const MarketData::ZCB& inMarketZCB ) :
        OneFactorAbstract( inNPath, inTerms, inMarketZCB.initialSpotRate(),
                           std::move( inuRandomPath ) ),
        mVol( inVol ),
        mKappa( inKappa ),
        mMeanFactor1( 1.0 / inKappa ),
        mMeanFactor2( inVol * inVol / ( 2.0 * inKappa * inKappa ) ),
        mMarketZCB( inMarketZCB )
    {
    }
};

class VasicekWithMarketBuilder : public OneFactorAbstractBuilder
{
private:
    double mVol, mKappa, mMean;
    std::unique_ptr<MarketData::ZCB> muMarketZCB;

public:
    VasicekWithMarketBuilder& setVol( double inVol )
    {
        mVol = inVol;
        return *this;
    }
    VasicekWithMarketBuilder& setKappa( double inKappa )
    {
        mKappa = inKappa;
        return *this;
    }
    VasicekWithMarketBuilder& setMarketZCB( const MarketData::ZCB& inMarketZCB )
    {
        muMarketZCB = std::make_unique<MarketData::ZCB>( inMarketZCB );
        return *this;
    }
    VasicekWithMarket build()
    {
        return VasicekWithMarket( mNPath, *muTerms, std::move( muRandomPath ),
                                  mVol, mKappa, *muMarketZCB );
    }
};

/**
 * @brief This is General One-Factor Gaussian Short Rate Model.
 */
class GSR : public OneFactorAbstract
{
private:
    Math::Interpolate1D::NewtonSpline mInterpVol, mInterpKappa, mInterpMean;
    double driftCoeff( std::size_t inIndPath, std::size_t inIndTerm,
                       const std::vector<std::vector<double>>& inSpots ) const;
    double volCoeff( std::size_t inIndPath, std::size_t inIndTerm,
                     const std::vector<std::vector<double>>& inSpots ) const;

public:
    GSR( std::size_t inNPath, const MarketData::Terms& inTerms,
         double inInitSpotRate,
         std::unique_ptr<Process::Random::PathAbstract> inuRandomPath,
         Math::Interpolate1D::NewtonSpline inInterpVol,
         Math::Interpolate1D::NewtonSpline inInterpKappa,
         Math::Interpolate1D::NewtonSpline inInterpMean ) :
        OneFactorAbstract( inNPath, inTerms, inInitSpotRate,
                           std::move( inuRandomPath ) ),
        mInterpVol( inInterpVol ),
        mInterpKappa( inInterpKappa ),
        mInterpMean( inInterpMean )
    {
    }
};

class GSRBuilder : public OneFactorAbstractBuilder
{
private:
    std::size_t mNDegVol, mNDegKappa, mNDegMean;
    std::unique_ptr<Math::Interpolate1D::NewtonSpline> muInterpVol,
        muInterpKappa, muInterpMean;

public:
    GSRBuilder( std::size_t inNDegVol = 3, std::size_t inNDegKappa = 3,
                std::size_t inNDegMean = 3 ) :
        mNDegVol( inNDegVol ),
        mNDegKappa( inNDegKappa ),
        mNDegMean( inNDegMean ),
        muInterpVol( nullptr ),
        muInterpKappa( nullptr ),
        muInterpMean( nullptr )
    {
    }
    GSRBuilder& setInterpVol( Math::Interpolate1D::NewtonSpline inInterpVol )
    {
        *muInterpVol = inInterpVol;
        return *this;
    }
    GSRBuilder& setInterpVol( const std::vector<double>& inVols )
    {
        muInterpVol = std::make_unique<Math::Interpolate1D::NewtonSpline>(
            muTerms->ptr(),
            std::make_shared<const std::vector<double>>( inVols ) );
        return *this;
    }
    GSRBuilder& setInterpVol( const std::vector<double>& inTerms,
                              const std::vector<double>& inVols )
    {
        muInterpVol = std::make_unique<Math::Interpolate1D::NewtonSpline>(
            std::make_shared<const std::vector<double>>( inTerms ),
            std::make_shared<const std::vector<double>>( inVols ) );
        return *this;
    }
    GSRBuilder& setInterpKappa(
        Math::Interpolate1D::NewtonSpline inInterpKappa )
    {
        *muInterpKappa = inInterpKappa;
        return *this;
    }
    GSRBuilder& setInterpKappa( const std::vector<double>& inKappas )
    {
        muInterpKappa = std::make_unique<Math::Interpolate1D::NewtonSpline>(
            muTerms->ptr(),
            std::make_shared<const std::vector<double>>( inKappas ) );
        return *this;
    }
    GSRBuilder& setInterpKappa( const std::vector<double>& inTerms,
                                const std::vector<double>& inKappas )
    {
        muInterpKappa = std::make_unique<Math::Interpolate1D::NewtonSpline>(
            std::make_shared<const std::vector<double>>( inTerms ),
            std::make_shared<const std::vector<double>>( inKappas ) );
        return *this;
    }
    GSRBuilder& setInterpMean( Math::Interpolate1D::NewtonSpline inInterpMean )
    {
        *muInterpMean = inInterpMean;
        return *this;
    }
    GSRBuilder& setInterpMean( const std::vector<double>& inMeans )
    {
        muInterpMean = std::make_unique<Math::Interpolate1D::NewtonSpline>(
            muTerms->ptr(),
            std::make_shared<const std::vector<double>>( inMeans ) );
        return *this;
    }
    GSRBuilder& setInterpMean( const std::vector<double>& inTerms,
                               const std::vector<double>& inMeans )
    {
        muInterpMean = std::make_unique<Math::Interpolate1D::NewtonSpline>(
            std::make_shared<const std::vector<double>>( inTerms ),
            std::make_shared<const std::vector<double>>( inMeans ) );
        return *this;
    }
    GSR build()
    {
        return GSR( mNPath, *muTerms, mInitSpotRate, std::move( muRandomPath ),
                    *muInterpVol, *muInterpKappa, *muInterpMean );
    }
};

/**
 * @brief This is General One-Factor Gaussian Short Rate Model.
 */
class GSRWithMarket : public OneFactorAbstract
{
private:
    Math::Interpolate1D::NewtonSpline mInterpVol, mInterpKappa;
    MarketData::ZCB mMarketZCB;
    double driftCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots ) const override
    {
        return 0.0;
    };
    double volCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots ) const override
    {
        return 0.0;
    };
    double driftCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots,
        const std::vector<std::vector<double>>& inFactors ) const;
    double volCoeff( std::size_t inIndPath, std::size_t inIndTerm,
                     const std::vector<std::vector<double>>& inSpots,
                     const std::vector<std::vector<double>>& inFactors ) const;
    double factorCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<double>>& inSpots,
        const std::vector<std::vector<double>>& inFactors ) const;

public:
    GSRWithMarket( std::size_t inNPath, const MarketData::Terms& inTerms,
                   std::unique_ptr<Process::Random::PathAbstract> inuRandomPath,
                   Math::Interpolate1D::NewtonSpline inInterpVol,
                   Math::Interpolate1D::NewtonSpline inInterpKappa,
                   MarketData::ZCB inMarketZCB ) :
        OneFactorAbstract( inNPath, inTerms, inMarketZCB.initialSpotRate(),
                           std::move( inuRandomPath ) ),
        mInterpVol( inInterpVol ),
        mInterpKappa( inInterpKappa ),
        mMarketZCB( inMarketZCB )
    {
    }
    MarketData::SpotRates calcSpotRates() const override;
};

class GSRWithMarketBuilder : public OneFactorAbstractBuilder
{
private:
    std::size_t mNDegVol, mNDegKappa;
    std::unique_ptr<Math::Interpolate1D::NewtonSpline> muInterpVol,
        muInterpKappa;
    std::unique_ptr<MarketData::ZCB> muMarketZCB;

public:
    GSRWithMarketBuilder( std::size_t inNDegVol   = 3,
                          std::size_t inNDegKappa = 3 ) :
        mNDegVol( inNDegVol ),
        mNDegKappa( inNDegKappa ),
        muInterpVol( nullptr ),
        muInterpKappa( nullptr )
    {
    }
    GSRWithMarketBuilder& setInterpVol(
        Math::Interpolate1D::NewtonSpline inInterpVol )
    {
        *muInterpVol = inInterpVol;
        return *this;
    }
    GSRWithMarketBuilder& setInterpVol( const std::vector<double>& inVols )
    {
        muInterpVol = std::make_unique<Math::Interpolate1D::NewtonSpline>(
            muTerms->ptr(),
            std::make_shared<const std::vector<double>>( inVols ) );
        return *this;
    }
    GSRWithMarketBuilder& setInterpVol( const std::vector<double>& inTerms,
                                        const std::vector<double>& inVols )
    {
        muInterpVol = std::make_unique<Math::Interpolate1D::NewtonSpline>(
            std::make_shared<const std::vector<double>>( inTerms ),
            std::make_shared<const std::vector<double>>( inVols ) );
        return *this;
    }
    GSRWithMarketBuilder& setInterpKappa(
        Math::Interpolate1D::NewtonSpline inInterpKappa )
    {
        *muInterpKappa = inInterpKappa;
        return *this;
    }
    GSRWithMarketBuilder& setInterpKappa( const std::vector<double>& inKappas )
    {
        muInterpKappa = std::make_unique<Math::Interpolate1D::NewtonSpline>(
            muTerms->ptr(),
            std::make_shared<const std::vector<double>>( inKappas ) );
        return *this;
    }
    GSRWithMarketBuilder& setInterpKappa( const std::vector<double>& inTerms,
                                          const std::vector<double>& inKappas )
    {
        muInterpKappa = std::make_unique<Math::Interpolate1D::NewtonSpline>(
            std::make_shared<const std::vector<double>>( inTerms ),
            std::make_shared<const std::vector<double>>( inKappas ) );
        return *this;
    }
    GSRWithMarketBuilder& setMarketZCB( const MarketData::ZCB& inMarketZCB )
    {
        muMarketZCB = std::make_unique<MarketData::ZCB>( inMarketZCB );
        return *this;
    }
    GSRWithMarket build()
    {
        return GSRWithMarket( mNPath, *muTerms, std::move( muRandomPath ),
                              *muInterpVol, *muInterpKappa, *muMarketZCB );
    }
};

}  // namespace ShortRateMCOne
}  // namespace Process

#endif