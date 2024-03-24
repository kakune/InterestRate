/**
 * @file Gauss.hpp
 * @brief This defines multi-factor Gaussian short-rate models. (see Chap.4 of
 * Brigo)
 * @author kakune
 * @date 3/24/2024
 */

#ifndef PROCESS_SHORT_RATE_MC_MULTI_GAUSS_HPP
#define PROCESS_SHORT_RATE_MC_MULTI_GAUSS_HPP

#include <iostream>
#include <memory>
#include <vector>

#include "process/market_data.hpp"
#include "process/short_rate_MC_multi/core.hpp"

namespace Process
{
namespace ShortRateMCMulti
{

/**
 * @brief This is multi-factor short-rate Gaussian model with constant coeff
 * without shift $\varphi$.
 */
class ConstantGauss : public MultiFactorAbstract
{
private:
    Math::Vec mDriftCoeff;  //! coefficient of dt
    Math::Mat
        mVolCoeff;  //! coefficient of dW which must be LowerTriangular Matrix

    // DriftCoeff * State * dt
    Math::Vec driftCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<Math::Vec>>& inSpots ) const override;

    // VolCoeff.dW
    Math::Vec volTerm( std::size_t inIndPath, std::size_t inIndTerm,
                       const std::vector<std::vector<Math::Vec>>& inSpots,
                       const Math::Vec& inRandomVec ) const override;

    // sum of components of state
    double transfStateToRate( const Math::Vec& inState,
                              double inTime ) const override;

public:
    ConstantGauss(
        std::size_t inNPath, const MarketData::Terms& inTerms,
        const Math::Vec& inInitState,
        std::unique_ptr<Process::RandomVec::PathAbstract> inuRandomPath,
        const Math::Vec& inDriftCoeff, const Math::Mat& inVolCoeff ) :
        MultiFactorAbstract( inNPath, inTerms, inInitState,
                             std::move( inuRandomPath ) ),
        mDriftCoeff( inDriftCoeff ),
        mVolCoeff( inVolCoeff )
    {
    }
    double analyticalPriceZCB( double inMaturityTime ) const;
};

class ConstantGaussBuilder : public MultiFactorAbstractBuilder
{
private:
    Math::Vec mDriftCoeff = Math::Vec( 0 );
    Math::Mat mVolCoeff   = Math::Mat( 0, 0 );

public:
    ConstantGaussBuilder& setDrift( const Math::Vec& inDriftCoeff )
    {
        mDriftCoeff = inDriftCoeff;
        return *this;
    }
    ConstantGaussBuilder& setVol( const Math::Mat& inVolCoeff )
    {
        mVolCoeff = inVolCoeff;
        return *this;
    }
    ConstantGauss build()
    {
        return ConstantGauss( mNPath, *muTerms, mInitState,
                              std::move( muRandomPath ), mDriftCoeff,
                              mVolCoeff );
    }
};

/**
 * @brief This is G2++ model with constant coeff and market data.
 */
class G2ppWithMarket : public MultiFactorAbstract
{
private:
    Math::Vec mDriftCoeff;  //! coefficient of dt
    Math::Mat
        mVolCoeff;  //! coefficient of dW which must be LowerTriangular Matrix
    MarketData::ZCB mMarketZCB;  //! the market data

    double mA, mB;
    double mFactorA, mFactorB, mFactorAB;  //! factor for transfStateToRate

    // (-a*x, -b*y)dt
    Math::Vec driftCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<Math::Vec>>& inSpots ) const override;
    // VolCoeff.(dW_1, dW_2)
    Math::Vec volTerm( std::size_t inIndPath, std::size_t inIndTerm,
                       const std::vector<std::vector<Math::Vec>>& inSpots,
                       const Math::Vec& inRandomVec ) const override;
    // r = x + y + \phi
    double transfStateToRate( const Math::Vec& inState,
                              double inTime ) const override;

public:
    G2ppWithMarket(
        std::size_t inNPath, const MarketData::Terms& inTerms,
        std::unique_ptr<Process::RandomVec::PathAbstract> inuRandomPath,
        const Math::Vec& inDriftCoeff, const Math::Mat& inVolCoeff,
        const MarketData::ZCB& inMarketZCB ) :
        MultiFactorAbstract( inNPath, inTerms, Math::Vec( 2, 0.0 ),
                             std::move( inuRandomPath ) ),
        mDriftCoeff( inDriftCoeff ),
        mVolCoeff( inVolCoeff ),
        mMarketZCB( inMarketZCB )
    {
        if ( inDriftCoeff.size() != 2 || inVolCoeff.sizeRow() != 2 ||
             inVolCoeff.sizeCol() != 2 )
        {
            std::cerr
                << "Error : "
                   "Process::ShortRateMCMulti::G2ppWithMarket::G2ppWithMarket()"
                << std::endl
                << "size must be 2." << std::endl;
            return;
        }
        mA            = -mDriftCoeff( 0 );
        mB            = -mDriftCoeff( 1 );
        double lSigma = mVolCoeff( 0, 0 );
        double lEta   = std::sqrt( mVolCoeff( 1, 0 ) * mVolCoeff( 1, 0 ) +
                                   mVolCoeff( 1, 1 ) * mVolCoeff( 1, 1 ) );
        double lRho   = mVolCoeff( 1, 0 ) / lEta;
        mFactorA      = 0.5 * lSigma * lSigma / ( mA * mA );
        mFactorB      = 0.5 * lEta * lEta / ( mB * mB );
        mFactorAB     = lRho * lSigma * lEta / ( mA * mB );
    }
};

class G2ppWithMarketBuilder : public MultiFactorAbstractBuilder
{
private:
    Math::Vec mDriftCoeff = Math::Vec( 0 );
    Math::Mat mVolCoeff   = Math::Mat( 0, 0 );
    std::unique_ptr<MarketData::ZCB> muMarketZCB;

public:
    G2ppWithMarketBuilder& setDrift( const Math::Vec& inDriftCoeff )
    {
        mDriftCoeff = inDriftCoeff;
        return *this;
    }
    G2ppWithMarketBuilder& setVol( const Math::Mat& inVolCoeff )
    {
        mVolCoeff = inVolCoeff;
        return *this;
    }
    G2ppWithMarketBuilder& setMarketZCB( const MarketData::ZCB& inMarketZCB )
    {
        muMarketZCB = std::make_unique<MarketData::ZCB>( inMarketZCB );
        return *this;
    }
    G2ppWithMarket build()
    {
        return G2ppWithMarket( mNPath, *muTerms, std::move( muRandomPath ),
                               mDriftCoeff, mVolCoeff, *muMarketZCB );
    }
    ModelAbstract& setInitState( const Math::Vec& inInitState ) = delete;
};

}  // namespace ShortRateMCMulti
}  // namespace Process

#endif