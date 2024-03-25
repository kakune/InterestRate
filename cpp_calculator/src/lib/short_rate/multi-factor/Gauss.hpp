/**
 * @file Gauss.hpp
 * @brief This defines multi-factor Gaussian short-rate models. (see Chap.4 of
 * Brigo)
 * @author kakune
 * @date 3/24/2024
 */

#ifndef SHORT_RATE_MULTI_FACTOR_GAUSS_HPP
#define SHORT_RATE_MULTI_FACTOR_GAUSS_HPP

#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include "process/market_data.hpp"
#include "short_rate/multi-factor/core.hpp"

namespace ShortRate
{
namespace MultiFactor
{

/**
 * @brief This is multi-factor short-rate Gaussian model with constant coeff
 * without shift $\varphi$.
 */
class ConstantGauss : public MultiFactorAbstract
{
private:
    const Math::Vec mDriftCoeff;  //! coefficient of dt
    const Math::Mat
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
                              std::size_t inIndTime ) const override;

public:
    ConstantGauss(
        std::size_t inNPath, const Process::MarketData::Terms& inTerms,
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
    const Math::Vec mDriftCoeff;  //! coefficient of dt
    const Math::Mat
        mVolCoeff;  //! coefficient of dW which must be LowerTriangular Matrix
    const Process::MarketData::ZCB mMarketZCB;  //! the market data

    double mA, mB;
    double mFactorA, mFactorB, mFactorAB;  //! factor for transfStateToRate
    std::vector<double>
        mInstantaneousFRs;  //! instantaneous forward rates at each time

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
                              std::size_t inIndTime ) const override;
    void prepareFactors();

public:
    G2ppWithMarket(
        std::size_t inNPath, const Process::MarketData::Terms& inTerms,
        std::unique_ptr<Process::RandomVec::PathAbstract> inuRandomPath,
        const Math::Vec& inDriftCoeff, const Math::Mat& inVolCoeff,
        const Process::MarketData::ZCB& inMarketZCB ) :
        MultiFactorAbstract( inNPath, inTerms, Math::Vec( 2, 0.0 ),
                             std::move( inuRandomPath ) ),
        mDriftCoeff( inDriftCoeff ),
        mVolCoeff( inVolCoeff ),
        mMarketZCB( inMarketZCB ),
        mInstantaneousFRs( inTerms.size() )
    {
        if ( inDriftCoeff.size() != 2 || inVolCoeff.sizeRow() != 2 ||
             inVolCoeff.sizeCol() != 2 )
        {
            throw std::invalid_argument(
                "Error : "
                "Process::ShortRateMCMulti::G2ppWithMarket::G2ppWithMarket()"
                "\nsize must be 2." );
        }
        prepareFactors();
    }
};

class G2ppWithMarketBuilder : public MultiFactorAbstractBuilder
{
private:
    Math::Vec mDriftCoeff = Math::Vec( 0 );
    Math::Mat mVolCoeff   = Math::Mat( 0, 0 );
    std::unique_ptr<Process::MarketData::ZCB> muMarketZCB;

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
    G2ppWithMarketBuilder& setMarketZCB(
        const Process::MarketData::ZCB& inMarketZCB )
    {
        muMarketZCB = std::make_unique<Process::MarketData::ZCB>( inMarketZCB );
        return *this;
    }
    G2ppWithMarket build()
    {
        return G2ppWithMarket( mNPath, *muTerms, std::move( muRandomPath ),
                               mDriftCoeff, mVolCoeff, *muMarketZCB );
    }
    ModelAbstract& setInitState( const Math::Vec& inInitState ) = delete;
};

}  // namespace MultiFactor
}  // namespace ShortRate

#endif