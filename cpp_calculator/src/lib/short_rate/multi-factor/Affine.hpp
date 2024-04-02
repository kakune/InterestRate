/**
 * @file Affine.hpp
 * @brief This defines multi-factor Affine short-rate models. (see Chap.4 of
 * Brigo)
 * @author kakune
 * @date 3/25/2024
 */

#ifndef SHORT_RATE_MULTI_FACTOR_AFFINE_HPP
#define SHORT_RATE_MULTI_FACTOR_AFFINE_HPP

#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include "math/matrix.hpp"
#include "process/market_data.hpp"
#include "short_rate/multi-factor/core.hpp"

namespace ShortRate
{
namespace MultiFactor
{

/**
 * @brief This is CIR++ model with constant coeff and market data.
 */
class CIR2ppWithMarket : public MultiFactorAbstract
{
private:
    const Math::Vec mConvSH, mMean;  //! coefficient of dt
    const Math::Vec
        mVol;  //! coefficient of dW which must be LowerTriangular Matrix
    const Process::MarketData::ZCB mMarketZCB;  //! the market data
    const std::vector<double>
        mShifts;  //! shift of short rate for transfStateToRate

    // ConvSH(Mean - x)dt
    Math::Vec driftCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<Math::Vec>>& inSpots ) const override;
    // VolCoeff.sqrt(x).(dW_1, dW_2)
    Math::Vec volTerm( std::size_t inIndPath, std::size_t inIndTerm,
                       const std::vector<std::vector<Math::Vec>>& inSpots,
                       const Math::Vec& inRandomVec ) const override;
    // r = x + y + \phi
    double transfStateToRate( const Math::Vec& inState,
                              std::size_t inIndTime ) const override;
    std::vector<double> calcShifts() const;

public:
    CIR2ppWithMarket(
        std::size_t inNPath, const Process::MarketData::Terms& inTerms,
        std::unique_ptr<Process::RandomVec::StdBrownAbstract> inuRandomPath,
        const Math::Vec& inConvSH, const Math::Vec& inMean,
        const Math::Vec& inVol, const Process::MarketData::ZCB& inMarketZCB ) :
        MultiFactorAbstract( inNPath, inTerms, Math::Vec( 2, 0.0 ),
                             std::move( inuRandomPath ) ),
        mConvSH( inConvSH ),
        mMean( inMean ),
        mVol( inVol ),
        mMarketZCB( inMarketZCB ),
        mShifts( calcShifts() )
    {
        if ( inConvSH.size() != 2 || inMean.size() != 2 || inVol.size() != 2 )
        {
            throw std::invalid_argument(
                "Error : "
                "Process::ShortRateMCMulti::CIR2ppWithMarket::CIR2ppWithMarket"
                "\nsize must be 2." );
        }
        if ( ( 2.0 * mConvSH * mMean - mVol * mVol ).min() <= 0.0 )
        {
            throw std::invalid_argument(
                "Error : "
                "Process::ShortRateMCMulti::CIR2ppWithMarket::CIR2ppWithMarket"
                "\n2 k theta must be larger than vol^2." );
        }
    }
};

class CIR2ppWithMarketBuilder : public MultiFactorAbstractBuilder
{
private:
    Math::Vec mConvSH = Math::Vec( 0 );
    Math::Vec mMean   = Math::Vec( 0 );
    Math::Vec mVol    = Math::Vec( 0 );
    std::unique_ptr<Process::MarketData::ZCB> muMarketZCB;

public:
    CIR2ppWithMarketBuilder& setDrift( const Math::Vec& inConvSH,
                                       const Math::Vec& inMean )
    {
        mConvSH = inConvSH;
        mMean   = inMean;
        return *this;
    }
    CIR2ppWithMarketBuilder& setVol( const Math::Vec& inVol )
    {
        mVol = inVol;
        return *this;
    }
    CIR2ppWithMarketBuilder& setMarketZCB(
        const Process::MarketData::ZCB& inMarketZCB )
    {
        muMarketZCB = std::make_unique<Process::MarketData::ZCB>( inMarketZCB );
        return *this;
    }
    CIR2ppWithMarket build()
    {
        return CIR2ppWithMarket( mNPath, *muTerms, std::move( muStdBrown ),
                                 mConvSH, mMean, mVol, *muMarketZCB );
    }
    ModelAbstractBuilder& setInitState( const Math::Vec& ) = delete;
};

}  // namespace MultiFactor
}  // namespace ShortRate

#endif