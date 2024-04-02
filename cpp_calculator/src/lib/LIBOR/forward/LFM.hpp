/**
 * @file LFM.hpp
 * @brief This defines class for lognormal forward model.
 * @author kakune
 * @date 3/29/2024
 */

#ifndef MARKET_MODEL_FORWARD_LFM_HPP
#define MARKET_MODEL_FORWARD_LFM_HPP

#include "LIBOR/forward/core.hpp"

namespace LIBOR
{
namespace Forward
{

class ConstantLFM : public ModelAbstract
{
private:
    Math::Vec mVols;
    Math::Vec mNegativeHalfVolVol;
    Math::Vec driftTerm(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<Math::Vec>>& inStates ) const override;
    Math::Vec volTerm( std::size_t inIndPath, std::size_t inIndTerm,
                       const std::vector<std::vector<Math::Vec>>& inStates,
                       const Math::Vec& inRandomVec ) const override;
    Math::Vec transfFRToState( const Math::Vec& inFR,
                               std::size_t inIndTime ) const override;
    Math::Vec transfStateToFR( const Math::Vec& inState,
                               std::size_t inIndTime ) const override;

public:
    ConstantLFM(
        std::size_t inNPath, const Process::MarketData::Terms& inTerms,
        const std::vector<std::size_t>& inIndTenor, const Math::Vec& inInitFRs,
        std::unique_ptr<Process::RandomVec::StdBrownAbstract> inuRandomPath,
        const Math::Mat& inCorr, const Math::Vec& inVols ) :
        ModelAbstract( inNPath, inTerms, inIndTenor, inInitFRs,
                       std::move( inuRandomPath ), inCorr ),
        mVols( inVols ),
        mNegativeHalfVolVol( -0.5 * inVols * inVols )
    {
    }
};

class ConstantLFMBuilder : public ModelAbstractBuilder
{
private:
    Math::Vec mVols = Math::Vec( 0 );

public:
    ConstantLFMBuilder& setVols( const Math::Vec& inVols )
    {
        mVols = inVols;
        return *this;
    }
    ConstantLFM build()
    {
        return ConstantLFM( mNPath, *muTerms, mIndTenor, mInitFRs,
                            std::move( muStdBrown ), mCorr, mVols );
    }
};

}  // namespace Forward
}  // namespace LIBOR

#endif