/**
 * @file vol_generator.hpp
 * @brief This defines volatility generator classes for forward model.
 * @author kakune
 * @date 4/15/2024
 */

#ifndef LIBOR_FORWARD_VOL_GENERATOR_HPP
#define LIBOR_FORWARD_VOL_GENERATOR_HPP

#include <cstddef>
#include <memory>
#include <vector>

#include "math/matrix.hpp"
#include "process/market_data.hpp"
#include "process/random.hpp"

namespace LIBOR::Forward::VolGen
{

class Constant;
template <Process::RandomVec::C_StdBrown StdBrownVecGenerator_> class SABR;

class Constant
{
private:
    const Math::Vec mVol;

public:
    Constant( const Math::Vec& inVol );
    Math::Vec operator()(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<Math::Vec>>& inForwardRates,
        const Math::Vec& inStdBrownForFR ) const;
};

template <Process::RandomVec::C_StdBrown StdBrownVecGenerator_> class SABR
{
private:
    const std::shared_ptr<std::vector<Math::Vec>> msVols;
    const double mExponent;
    const double mVolVol;
    const Math::Vec mCorrSV, mCorrSVInv;
    mutable std::size_t mTmpIndTerm;
    mutable double mVolVolSqrtDt;
    const Process::MarketData::Terms mTerms;
    const Process::MarketData::Tenor mTenor;
    const std::shared_ptr<StdBrownVecGenerator_> msStdBrownGen;

public:
    SABR( const Math::Vec& inInitVol, double inExponent, double inVolVol,
          const Math::Vec& inCorrSV, std::size_t inNPath,
          const Process::MarketData::Terms& inTerms,
          const Process::MarketData::Tenor& inTenor );
    Math::Vec operator()(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<Math::Vec>>& inForwardRates,
        const Math::Vec& inStdBrownForFR ) const;
};

}  // namespace LIBOR::Forward::VolGen

#ifndef NINCLUDE_TPP
#include "LIBOR/forward/vol_generator.tpp"
#endif

#endif