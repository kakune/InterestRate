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

}  // namespace ShortRate
}  // namespace Process

#endif