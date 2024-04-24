#ifndef APP_LIB_SHORT_RATE_CREATOR_HPP
#define APP_LIB_SHORT_RATE_CREATOR_HPP

#include "process/market_data.hpp"
#include "short_rate/multi-factor.hpp"
#include "short_rate/one-factor.hpp"
#include "utils/parameters.hpp"

namespace APP
{

ShortRate::SpotRates createSpotRateFromParam(
    const Utils::Parameters& inParams,
    const Process::MarketData::Terms& inTerms );

ShortRate::SpotRates createSpotRateFromParam(
    const Utils::Parameters& inParams );

ShortRate::SpotRates createSpotRateFromMarket(
    const Utils::Parameters& inParams, const Process::MarketData::Terms inTerms,
    const Process::MarketData::ZCB inMarketZCB );

ShortRate::SpotRates createSpotRateFromMarket(
    const Utils::Parameters& inParams,
    const Process::MarketData::ZCB inMarketZCB );

}  // namespace APP

#endif