#ifndef APP_LIB_TERMS_CREATOR_HPP
#define APP_LIB_TERMS_CREATOR_HPP

#include "process/market_data.hpp"
#include "utils/parameters.hpp"

namespace APP
{

Process::MarketData::Terms prepareTerms( const Utils::Parameters& inParams );

Process::MarketData::Tenor prepareTenor(
    const Utils::Parameters& inParams,
    const Process::MarketData::Terms& inTerms );

}  // namespace APP

#endif