#ifndef APP_LIB_SHORT_RATE_CREATOR_HPP
#define APP_LIB_SHORT_RATE_CREATOR_HPP

#include "LIBOR/forward.hpp"
#include "process/market_data.hpp"
#include "utils/parameters.hpp"

namespace APP
{

LIBOR::Forward::Data::TerminalMeas createForwardTerminalFromParam(
    const Utils::Parameters& inParams,
    const Process::MarketData::Terms& inTerms,
    const Process::MarketData::Tenor& inTenor );

LIBOR::Forward::Data::TerminalMeas createForwardTerminalFromParam(
    const Utils::Parameters& inParams );

}  // namespace APP

#endif