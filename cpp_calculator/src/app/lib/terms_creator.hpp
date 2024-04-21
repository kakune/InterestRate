#ifndef APP_LIB_TERMS_CREATOR_HPP
#define APP_LIB_TERMS_CREATOR_HPP

#include <algorithm>
#include <vector>

#include "process/market_data.hpp"
#include "utils/parameters.hpp"

namespace APP
{

Process::MarketData::Terms prepareTerms( const Utils::Parameters& inParams )
{
    std::size_t lNTerms = inParams.operator()<int>( "NTerms" );
    std::vector<double> lTerms( lNTerms, 0 );

    double lDt = inParams( "TimeMaturity" ) / double( lNTerms - 1 );
    for ( std::size_t iTerm = 1; iTerm < lNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }

    return Process::MarketData::Terms( lTerms );
}

Process::MarketData::Tenor prepareTenor(
    const Utils::Parameters& inParams,
    const Process::MarketData::Terms& inTerms )
{
    std::vector<int> lInds =
        inParams.operator()<std::vector<int>>( "IndTenor" );
    std::vector<std::size_t> lIndTenor( lInds.begin(), lInds.end() );
    return Process::MarketData::Tenor( inTerms, lIndTenor );
}

}  // namespace APP

#endif