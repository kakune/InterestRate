/**
 * @file model_data.tpp
 * @brief This implements model data.
 * @author kakune
 * @date 4/21/2024
 */

#include "LIBOR/forward/model_data.hpp"

namespace LIBOR::Forward
{
namespace Data
{

Process::MarketData::ZCB createZCBFromForwardRates(
    const Process::MarketData::Tenor& inTenor,
    const std::vector<std::vector<Math::Vec>>& inFRs, std::size_t inDeg )
{
    std::vector<double> lTmpTerms( inTenor.size() + 1, 0.0 );
    std::vector<double> lTmpZCBs( inTenor.size() + 1, 1.0 );
    for ( std::size_t i = 1; i < inTenor.size() + 1; ++i )
    {
        lTmpTerms[i] = inTenor.term( i );
        lTmpZCBs[i]  = lTmpZCBs[i - 1] /
                      ( inTenor.tau( i - 1 ) * inFRs[0][0][i - 1] + 1.0 );
    }
    return Process::MarketData::ZCB( Process::MarketData::Terms( lTmpTerms ),
                                     lTmpZCBs,
                                     std::min( inDeg, inTenor.size() ) );
}

}  // namespace Data
}  // namespace LIBOR::Forward