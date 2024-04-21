/**
 * @file model_data.cpp
 * @brief This implements classes for model data.
 * @author kakune
 * @date 3/30/2024
 */

#include <cmath>
#include <memory>
#include <vector>

#include "analytical/Black76.hpp"
#include "math/findroot_1d.hpp"
#include "process/model_data.hpp"

namespace Process
{
namespace ModelData
{

MarketData::ZCB createZCBFromForwardRates(
    const MarketData::Tenor& inTenor,
    const std::vector<std::vector<Math::Vec>>& inFRs, std::size_t inDeg )
{
    std::vector<double> lTmpTerms( inTenor.size() + 1, 0.0 );
    std::vector<double> lTmpZCBs( inTenor.size() + 1, 1.0 );
    for ( std::size_t i = 1; i < inTenor.size() + 1; ++i )
    {
        lTmpTerms[i] = inTenor.term( i );
        lTmpZCBs[i]  = lTmpZCBs[i - 1] /
                      ( inTenor.tau( i - 1 ) * inFRs[0][0]( i - 1 ) + 1.0 );
    }
    return MarketData::ZCB( MarketData::Terms( lTmpTerms ), lTmpZCBs,
                            std::min( inDeg, inTenor.size() ) );
}

}  // namespace ModelData
}  // namespace Process
