/**
 * @file core.cpp
 * @brief This implements class for forward model.
 * @author kakune
 * @date 3/29/2024
 */

#include "LIBOR/forward/core.hpp"

namespace LIBOR
{
namespace Forward
{

Math::Vec ModelAbstract::transfFRToState( const Math::Vec& inFR,
                                          std::size_t inIndTime ) const
{
    return inFR;
}
Math::Vec ModelAbstract::transfStateToFR( const Math::Vec& inState,
                                          std::size_t inIndTime ) const
{
    return inState;
}

Process::ModelData::ForwardRates ModelAbstract::createForwardRates() const
{
    std::vector<std::vector<Math::Vec>> lStates(
        mNPath, std::vector<Math::Vec>( mTerms.size(),
                                        transfFRToState( mInitFRs, 0 ) ) );
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
    {
        muStdBrown->initialize();
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lStates[iPath][iTerm] =
                lStates[iPath][iTerm - 1] + driftTerm( iPath, iTerm, lStates ) +
                volTerm( iPath, iTerm, lStates,
                         ( *muStdBrown )() * mTerms.sqrtDifTime( iTerm ) );
        }
    }
    std::vector<std::vector<Math::Vec>> lFRs(
        mNPath, std::vector<Math::Vec>( mTerms.size(), Math::Vec( mNFR ) ) );
    for ( std::size_t iTerm = 0; iTerm < mTerms.size(); ++iTerm )
    {
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lFRs[iPath][iTerm] =
                transfStateToFR( lStates[iPath][iTerm], iTerm );
        }
    }
    return Process::ModelData::ForwardRates( mTerms, mIndTenor, lFRs );
}

std::vector<double> ModelAbstract::calcTauTenor(
    const std::vector<std::size_t>& inIndTenor,
    const Process::MarketData::Terms& inTerms ) const
{
    std::vector<double> lResult( inIndTenor.size() - 1 );
    for ( std::size_t iTenor = 0; iTenor < inIndTenor.size(); ++iTenor )
    {
        lResult[iTenor] =
            inTerms[inIndTenor[iTenor + 1]] - inTerms[inIndTenor[iTenor]];
    }
    return lResult;
}

}  // namespace Forward
}  // namespace LIBOR