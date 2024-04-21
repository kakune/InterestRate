/**
 * @file core.cpp
 * @brief This implements class for forward model.
 * @author kakune
 * @date 3/29/2024
 */

#ifndef LIBOR_FORWARD_CORE_TPP
#define LIBOR_FORWARD_CORE_TPP

#include <type_traits>

#include "LIBOR/forward/core.hpp"

namespace LIBOR::Forward
{

template <C_StepCalc StepCalculator_>
template <Process::RandomVec::C_StdBrown StdBrownVecGenerator_>
typename StepCalculator_::ForwardRatesType
Factory<StepCalculator_>::createForwardRates() const
{
    StdBrownVecGenerator_ lStdBrownGen( mInitFR.size() );
    typename StepCalculator_::StatesType lStates( mNPath, mTerms.size(),
                                                  mInitFR );
    for ( std::size_t iTerm = 1; iTerm < mTerms.size(); ++iTerm )
    {
        lStdBrownGen.initialize();
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lStates.setStateElement(
                iPath, iTerm,
                lStates.getStates()[iPath][iTerm - 1] +
                    mStepCalc( iPath, iTerm, lStates, lStdBrownGen() ) );
        }
    }
    return typename StepCalculator_::ForwardRatesType( mTerms, mTenor,
                                                       lStates.getValues() );
}

}  // namespace LIBOR::Forward

#endif