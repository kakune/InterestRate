/**
 * @file core.cpp
 * @brief This implements class for forward model.
 * @author kakune
 * @date 3/29/2024
 */

#ifdef NINCLUDE_TPP
#include "LIBOR/forward/core.hpp"
#endif

namespace LIBOR::Forward
{

template <C_StepCalc StepCalculator_>
template <Process::RandomVec::C_StdBrown StdBrownVecGenerator_>
typename StepCalculator_::ForwardRatesType
Factory<StepCalculator_>::createForwardRates() const
{
    StdBrownVecGenerator_ lStdBrownGen( mInitFR.size() );
    typename StepCalculator_::StatesType lStates(
        mNPath, mTenor[mTenor.size() - 1] + 1, mInitFR );
    for ( std::size_t iTerm = 1; iTerm <= mTenor[mTenor.size() - 1]; ++iTerm )
    {
        lStdBrownGen.initialize( mTenor.minIndex( iTerm ) );
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lStates.setStateElement(
                iPath, iTerm,
                lStates.getStates()[iPath][iTerm - 1] +
                    mStepCalc( iPath, iTerm, lStates, lStdBrownGen() ) );
        }
    }
    return typename StepCalculator_::ForwardRatesType( mTerms, mTenor,
                                                       lStates.getForwards() );
}

}  // namespace LIBOR::Forward