/**
 * @file core.hpp
 * @brief This defines class for forward model.
 * @author kakune
 * @date 3/29/2024
 */

#ifndef LIBOR_FORWARD_CORE_HPP
#define LIBOR_FORWARD_CORE_HPP

#include <memory>
#include <vector>

#include "LIBOR/forward/model_data.hpp"
#include "math/matrix.hpp"
#include "process/market_data.hpp"
#include "process/random.hpp"

namespace LIBOR::Forward
{

/**
 * @brief concept step calculators must satisfy.
 */
template <class T_>
concept C_StepCalc =
    requires( T_ inObj, std::size_t inIndPath, std::size_t inIndTerm,
              const typename T_::StatesType& inStates,
              const Math::Vec& inStdBrownianVec ) {
        typename T_::ForwardRatesType;
        typename T_::StatesType;
        {
            inObj( inIndPath, inIndTerm, inStates, inStdBrownianVec )
        } -> std::convertible_to<Math::Vec>;
    };
/**
 * @brief concept volatility generators must satisfy.
 */
template <class T_>
concept C_VolGen =
    requires( T_ inObj, std::size_t inIndPath, std::size_t inIndTerm,
              const std::vector<std::vector<Math::Vec>>& inForwardRates,
              const Math::Vec& inStdBrownForFR ) {
        {
            inObj( inIndPath, inIndTerm, inForwardRates, inStdBrownForFR )
        } -> std::convertible_to<Math::Vec>;
    };

/**
 * @brief This is the factory class for market forward rate models.
 */
template <C_StepCalc StepCalculator_> class Factory
{
private:
    const std::size_t mNPath;                 //! the number of Path
    const Process::MarketData::Terms mTerms;  //! terms
    const Process::MarketData::Tenor mTenor;  //! tenor
    const Math::Vec mInitFR;                  //! initial forward rate
    const StepCalculator_ mStepCalc;          //! step calculator

public:
    /**
     * @brief This constructs a new Factory.
     *
     * @param inNPath the number of Monte Carlo path.
     * @param inTerms term structure
     * @param inTenor tenor structure
     * @param inInitFR initial forward rate
     * @param inStepCalc the calculator for each steps
     */
    Factory( std::size_t inNPath, const Process::MarketData::Terms& inTerms,
             const Process::MarketData::Tenor& inTenor,
             const Math::Vec& inInitFR, const StepCalculator_& inStepCalc ) :
        mNPath( inNPath ),
        mTerms( inTerms ),
        mTenor( inTenor ),
        mInitFR( inInitFR ),
        mStepCalc( inStepCalc )
    {
    }
    template <Process::RandomVec::C_StdBrown StdBrownVecGenerator_>
    typename StepCalculator_::ForwardRatesType createForwardRates() const;
};

}  // namespace LIBOR::Forward
#include "LIBOR/forward/core.tpp"
#endif
