/**
 * @file model_data.hpp
 * @brief This defines classes for model data. Unless there is a special
 * reason, calculation results should be stored in these classes. Since most
 * data is stored by smart pointers, there is no need to worry about the cost of
 * copying instances.
 * @author kakune
 * @date 4/21/2024
 */

#ifndef LIBOR_FORWARD_MODEL_DATA_HPP
#define LIBOR_FORWARD_MODEL_DATA_HPP

#include <type_traits>

#include "LIBOR/forward/payoff.hpp"
#include "analytical/Black76.hpp"
#include "math/findroot_1d.hpp"
#include "process/market_data.hpp"

namespace LIBOR::Forward
{

Process::MarketData::ZCB createZCBFromForwardRates(
    const Process::MarketData::Tenor& inTenor,
    const std::vector<std::vector<Math::Vec>>& inFRs, std::size_t inDeg );

/**
 * @brief This stores forward rate data at each path and each term.
 */
template <class DerivedForwardRates> class ForwardRatesAbstract
{
protected:
    const std::size_t mNPath;                 //! the number of path
    const Process::MarketData::Terms mTerms;  //! terms
    const Process::MarketData::Tenor mTenor;  //! tenor
    const std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
        msDataForwardRates;  //! forward rates
    const std::shared_ptr<const Process::MarketData::ZCB> msZCB;

public:
    /**
     * @brief This constructs a new Forward Rates.
     * @param inTerms Term structure
     * @param inIndTenor indices of tenor of FR
     * @param insDataForwardRates forward rate
     * @param inDegZCB the degree for the spline of ZCB
     */
    ForwardRatesAbstract(
        const Process::MarketData::Terms& inTerms,
        const Process::MarketData::Tenor& inTenor,
        std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
            insDataForwardRates,
        std::size_t inDegZCB = 3 );
    ForwardRatesAbstract(
        const Process::MarketData::Terms& inTerms,
        const Process::MarketData::Tenor& inTenor,
        const std::vector<std::vector<Math::Vec>>& inDataForwardRate,
        std::size_t inDegZCB = 3 );
    const std::vector<Math::Vec>& operator[]( std::size_t inIndex ) const
    {
        return ( *msDataForwardRates )[inIndex];
    }

    double calcCaplet( double inStrike, std::size_t inIndTenor ) const;
    double calcFloorlet( double inStrike, std::size_t inIndTenor ) const;
    double calcBlackImpVol( double inStrike, std::size_t inIndTenor,
                            bool inIsUseCaplet = true ) const;
};

class ForwardRatesTerminalMeas
    : public ForwardRatesAbstract<ForwardRatesTerminalMeas>
{
public:
    ForwardRatesTerminalMeas(
        const Process::MarketData::Terms& inTerms,
        const Process::MarketData::Tenor& inTenor,
        std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
            insDataForwardRates,
        std::size_t inDegZCB = 3 ) :
        ForwardRatesAbstract( inTerms, inTenor, insDataForwardRates, inDegZCB )
    {
    }
    ForwardRatesTerminalMeas(
        const Process::MarketData::Terms& inTerms,
        const Process::MarketData::Tenor& inTenor,
        const std::vector<std::vector<Math::Vec>>& inDataForwardRate,
        std::size_t inDegZCB = 3 ) :
        ForwardRatesAbstract( inTerms, inTenor, inDataForwardRate, inDegZCB )
    {
    }
    template <class PayoffObject_>
    double calcExpectation( PayoffObject_ inPayoff ) const;
};

template <typename T_, typename ElementState_>
concept C_State = requires( T_ inObj, std::size_t inIndPath,
                            std::size_t inIndTerms, ElementState_ inVal ) {
    inObj.setStateElement( inIndPath, inIndTerms, inVal );
    inObj.setValueElement( inIndPath, inIndTerms, inVal );
    {
        inObj.getStates()
    } -> std::same_as<const std::vector<std::vector<ElementState_>>&>;
    {
        inObj.getValues()
    } -> std::same_as<const std::vector<std::vector<ElementState_>>&>;
};
template <typename T_>
concept C_LogExpElement = requires( T_ inObj ) {
    log( inObj );
    exp( inObj );
};

template <typename ElementState_> class StatesPlain;
template <C_LogExpElement ElementState_> class StatesLog;

}  // namespace LIBOR::Forward

#include "LIBOR/forward/model_data.tpp"

#endif