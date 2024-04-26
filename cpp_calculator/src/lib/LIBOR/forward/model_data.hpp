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

#include <concepts>

#include "LIBOR/forward/payoff.hpp"
#include "analytical/Black76.hpp"
#include "process/market_data.hpp"

namespace LIBOR::Forward
{
namespace States
{

template <typename T_, typename ElementState_>
concept C_States = requires( T_ inObj, std::size_t inIndPath,
                             std::size_t inIndTerms, ElementState_ inVal ) {
    inObj.setStateElement( inIndPath, inIndTerms, inVal );
    inObj.setForwardElement( inIndPath, inIndTerms, inVal );
    {
        inObj.getStates()
    } -> std::same_as<const std::vector<std::vector<ElementState_>>&>;
    {
        inObj.getForwards()
    } -> std::same_as<const std::vector<std::vector<ElementState_>>&>;
};
template <typename T_>
concept C_LogExpElement = requires( T_ inObj ) {
    log( inObj );
    exp( inObj );
};

/**
 * @brief State = forward rate
 */
template <typename ElementState_> class Plain;
/**
 * @brief State = log( forward rate )
 */
template <C_LogExpElement ElementState_> class Log;

}  // namespace States

namespace Data
{

/**
 * @brief member function of Black76 double(double, size_t)
 */
template <auto Func_>
concept C_Black76OneTermFunction =
    requires( const Analytical::Black76::Model& inObj, double inStrike,
              std::size_t inIndex ) {
        { ( inObj.*Func_ )( inStrike, inIndex ) } -> std::same_as<double>;
    };

/**
 * @brief member function of Black76 double(double, size_t, size_t)
 */
template <auto Func_>
concept C_Black76MultiTermFunction =
    requires( const Analytical::Black76::Model& inObj, double inStrike,
              std::size_t inIndStart, std::size_t inIndLast ) {
        {
            ( inObj.*Func_ )( inStrike, inIndStart, inIndLast )
        } -> std::same_as<double>;
    };

Process::MarketData::ZCB createZCBFromForwardRates(
    const Process::MarketData::Tenor& inTenor,
    const std::vector<std::vector<Math::Vec>>& inFRs, std::size_t inDeg );

/**
 * @brief This stores forward rate data at each path and each term.
 */
template <class Derived> class Abstract
{
protected:
    const std::size_t mNPath;                 //! the number of path
    const Process::MarketData::Terms mTerms;  //! terms
    const Process::MarketData::Tenor mTenor;  //! tenor
    const std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
        msDataForwardRates;  //! forward rates
    const std::shared_ptr<const Process::MarketData::ZCB>
        msZCB;  //! initial ZCB

    template <auto BlackPriceFunc_>
        requires C_Black76OneTermFunction<BlackPriceFunc_>
    double calcBlackImpVolByOneTerm( double inStrike, std::size_t inIndTenor,
                                     double inModelPrice ) const;
    template <auto BlackPriceFunc_>
        requires C_Black76MultiTermFunction<BlackPriceFunc_>
    double calcBlackImpVolByMultiTerm( double inStrike,
                                       std::size_t inIndTenorStart,
                                       std::size_t inIndTenorLast,
                                       double inModelPrice ) const;

public:
    /**
     * @brief This constructs a new Forward Rates.
     * @param inTerms Term structure
     * @param inIndTenor indices of tenor of FR
     * @param insDataForwardRates forward rate
     * @param inDegZCB the degree for the spline of ZCB
     */
    Abstract( const Process::MarketData::Terms& inTerms,
              const Process::MarketData::Tenor& inTenor,
              std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
                  insDataForwardRates,
              std::size_t inDegZCB = 3 );
    Abstract( const Process::MarketData::Terms& inTerms,
              const Process::MarketData::Tenor& inTenor,
              const std::vector<std::vector<Math::Vec>>& inDataForwardRate,
              std::size_t inDegZCB = 3 );
    const std::vector<Math::Vec>& operator[]( std::size_t inIndex ) const
    {
        return ( *msDataForwardRates )[inIndex];
    }

    double calcCaplet( double inStrike, std::size_t inIndTenor ) const;
    double calcFloorlet( double inStrike, std::size_t inIndTenor ) const;
    double calcPayerSwaption( double inStrike, std::size_t inIndTenorStart,
                              std::size_t inIndTenorLast ) const;
    double calcReceiverSwaption( double inStrike, std::size_t inIndTenorStart,
                                 std::size_t inIndTenorLast ) const;
    double calcBlackImpVolByCaplet( double inStrike,
                                    std::size_t inIndTenor ) const;
    double calcBlackImpVolByFloorlet( double inStrike,
                                      std::size_t inIndTenor ) const;
    double calcBlackImpVolByPayerSwaption( double inStrike,
                                           std::size_t inIndTenorStart,
                                           std::size_t inIndTenorLast ) const;
    double calcBlackImpVolByReceiverSwaption(
        double inStrike, std::size_t inIndTenorStart,
        std::size_t inIndTenorLast ) const;
};

class TerminalMeas : public Abstract<TerminalMeas>
{
public:
    TerminalMeas( const Process::MarketData::Terms& inTerms,
                  const Process::MarketData::Tenor& inTenor,
                  std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
                      insDataForwardRates,
                  std::size_t inDegZCB = 3 ) :
        Abstract( inTerms, inTenor, insDataForwardRates, inDegZCB )
    {
    }
    TerminalMeas( const Process::MarketData::Terms& inTerms,
                  const Process::MarketData::Tenor& inTenor,
                  const std::vector<std::vector<Math::Vec>>& inDataForwardRates,
                  std::size_t inDegZCB = 3 ) :
        Abstract( inTerms, inTenor, inDataForwardRates, inDegZCB )
    {
    }
    template <LIBOR::Forward::Payoff::C_OneTerm PayoffObject_>
    double calcExpectation( PayoffObject_ inPayoff ) const;
    template <LIBOR::Forward::Payoff::C_MultiTerm PayoffObject_>
    double calcExpectation( PayoffObject_ inPayoff ) const;
};

class SpotMeas : public Abstract<SpotMeas>
{
public:
    SpotMeas( const Process::MarketData::Terms& inTerms,
              const Process::MarketData::Tenor& inTenor,
              std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
                  insDataForwardRates,
              std::size_t inDegZCB = 3 ) :
        Abstract( inTerms, inTenor, insDataForwardRates, inDegZCB )
    {
    }
    SpotMeas( const Process::MarketData::Terms& inTerms,
              const Process::MarketData::Tenor& inTenor,
              const std::vector<std::vector<Math::Vec>>& inDataForwardRates,
              std::size_t inDegZCB = 3 ) :
        Abstract( inTerms, inTenor, inDataForwardRates, inDegZCB )
    {
    }
    template <LIBOR::Forward::Payoff::C_OneTerm PayoffObject_>
    double calcExpectation( PayoffObject_ inPayoff ) const;
    template <LIBOR::Forward::Payoff::C_MultiTerm PayoffObject_>
    double calcExpectation( PayoffObject_ inPayoff ) const;
};

}  // namespace Data

}  // namespace LIBOR::Forward

#ifndef NINCLUDE_TPP
#include "LIBOR/forward/model_data.tpp"
#endif

#endif