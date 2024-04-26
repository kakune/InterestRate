/**
 * @file payoff.hpp
 * @brief This defines payoff function objects.
 * @author kakune
 * @date 4/21/2024
 */

#ifndef LIBOR_FORWARD_PAYOFF_HPP
#define LIBOR_FORWARD_PAYOFF_HPP

#include <memory>
#include <vector>

#include "math/matrix.hpp"
#include "process/market_data.hpp"

namespace LIBOR::Forward::Payoff
{

template <typename T_>
concept C_OneTerm = requires( T_ inObj, std::size_t inIndPath ) {
    { inObj( inIndPath ) } -> std::same_as<double>;
};

template <typename T_>
concept C_MultiTerm = requires( T_ inObj, std::size_t inIndPath ) {
    { inObj( inIndPath ) } -> std::same_as<Math::Vec>;
};

class CapletFloorlet
{
private:
    const std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
        msDataForwardRates;
    const std::size_t mIndTenorPay;
    const std::size_t mIndTenorFR;
    const std::size_t mIndTermsFR;
    const double mTau;
    const double mStrike;
    const bool mIsCaplet;

public:
    CapletFloorlet( const Process::MarketData::Tenor& inTenor,
                    std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
                        insDataForwardRates,
                    std::size_t inIndTenorPay, std::size_t inIndTenorFR,
                    double inStrike, bool inIsCaplet );
    std::size_t getIndexTenorPay() const;
    double operator()( std::size_t inIndPath ) const;
};

/**
 * @brief this calculates the swaption price AT START TIME.
 */
class Swaption
{
private:
    const Process::MarketData::Tenor mTenor;
    const std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
        msDataForwardRates;
    const std::size_t mIndTenorStart;
    const std::size_t mIndTermsStart;
    const std::size_t mIndTenorLast;
    const std::size_t mIndTermsLast;
    const Math::Vec mTaus;
    const double mStrike;
    const bool mIsPayer;

public:
    Swaption( const Process::MarketData::Tenor& inTenor,
              std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
                  insDataForwardRates,
              std::size_t inIndTenorStart, std::size_t inIndTenorLast,
              double inStrike, bool inIsPayer );
    std::size_t getIndexTenorFirstPay() const;
    std::size_t getIndexTenorLastPay() const;
    Math::Vec operator()( std::size_t inIndPath ) const;
};

// class Swaption
// {
// private:
//     const std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
//         msDataForwardRates;
//     const double mStrike;
//     const std::size_t mIndTermsFR;
//     const std::size_t mIndTenorFR;
//     const std::size_t mIndTenorPay;
//     const double mTau;
//     bool mIsCaplet;

// public:
//     Swaption( const Process::MarketData::Tenor& inTenor,
//               std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
//                   insDataForwardRates,
//               double inStrike, std::size_t inIndTenorFR,
//               std::size_t inIndTenorPay, bool inIsCaplet );
//     std::size_t getIndexTenorPay() const;
//     double operator()( std::size_t inIndPath ) const;
// };

}  // namespace LIBOR::Forward::Payoff

#endif