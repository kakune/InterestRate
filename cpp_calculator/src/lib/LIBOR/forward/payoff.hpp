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

namespace LIBOR::Forward
{

class CapletFloorletPayoff
{
private:
    const std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
        msDataForwardRates;
    const double mStrike;
    const std::size_t mIndTermsFR;
    const std::size_t mIndTenorFR;
    const std::size_t mIndTenorPay;
    const double mTau;
    bool mIsCaplet;

public:
    CapletFloorletPayoff(
        const Process::MarketData::Tenor& inTenor,
        std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
            insDataForwardRates,
        double inStrike, std::size_t inIndTenorFR, std::size_t inIndTenorPay,
        bool inIsCaplet );
    std::size_t getIndexTenorPay() const;
    double operator()( std::size_t inIndPath ) const;
};

}  // namespace LIBOR::Forward

#endif