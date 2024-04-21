/**
 * @file step_calculator.hpp
 * @brief This defines step calculator classes for forward model.
 * @author kakune
 * @date 4/15/2024
 */

#ifndef LIBOR_FORWARD_STEP_CALCULATOR_HPP
#define LIBOR_FORWARD_STEP_CALCULATOR_HPP

#include "LIBOR/forward/vol_generator.hpp"

namespace LIBOR::Forward::StepCalc
{

Math::Mat calcUpperTriangularRhoWithOutDiag( const Math::Mat& inCorr );

template <C_VolGen VolatilityGenerator_> class LogNormalTerminalMeas;
template <C_VolGen VolatilityGenerator_> class LogNormalTerminalMeasWithLog;
template <C_VolGen VolatilityGenerator_> class NormalTerminalMeas;

template <C_VolGen VolatilityGenerator_> class LogNormalTerminalMeas
{
private:
    const Process::MarketData::Terms mTerms;
    const Process::MarketData::Tenor mTenor;
    const Math::Mat mCorr;
    const Math::Mat mRho;
    const Math::Mat mUpperRho;
    const VolatilityGenerator_ mVolGen;

public:
    using ForwardRatesType = LIBOR::Forward::Data::TerminalMeas;
    using StatesType       = LIBOR::Forward::States::Plain<Math::Vec>;
    LogNormalTerminalMeas( const Process::MarketData::Terms& inTerms,
                           const Process::MarketData::Tenor& inTenor,
                           const Math::Mat& inCorr,
                           const VolatilityGenerator_& inVol );
    Math::Vec operator()( std::size_t inIndPath, std::size_t inIndTerm,
                          const StatesType& inStates,
                          const Math::Vec& inStdBrownianVec ) const;
};

template <C_VolGen VolatilityGenerator_> class LogNormalTerminalMeasWithLog
{
private:
    const Process::MarketData::Terms mTerms;
    const Process::MarketData::Tenor mTenor;
    const Math::Mat mCorr;
    const Math::Mat mRho;
    const Math::Mat mUpperRho;
    const VolatilityGenerator_ mVolGen;

public:
    using ForwardRatesType = LIBOR::Forward::Data::TerminalMeas;
    using StatesType       = LIBOR::Forward::States::Log<Math::Vec>;
    LogNormalTerminalMeasWithLog( const Process::MarketData::Terms& inTerms,
                                  const Process::MarketData::Tenor& inTenor,
                                  const Math::Mat& inCorr,
                                  const VolatilityGenerator_& inVol );
    Math::Vec operator()( std::size_t inIndPath, std::size_t inIndTerm,
                          const StatesType& inStates,
                          const Math::Vec& inStdBrownianVec ) const;
};

template <C_VolGen VolatilityGenerator_> class NormalTerminalMeas
{
private:
    const Process::MarketData::Terms mTerms;
    const Process::MarketData::Tenor mTenor;
    const Math::Mat mCorr;
    const Math::Mat mRho;
    const Math::Mat mUpperRho;
    const VolatilityGenerator_ mVolGen;

public:
    using ForwardRatesType = LIBOR::Forward::Data::TerminalMeas;
    using StatesType       = LIBOR::Forward::States::Plain<Math::Vec>;
    NormalTerminalMeas( const Process::MarketData::Terms& inTerms,
                        const Process::MarketData::Tenor& inTenor,
                        const Math::Mat& inCorr,
                        const VolatilityGenerator_& inVol );
    Math::Vec operator()( std::size_t inIndPath, std::size_t inIndTerm,
                          const StatesType& inStates,
                          const Math::Vec& inStdBrownianVec ) const;
};

}  // namespace LIBOR::Forward::StepCalc

#include "LIBOR/forward/step_calculator.tpp"

#endif
