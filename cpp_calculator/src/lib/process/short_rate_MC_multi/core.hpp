/**
 * @file core.hpp
 * @brief This defines classes to calculate multi-factor short rate paths.
 * @author kakune
 * @date 3/23/2024
 */

#ifndef PROCESS_SHORT_RATE_MC_MULTI_CORE_HPP
#define PROCESS_SHORT_RATE_MC_MULTI_CORE_HPP

#include <memory>
#include <vector>

#include "math/interpolate_1d.hpp"
#include "math/matrix.hpp"
#include "process/market_data.hpp"
#include "process/random.hpp"

namespace Process
{
namespace ShortRateMCMulti
{

/**
 * @brief This is the abstract class for short rate models.
 */
class ModelAbstract
{
protected:
    const std::size_t mDim;          //! the dimension of state
    const std::size_t mNPath;        //! the number of Path
    const MarketData::Terms mTerms;  //! term structure
    const Math::Vec mInitState;      //! initial spot state
    /**
     * @brief The coefficient of dt in SDE of r[inIndPath][inIndTerm]
     * @param inIndPath the index of path
     * @param inIndTerm the index of term
     * @return double the coefficient
     */
    virtual Math::Vec driftCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<Math::Vec>>& inSpots ) const = 0;

    virtual double transfStateToRate( const Math::Vec& inState,
                                      double inTime ) const;

public:
    /**
     * @brief This constructs a ModelAbstract.
     * @param inNPath the number of Path
     * @param insTerms term structure
     * @param inInitState the initial state
     */
    ModelAbstract( std::size_t inNPath, const MarketData::Terms& inTerms,
                   const Math::Vec& inInitState ) :
        mDim( inInitState.size() ),
        mNPath( inNPath ),
        mTerms( inTerms ),
        mInitState( inInitState )
    {
    }
    virtual ~ModelAbstract() = default;
    /**
     * @brief This calcurate spot rates and Disconunt Factors.
     */
    virtual MarketData::SpotRates calcSpotRates() const;
};

/**
 * @brief This build the object of ModelAbstract.
 */
class ModelAbstractBuilder
{
protected:
    Math::Vec mInitState = Math::Vec( 0 );
    std::size_t mNPath;                          //! the number of Path
    std::unique_ptr<MarketData::Terms> muTerms;  //! term structure

public:
    ModelAbstractBuilder& setNPath( std::size_t inNPath )
    {
        mNPath = inNPath;
        return *this;
    }
    ModelAbstractBuilder& setTerms(
        std::shared_ptr<const std::vector<double>> insTerms )
    {
        muTerms = std::make_unique<MarketData::Terms>( insTerms );
        return *this;
    }
    ModelAbstractBuilder& setTerms( const MarketData::Terms& inTerms )
    {
        muTerms = std::make_unique<MarketData::Terms>( inTerms );
        return *this;
    }
    ModelAbstractBuilder& setInitState( Math::Vec inInitState )
    {
        mInitState = inInitState;
        return *this;
    }
    virtual ~ModelAbstractBuilder() = default;
};

/**
 * @brief This is (test) class of constant short-rate model.
 */
class ConstantRate : public ModelAbstract
{
private:
    Math::Vec driftCoeff(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<Math::Vec>>& inSpots ) const override;

public:
    ConstantRate( const MarketData::Terms& inTerms,
                  const Math::Vec& inInitState ) :
        ModelAbstract( 1, inTerms, inInitState )
    {
    }
};

/**
 * @brief This is abstract class for short-rate models generated by
 * one-dimensional Brownian motion.
 */
class MultiFactorAbstract : public ModelAbstract
{
protected:
    std::unique_ptr<Process::RandomVec::PathAbstract>
        muRandomPath;  //! random path
    /**
     * @brief The coefficient of dW in SDE of r[inIndPath][inIndTerm]
     * @param inIndPath the index of path
     * @param inIndTerm the index of term
     * @return double the coefficient
     */
    virtual Math::Vec volTerm(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<Math::Vec>>& inSpots,
        const Math::Vec& inRandomVec ) const = 0;

public:
    MultiFactorAbstract(
        std::size_t inNPath, const MarketData::Terms& inTerms,
        const Math::Vec& inInitState,
        std::unique_ptr<Process::RandomVec::PathAbstract> inuRandomPath ) :
        ModelAbstract( inNPath, inTerms, inInitState ),
        muRandomPath( std::move( inuRandomPath ) )
    {
    }
    virtual MarketData::SpotRates calcSpotRates() const override;
};

class MultiFactorAbstractBuilder : public ModelAbstractBuilder
{
protected:
    std::unique_ptr<Process::RandomVec::PathAbstract>
        muRandomPath;  //! random path
public:
    MultiFactorAbstractBuilder& setRandom(
        std::unique_ptr<Process::RandomVec::PathAbstract> inuRandomPath )
    {
        muRandomPath = std::move( inuRandomPath );
        return *this;
    }
};

}  // namespace ShortRateMCMulti
}  // namespace Process

#endif