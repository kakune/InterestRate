/**
 * @file asset.hpp
 * @brief This defines calculator of assets.
 * @author kakune
 * @date 1/30/2024
 */

#ifndef PROCESS_ASSET_HPP
#define PROCESS_ASSET_HPP

#include <memory>
#include <vector>

#include "math/interpolate_1d.hpp"
#include "process/random.hpp"

namespace Process
{
namespace Asset
{

/**
 * @brief This is the abstract class for SDE Model of forward price.
 */
class ModelForwardAbstract
{
protected:
    std::size_t mNPath;  //! the number of Path
    std::shared_ptr< const std::vector< double > > msTerms;  //! term structure
    double mInitPrice;
    std::vector< std::vector< double > >
        mForwardPrice;  //! calcurated forward price
    /**
     * @brief This calcurate interest rates of each random path.
     */
    virtual void calcEachForwardPrice() = 0;

public:
    /**
     * @brief This constructs a ModelAbstract.
     * @param inNPath the number of Path
     * @param insTerms term structure
     * @param inInitPrice the initial price of forward price
     */
    ModelForwardAbstract(
        std::size_t inNPath,
        std::shared_ptr< const std::vector< double > > insTerms,
        double inInitPrice ) :
        mNPath( inNPath ),
        msTerms( insTerms ),
        mInitPrice( inInitPrice ),
        mForwardPrice( mNPath,
                       std::vector< double >( insTerms->size(), inInitPrice ) )
    {
    }
    virtual ~ModelForwardAbstract() = default;
    /**
     * @brief This computes the price of call option WITHOUT discount factor.
     * @param inStrike the price of strike
     * @param inIndTime the index of time of maturity
     * @return double price of call option
     */
    double priceCallOption( double inStrike, std::size_t inIndTime );
    /**
     * @brief This computes the implied volatility.
     * @param inStrike the price of strike
     * @param inIndTime the index of time of maturity
     * @return double implied volatility
     */
    double impliedVolatility( double inStrike, std::size_t inIndTime );
};

/**
 * @brief This build the object of ModelAbstract
 */
class ModelForwardAbstractBuilder
{
protected:
    std::size_t mNPath;  //! the number of Path
    std::shared_ptr< const std::vector< double > > msTerms;  //! term structure
    double mInitPrice;  //! initial forward price

public:
    ModelForwardAbstractBuilder& setNPath( std::size_t inNPath )
    {
        mNPath = inNPath;
        return *this;
    }
    ModelForwardAbstractBuilder& setTerms(
        std::shared_ptr< const std::vector< double > > insTerms )
    {
        msTerms = insTerms;
        return *this;
    }
    ModelForwardAbstractBuilder& setInitPrice( double inInitPrice )
    {
        mInitPrice = inInitPrice;
        return *this;
    }
    virtual ~ModelForwardAbstractBuilder() = default;
};

/**
 * @brief This is the abstract class for Local Volatility Model of forward
 * price.
 */
class LocalVolatilityForwardAbstract : public ModelForwardAbstract
{
protected:
    std::unique_ptr< Process::Random::PathAbstract >
        muRandomPath;  //! random path
    /**
     * @brief This is deterministic function of local volatility.
     * @param inPrice current forward price
     * @param inTime current time
     * @return double
     */
    virtual double localVolatility( double inPrice, double inTime ) = 0;
    void calcEachForwardPrice() override;

public:
    LocalVolatilityForwardAbstract(
        std::size_t inNPath,
        std::shared_ptr< const std::vector< double > > insTerms,
        double inInitPrice,
        std::unique_ptr< Process::Random::PathAbstract > inuRandomPath ) :
        ModelForwardAbstract( inNPath, insTerms, inInitPrice ),
        muRandomPath( std::move( inuRandomPath ) )
    {
    }
    virtual ~LocalVolatilityForwardAbstract() = default;
};

/**
 * @brief This build the object of LocalVolatilityForward
 */
class LocalVolatilityForwardAbstractBuilder : public ModelForwardAbstractBuilder
{
protected:
    std::unique_ptr< Process::Random::PathAbstract >
        muRandomPath;  //! random path
public:
    LocalVolatilityForwardAbstractBuilder& setRandomPath(
        std::unique_ptr< Process::Random::PathAbstract > inuRandomPath )
    {
        muRandomPath = std::move( inuRandomPath );
        return *this;
    }
};

/**
 * @brief This is SDE calculator for Black-Sholes model.
 */
class BlackSholesForward : public LocalVolatilityForwardAbstract
{
private:
    double mVol;  //! constant volatility
    double localVolatility( double inPrice, double inTime = 0.0 ) override;

public:
    BlackSholesForward(
        std::size_t inNPath,
        std::shared_ptr< const std::vector< double > > insTerms,
        double inInitPrice,
        std::unique_ptr< Process::Random::PathAbstract > inuRandomPath,
        double inVol ) :
        LocalVolatilityForwardAbstract( inNPath, insTerms, inInitPrice,
                                        std::move( inuRandomPath ) ),
        mVol( inVol )
    {
        calcEachForwardPrice();
    }
};

/**
 * @brief This build the object of BlackSholesForward
 */
class BlackSholesForwardBuilder : public LocalVolatilityForwardAbstractBuilder
{
private:
    double mVol;  //! constant volatility

public:
    BlackSholesForwardBuilder& setVol( double inVol )
    {
        mVol = inVol;
        return *this;
    }
    BlackSholesForward build()
    {
        return BlackSholesForward( mNPath, msTerms, mInitPrice,
                                   std::move( muRandomPath ), mVol );
    }
};

}  // namespace Asset
}  // namespace Process

#endif