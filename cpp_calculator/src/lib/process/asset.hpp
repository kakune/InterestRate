/**
 * @file asset.hpp
 * @brief This defines calculator of assets.
 * @author kakune
 * @date 1/30/2024
 */

#ifndef PROCESS_ASSET_HPP
#define PROCESS_ASSET_HPP

#include <cmath>
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
     * @brief This computes the price of put option WITHOUT discount factor.
     * @param inStrike the price of strike
     * @param inIndTime the index of time of maturity
     * @return double price of call option
     */
    double pricePutOption( double inStrike, std::size_t inIndTime );
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
    virtual void calcEachForwardPrice() override;

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
    std::unique_ptr< Process::Random::PathAbstract > muRandomPath =
        nullptr;  //! random path
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
        if ( muRandomPath == nullptr )
        {
            muRandomPath =
                std::make_unique< Process::Random::PathBrownAntithetic >(
                    mNPath, msTerms );
        }
        return BlackSholesForward( mNPath, msTerms, mInitPrice,
                                   std::move( muRandomPath ), mVol );
    }
};

/**
 * @brief This is the abstract class for Local Volatility Model of forward
 * price.
 */
class StochasticVolatilityForwardAbstract : public ModelForwardAbstract
{
protected:
    double mInitVol;  //! initial volatility
    double mCorr;     //! correlation between brownians
    double mAuxCorr;  //! sqrt{1 - mCorr^2}
    std::vector< std::vector< double > >
        mVolatility;  //! calcurated forward price
    std::unique_ptr< Process::Random::PathAbstract >
        muRandomPath;  //! random path
    std::unique_ptr< Process::Random::PathAbstract >
        muRandomVol;  //! random path for vol

    /**
     * @brief This is the coefficient of Brownian motion (dW) in forward price
     * SDE.
     * @param inPrice current forward price
     * @param inVol current volatility
     * @param inTime current time
     * @return double
     */
    virtual double volForward( double inPrice, double inVol,
                               double inTime ) = 0;
    /**
     * @brief This is the coefficient of Brownian motion (dW) in volatility SDE.
     * @param inVol current volatility
     * @param inTime current time
     * @return double
     */
    virtual double volVolatility( double inVol, double inTime ) = 0;
    virtual void calcEachForwardPrice() override;

public:
    StochasticVolatilityForwardAbstract(
        std::size_t inNPath,
        std::shared_ptr< const std::vector< double > > insTerms,
        double inInitPrice,
        std::unique_ptr< Process::Random::PathAbstract > inuRandomPath,
        double inInitVol,
        std::unique_ptr< Process::Random::PathAbstract > inuRandomVol,
        double inCorr ) :
        ModelForwardAbstract( inNPath, insTerms, inInitPrice ),
        mInitVol( inInitVol ),
        mVolatility( mNPath,
                     std::vector< double >( insTerms->size(), inInitVol ) ),
        muRandomPath( std::move( inuRandomPath ) ),
        muRandomVol( std::move( inuRandomVol ) ),
        mCorr( inCorr ),
        mAuxCorr( std::sqrt( 1.0 - inCorr * inCorr ) )
    {
    }
    virtual ~StochasticVolatilityForwardAbstract() = default;
};

class StochasticVolatilityWithLogForwardAbstract
    : public StochasticVolatilityForwardAbstract
{
protected:
    virtual void calcEachForwardPrice() override;
    /**
     * @brief This is the coefficient of drift (dt) in forward price SDE.
     * @param inPrice current forward price
     * @param inVol current volatility
     * @param inTime current time
     * @return double
     */
    virtual double driftForward( double inPrice, double inVol,
                                 double inTime ) = 0;
    /**
     * @brief This is the coefficient of drift (dt) in forward price SDE.
     * @param inVol current volatility
     * @param inTime current time
     * @return double
     */
    virtual double driftVolatility( double inVol, double inTime ) = 0;

public:
    using StochasticVolatilityForwardAbstract::
        StochasticVolatilityForwardAbstract;
};

/**
 * @brief This build the object of StochasticVolatilityForward
 */
class StochasticVolatilityForwardAbstractBuilder
    : public ModelForwardAbstractBuilder
{
protected:
    double mInitVol;  //! initial volatility
    double mCorr;     //! correlation between brownians
    std::unique_ptr< Process::Random::PathAbstract > muRandomPath =
        nullptr;  //! random path
    std::unique_ptr< Process::Random::PathAbstract > muRandomVol =
        nullptr;  //! random path for vol
public:
    StochasticVolatilityForwardAbstractBuilder& setInitVol( double inInitVol )
    {
        mInitVol = inInitVol;
        return *this;
    }
    StochasticVolatilityForwardAbstractBuilder& setCorr( double inCorr )
    {
        mCorr = inCorr;
        return *this;
    }
    StochasticVolatilityForwardAbstractBuilder& setRandomPath(
        std::unique_ptr< Process::Random::PathAbstract > inuRandomPath )
    {
        muRandomPath = std::move( inuRandomPath );
        return *this;
    }
    StochasticVolatilityForwardAbstractBuilder& setRandomVol(
        std::unique_ptr< Process::Random::PathAbstract > inuRandomVol )
    {
        muRandomVol = std::move( inuRandomVol );
        return *this;
    }
};

/**
 * @brief This is calculator for SABR.
 */
class SABRForward : public StochasticVolatilityForwardAbstract
{
private:
    double mExponent;  //! exponent of forward price in SDE
    double mVolvol;    //! volvol
    double volForward( double inPrice, double inVol,
                       double inTime = 0.0 ) override;
    double volVolatility( double inVol, double inTime = 0.0 ) override;

public:
    SABRForward( std::size_t inNPath,
                 std::shared_ptr< const std::vector< double > > insTerms,
                 double inInitPrice,
                 std::unique_ptr< Process::Random::PathAbstract > inuRandomPath,
                 double inInitVol,
                 std::unique_ptr< Process::Random::PathAbstract > inuRandomVol,
                 double inCorr, double inExponent, double inVolvol ) :
        StochasticVolatilityForwardAbstract(
            inNPath, insTerms, inInitPrice, std::move( inuRandomPath ),
            inInitVol, std::move( inuRandomVol ), inCorr ),
        mExponent( inExponent ),
        mVolvol( inVolvol )
    {
        calcEachForwardPrice();
    }
};

/**
 * @brief This build the object of SABRForward
 */
class SABRForwardBuilder : public StochasticVolatilityForwardAbstractBuilder
{
private:
    double mExponent;  //! exponent of forward price in SDE
    double mVolvol;    //! volvol

public:
    SABRForwardBuilder& setExponent( double inExponent )
    {
        mExponent = inExponent;
        return *this;
    }
    SABRForwardBuilder& setVolvol( double inVolvol )
    {
        mVolvol = inVolvol;
        return *this;
    }
    SABRForward build()
    {
        if ( muRandomPath == nullptr )
        {
            muRandomPath =
                std::make_unique< Process::Random::PathBrownAntithetic >(
                    mNPath, msTerms );
        }
        if ( muRandomVol == nullptr )
        {
            muRandomVol =
                std::make_unique< Process::Random::PathBrownAntithetic >(
                    mNPath, msTerms );
        }
        return SABRForward(
            mNPath, msTerms, mInitPrice, std::move( muRandomPath ), mInitVol,
            std::move( muRandomVol ), mCorr, mExponent, mVolvol );
    }
};

/**
 * @brief This is calculator for SABR using SDE of log(forward price) and
 * log(volatility).
 */
class SABRWithLogForward : public StochasticVolatilityWithLogForwardAbstract
{
private:
    double mExponent;  //! exponent of forward price in SDE
    double mVolvol;    //! volvol
    double volForward( double inPrice, double inVol,
                       double inTime = 0.0 ) override;
    double volVolatility( double inVol, double inTime = 0.0 ) override;
    double driftForward( double inPrice, double inVol,
                         double inTime = 0.0 ) override;
    double driftVolatility( double inVol, double inTime = 0.0 ) override;

public:
    SABRWithLogForward(
        std::size_t inNPath,
        std::shared_ptr< const std::vector< double > > insTerms,
        double inInitPrice,
        std::unique_ptr< Process::Random::PathAbstract > inuRandomPath,
        double inInitVol,
        std::unique_ptr< Process::Random::PathAbstract > inuRandomVol,
        double inCorr, double inExponent, double inVolvol ) :
        StochasticVolatilityWithLogForwardAbstract(
            inNPath, insTerms, inInitPrice, std::move( inuRandomPath ),
            inInitVol, std::move( inuRandomVol ), inCorr ),
        mExponent( inExponent ),
        mVolvol( inVolvol )
    {
        calcEachForwardPrice();
    }
};

/**
 * @brief This build the object of SABRWithLogForward
 */
class SABRWithLogForwardBuilder
    : public StochasticVolatilityForwardAbstractBuilder
{
private:
    double mExponent;  //! exponent of forward price in SDE
    double mVolvol;    //! volvol

public:
    SABRWithLogForwardBuilder& setExponent( double inExponent )
    {
        mExponent = inExponent;
        return *this;
    }
    SABRWithLogForwardBuilder& setVolvol( double inVolvol )
    {
        mVolvol = inVolvol;
        return *this;
    }
    SABRWithLogForward build()
    {
        if ( muRandomPath == nullptr )
        {
            muRandomPath =
                std::make_unique< Process::Random::PathBrownAntithetic >(
                    mNPath, msTerms );
        }
        if ( muRandomVol == nullptr )
        {
            muRandomVol =
                std::make_unique< Process::Random::PathBrownAntithetic >(
                    mNPath, msTerms );
        }
        return SABRWithLogForward(
            mNPath, msTerms, mInitPrice, std::move( muRandomPath ), mInitVol,
            std::move( muRandomVol ), mCorr, mExponent, mVolvol );
    }
};

}  // namespace Asset
}  // namespace Process

#endif