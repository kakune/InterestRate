/**
 * @file core.hpp
 * @brief This defines short rate path classes.
 * @author kakune
 * @date 1/29/2024
 */

#ifndef PROCESS_SHORT_RATE_CORE_HPP
#define PROCESS_SHORT_RATE_CORE_HPP

#include <memory>
#include <vector>

#include "math/interpolate_1d.hpp"
#include "process/market.hpp"
#include "process/random.hpp"

namespace Process
{
namespace ShortRate
{

/**
 * @brief This is the abstract class for short rate models.
 */
class ModelAbstract
{
protected:
    std::size_t mNPath;                                   //! the number of Path
    std::shared_ptr<const std::vector<double> > msTerms;  //! term structure
    std::vector<std::vector<double> > mSpotRates;  //! calcurated interest rate
    std::vector<std::vector<double> > mDFs;  //! calcurated discount factor
    std::vector<double>
        mExpectedSpotRates;     //! calcurated expectation value of spot rate
    std::vector<double> mZCBs;  //! price of zero-coupon bond
    Math::Interpolate1d::NewtonSpline
        mInterpZCB;        //! interpolated function of ZCB
    double mInitSpotRate;  //! initial spot rate
    std::shared_ptr<const Market::Data> msMarketData;  //! the data of market

    /**
     * @brief The coefficient of dt in SDE of r[inIndPath][inIndTerm]
     * @param inIndPath the index of path
     * @param inIndTerm the index of term
     * @return double the coefficient
     */
    virtual double driftCoeff( std::size_t inIndPath,
                               std::size_t inIndTerm ) const = 0;
    virtual ~ModelAbstract()                                 = default;
    /**
     * @brief This calcurate spot rates and Disconunt Factors.
     */
    virtual void buildSpot();
    /**
     * @brief This calcurate ExpectedSpotRate and ZCB.
     */
    virtual void buildZCB();

public:
    /**
     * @brief This constructs a ModelAbstract.
     * @param inNPath the number of Path
     * @param insTerms term structure
     */
    ModelAbstract( std::size_t inNPath,
                   std::shared_ptr<const std::vector<double> > insTerms,
                   std::shared_ptr<const Market::Data> insMarketData,
                   double inInitSpotRate ) :
        mNPath( inNPath ),
        msTerms( insTerms ),
        mSpotRates( mNPath, std::vector<double>( insTerms->size(), 0 ) ),
        mDFs( mNPath, std::vector<double>( insTerms->size(), 1.0 ) ),
        mExpectedSpotRates( insTerms->size(), 0.0 ),
        mZCBs( insTerms->size(), 1.0 ),
        mInterpZCB( 3 ),
        mInitSpotRate(
            insMarketData == nullptr
                ? inInitSpotRate
                : insMarketData->mInterpInstantaneousForwardRate(
                      0.95 * insTerms->at( 0 ) + 0.05 * insTerms->at( 1 ) ) ),
        msMarketData( insMarketData )
    {
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            mSpotRates.at( iPath ).at( 0 ) = mInitSpotRate;
        }
    }
    /**
     * @brief This calcurate interest rates of each random path.
     */
    virtual void build();

    /**
     * @brief This calculates ZCB starting at
     * Terms[0].
     * @param inMaturityTime maturity time of ZCB
     * @return double P(inStartTime, inMaturityTime)
     */
    virtual double priceZCB( double inMaturityTime ) const;
    /**
     * @brief This calculates ZCB price within arbitrary interval observed at
     * Terms[0].
     * @param inStartTime start time of ZCB
     * @param inMaturityTime maturity time of ZCB
     * @return double P(inStartTime, inMaturityTime)
     */
    double priceZCB( double inStartTime, double inMaturityTime ) const;
    double priceZCB( std::size_t inIndStartTime,
                     std::size_t inIndMaturityTime ) const;

    /**
     * @brief This calculates forward rate within arbitary interval observed at
     * Terms[0].
     * @param inStartTime start time of forward rate
     * @param inTerminalTime terminal time of forward rate
     * @return double f(inStartTime, inTerminalTime)
     */
    virtual double forwardRate( double inStartTime,
                                double inTerminalTime ) const;
    double forwardRate( std::size_t inIndStartTime,
                        std::size_t inIndTerminalTime ) const;

    /**
     * @brief This calculates instantaneous forward rate at arbitary time
     * observed at Terms[0].
     * @param inTime time of forward rate
     * @return double f(inTime)
     */
    virtual double instantaneousForwardRate( double inTime ) const;
    double instantaneousForwardRate( std::size_t inIndTime ) const;
};

/**
 * @brief This build the object of ModelAbstract.
 */
class ModelAbstractBuilder
{
protected:
    double mInitSpotRate;
    std::size_t mNPath;                                   //! the number of Path
    std::shared_ptr<const std::vector<double> > msTerms;  //! term structure
    std::shared_ptr<const Market::Data> msMarketData = nullptr;

public:
    ModelAbstractBuilder& setNPath( std::size_t inNPath )
    {
        mNPath = inNPath;
        return *this;
    }
    ModelAbstractBuilder& setTerms(
        std::shared_ptr<const std::vector<double> > insTerms )
    {
        msTerms = insTerms;
        return *this;
    }
    ModelAbstractBuilder& setMarketData(
        std::shared_ptr<const Market::Data> insMarketData )
    {
        msMarketData = insMarketData;
        return *this;
    }
    ModelAbstractBuilder& setInitSpotRate( double inInitSpotRate )
    {
        mInitSpotRate = inInitSpotRate;
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
    double driftCoeff( std::size_t inIndPath,
                       std::size_t inIndTerm ) const override;

public:
    ConstantRate( std::shared_ptr<const std::vector<double> > insTerms,
                  std::shared_ptr<const Market::Data> insMarketData,
                  double inInitSpotRate ) :
        ModelAbstract( 1, insTerms, insMarketData, inInitSpotRate )
    {
    }
};

/**
 * @brief This is abstract class for short-rate models generated by
 * one-dimensional Brownian motion.
 */
class OneFactorAbstract : public ModelAbstract
{
protected:
    std::unique_ptr<Process::Random::PathAbstract>
        muRandomPath;  //! random path
    /**
     * @brief The coefficient of dW in SDE of r[inIndPath][inIndTerm]
     * @param inIndPath the index of path
     * @param inIndTerm the index of term
     * @return double the coefficient
     */
    virtual double volCoeff( std::size_t inIndPath,
                             std::size_t inIndTerm ) const = 0;
    virtual void buildSpot() override;

public:
    OneFactorAbstract(
        std::size_t inNPath,
        std::shared_ptr<const std::vector<double> > insTerms,
        std::shared_ptr<const Market::Data> insMarketData,
        double inInitSpotRate,
        std::unique_ptr<Process::Random::PathAbstract> inuRandomPath ) :
        ModelAbstract( inNPath, insTerms, insMarketData, inInitSpotRate ),
        muRandomPath( std::move( inuRandomPath ) )
    {
    }
};

class OneFactorAbstractBuilder : public ModelAbstractBuilder
{
protected:
    std::unique_ptr<Process::Random::PathAbstract>
        muRandomPath;  //! random path
public:
    OneFactorAbstractBuilder& setRandom(
        std::unique_ptr<Process::Random::PathAbstract> inuRandomPath )
    {
        muRandomPath = std::move( inuRandomPath );
        return *this;
    }
};

}  // namespace ShortRate
}  // namespace Process

#endif