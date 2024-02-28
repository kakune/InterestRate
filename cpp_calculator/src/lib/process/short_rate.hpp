/**
 * @file short_rate.hpp
 * @brief This defines short rate path classes.
 * @author kakune
 * @date 1/29/2024
 */

#ifndef PROCESS_SHORT_RATE_HPP
#define PROCESS_SHORT_RATE_HPP

#include <memory>
#include <vector>

#include "math/interpolate_1d.hpp"
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

    virtual double driftCoeff( std::size_t inIndPath,
                               std::size_t inIndTerm ) = 0;
    virtual ~ModelAbstract()                           = default;

public:
    /**
     * @brief This constructs a ModelAbstract.
     * @param inNPath the number of Path
     * @param insTerms term structure
     */
    ModelAbstract( std::size_t inNPath,
                   std::shared_ptr<const std::vector<double> > insTerms,
                   double inInitSpotRate ) :
        mNPath( inNPath ),
        msTerms( insTerms ),
        mSpotRates( mNPath, std::vector<double>( insTerms->size(), 0 ) ),
        mDFs( mNPath, std::vector<double>( insTerms->size(), 1.0 ) ),
        mExpectedSpotRates( insTerms->size(), 0.0 ),
        mZCBs( insTerms->size(), 1.0 ),
        mInterpZCB( 3 ),
        mInitSpotRate( inInitSpotRate )
    {
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            mSpotRates.at( 0 ).at( iPath ) = inInitSpotRate;
        }
    }
    /**
     * @brief This calcurate interest rates of each random path.
     */
    virtual void build();
    /**
     * @brief This calculates ZCB price within arbitrary interval observed at
     * Terms[0].
     * @param inStartTime start time of ZCB
     * @param inMaturityTime maturity time of ZCB
     * @return double P(inStartTime, inMaturityTime)
     */
    virtual double priceZCB( double inStartTime, double inMaturityTime );

    /**
     * @brief This calculates forward rate within arbitary interval observed at
     * Terms[0].
     * @param inStartTime start time of forward rate
     * @param inTerminalTime terminal time of forward rate
     * @return double f(inStartTime, inTerminalTime)
     */
    virtual double forwardRate( double inStartTime, double inTerminalTime );

    /**
     * @brief This calculates instantaneous forward rate at arbitary time
     * observed at Terms[0].
     * @param inTime time of forward rate
     * @return double f(inTime)
     */
    virtual double instantaneousForwardRate( double inTime );
};

/**
 * @brief This build the object of ModelAbstract
 */
class ModelAbstractBuilder
{
protected:
    double mInitSpotRate;
    std::size_t mNPath;                                   //! the number of Path
    std::shared_ptr<const std::vector<double> > msTerms;  //! term structure
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
    ModelAbstractBuilder& setInitSpotRate( double inInitSpotRate )
    {
        mInitSpotRate = inInitSpotRate;
        return *this;
    }
    virtual ~ModelAbstractBuilder() = default;
};

class ConstantRate : public ModelAbstract
{
private:
    double driftCoeff( std::size_t inIndPath, std::size_t inIndTerm ) override;

public:
    ConstantRate( std::shared_ptr<const std::vector<double> > insTerms,
                  double inInitSpotRate ) :
        ModelAbstract( 1, insTerms, inInitSpotRate )
    {
    }
};

}  // namespace ShortRate
}  // namespace Process

#endif