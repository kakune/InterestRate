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
    std::size_t mNPath;  //! the number of Path
    std::shared_ptr< const std::vector< double > > msTerms;  //! term structure
    std::unique_ptr< Process::Random::PathAbstract >
        muRandomPath;  //! random path
    std::vector< std::vector< double > >
        mSpotRates;              //! calcurated interest rate
    std::vector< double > mZCB;  //! price of zero-coupon bond
    Math::Interpolate1d::NewtonSpline
        mInterpZCB;  //! interpolated function of ZCB
    /**
     * @brief This calcurate interest rates of each random path.
     */
    virtual void calcEachRates() = 0;
    virtual ~ModelAbstract()     = default;

public:
    /**
     * @brief This constructs a ModelAbstract.
     * @param inNPath the number of Path
     * @param insTerms term structure
     * @param inuRandomPath unique_ptr of class deriving
     * Process::Random::PathAbstract
     */
    ModelAbstract(
        std::size_t inNPath,
        std::shared_ptr< const std::vector< double > > insTerms,
        std::unique_ptr< Process::Random::PathAbstract > inuRandomPath ) :
        mNPath( inNPath ),
        msTerms( insTerms ),
        muRandomPath( std::move( inuRandomPath ) ),
        mSpotRates( mNPath, std::vector< double >( insTerms->size(), 0 ) ),
        mZCB( insTerms->size(), 1.0 ),
        mInterpZCB( 3 )
    {
    }
    /**
     * @brief This calculates ZCB price.
     */
    void calcZCB();
    /**
     * @brief This calculates ZCB price within arbitrary interval observed at
     * Terms[0].
     * @param inStartTime start time of ZCB
     * @param inMaturityTime maturity time of ZCB
     * @return double P(inStartTime, inMaturityTime)
     */
    double priceZCB( double inStartTime, double inMaturityTime );

    /**
     * @brief This calculates forward rate within arbitary interval observed at
     * Terms[0].
     *
     * @param inStartTime start time of forward rate
     * @param inTerminalTime terminal time of forward rate
     * @return double f(inStartTime, inTerminalTime)
     */
    double forwardRate( double inStartTime, double inTerminalTime );
};

/**
 * @brief This build the object of ModelAbstract
 */
class ModelAbstractBuilder
{
protected:
    std::size_t mNPath;  //! the number of Path
    std::shared_ptr< const std::vector< double > > msTerms;  //! term structure
    std::unique_ptr< Process::Random::PathAbstract >
        muRandomPath;  //! random path
public:
    ModelAbstractBuilder& setNPath( std::size_t inNPath )
    {
        mNPath = inNPath;
        return *this;
    }
    ModelAbstractBuilder& setTerms(
        std::shared_ptr< const std::vector< double > > insTerms )
    {
        msTerms = insTerms;
        return *this;
    }
    ModelAbstractBuilder& setRandomPath(
        std::unique_ptr< Process::Random::PathAbstract > inuRandomPath )
    {
        muRandomPath = std::move( inuRandomPath );
        return *this;
    }
    virtual ~ModelAbstractBuilder() = default;
};

}  // namespace ShortRate
}  // namespace Process

#endif