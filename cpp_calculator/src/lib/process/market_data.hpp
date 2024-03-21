/**
 * @file interpolate_1d.cpp
 * @brief This defines classes for market data. Unless there is a special
 * reason, calculation results and actual market data should be exchanged
 * through these classes. Since most data is stored by smart pointers, there is
 * no need to worry about the cost of copying instances.
 * @author kakune
 * @date 3/21/2024
 */

#ifndef PROCESS_MARKET_DATA_HPP
#define PROCESS_MARKET_DATA_HPP

#include <cmath>
#include <memory>
#include <vector>

#include "math/interpolate_1d.hpp"

namespace Process
{
namespace MarketData
{

/**
 * @brief This stores the tenor.
 */
class Terms
{
private:
    const std::shared_ptr<const std::vector<double>> msData;

public:
    /**
     * @brief This constructs a new Terms.
     * @param insData The term structure. It must be sorted.
     */
    Terms( std::shared_ptr<const std::vector<double>> insData );
    Terms( std::vector<double> inData );
    double operator[]( std::size_t inIndex ) const;
    double at( std::size_t inIndex ) const;
    std::size_t size() const;
    const std::shared_ptr<const std::vector<double>> ptr() const;
    double front() const;
    double back() const;
};

/**
 * @brief This stores spot rate data at each path and each term. It
 * automatically calculate Discount factors.
 */
class SpotRates
{
private:
    std::size_t mNPath;  //! the number of path
    const Terms mTerms;  //! tenor
    const std::shared_ptr<const std::vector<std::vector<double>>>
        msDataSpotRate;  //! spot rates
    const std::shared_ptr<const std::vector<std::vector<double>>>
        msDataDF;  //! discount factor

    /**
     * @brief This calculates DF from spot rates.
     * @param inTerms
     * @param insDataSpotRate
     * @return std::vector<std::vector<double>> DFs
     */
    std::vector<std::vector<double>> calcDFFromSpotRate(
        const Terms& inTerms,
        std::shared_ptr<const std::vector<std::vector<double>>>
            insDataSpotRate );

public:
    /**
     * @brief This constructs a new Spot Rates.
     * @param inTerms Term structure
     * @param insDataSpotRate Spot rate
     */
    SpotRates( const Terms& inTerms,
               std::shared_ptr<const std::vector<std::vector<double>>>
                   insDataSpotRate );
    SpotRates( const Terms& inTerms,
               std::vector<std::vector<double>> inDataSpotRate );
    const std::vector<double>& operator[]( std::size_t inIndex ) const;
    double term( std::size_t inIndex ) const;
    const Terms& getTerms() const;
    const std::vector<std::vector<double>>& getDF() const;
    std::size_t sizeTerms() const;
    std::size_t sizePath() const;
};

/**
 * @brief This stores prices of ZCB at each term.
 */
class ZCB
{
private:
    const Terms mTerms;
    const std::shared_ptr<const std::vector<double>> msData;
    const std::shared_ptr<Math::Interpolate1D::NewtonSpline> msSpline;

    /**
     * @brief This calculates ZCB from (Monte-Carlo) spot rates.
     * @param inSpotRates
     * @return std::vector<double>
     */
    std::vector<double> calcZCBFromSpotRates(
        const SpotRates& inSpotRates ) const;

public:
    ZCB( const Terms& inTerms,
         std::shared_ptr<const std::vector<double>> insData,
         std::size_t inDeg = 3 );
    ZCB( const Terms& inTerms, std::vector<double> inData,
         std::size_t inDeg = 3 );
    ZCB( const SpotRates& inSpotRate, std::size_t inDeg = 3 );
    double operator[]( std::size_t inIndex ) const;
    double term( std::size_t inIndex ) const;
    const Terms& getTerms() const;
    std::size_t sizeTerms() const;
    /**
     * @brief This calculates ZCB starting at
     * Terms[0].
     * @param inMaturityTime maturity time of ZCB
     * @return double P(inStartTime, inMaturityTime)
     */
    double operator()( double inTime ) const;
    /**
     * @brief This calculates ZCB price within arbitrary interval observed at
     * Terms[0].
     * @param inStartTime start time of ZCB
     * @param inMaturityTime maturity time of ZCB
     * @return double P(Terms[0], inMaturityTime)
     */
    double operator()( double inStartTime, double inMaturityTime ) const;
    /**
     * @brief This calculates forward rate within arbitary interval observed at
     * Terms[0].
     * @param inStartTime start time of forward rate
     * @param inTerminalTime terminal time of forward rate
     * @return double f(inStartTime, inTerminalTime)
     */
    double forwardRate( double inStartTime, double inTerminalTime ) const;
    /**
     * @brief This calculates instantaneous forward rate at arbitary time
     * observed at Terms[0].
     * @param inTime time of forward rate
     * @return double f(inTime)
     */
    double instantaneousForwardRate( double inTime ) const;
    /**
     * @brief This calculates 1st derivative of instantaneous forward rate at
     * arbitary time observed at Terms[0].
     * @param inTime time of forward rate
     * @return double df(inTime)/dT
     */
    double derivInstantaneousForwardRate( double inTime ) const;
    /**
     * @brief This gives the spot rate at term[0].
     * @return double r(Term[0])
     */
    double initSpotRate() const;
};

;
}  // namespace MarketData
}  // namespace Process

#endif