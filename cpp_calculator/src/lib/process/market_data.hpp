/**
 * @file market_data.hpp
 * @brief This defines classes for market data. Unless there is a special
 * reason, actual or model market data should be exchanged
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

#include "analytical/Black76.hpp"
#include "math/interpolate_1d.hpp"
#include "math/matrix.hpp"

namespace Process
{
namespace MarketData
{

/**
 * @brief This stores the term structure {t_0, t_1, ... t_{N-1}}.
 */
class Terms
{
private:
    const std::shared_ptr<const std::vector<double>> msData;  //! term structure
    const std::shared_ptr<const std::vector<double>> msDifTime, msSqrtDifTime;

public:
    /**
     * @brief This constructs a new Terms.
     * @param insData The term structure. It must be sorted.
     */
    Terms( std::shared_ptr<const std::vector<double>> insData );
    Terms( std::vector<double> inData );
    /**
     * @brief This returns t_{inIndex}.
     */
    double operator[]( std::size_t inIndex ) const;
    double at( std::size_t inIndex ) const;

    /**
     * @brief This returns t_{inIndex} - t_{inIndex-1}.
     */
    double difTime( std::size_t inIndex ) const;
    /**
     * @brief This returns sqrt(t_{inIndex} - t_{inIndex-1}).
     */
    double sqrtDifTime( std::size_t inIndex ) const;
    /**
     * @brief This returns the number of terms N.
     */
    std::size_t size() const;
    const std::shared_ptr<const std::vector<double>> ptr() const;
    double front() const;
    double back() const;
};

/**
 * @brief This stores Tenor structure {T_0, T_1, ... T_N} and {tau_0, tau_1, ...
 * tau_{N-1}}.
 *
 */
class Tenor
{
private:
    const std::size_t mNSize;  //! the size of tenor N.
    const Terms mTerms;
    const std::shared_ptr<const std::vector<std::size_t>> msData;
    const std::shared_ptr<const Math::Vec> msTau;
    const std::shared_ptr<const std::vector<std::size_t>> msMinIndAtEachTime;

public:
    Tenor( const Terms& inTerms,
           const std::shared_ptr<const std::vector<std::size_t>> insData );
    Tenor( const Terms& inTerms, const std::vector<std::size_t>& inData );
    /**
     * @brief this returns the index of T_{inIndex} in Terms.
     */
    std::size_t operator[]( std::size_t inIndex ) const;
    /**
     * @brief This returns T_{inIndex}.
     */
    double term( std::size_t inIndex ) const;
    /**
     * @brief This returns tau_{inIndex} = T_{inIndex+1} - T_{inIndex}.
     */
    double tau( std::size_t inIndex ) const;
    /**
     * @brief This gets the Tau Vec [tau_0, tau_1, ..., tau_{N-1}].
     */
    const Math::Vec& getTauVec() const;
    /**
     * @brief This returns the minimum tenor index i which satisfies Tenor[i] >=
     * Terms[inIndex].
     */
    std::size_t minIndex( std::size_t inIndex ) const;
    std::size_t size() const;
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

public:
    ZCB( const Terms& inTerms,
         std::shared_ptr<const std::vector<double>> insData,
         std::size_t inDeg = 3 );
    ZCB( const Terms& inTerms, std::vector<double> inData,
         std::size_t inDeg = 3 );
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
    double initialSpotRate() const;
};

/**
 * @brief This stores prices of Caplet at each term.
 */
class Caplets
{
private:
    const Terms mTerms;
    const std::shared_ptr<const std::vector<double>> msStrikes;
    const std::shared_ptr<const std::vector<std::vector<double>>> msCaplets;
    const std::shared_ptr<const std::vector<double>> msZCBs;
    const std::shared_ptr<Analytical::Black76::Model> msBlack;

public:
    /**
     * @brief This constructs a new Caplets.
     * @param inTerms
     * @param insStrikes
     * @param insCaplets Caplet[i][j] must be the price of caplet at strike[i]
     * in Term[j]~Term[j+1].
     * @param inZCB zero coupon bond to calc DF.
     */
    Caplets( const Terms& inTerms,
             std::shared_ptr<const std::vector<double>> insStrikes,
             std::shared_ptr<const std::vector<std::vector<double>>> insCaplets,
             const ZCB& inZCB );
    Caplets( const Terms& inTerms, const std::vector<double>& inStrikes,
             const std::vector<std::vector<double>>& inCaplets,
             const ZCB& inZCB );
    double strike( std::size_t inIndex ) const;
    const std::vector<double>& operator[]( std::size_t inIndex ) const;
    double impliedBlackVol( std::size_t inIndexStrike,
                            std::size_t inIndexTerm ) const;
};

}  // namespace MarketData
}  // namespace Process

#endif