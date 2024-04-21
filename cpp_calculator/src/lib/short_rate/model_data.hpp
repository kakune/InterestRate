/**
 * @file model_data.hpp
 * @brief This defines classes for model data. Unless there is a special
 * reason, calculation results should be stored in these classes. Since most
 * data is stored by smart pointers, there is no need to worry about the cost of
 * copying instances.
 * @author kakune
 * @date 4/21/2024
 */

#ifndef SHORT_RATE_MODEL_DATA_HPP
#define SHORT_RATE_MODEL_DATA_HPP
#include "process/market_data.hpp"

namespace ShortRate
{

Process::MarketData::ZCB createZCBFromSpotRates(
    const Process::MarketData::Terms& inTerms,
    const std::vector<std::vector<double>>& inDFs, std::size_t inDeg );
/**
 * @brief This calculates DF from spot rates.
 * @param inTerms
 * @param insDataSpotRate
 * @return std::vector<std::vector<double>> DFs
 */
std::vector<std::vector<double>> calcDFFromSpotRates(
    const Process::MarketData::Terms& inTerms,
    std::shared_ptr<const std::vector<std::vector<double>>> insDataSpotRate );
/**
 * @brief This stores spot rate data at each path and each term UNDER RISK
 * NEUTRAL MEASURE. It automatically calculate Discount factors and ZCB.
 */
class SpotRates
{
private:
    const std::size_t mNPath;                 //! the number of path
    const Process::MarketData::Terms mTerms;  //! terms
    const std::shared_ptr<const std::vector<std::vector<double>>>
        msDataSpotRate;  //! spot rates
    const std::shared_ptr<const std::vector<std::vector<double>>>
        msDataDF;  //! discount factor
    const std::shared_ptr<const Process::MarketData::ZCB> msZCB;

public:
    /**
     * @brief This constructs a new Spot Rates.
     * @param inTerms Term structure
     * @param insDataSpotRate Spot rate
     */
    SpotRates(
        const Process::MarketData::Terms& inTerms,
        std::shared_ptr<const std::vector<std::vector<double>>> insDataSpotRate,
        std::size_t inDegZCB = 3 );
    SpotRates( const Process::MarketData::Terms& inTerms,
               std::vector<std::vector<double>> inDataSpotRate,
               std::size_t inDegZCB = 3 );
    const std::vector<double>& operator[]( std::size_t inIndex ) const;
    double term( std::size_t inIndex ) const;
    const Process::MarketData::Terms& getTerms() const;
    const std::vector<std::vector<double>>& getDF() const;
    const Process::MarketData::ZCB& getZCB() const;
    std::size_t sizeTerms() const;
    std::size_t sizePath() const;

    double calcMeanMMA( std::size_t inIndStart, std::size_t inIndEnd ) const;
    double calcCaplet( double inStrike, std::size_t inIndStart,
                       std::size_t inIndEnd ) const;
    double calcFloorlet( double inStrike, std::size_t inIndStart,
                         std::size_t inIndEnd ) const;
    double calcCap( double inStrike, std::vector<std::size_t> inTenor ) const;
    double calcFloor( double inStrike, std::vector<std::size_t> inTenor ) const;
};

}  // namespace ShortRate

#endif