/**
 * @file model_data.hpp
 * @brief This defines classes for model data. Unless there is a special
 * reason, calculation results should be stored in these classes. Since most
 * data is stored by smart pointers, there is no need to worry about the cost of
 * copying instances.
 * @author kakune
 * @date 3/30/2024
 */

#ifndef PROCESS_MODEL_DATA_HPP
#define PROCESS_MODEL_DATA_HPP

#include "process/market_data.hpp"

namespace Process
{
namespace ModelData
{

/**
 * @brief This stores spot rate data at each path and each term. It
 * automatically calculate Discount factors.
 */
class SpotRates
{
private:
    const std::size_t mNPath;        //! the number of path
    const MarketData::Terms mTerms;  //! terms
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
        const MarketData::Terms& inTerms,
        std::shared_ptr<const std::vector<std::vector<double>>>
            insDataSpotRate );

public:
    /**
     * @brief This constructs a new Spot Rates.
     * @param inTerms Term structure
     * @param insDataSpotRate Spot rate
     */
    SpotRates( const MarketData::Terms& inTerms,
               std::shared_ptr<const std::vector<std::vector<double>>>
                   insDataSpotRate );
    SpotRates( const MarketData::Terms& inTerms,
               std::vector<std::vector<double>> inDataSpotRate );
    const std::vector<double>& operator[]( std::size_t inIndex ) const;
    double term( std::size_t inIndex ) const;
    const MarketData::Terms& getTerms() const;
    const std::vector<std::vector<double>>& getDF() const;
    std::size_t sizeTerms() const;
    std::size_t sizePath() const;

    MarketData::ZCB createZCB( std::size_t inDeg = 3 ) const;
};

/**
 * @brief This stores forward rate data at each path and each term.
 */
class ForwardRates
{
private:
    const std::size_t mNPath;                  //! the number of path
    const MarketData::Terms mTerms;            //! terms
    const std::vector<std::size_t> mIndTenor;  //! indices of tenor of FR
    const std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
        msDataForwardRate;  //! forward rates

public:
    /**
     * @brief This constructs a new Forward Rates.
     * @param inTerms Term structure
     * @param inIndTenor indices of tenor of FR
     * @param insDataForwardRate forward rate
     */
    ForwardRates( const MarketData::Terms& inTerms,
                  const std::vector<std::size_t>& inIndTenor,
                  std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
                      insDataForwardRate );
    ForwardRates( const MarketData::Terms& inTerms,
                  const std::vector<std::size_t>& inIndTenor,
                  std::vector<std::vector<Math::Vec>> inDataForwardRate );
    const std::vector<Math::Vec>& operator[]( std::size_t inIndex ) const;
    double term( std::size_t inIndex ) const;
    const MarketData::Terms& getTerms() const;
    std::size_t sizeTerms() const;
    std::size_t sizePath() const;
};

}  // namespace ModelData
}  // namespace Process

#endif