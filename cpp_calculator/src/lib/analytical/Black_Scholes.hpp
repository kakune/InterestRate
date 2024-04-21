/**
 * @file Black_Scholes.hpp
 * @brief This defines analytical calculators relating to Black Scholes Model.
 * @author kakune
 * @date 1/29/2024
 */

#ifndef ANALYTICAL_BLACK_SCHOLES_HPP
#define ANALYTICAL_BLACK_SCHOLES_HPP

#include <cmath>

#include "math/special_functions.hpp"

namespace Analytical
{
namespace BlackScholes
{

/**
 * @brief This stores BS parameters and has calculators.
 */
struct Model;

/**
 * @brief This calculates european Call Option Price of BS Model.
 * @param inInitS initial spot price
 * @param inStrike strike price of the option
 * @param inRate risk-free interest rate
 * @param inVol volatility
 * @param inTMat time until maturity
 * @return double price
 */
double priceEuropeanCallOption( double inInitS, double inStrike, double inRate,
                                double inVol, double inTMat );
/**
 * @brief This calculates european Put Option Price of BS Model.
 * @param inInitS initial spot price
 * @param inStrike strike price of the option
 * @param inRate risk-free interest rate
 * @param inVol volatility
 * @param inTMat time until maturity
 * @return double price
 */
double priceEuropeanPutOption( double inInitS, double inStrike, double inRate,
                               double inVol, double inTMat );

struct Model
{
    double mInitS;   //! initial spot price
    double mStrike;  //! strike price of options
    double mRate;    //! risk-free interest rate
    double mVol;     //! volatility
    double mTMat;    //! time until maturity

    /**
     * @brief This calculates european Call Option Price of BS Model.
     * @return double price
     */
    double priceEuropeanCallOption();
    /**
     * @brief This calculates european Put Option Price of BS Model.
     * @return double price
     */
    double priceEuropeanPutOption();
};

}  // namespace BlackScholes
}  // namespace Analytical

#endif