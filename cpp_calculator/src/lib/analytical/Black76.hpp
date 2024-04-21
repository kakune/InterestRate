/**
 * @file Black76.hpp
 * @brief This defines analytical calculators relating to Black-76 Model.
 * @author kakune
 * @date 4/13/2024
 */

#ifndef ANALYTICAL_BLACK76_HPP
#define ANALYTICAL_BLACK76_HPP

#include <cmath>

#include "math/special_functions.hpp"

namespace Analytical
{
namespace Black76
{

/**
 * @brief This stores Black parameters and has calculators.
 */
struct Model;

/**
 * @brief This calculates caplet Price of Black Model.
 * @param inCurrentFR current forward rate
 * @param inStrike strike price of caplet
 * @param inTimeStart  start time of forward rate
 * @param inTau the length of forward rate
 * @param inVol volatility
 * @param inZCB price of ZCB whose maturity is inTimeStart + inTau
 * @return double price
 */
double priceCaplet( double inCurrentFR, double inStrike, double inTimeStart,
                    double inTau, double inVol, double inZCB );
/**
 * @brief This calculates floorlet Price of Black Model.
 * @param inCurrentFR current forward rate
 * @param inStrike strike price of caplet
 * @param inTimeStart  start time of forward rate
 * @param inTau the length of forward rate
 * @param inVol volatility
 * @param inZCB price of ZCB whose maturity is inTimeStart + inTau
 * @return double price
 */
double priceFloorlet( double inCurrentFR, double inStrike, double inTimeStart,
                      double inTau, double inVol, double inZCB );

struct Model
{
    double mCurrentFR;  //! initial spot price
    double mStrike;     //! strike price of options
    double mTimeStart;  //! risk-free interest rate
    double mTau;        //! volatility
    double mVol;        //! time until maturity
    double mZCB;        //! price of ZCB whose maturity is mTimeStart+mTau

    double operator()( bool inIsCaplet );
    /**
     * @brief This calculates caplet Price of Black Model.
     * @return double price
     */
    double priceCaplet();
    /**
     * @brief This calculates floorlet Price of Black Model.
     * @return double price
     */
    double priceFloorlet();
};

}  // namespace Black76
}  // namespace Analytical

#endif