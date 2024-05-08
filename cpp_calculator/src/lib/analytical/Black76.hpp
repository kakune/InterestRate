/**
 * @file Black76.hpp
 * @brief This defines analytical calculators relating to Black-76 Model.
 * @author kakune
 * @date 4/13/2024
 */

#ifndef ANALYTICAL_BLACK76_HPP
#define ANALYTICAL_BLACK76_HPP

#include <cmath>
#include <vector>

namespace Analytical
{
namespace Black76
{

double funcBlackPositive( double inStrike, double inPrice, double inVol );
double funcBlackNegative( double inStrike, double inPrice, double inVol );

class Model
{
private:
    const std::vector<double> mTerms;  //! term structure
    std::vector<double> mForwardRate;  //! forward rate at time 0
    std::vector<double> mZCB;          //! ZCB at time 0
    std::vector<double> mSumTauZCB;
    double mVol;  //! volatility

    std::vector<double> calcZCBFromForwardRate(
        const std::vector<double>& inForwardRate ) const;
    std::vector<double> calcForwardRateFromZCB(
        const std::vector<double>& inZCB ) const;
    std::vector<double> calcSumTauZCBFromZCB(
        const std::vector<double>& inZCB ) const;

public:
    Model( const std::vector<double>& inTerms );
    void setInitForwardRate( const std::vector<double>& inForwardRate );
    void setInitZCB( const std::vector<double>& inZCB );
    void setVol( double inVol );
    double priceCaplet( double inStrike, std::size_t inIndex ) const;
    double priceFloorlet( double inStrike, std::size_t inIndex ) const;
    double priceCap( double inStrike, std::size_t inIndStart,
                     std::size_t inIndLast ) const;
    double priceFloor( double inStrike, std::size_t inIndStart,
                       std::size_t inIndLast ) const;
    double priceSwapRate( std::size_t inIndStart, std::size_t inIndLast ) const;
    double pricePayerSwaption( double inStrike, std::size_t inIndStart,
                               std::size_t inIndLast ) const;
    double priceReceiverSwaption( double inStrike, std::size_t inIndStart,
                                  std::size_t inIndLast ) const;
};

}  // namespace Black76
}  // namespace Analytical

#endif