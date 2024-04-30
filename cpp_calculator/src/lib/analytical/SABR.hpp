/**
 * @file SABR.hpp
 * @brief This defines analytical calculators relating to SABR Model. See also
 * "SABR and SABR LIBOR Market Models in Practice."
 * @author kakune
 * @date 4/30/2024
 */

#ifndef ANALYTICAL_SABR_HPP
#define ANALYTICAL_SABR_HPP

#include <cmath>
#include <vector>

namespace Analytical::SABR
{

constexpr double gInv24   = 1.0 / 24.0;
constexpr double gInv1920 = 1.0 / 1920.0;

class OneTerm
{
private:
    const double mTime;  //! time to start of forward rate
    double mInitPrice;   //! F
    double mInitVol;     //! alpha
    double mExponent;    //! beta
    double mCorr;        //! rho
    double mVolVol;      //! nu
    double mStrike;      //! K

    double approxBlackImpVolByHaganATM() const;

public:
    OneTerm( double inTime );
    OneTerm& setInitPrice( double inInitPrice );
    OneTerm& setInitVol( double inInitVol );
    OneTerm& setExponent( double inExponent );
    OneTerm& setCorr( double inCorr );
    OneTerm& setVolVol( double inVolVol );
    OneTerm& setStrike( double inStrike );
    double approxBlackImpVolByHagan() const;
    double approxNormalImpVolByHagan() const;
};

}  // namespace Analytical::SABR

#endif