#include <gtest/gtest.h>

#include "analytical/SABR.hpp"

double testBlackImpVolOfSABR( double inTime, double inInitPrice,
                              double inInitVol, double inExponent,
                              double inCorr, double inVolVol, double inStrike )
{
    Analytical::SABR::OneTerm lObj( inTime );
    lObj.setInitPrice( inInitPrice )
        .setInitVol( inInitVol )
        .setExponent( inExponent )
        .setCorr( inCorr )
        .setVolVol( inVolVol )
        .setStrike( inStrike );
    for ( double lStrike = 0.75; lStrike <= 1.25; lStrike += 0.02 )
    {
        lObj.setStrike( inInitPrice * lStrike );
        std::cout << lObj.approxBlackImpVolByHagan() << std::endl;
    }
    lObj.setStrike( inStrike );
    return lObj.approxBlackImpVolByHagan();
}

TEST( SABRTest, BlackImpVol )
{
    EXPECT_NEAR( 0.2003975,
                 testBlackImpVolOfSABR( 0.25, 0.05, 0.05, 0.8, 0.1, 0.8, 0.05 ),
                 1e-6 );
}
