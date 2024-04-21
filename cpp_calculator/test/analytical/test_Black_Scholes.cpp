#include <gtest/gtest.h>

#include "analytical/Black_Scholes.hpp"

TEST( BlackScholesTest, EuropeanCall )
{
    EXPECT_NEAR( 0.1098954915262599,
                 Analytical::BlackScholes::priceEuropeanCallOption(
                     1.0, 1.0, 0.06, 0.2, 1.0 ),
                 1e-6 );
    EXPECT_NEAR( 0.00004773739958771617,
                 Analytical::BlackScholes::priceEuropeanCallOption(
                     1.0, 1.2, 0.01, 0.01, 10.0 ),
                 1e-10 );
}

TEST( BlackScholesTest, EuropeanPut )
{
    EXPECT_NEAR( 0.0516600251105086,
                 Analytical::BlackScholes::priceEuropeanPutOption(
                     1.0, 1.0, 0.06, 0.2, 1.0 ),
                 1e-6 );
    EXPECT_NEAR( 0.0858526390427392,
                 Analytical::BlackScholes::priceEuropeanPutOption(
                     1.0, 1.2, 0.01, 0.01, 10.0 ),
                 1e-10 );
}