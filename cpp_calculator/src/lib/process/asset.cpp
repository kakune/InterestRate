/**
 * @file asset.hpp
 * @brief This implements calculator of assets.
 * @author kakune
 * @date 1/30/2024
 */

#include "process/asset.hpp"

#include <memory>
#include <vector>

#include "analytical/Black_Scholes.hpp"
#include "math/findroot_1d.hpp"

namespace Process
{
namespace Asset
{

double ModelForwardAbstract::priceCallOption( double inStrike,
                                              std::size_t inIndTime )
{
    double lResult = 0;
    for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
    {
        lResult += std::max(
            0.0, mForwardPrice.at( iPath ).at( inIndTime ) - inStrike );
    }
    return ( lResult / double( mNPath ) );
}
double ModelForwardAbstract::pricePutOption( double inStrike,
                                             std::size_t inIndTime )
{
    double lResult = 0;
    for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
    {
        lResult += std::max(
            0.0, inStrike - mForwardPrice.at( iPath ).at( inIndTime ) );
    }
    return ( lResult / double( mNPath ) );
}
double ModelForwardAbstract::impliedVolatility( double inStrike,
                                                std::size_t inIndTime )
{
    // if ( inStrike > mInitPrice )
    {
        double lCallPrice    = priceCallOption( inStrike, inIndTime );
        double lTimeMaturity = msTerms->at( inIndTime ) - msTerms->front();
        auto lFunction       = [this, inStrike, lCallPrice,
                          lTimeMaturity]( double inVol ) -> double
        {
            return lCallPrice -
                   Analytical::BlackScholes::europeanCallOptionPrice(
                       this->mInitPrice, inStrike, 0.0, inVol, lTimeMaturity );
        };
        return Math::FindRoot1D::Brent( lFunction, 0.000000001, 10.0 );
    }
    double lPutPrice     = pricePutOption( inStrike, inIndTime );
    double lTimeMaturity = msTerms->at( inIndTime ) - msTerms->front();
    auto lFunction       = [this, inStrike, lPutPrice,
                      lTimeMaturity]( double inVol ) -> double
    {
        return lPutPrice -
               Analytical::BlackScholes::europeanPutOptionPrice(
                   this->mInitPrice, inStrike, 0.0, inVol, lTimeMaturity );
    };
    return Math::FindRoot1D::Brent( lFunction, 0.000000001, 10.0 );
}

void LocalVolatilityForwardAbstract::calcEachForwardPrice()
{
    for ( std::size_t iTime = 1; iTime < msTerms->size(); ++iTime )
    {
        muRandomPath->setIndexTime( iTime );
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            mForwardPrice.at( iPath ).at( iTime ) =
                mForwardPrice.at( iPath ).at( iTime - 1 ) +
                muRandomPath->generateRandomVal() *
                    localVolatility( mForwardPrice.at( iPath ).at( iTime - 1 ),
                                     msTerms->at( iTime - 1 ) );
        }
    }
}

double BlackSholesForward::localVolatility( double inPrice, double inTime )
{
    return mVol * inPrice;
}

void StochasticVolatilityForwardAbstract::calcEachForwardPrice()
{
    for ( std::size_t iTime = 1; iTime < msTerms->size(); ++iTime )
    {
        muRandomPath->setIndexTime( iTime );
        muRandomVol->setIndexTime( iTime );
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            double lW1 = muRandomPath->generateRandomVal();
            double lW2 = muRandomVol->generateRandomVal();

            mVolatility.at( iPath ).at( iTime ) =
                mVolatility.at( iPath ).at( iTime - 1 ) +
                lW2 * volVolatility( mVolatility.at( iPath ).at( iTime - 1 ),
                                     msTerms->at( iTime - 1 ) );

            mForwardPrice.at( iPath ).at( iTime ) =
                mForwardPrice.at( iPath ).at( iTime - 1 ) +
                ( lW2 * mCorr + lW1 * mAuxCorr ) *
                    volForward( mForwardPrice.at( iPath ).at( iTime - 1 ),
                                mVolatility.at( iPath ).at( iTime - 1 ),
                                msTerms->at( iTime - 1 ) );
            if ( mVolatility.at( iPath ).at( iTime ) < 0.0 )
            {
                mVolatility.at( iPath ).at( iTime ) = 0.0001;
            }
            if ( mForwardPrice.at( iPath ).at( iTime ) < 0.0 )
            {
                mForwardPrice.at( iPath ).at( iTime ) = 0.0001;
            }
        }
    }
}

double SABRForward::volForward( double inPrice, double inVol, double inTime )
{
    return inVol * std::pow( inPrice, mExponent );
}
double SABRForward::volVolatility( double inVol, double inTime )
{
    return mVolvol * inVol;
}

double SABRWithLogForward::volForward( double inPrice, double inVol,
                                       double inTime )
{
    return inVol * std::pow( inPrice, mExponent - 1.0 );
}
double SABRWithLogForward::volVolatility( double inVol, double inTime )
{
    return mVolvol;
}
double SABRWithLogForward::driftForward( double inPrice, double inVol,
                                         double inTime )
{
    double lTmp = inVol * std::pow( inPrice, mExponent - 1.0 );
    return -0.5 * lTmp * lTmp;
}
double SABRWithLogForward::driftVolatility( double inVol, double inTime )
{
    return -0.5 * mVolvol * mVolvol;
}

void StochasticVolatilityWithLogForwardAbstract::calcEachForwardPrice()
{
    for ( std::size_t iTime = 1; iTime < msTerms->size(); ++iTime )
    {
        muRandomPath->setIndexTime( iTime );
        muRandomVol->setIndexTime( iTime );
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            double lW1 = muRandomPath->generateRandomVal();
            double lW2 = muRandomVol->generateRandomVal();

            double lTerm = msTerms->at( iTime ) - msTerms->at( iTime - 1 );
            double lLogVol =
                std::log( mVolatility.at( iPath ).at( iTime - 1 ) );
            lLogVol +=
                lW2 * volVolatility( mVolatility.at( iPath ).at( iTime - 1 ),
                                     msTerms->at( iTime - 1 ) );
            lLogVol += lTerm *
                       driftVolatility( mVolatility.at( iPath ).at( iTime - 1 ),
                                        msTerms->at( iTime - 1 ) );
            mVolatility.at( iPath ).at( iTime ) = std::exp( lLogVol );

            double lLogForward =
                std::log( mForwardPrice.at( iPath ).at( iTime - 1 ) );
            lLogForward +=
                ( lW2 * mCorr + lW1 * mAuxCorr ) *
                volForward( mForwardPrice.at( iPath ).at( iTime - 1 ),
                            mVolatility.at( iPath ).at( iTime - 1 ),
                            msTerms->at( iTime - 1 ) );
            lLogForward +=
                lTerm * driftForward( mForwardPrice.at( iPath ).at( iTime - 1 ),
                                      mVolatility.at( iPath ).at( iTime - 1 ),
                                      msTerms->at( iTime - 1 ) );
            mForwardPrice.at( iPath ).at( iTime ) = std::exp( lLogForward );
        }
    }
}

}  // namespace Asset
}  // namespace Process