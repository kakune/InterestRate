#ifndef APP_SHORT_RATE_PREPARE_MODEL_HPP
#define APP_SHORT_RATE_PREPARE_MODEL_HPP

#include <memory>
#include <string>
#include <vector>

#include "process/market_data.hpp"
#include "short_rate/multi-factor.hpp"
#include "short_rate/one-factor.hpp"
#include "utils/parameters.hpp"

namespace APP
{

Process::MarketData::Terms prepareTerms( const Utils::Parameters& inParams )
{
    std::size_t lNTerms = std::size_t( inParams( "NTerms" ) );
    std::vector<double> lTerms( lNTerms, 0 );

    double lDt = inParams( "TimeMaturity" ) / double( lNTerms - 1 );
    for ( std::size_t iTerm = 1; iTerm < lNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }

    return Process::MarketData::Terms( lTerms );
}

Process::MarketData::SpotRates calcSpotRateFromParam(
    std::string inNameModel, const Utils::Parameters& inParams,
    const Process::MarketData::Terms& inTerms )
{
    if ( inNameModel == "Vasicek" )
    {
        auto luRandom = std::make_unique<Process::Random::PathBrownAntithetic>(
            inParams( "NPath" ), inTerms );
        ShortRate::OneFactor::VasicekBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setInitSpotRate( inParams( "InitRate" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        lBuilder.setKappa( inParams( "Kappa" ) );
        lBuilder.setMean( inParams( "Mean" ) );
        return Process::MarketData::SpotRates(
            lBuilder.build().calcSpotRates() );
    }
    if ( inNameModel == "HoLee" )
    {
        auto luRandom = std::make_unique<Process::Random::PathBrownAntithetic>(
            inParams( "NPath" ), inTerms );
        ShortRate::OneFactor::HoLeeBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setInitSpotRate( inParams( "InitRate" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        return Process::MarketData::SpotRates(
            lBuilder.build().calcSpotRates() );
    }
    if ( inNameModel == "ConstantAffine" )
    {
        auto luRandom = std::make_unique<Process::Random::PathBrownAntithetic>(
            inParams( "NPath" ), inTerms );
        ShortRate::OneFactor::ConstantAffineBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setInitSpotRate( inParams( "InitRate" ) );
        lBuilder.setDrift( inParams( "Lambda" ), inParams( "Eta" ) );
        lBuilder.setVol( inParams( "Gamma" ), inParams( "Delta" ) );
        return Process::MarketData::SpotRates(
            lBuilder.build().calcSpotRates() );
    }
    return Process::MarketData::SpotRates(
        inTerms,
        std::vector<std::vector<double> >(
            inParams( "NPath" ), std::vector<double>( inTerms.size() ) ) );
}

Process::MarketData::SpotRates calcSpotRateFromParam(
    std::string inNameModel, const Utils::Parameters& inParams )
{
    return ( calcSpotRateFromParam( inNameModel, inParams,
                                    prepareTerms( inParams ) ) );
}

Process::MarketData::SpotRates calcSpotRateFromMarket(
    std::string inNameModel, const Utils::Parameters& inParams,
    const Process::MarketData::Terms inTerms,
    const Process::MarketData::ZCB inMarketZCB )
{
    if ( inNameModel == "Vasicek" )
    {
        auto luRandom = std::make_unique<Process::Random::PathBrownAntithetic>(
            inParams( "NPath" ), inTerms );
        ShortRate::OneFactor::VasicekWithMarketBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        lBuilder.setMarketZCB( inMarketZCB );
        return Process::MarketData::SpotRates(
            lBuilder.build().calcSpotRates() );
    }
    if ( inNameModel == "HoLee" )
    {
        auto luRandom = std::make_unique<Process::Random::PathBrownAntithetic>(
            inParams( "NPath" ), inTerms );
        ShortRate::OneFactor::HoLeeWithMarketBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        lBuilder.setMarketZCB( inMarketZCB );
        return Process::MarketData::SpotRates(
            lBuilder.build().calcSpotRates() );
    }
    if ( inNameModel == "G2pp" )
    {
        auto luRandom =
            std::make_unique<Process::RandomVec::PathBrownAntithetic>(
                inParams( "NPath" ), inTerms, 2 );
        ShortRate::MultiFactor::G2ppWithMarketBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        Math::Mat lVol( 2, 2, 0.0 );
        lVol( 0, 0 ) = inParams( "Sigma" );
        lVol( 1, 0 ) = inParams( "Eta" ) * inParams( "Rho" );
        lVol( 1, 1 ) = inParams( "Eta" ) *
                       std::sqrt( 1.0 - inParams( "Rho" ) * inParams( "Rho" ) );
        lBuilder.setDrift( { -inParams( "A" ), -inParams( "B" ) } );
        lBuilder.setVol( lVol );
        lBuilder.setMarketZCB( inMarketZCB );
        return Process::MarketData::SpotRates(
            lBuilder.build().calcSpotRates() );
    }
    return Process::MarketData::SpotRates(
        inTerms,
        std::vector<std::vector<double> >(
            inParams( "NPath" ), std::vector<double>( inTerms.size() ) ) );
}
Process::MarketData::SpotRates calcSpotRateFromMarket(
    std::string inNameModel, const Utils::Parameters& inParams,
    const Process::MarketData::ZCB inMarketZCB )
{
    return calcSpotRateFromMarket( inNameModel, inParams,
                                   prepareTerms( inParams ), inMarketZCB );
}

}  // namespace APP

#endif