#include "short_rate_creator.hpp"

#include <string>
#include <vector>

#include "process/market_data.hpp"
#include "short_rate/multi-factor.hpp"
#include "short_rate/one-factor.hpp"
#include "terms_creator.hpp"
#include "utils/parameters.hpp"

namespace APP
{

ShortRate::SpotRates createSpotRateFromParam(
    const Utils::Parameters& inParams,
    const Process::MarketData::Terms& inTerms )
{
    if ( inParams.operator()<std::string>( "NameModel" ) == "Vasicek" )
    {
        auto luRandom = std::make_unique<Process::Random::StdBrownAntithetic>();
        ShortRate::OneFactor::VasicekBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setInitSpotRate( inParams( "InitRate" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        lBuilder.setKappa( inParams( "Kappa" ) );
        lBuilder.setMean( inParams( "Mean" ) );
        return ShortRate::SpotRates( lBuilder.build().createSpotRates() );
    }
    if ( inParams.operator()<std::string>( "NameModel" ) == "HoLee" )
    {
        auto luRandom = std::make_unique<Process::Random::StdBrownAntithetic>();
        ShortRate::OneFactor::HoLeeBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setInitSpotRate( inParams( "InitRate" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        return ShortRate::SpotRates( lBuilder.build().createSpotRates() );
    }
    if ( inParams.operator()<std::string>( "NameModel" ) == "ConstantAffine" )
    {
        auto luRandom = std::make_unique<Process::Random::StdBrownAntithetic>();
        ShortRate::OneFactor::ConstantAffineBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setInitSpotRate( inParams( "InitRate" ) );
        lBuilder.setDrift( inParams( "Lambda" ), inParams( "Eta" ) );
        lBuilder.setVol( inParams( "Gamma" ), inParams( "Delta" ) );
        return ShortRate::SpotRates( lBuilder.build().createSpotRates() );
    }
    return ShortRate::SpotRates(
        inTerms,
        std::vector<std::vector<double> >(
            inParams( "NPath" ), std::vector<double>( inTerms.size() ) ) );
}

ShortRate::SpotRates createSpotRateFromParam(
    const Utils::Parameters& inParams )
{
    return ( createSpotRateFromParam( inParams, prepareTerms( inParams ) ) );
}

ShortRate::SpotRates createSpotRateFromMarket(
    const Utils::Parameters& inParams, const Process::MarketData::Terms inTerms,
    const Process::MarketData::ZCB inMarketZCB )
{
    if ( inParams.operator()<std::string>( "NameModel" ) == "Vasicek" )
    {
        auto luRandom = std::make_unique<Process::Random::StdBrownAntithetic>();
        ShortRate::OneFactor::VasicekWithMarketBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        lBuilder.setMarketZCB( inMarketZCB );
        return ShortRate::SpotRates( lBuilder.build().createSpotRates() );
    }
    if ( inParams.operator()<std::string>( "NameModel" ) == "HoLee" )
    {
        auto luRandom = std::make_unique<Process::Random::StdBrownAntithetic>();
        ShortRate::OneFactor::HoLeeWithMarketBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        lBuilder.setMarketZCB( inMarketZCB );
        return ShortRate::SpotRates( lBuilder.build().createSpotRates() );
    }
    if ( inParams.operator()<std::string>( "NameModel" ) == "G2pp" )
    {
        auto luRandom =
            std::make_unique<Process::RandomVec::StdBrownAntithetic>( 2 );
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
        return ShortRate::SpotRates( lBuilder.build().createSpotRates() );
    }
    return ShortRate::SpotRates(
        inTerms,
        std::vector<std::vector<double> >(
            inParams( "NPath" ), std::vector<double>( inTerms.size() ) ) );
}
ShortRate::SpotRates createSpotRateFromMarket(
    const Utils::Parameters& inParams,
    const Process::MarketData::ZCB inMarketZCB )
{
    return createSpotRateFromMarket( inParams, prepareTerms( inParams ),
                                     inMarketZCB );
}

}  // namespace APP