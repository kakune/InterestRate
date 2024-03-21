#ifndef APP_SHORT_RATE_PREPARE_MODEL_HPP
#define APP_SHORT_RATE_PREPARE_MODEL_HPP

#include <memory>
#include <string>
#include <vector>

#include "process/market_data.hpp"
#include "process/short_rate_MC.hpp"
#include "utils/parameters.hpp"

namespace APP
{
namespace ShortRate
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

std::unique_ptr<Process::ShortRateMC::ModelAbstract> prepareModelFromParam(
    std::string inNameModel, const Utils::Parameters& inParams,
    const Process::MarketData::Terms& inTerms )
{
    auto luRandom = std::make_unique<Process::Random::PathBrownAntithetic>(
        inParams( "NPath" ), inTerms );
    if ( inNameModel == "Vasicek" )
    {
        Process::ShortRateMC::VasicekBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setInitSpotRate( inParams( "InitRate" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        lBuilder.setKappa( inParams( "Kappa" ) );
        lBuilder.setMean( inParams( "Mean" ) );
        return std::move( std::make_unique<Process::ShortRateMC::Vasicek>(
            lBuilder.build() ) );
    }
    if ( inNameModel == "HoLee" )
    {
        Process::ShortRateMC::HoLeeBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setInitSpotRate( inParams( "InitRate" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        return std::move(
            std::make_unique<Process::ShortRateMC::HoLee>( lBuilder.build() ) );
    }
    if ( inNameModel == "ConstantAffine" )
    {
        Process::ShortRateMC::ConstantAffineBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setInitSpotRate( inParams( "InitRate" ) );
        lBuilder.setDrift( inParams( "Lambda" ), inParams( "Eta" ) );
        lBuilder.setVol( inParams( "Gamma" ), inParams( "Delta" ) );
        return std::move(
            std::make_unique<Process::ShortRateMC::ConstantAffine>(
                lBuilder.build() ) );
    }
    return nullptr;
}

std::unique_ptr<Process::ShortRateMC::ModelAbstract> prepareModelFromParam(
    std::string inNameModel, const Utils::Parameters& inParams )
{
    return ( prepareModelFromParam( inNameModel, inParams,
                                    prepareTerms( inParams ) ) );
}

std::unique_ptr<Process::ShortRateMC::ModelAbstract> prepareModelFromMarket(
    std::string inNameModel, const Utils::Parameters& inParams,
    const Process::MarketData::Terms inTerms,
    const Process::MarketData::ZCB inMarketZCB )
{
    auto luRandom = std::make_unique<Process::Random::PathBrownAntithetic>(
        inParams( "NPath" ), inTerms );
    if ( inNameModel == "Vasicek" )
    {
        Process::ShortRateMC::VasicekWithMarketBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        lBuilder.setMarketZCB( inMarketZCB );
        return std::move(
            std::make_unique<Process::ShortRateMC::VasicekWithMarket>(
                lBuilder.build() ) );
    }
    if ( inNameModel == "HoLee" )
    {
        Process::ShortRateMC::HoLeeWithMarketBuilder lBuilder;
        lBuilder.setTerms( inTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        lBuilder.setMarketZCB( inMarketZCB );
        return std::move(
            std::make_unique<Process::ShortRateMC::HoLeeWithMarket>(
                lBuilder.build() ) );
    }
    return nullptr;
}
std::unique_ptr<Process::ShortRateMC::ModelAbstract> prepareModelFromMarket(
    std::string inNameModel, const Utils::Parameters& inParams,
    const Process::MarketData::ZCB inMarketZCB )
{
    return prepareModelFromMarket( inNameModel, inParams,
                                   prepareTerms( inParams ), inMarketZCB );
}

}  // namespace ShortRate
}  // namespace APP

#endif