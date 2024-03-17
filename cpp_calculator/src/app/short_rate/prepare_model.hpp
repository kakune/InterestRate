#ifndef APP_SHORT_RATE_PREPARE_MODEL_HPP
#define APP_SHORT_RATE_PREPARE_MODEL_HPP

#include <memory>
#include <string>
#include <vector>

#include "process/short_rate.hpp"
#include "utils/parameters.hpp"

namespace APP
{
namespace ShortRate
{

std::shared_ptr<std::vector<double> > prepareTerms(
    const Utils::Parameters& inParams )
{
    std::size_t lNTerms = std::size_t( inParams( "NTerms" ) );
    std::vector<double> lTerms( lNTerms, 0 );

    double lDt = inParams( "TimeMaturity" ) / double( lNTerms - 1 );
    for ( std::size_t iTerm = 1; iTerm < lNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }

    return std::make_shared<std::vector<double> >( lTerms );
}

std::unique_ptr<Process::ShortRate::ModelAbstract> prepareModelFromParam(
    std::string inNameModel, const Utils::Parameters& inParams,
    std::shared_ptr<std::vector<double> > insTerms )
{
    auto luRandom = std::make_unique<Process::Random::PathBrownAntithetic>(
        inParams( "NPath" ), insTerms );
    if ( inNameModel == "Vasicek" )
    {
        Process::ShortRate::VasicekBuilder lBuilder;
        lBuilder.setTerms( insTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setInitSpotRate( inParams( "InitRate" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        lBuilder.setKappa( inParams( "Kappa" ) );
        lBuilder.setMean( inParams( "Mean" ) );
        return std::move(
            std::make_unique<Process::ShortRate::Vasicek>( lBuilder.build() ) );
    }
    if ( inNameModel == "HoLee" )
    {
        Process::ShortRate::HoLeeBuilder lBuilder;
        lBuilder.setTerms( insTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setInitSpotRate( inParams( "InitRate" ) );
        lBuilder.setVol( inParams( "Vol" ) );
        return std::move(
            std::make_unique<Process::ShortRate::HoLee>( lBuilder.build() ) );
    }
    if ( inNameModel == "ConstantAffine" )
    {
        Process::ShortRate::ConstantAffineBuilder lBuilder;
        lBuilder.setTerms( insTerms );
        lBuilder.setRandom( std::move( luRandom ) );
        lBuilder.setNPath( inParams( "NPath" ) );
        lBuilder.setInitSpotRate( inParams( "InitRate" ) );
        lBuilder.setDrift( inParams( "Lambda" ), inParams( "Eta" ) );
        lBuilder.setVol( inParams( "Gamma" ), inParams( "Delta" ) );
        return std::move( std::make_unique<Process::ShortRate::ConstantAffine>(
            lBuilder.build() ) );
    }
    return nullptr;
}

std::unique_ptr<Process::ShortRate::ModelAbstract> prepareModelFromParam(
    std::string inNameModel, const Utils::Parameters& inParams )
{
    return ( prepareModelFromParam( inNameModel, inParams,
                                    prepareTerms( inParams ) ) );
}

}  // namespace ShortRate
}  // namespace APP

#endif