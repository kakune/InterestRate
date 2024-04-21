#ifndef APP_LIB_SHORT_RATE_CREATOR_HPP
#define APP_LIB_SHORT_RATE_CREATOR_HPP

#include <string>
#include <vector>

#include "LIBOR/forward.hpp"
#include "math/matrix.hpp"
#include "process/market_data.hpp"
#include "terms_creator.hpp"
#include "utils/parameters.hpp"

namespace APP
{

static Math::Vec getVecFromParam( const std::string& inName,
                                  const Utils::Parameters& inParams )
{
    return Math::Vec( inParams.operator()<std::vector<double>>( inName ) );
}
static Math::Mat getMatFromParam( const std::string& inName,
                                  const Utils::Parameters& inParams )
{
    return Math::Mat(
        inParams.operator()<std::vector<std::vector<double>>>( inName ) );
}

LIBOR::Forward::Data::TerminalMeas createForwardTerminalFromParam(
    const Utils::Parameters& inParams,
    const Process::MarketData::Terms& inTerms,
    const Process::MarketData::Tenor& inTenor )
{
    std::size_t lNPath = inParams.operator()<int>( "NPath" );
    if ( inParams.operator()<std::string>( "NameModel" ) == "ConstantVol" )
    {
        LIBOR::Forward::VolGen::Constant lVolGen(
            getVecFromParam( "Vol", inParams ) );
        LIBOR::Forward::StepCalc::LogNormalTerminalMeasWithLog lStep(
            inTerms, inTenor, getMatFromParam( "Corr", inParams ), lVolGen );
        return LIBOR::Forward::Factory( lNPath, inTerms, inTenor,
                                        getVecFromParam( "InitFR", inParams ),
                                        lStep )
            .template createForwardRates<
                Process::RandomVec::StdBrownAntithetic>();
    }
    if ( inParams.operator()<std::string>( "NameModel" ) == "SABR" )
    {
        LIBOR::Forward::VolGen::SABR<Process::RandomVec::StdBrownAntithetic>
            lVolGen( getVecFromParam( "InitVol", inParams ),
                     inParams( "Exponent" ), inParams( "VolVol" ),
                     getVecFromParam( "CorrSV", inParams ), lNPath, inTerms );
        LIBOR::Forward::StepCalc::NormalTerminalMeas lStep(
            inTerms, inTenor, getMatFromParam( "Corr", inParams ), lVolGen );
        return LIBOR::Forward::Factory( lNPath, inTerms, inTenor,
                                        getVecFromParam( "InitFR", inParams ),
                                        lStep )
            .template createForwardRates<
                Process::RandomVec::StdBrownAntithetic>();
    }
    throw std::invalid_argument(
        std::string( "Error: APP::createForwardTerminalFromParam()\n" ) +
        std::string( "No Corresponding name of model." ) );
}

LIBOR::Forward::Data::TerminalMeas createForwardTerminalFromParam(
    const Utils::Parameters& inParams )
{
    auto lTerms = prepareTerms( inParams );
    return createForwardTerminalFromParam( inParams, lTerms,
                                           prepareTenor( inParams, lTerms ) );
}

}  // namespace APP

#endif