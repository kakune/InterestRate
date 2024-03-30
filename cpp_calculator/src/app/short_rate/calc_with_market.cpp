#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "process/market_data.hpp"
#include "process/model_data.hpp"
#include "spot_calculator.hpp"
#include "utils/csv.hpp"
#include "utils/parameters.hpp"

int main( int argc, char* argv[] )
{
    std::string lPathOutput        = argv[1];
    std::string lNameModel         = argv[2];
    std::string lNameSection       = argv[3];
    std::string lPathParam         = argv[4];
    std::string lNameSectionMarket = argv[5];
    std::string lPathParamMarket   = argv[6];
    std::string lPathCSVMarket     = argv[7];

    Utils::Parameters lParams;
    lParams.readParameters( lPathParam );
    lParams.setNameCommonSection( "COMMON" );
    lParams.setNameCurrentSection( lNameSection );
    Utils::Parameters lParamsMarket;
    lParamsMarket.readParameters( lPathParamMarket );
    lParamsMarket.setNameCommonSection( "COMMON" );
    lParamsMarket.setNameCurrentSection( lNameSectionMarket );

    auto lMapMarket = Utils::CSV::readFile( lPathCSVMarket );

    for ( std::size_t i = 0; i < lMapMarket["Maturity"].size(); ++i )
    {
        lMapMarket["ZCB"].push_back(
            pow( 1.0 + 0.01 * lMapMarket["ZCBRate"][i],
                 -lMapMarket["Maturity"][i] / lParamsMarket( "Compound" ) ) );
    }

    Process::MarketData::ZCB lMarketZCB(
        Process::MarketData::Terms( lMapMarket["Maturity"] ),
        lMapMarket["ZCB"] );

    Process::MarketData::Terms lTerms = APP::prepareTerms( lParams );
    Process::ModelData::SpotRates lSpots =
        APP::calcSpotRateFromMarket( lNameModel, lParams, lTerms, lMarketZCB );
    Process::MarketData::ZCB lZCB = lSpots.createZCB();

    std::ofstream lFileOutput( lPathOutput );
    if ( lFileOutput.is_open() )
    {
        lFileOutput
            << "Start,Maturity,PriceZCB,ForwardRate,InstantaneousForwardRate"
            << std::endl;
        for ( std::size_t iStart = 0; iStart < lTerms.size(); ++iStart )
        {
            double lTmpStartTime = lTerms[iStart];
            for ( std::size_t iMaturity = iStart + 1; iMaturity < lTerms.size();
                  ++iMaturity )
            {
                double lTmpMaturityTime = lTerms[iMaturity];
                lFileOutput
                    << std::setprecision( 12 ) << lTmpStartTime << ","
                    << lTmpMaturityTime << ","
                    << lZCB( lTmpStartTime, lTmpMaturityTime ) << ","
                    << lZCB.forwardRate( lTmpStartTime, lTmpMaturityTime )
                    << "," << lZCB.instantaneousForwardRate( lTmpMaturityTime )
                    << std::endl;
            }
        }
        lFileOutput.close();
    }
    else
    {
        std::cerr << "Could not open parameter file: " << lPathOutput
                  << std::endl;
    }
}
