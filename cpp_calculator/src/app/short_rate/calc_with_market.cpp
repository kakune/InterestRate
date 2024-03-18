#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "prepare_model.hpp"
#include "process/market.hpp"
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

    Process::Market::Data lMarketData(
        std::make_shared<std::vector<double> >( lMapMarket["Maturity"] ), 100,
        2 );
    lMarketData.setZCB( lMapMarket["ZCB"] );

    auto lsTerms = APP::ShortRate::prepareTerms( lParams );
    std::unique_ptr<Process::ShortRate::ModelAbstract> luModel =
        APP::ShortRate::prepareModelFromMarket( lNameModel, lParams, lsTerms,
                                                lMarketData );
    luModel->build();

    std::ofstream lFileOutput( lPathOutput );
    if ( lFileOutput.is_open() )
    {
        lFileOutput
            << "Start,Maturity,PriceZCB,ForwardRate,InstantaneousForwardRate"
            << std::endl;
        for ( std::size_t iStart = 0; iStart < lsTerms->size(); ++iStart )
        {
            double lTmpStartTime = lsTerms->operator[]( iStart );
            for ( std::size_t iMaturity = iStart + 1;
                  iMaturity < lsTerms->size(); ++iMaturity )
            {
                double lTmpMaturityTime = lsTerms->operator[]( iMaturity );
                lFileOutput
                    << std::setprecision( 12 ) << lTmpStartTime << ","
                    << lTmpMaturityTime << ","
                    << luModel->priceZCB( lTmpStartTime, lTmpMaturityTime )
                    << ","
                    << luModel->forwardRate( lTmpStartTime, lTmpMaturityTime )
                    << ","
                    << luModel->instantaneousForwardRate( lTmpMaturityTime )
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
