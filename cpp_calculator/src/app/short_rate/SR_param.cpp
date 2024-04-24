#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "short_rate_creator.hpp"
#include "terms_creator.hpp"
#include "utils/parameters.hpp"

int main( int argc, char* argv[] )
{
    std::string lPathOutput  = argv[1];
    std::string lNameSection = argv[2];
    std::string lPathParam   = argv[3];

    Utils::Parameters lParams;
    lParams.readParameters( lPathParam );
    lParams.setNameCommonSection( "COMMON" );
    lParams.setNameCurrentSection( lNameSection );

    Process::MarketData::Terms lTerms = APP::prepareTerms( lParams );
    ShortRate::SpotRates lSpots =
        APP::createSpotRateFromParam( lParams, lTerms );
    Process::MarketData::ZCB lZCB = lSpots.getZCB();

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
