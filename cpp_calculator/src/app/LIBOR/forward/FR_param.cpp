#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "LIBOR_forward_creator.hpp"
#include "process/market_data.hpp"
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
    Process::MarketData::Tenor lTenor = APP::prepareTenor( lParams, lTerms );
    auto lFR = APP::createForwardTerminalFromParam( lParams, lTerms, lTenor );

    std::ofstream lFileOutput( lPathOutput );
    lFileOutput << std::fixed << std::setprecision( 12 );
    std::vector<double> lStrikes =
        lParams.operator()<std::vector<double>>( "Strikes" );
    if ( lFileOutput.is_open() )
    {
        lFileOutput << "IndTenor,Start,Maturity,Strike,ImpVol" << std::endl;
        for ( std::size_t i = 1; i < lTenor.size(); ++i )
        {
            for ( double lStrike : lStrikes )
            {
                lFileOutput << i;
                lFileOutput << "," << lTenor.term( i );
                lFileOutput << "," << lTenor.term( i + 1 );
                lFileOutput << "," << lStrike;
                lFileOutput << "," << lFR.calcBlackImpVolByCaplet( lStrike, i );
                lFileOutput << std::endl;
            }
        }
        lFileOutput.close();
    }
    else
    {
        std::cerr << "Could not open parameter file: " << lPathOutput
                  << std::endl;
    }
    return 0;
}
