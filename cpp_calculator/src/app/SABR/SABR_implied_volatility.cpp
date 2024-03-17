#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "process/asset.hpp"
#include "utils/parameters.hpp"

int main( int argc, char* argv[] )
{
    std::string lPathParam   = argv[1];
    std::string lSectionName = argv[2];
    std::string lPathOutput  = argv[3];

    Utils::Parameters lParams;
    lParams.readParameters( lPathParam );
    lParams.setNameCommonSection( "COMMON" );
    lParams.setNameCurrentSection( lSectionName );

    std::size_t lNTerms = std::size_t( lParams( "NTerms" ) );
    std::vector<double> lTerms( lNTerms, 0 );

    double lDt = lParams( "TimeMaturity" ) / double( lNTerms - 1 );
    for ( std::size_t iTerm = 1; iTerm < lNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }

    auto lsTerms = std::make_shared<std::vector<double> >( lTerms );

    Process::Asset::SABRWithLogForwardBuilder lSABRBuilder;
    lSABRBuilder.setNPath( lParams( "NPath" ) );
    lSABRBuilder.setTerms( lsTerms );
    lSABRBuilder.setInitPrice( lParams( "InitPrice" ) );
    lSABRBuilder.setInitVol( lParams( "InitVol" ) );
    lSABRBuilder.setCorr( lParams( "Corr" ) );
    lSABRBuilder.setExponent( lParams( "Exponent" ) );
    lSABRBuilder.setVolvol( lParams( "Volvol" ) );

    auto lSABRObj = lSABRBuilder.build();

    std::ofstream lFileOutput( lPathOutput );

    if ( lFileOutput.is_open() )
    {
        lFileOutput << "Strike,ImpVol" << std::endl;
        for ( double lStrike = lParams( "MinStrike" );
              lStrike < lParams( "MaxStrike" );
              lStrike += lParams( "DStrike" ) )
        {
            lFileOutput << lStrike << "," << std::setprecision( 12 )
                        << lSABRObj.impliedVolatility( lStrike, lNTerms - 1 )
                        << std::endl;
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