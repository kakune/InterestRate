#include <iomanip>
#include <iostream>
#include <vector>

#include "process/asset.hpp"
#include "utils/parameters.hpp"

int main( int argc, char* argv[] )
{
    std::string lFilePath    = argv[1];
    std::string lSectionName = argv[2];
    Utils::Parameters lParams;
    lParams.readParameters( lFilePath );
    lParams.setCommonSectionName( "COMMON" );
    lParams.setCurrentSectionName( lSectionName );
    std::size_t lNTerms = lParams( "NTerms" );
    std::size_t lNPath  = 1000000;
    double lDt          = 1.0 / double( lNTerms );
    std::vector< double > lTerms( lNTerms, 0 );
    for ( std::size_t iTerm = 1; iTerm < lNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    auto lsTerms = std::make_shared< std::vector< double > >( lTerms );

    Process::Asset::SABRWithLogForwardBuilder lSABRBuilder;
    lSABRBuilder.setNPath( lParams( "NPath" ) );
    lSABRBuilder.setTerms( lsTerms );
    lSABRBuilder.setInitVol( lParams( "InitVol" ) );
    lSABRBuilder.setInitPrice( lParams( "InitPrice" ) );
    lSABRBuilder.setVolvol( lParams( "Volvol" ) );
    lSABRBuilder.setExponent( lParams( "Exponent" ) );
    lSABRBuilder.setCorr( lParams( "Corr" ) );

    auto lSABRObj = lSABRBuilder.build();

    for ( double lStrike = 40.0; lStrike < 160.0; lStrike += 1.0 )
    {
        std::cout << lStrike << "," << std::setprecision( 12 )
                  << lSABRObj.impliedVolatility( lStrike, lNTerms - 1 )
                  << std::endl;
    }
    std::cout << std::endl;
    return 0;
}