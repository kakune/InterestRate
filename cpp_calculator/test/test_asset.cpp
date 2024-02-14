#include <iomanip>
#include <iostream>
#include <vector>

#include "process/asset.hpp"

int main( int argc, char* argv[] )
{
    std::size_t lNTerm = 10;
    std::size_t lNPath = 100000000;
    double lDt         = 0.1;
    double lVol        = 0.1;
    double lInitPrice  = 100.0;
    std::vector<double> lTerms( lNTerm, 0 );
    for ( std::size_t iTerm = 1; iTerm < lNTerm; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }
    auto lsTerms = std::make_shared<std::vector<double> >( lTerms );

    Process::Asset::BlackSholesForwardBuilder lBSBuilder;
    lBSBuilder.setNPath( lNPath );
    lBSBuilder.setTerms( lsTerms );
    lBSBuilder.setVol( lVol );
    lBSBuilder.setInitPrice( lInitPrice );

    auto lBSObj = lBSBuilder.build();

    std::cout << std::endl;
    std::cout << lBSObj.impliedVolatility( 100.0, lNTerm - 1 ) << std::endl;
    // for ( double lStrike = 40.0; lStrike < 160.0; lStrike += 1.0 )
    // {
    //     std::cout << lStrike << " " << std::setprecision( 12 )
    //               << lBSObj.impliedVolatility( lStrike, lNTerms - 1 )
    //               << std::endl;
    // }
    // std::cout << std::endl;

    // double lVolvol   = 0.8;
    // double lExponent = 0.5;
    // double lCorr     = -0.0;
    // Process::Asset::SABRWithLogForwardBuilder lSABRBuilder;
    // lSABRBuilder.setNPath( lNPath );
    // lSABRBuilder.setTerms( lsTerms );
    // lSABRBuilder.setInitVol( lVol );
    // lSABRBuilder.setInitPrice( lInitPrice );
    // lSABRBuilder.setVolvol( lVolvol );
    // lSABRBuilder.setExponent( lExponent );
    // lSABRBuilder.setCorr( lCorr );

    // auto lSABRObj = lSABRBuilder.build();

    // std::cout << std::endl;
    // for ( double lStrike = 40.0; lStrike < 160.0; lStrike += 1.0 )
    // {
    //     std::cout << lStrike << " " << std::setprecision( 12 )
    //               << lSABRObj.impliedVolatility( lStrike, lNTerms - 1 )
    //               << std::endl;
    // }
    // std::cout << std::endl;
    return 0;
}