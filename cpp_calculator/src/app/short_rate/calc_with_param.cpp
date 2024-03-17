#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "prepare_model.hpp"
#include "utils/parameters.hpp"

int main( int argc, char* argv[] )
{
    std::string lPathParam   = argv[1];
    std::string lNameSection = argv[2];
    std::string lNameModel   = argv[3];
    std::string lPathOutput  = argv[4];

    Utils::Parameters lParams;
    lParams.readParameters( lPathParam );
    lParams.setNameCommonSection( "COMMON" );
    lParams.setNameCurrentSection( lNameSection );

    auto lsTerms = APP::ShortRate::prepareTerms( lParams );
    std::unique_ptr<Process::ShortRate::ModelAbstract> luModel =
        APP::ShortRate::prepareModelFromParam( lNameModel, lParams, lsTerms );
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
