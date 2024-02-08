#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "process/asset.hpp"
#include "utils/parameters.hpp"

namespace Process
{
namespace Asset
{
class ExerciseWithLogForward : public StochasticVolatilityWithLogForwardAbstract
{
private:
    double mExponent;
    double mShift;
    double mStandard;
    double mShiftExp;
    double mVolvol;
    double volForward( double inPrice, double inVol,
                       double inTime = 0.0 ) override;
    double volVolatility( double inVol, double inTime = 0.0 ) override;
    double driftForward( double inPrice, double inVol,
                         double inTime = 0.0 ) override;
    double driftVolatility( double inVol, double inTime = 0.0 ) override;

public:
    ExerciseWithLogForward(
        std::size_t inNPath,
        std::shared_ptr< const std::vector< double > > insTerms,
        double inInitPrice,
        std::unique_ptr< Process::Random::PathAbstract > inuRandomPath,
        double inInitVol,
        std::unique_ptr< Process::Random::PathAbstract > inuRandomVol,
        double inCorr, double inExponent, double inShift, double inStandard,
        double inVolvol ) :
        StochasticVolatilityWithLogForwardAbstract(
            inNPath, insTerms, inInitPrice, std::move( inuRandomPath ),
            inInitVol, std::move( inuRandomVol ), inCorr ),
        mExponent( inExponent ),
        mShift( inShift ),
        mStandard( inStandard ),
        mShiftExp( std::pow( inShift, mExponent ) ),
        mVolvol( inVolvol )
    {
        calcEachForwardPrice();
    }
};

class ExerciseWithLogForwardBuilder
    : public StochasticVolatilityForwardAbstractBuilder
{
private:
    double mExponent;
    double mShift;
    double mStandard;
    double mVolvol;

public:
    ExerciseWithLogForwardBuilder& setExponent( double inExponent )
    {
        mExponent = inExponent;
        return *this;
    }
    ExerciseWithLogForwardBuilder& setShift( double inShift )
    {
        mShift = inShift;
        return *this;
    }
    ExerciseWithLogForwardBuilder& setStandard( double inStandard )
    {
        mStandard = inStandard;
        return *this;
    }
    ExerciseWithLogForwardBuilder& setVolvol( double inVolvol )
    {
        mVolvol = inVolvol;
        return *this;
    }
    ExerciseWithLogForward build()
    {
        if ( muRandomPath == nullptr )
        {
            muRandomPath =
                std::make_unique< Process::Random::PathBrownAntithetic >(
                    mNPath, msTerms );
        }
        if ( muRandomVol == nullptr )
        {
            muRandomVol =
                std::make_unique< Process::Random::PathBrownAntithetic >(
                    mNPath, msTerms );
        }
        return ExerciseWithLogForward( mNPath, msTerms, mInitPrice,
                                       std::move( muRandomPath ), mInitVol,
                                       std::move( muRandomVol ), mCorr,
                                       mExponent, mShift, mStandard, mVolvol );
    }
};
double ExerciseWithLogForward::volForward( double inPrice, double inVol,
                                           double inTime )
{
    double lX = std::pow( inPrice, mExponent );
    return inVol * mStandard * lX / ( ( lX + mShiftExp ) * inPrice );
}
double ExerciseWithLogForward::volVolatility( double inVol, double inTime )
{
    return mVolvol;
}
double ExerciseWithLogForward::driftForward( double inPrice, double inVol,
                                             double inTime )
{
    double lX   = std::pow( inPrice, mExponent );
    double lTmp = inVol * mStandard * lX / ( ( lX + mShiftExp ) * inPrice );
    return -0.5 * lTmp * lTmp;
}
double ExerciseWithLogForward::driftVolatility( double inVol, double inTime )
{
    return -0.5 * mVolvol * mVolvol;
}

}  // namespace Asset
}  // namespace Process

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
    std::vector< double > lTerms( lNTerms, 0 );

    double lDt = lParams( "TimeMaturity" ) / double( lNTerms );
    for ( std::size_t iTerm = 1; iTerm < lNTerms; ++iTerm )
    {
        lTerms[iTerm] = lTerms[iTerm - 1] + lDt;
    }

    auto lsTerms = std::make_shared< std::vector< double > >( lTerms );

    Process::Asset::ExerciseWithLogForwardBuilder lExerciseBuilder;
    lExerciseBuilder.setNPath( lParams( "NPath" ) );
    lExerciseBuilder.setTerms( lsTerms );
    lExerciseBuilder.setInitPrice( lParams( "InitPrice" ) );
    lExerciseBuilder.setInitVol( lParams( "InitVol" ) );
    lExerciseBuilder.setCorr( lParams( "Corr" ) );
    lExerciseBuilder.setExponent( lParams( "Exponent" ) );
    lExerciseBuilder.setShift( lParams( "Shift" ) );
    lExerciseBuilder.setStandard( lParams( "Standard" ) );
    lExerciseBuilder.setVolvol( lParams( "Volvol" ) );

    auto lExerciseObj = lExerciseBuilder.build();

    std::ofstream lFileOutput( lPathOutput );

    if ( lFileOutput.is_open() )
    {
        lFileOutput << "Strike,ImpVol" << std::endl;
        for ( double lStrike = lParams( "MinStrike" );
              lStrike < lParams( "MaxStrike" );
              lStrike += lParams( "DStrike" ) )
        {
            lFileOutput << lStrike << "," << std::setprecision( 12 )
                        << lExerciseObj.impliedVolatility( lStrike,
                                                           lNTerms - 1 )
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