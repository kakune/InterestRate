/**
 * @file model_data.tpp
 * @brief This implements template classes of model data.
 * @author kakune
 * @date 4/21/2024
 */

#ifdef NINCLUDE_TPP
#include "LIBOR/forward/model_data.hpp"
#endif

#include "LIBOR/forward/payoff.hpp"
#include "math/findroot_1d.hpp"

namespace LIBOR::Forward
{

namespace States
{
template <typename ElementState_> class Plain
{
private:
    std::vector<std::vector<ElementState_>> mValues;

public:
    Plain( std::size_t inNPath, std::size_t inNTerms,
           ElementState_ inInitValue ) :
        mValues( inNPath, std::vector<ElementState_>( inNTerms, inInitValue ) )
    {
    }
    void setStateElement( std::size_t inIndPath, std::size_t inIndTerms,
                          ElementState_ inState )
    {
        mValues[inIndPath][inIndTerms] = inState;
    }
    void setForwardElement( std::size_t inIndPath, std::size_t inIndTerms,
                            ElementState_ inValue )
    {
        mValues[inIndPath][inIndTerms] = inValue;
    }
    const std::vector<std::vector<ElementState_>>& getStates() const
    {
        return mValues;
    }
    const std::vector<std::vector<ElementState_>>& getForwards() const
    {
        return mValues;
    }
};

template <C_LogExpElement ElementState_> class Log
{
private:
    std::vector<std::vector<ElementState_>> mValues;
    std::vector<std::vector<ElementState_>> mLogValues;

public:
    Log( std::size_t inNPath, std::size_t inNTerms,
         ElementState_ inInitValue ) :
        mValues( inNPath, std::vector<ElementState_>( inNTerms, inInitValue ) ),
        mLogValues( inNPath,
                    std::vector<ElementState_>( inNTerms, log( inInitValue ) ) )
    {
    }
    void setStateElement( std::size_t inIndPath, std::size_t inIndTerms,
                          ElementState_ inState )
    {
        mLogValues[inIndPath][inIndTerms] = inState;
        mValues[inIndPath][inIndTerms]    = exp( inState );
    }
    void setForwardElement( std::size_t inIndPath, std::size_t inIndTerms,
                            ElementState_ inValue )
    {
        mValues[inIndPath][inIndTerms]    = inValue;
        mLogValues[inIndPath][inIndTerms] = log( inValue );
    }
    const std::vector<std::vector<ElementState_>>& getStates() const
    {
        return mLogValues;
    }
    const std::vector<std::vector<ElementState_>>& getForwards() const
    {
        return mValues;
    }
};
}  // namespace States

namespace Data
{
template <class Derived>
Abstract<Derived>::Abstract(
    const Process::MarketData::Terms& inTerms,
    const Process::MarketData::Tenor& inTenor,
    std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
        insDataForwardRates,
    std::size_t inDegZCB ) :
    mNPath( insDataForwardRates->size() ),
    mTerms( inTerms ),
    mTenor( inTenor ),
    msDataForwardRates( insDataForwardRates ),
    msZCB( std::make_shared<const Process::MarketData::ZCB>(
        createZCBFromForwardRates( inTenor, *insDataForwardRates, inDegZCB ) ) )
{
}

template <class Derived>
Abstract<Derived>::Abstract(
    const Process::MarketData::Terms& inTerms,
    const Process::MarketData::Tenor& inTenor,
    const std::vector<std::vector<Math::Vec>>& inDataForwardRate,
    std::size_t inDegZCB ) :
    Abstract( inTerms, inTenor,
              std::make_shared<const std::vector<std::vector<Math::Vec>>>(
                  inDataForwardRate ),
              inDegZCB )
{
}

template <class Derived>
template <auto BlackPriceFunc_>
    requires C_Black76OneTermFunction<BlackPriceFunc_>
double Abstract<Derived>::calcBlackImpVolByOneTerm( double inStrike,
                                                    std::size_t inIndTenor,
                                                    double inModelPrice ) const
{
    Analytical::Black76::Model lBlackModel(
        { 0.0, mTenor.term( inIndTenor ), mTenor.term( inIndTenor + 1 ) } );
    lBlackModel.setInitZCB( { ( *msZCB )( mTenor.term( inIndTenor ) ),
                              ( *msZCB )( mTenor.term( inIndTenor + 1 ) ) } );
    auto lFuncDif = [&lBlackModel, inModelPrice,
                     inStrike]( double inVol ) -> double
    {
        lBlackModel.setVol( inVol );
        return inModelPrice - ( lBlackModel.*BlackPriceFunc_ )( inStrike, 1 );
    };
    return Math::FindRoot1D::Brent( lFuncDif, 1e-10, 1e2 );
}

template <class Derived>
double Abstract<Derived>::calcCaplet( double inStrike,
                                      std::size_t inIndTenor ) const
{
    CapletFloorletPayoff lCaplet( mTenor, msDataForwardRates, inStrike,
                                  inIndTenor, inIndTenor + 1, true );
    return static_cast<const Derived*>( this )->template calcExpectation(
        lCaplet );
}
template <class Derived>
double Abstract<Derived>::calcFloorlet( double inStrike,
                                        std::size_t inIndTenor ) const
{
    CapletFloorletPayoff lCaplet( mTenor, msDataForwardRates, inStrike,
                                  inIndTenor, inIndTenor + 1, false );
    return static_cast<const Derived*>( this )->template calcExpectation(
        lCaplet );
}
template <class Derived>
double Abstract<Derived>::calcBlackImpVolByCaplet(
    double inStrike, std::size_t inIndTenor ) const
{
    return this
        ->calcBlackImpVolByOneTerm<&Analytical::Black76::Model::priceCaplet>(
            inStrike, inIndTenor, calcCaplet( inStrike, inIndTenor ) );
}
template <class Derived>
double Abstract<Derived>::calcBlackImpVolByFloorlet(
    double inStrike, std::size_t inIndTenor ) const
{
    return this
        ->calcBlackImpVolByOneTerm<&Analytical::Black76::Model::priceFloorlet>(
            inStrike, inIndTenor, calcFloorlet( inStrike, inIndTenor ) );
}
template <class PayoffObject_>
double TerminalMeas::calcExpectation( PayoffObject_ inPayoff ) const
{
    double lResult           = 0.0;
    std::size_t lIndTenorPay = inPayoff.getIndexTenorPay();
    std::size_t lIndTermsPay = mTerms[lIndTenorPay];
    bool lIsTerminal         = ( lIndTenorPay == mTenor.size() );
    if ( lIsTerminal )
    {
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            lResult += inPayoff( iPath );
        }
        // discount by initial ZCB
        lResult *= ( *msZCB )( mTenor.term( mTenor.size() ) ) / mNPath;
        return lResult;
    }

    for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
    {
        double lValuePayoff = inPayoff( iPath );

        if ( lValuePayoff == 0.0 ) { continue; }

        // calculate using E_t[Payoff / P_{PayTime}(TerminalTime)]
        Math::Vec lZCB =
            mTenor.getTauVec() * ( *msDataForwardRates )[iPath][lIndTermsPay] +
            1.0;
        for ( std::size_t iZCB = lIndTenorPay; iZCB < mTenor.size(); ++iZCB )
        {
            lValuePayoff *= lZCB( iZCB );
        }
        lResult += lValuePayoff;
    }
    // discount by initial ZCB
    lResult *= ( *msZCB )( mTenor.term( mTenor.size() ) ) / mNPath;
    return lResult;
}

template <class PayoffObject_>
double SpotMeas::calcExpectation( PayoffObject_ inPayoff ) const
{
    double lResult           = 0.0;
    std::size_t lIndTenorPay = inPayoff.getIndexTenorPay();
    std::size_t lIndTermsPay = mTerms[lIndTenorPay];

    for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
    {
        double lValuePayoff = inPayoff( iPath );
        if ( lValuePayoff == 0.0 ) { continue; }
        double lFactor = 1.0;

        // calculate using E_t[Payoff / \Prod_{n=0}^{IndTenorPay-1} (1 +
        // \tau_n L_{T_n}(T_{n+1})]
        Math::Vec lRate =
            mTenor.getTauVec() * ( *msDataForwardRates )[iPath][lIndTermsPay] +
            1.0;
        for ( std::size_t iRate = 0; iRate < lIndTenorPay; ++iRate )
        {
            lFactor *= lRate( iRate );
        }
        lResult += lValuePayoff / lFactor;
    }
    return lResult / mNPath;
}

}  // namespace Data

}  // namespace LIBOR::Forward