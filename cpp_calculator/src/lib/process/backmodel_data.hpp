/**
 * @file model_data.hpp
 * @brief This defines classes for model data. Unless there is a special
 * reason, calculation results should be stored in these classes. Since most
 * data is stored by smart pointers, there is no need to worry about the cost of
 * copying instances.
 * @author kakune
 * @date 3/30/2024
 */

#ifndef PROCESS_MODEL_DATA_HPP
#define PROCESS_MODEL_DATA_HPP

#include <type_traits>

#include "analytical/Black76.hpp"
#include "math/findroot_1d.hpp"
#include "process/market_data.hpp"

namespace Process
{
namespace ModelData
{

class CapletFloorletPayoff
{
private:
    const std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
        msDataForwardRates;
    const double mStrike;
    const std::size_t mIndTermsFR;
    const std::size_t mIndTenorFR;
    const std::size_t mIndTenorPay;
    const double mTau;
    bool mIsCaplet;

public:
    CapletFloorletPayoff(
        const Process::MarketData::Tenor& inTenor,
        std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
            insDataForwardRates,
        double inStrike, std::size_t inIndTenorFR, std::size_t inIndTenorPay,
        bool inIsCaplet ) :
        msDataForwardRates( insDataForwardRates ),
        mStrike( inStrike ),
        mIndTermsFR( inTenor[inIndTenorFR] ),
        mIndTenorFR( inIndTenorFR ),
        mIndTenorPay( inIndTenorPay ),
        mTau( inTenor.tau( inIndTenorFR ) ),
        mIsCaplet( inIsCaplet )
    {
    }
    std::size_t getIndexTenorPay() const { return mIndTenorPay; }
    double operator()( std::size_t inIndPath ) const
    {
        if ( mIsCaplet )
        {
            return std::max(
                0.0, mTau * ( ( *msDataForwardRates )[inIndPath][mIndTermsFR](
                                  mIndTenorFR ) -
                              mStrike ) );
        }
        return std::max(
            0.0,
            mTau * ( mStrike - ( *msDataForwardRates )[inIndPath][mIndTermsFR](
                                   mIndTenorFR ) ) );
    }
};

MarketData::ZCB createZCBFromForwardRates(
    const MarketData::Tenor& inTenor,
    const std::vector<std::vector<Math::Vec>>& inFRs, std::size_t inDeg );

/**
 * @brief This stores forward rate data at each path and each term.
 */
template <class Derived> class ForwardRatesAbstract
{
protected:
    const std::size_t mNPath;        //! the number of path
    const MarketData::Terms mTerms;  //! terms
    const MarketData::Tenor mTenor;
    const std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
        msDataForwardRates;  //! forward rates
    const std::shared_ptr<const MarketData::ZCB> msZCB;

public:
    /**
     * @brief This constructs a new Forward Rates.
     * @param inTerms Term structure
     * @param inIndTenor indices of tenor of FR
     * @param insDataForwardRates forward rate
     */
    ForwardRatesAbstract(
        const MarketData::Terms& inTerms, const MarketData::Tenor& inTenor,
        std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
            insDataForwardRates,
        std::size_t inDegZCB = 3 ) :
        mNPath( insDataForwardRates->size() ),
        mTerms( inTerms ),
        mTenor( inTenor ),
        msDataForwardRates( insDataForwardRates ),
        msZCB(
            std::make_shared<const MarketData::ZCB>( createZCBFromForwardRates(
                inTenor, *insDataForwardRates, inDegZCB ) ) )
    {
    }
    ForwardRatesAbstract(
        const MarketData::Terms& inTerms, const MarketData::Tenor& inTenor,
        const std::vector<std::vector<Math::Vec>>& inDataForwardRate,
        std::size_t inDegZCB = 3 ) :
        ForwardRatesAbstract(
            inTerms, inTenor,
            std::make_shared<const std::vector<std::vector<Math::Vec>>>(
                inDataForwardRate ),
            inDegZCB )
    {
    }
    const std::vector<Math::Vec>& operator[]( std::size_t inIndex ) const
    {
        return ( *msDataForwardRates )[inIndex];
    }
    double term( std::size_t inIndex ) const { return mTerms[inIndex]; }
    const MarketData::Terms& getTerms() const { return mTerms; }
    const MarketData::ZCB& getZCB() const { return *msZCB; }
    std::size_t sizeTerms() const { return mTerms.size(); }
    std::size_t sizePath() const { return mNPath; }

    double calcCaplet( double inStrike, std::size_t inIndTenor ) const
    {
        CapletFloorletPayoff lCaplet( mTenor, msDataForwardRates, inStrike,
                                      inIndTenor, inIndTenor + 1, true );
        return static_cast<const Derived*>( this )->template calcExpectation(
            lCaplet );
    }
    double calcFloorlet( double inStrike, std::size_t inIndTenor ) const
    {
        CapletFloorletPayoff lCaplet( mTenor, msDataForwardRates, inStrike,
                                      inIndTenor, inIndTenor + 1, false );
        return static_cast<const Derived*>( this )->template calcExpectation(
            lCaplet );
    }
    double calcBlackImpVol( double inStrike, std::size_t inIndTenor,
                            bool inIsUseCaplet = true )
    {
        double lInitFR    = ( *msDataForwardRates )[0][0]( inIndTenor );
        double lTimeStart = mTenor.term( inIndTenor );
        double lTau       = mTenor.tau( inIndTenor );
        double lZCB       = ( *msZCB )( mTenor.term( inIndTenor + 1 ) );
        double lPrice     = ( inIsUseCaplet )
                                ? calcCaplet( inStrike, inIndTenor )
                                : calcFloorlet( inStrike, inIndTenor );
        Analytical::Black76::Model lBlackModel{ lInitFR, inStrike, lTimeStart,
                                                lTau,    0.0,      lZCB };
        auto lFuncDif = [this, &lBlackModel, lPrice,
                         inIsUseCaplet]( double inVol ) -> double
        {
            lBlackModel.mVol = inVol;
            return lPrice - lBlackModel( inIsUseCaplet );
        };
        return Math::FindRoot1D::Brent( lFuncDif, 1e-10, 1e2 );
    }
};

class ForwardRatesTerminalMeas
    : public ForwardRatesAbstract<ForwardRatesTerminalMeas>
{
public:
    ForwardRatesTerminalMeas(
        const MarketData::Terms& inTerms, const MarketData::Tenor& inTenor,
        std::shared_ptr<const std::vector<std::vector<Math::Vec>>>
            insDataForwardRates,
        std::size_t inDegZCB = 3 ) :
        ForwardRatesAbstract( inTerms, inTenor, insDataForwardRates, inDegZCB )
    {
    }
    ForwardRatesTerminalMeas(
        const MarketData::Terms& inTerms, const MarketData::Tenor& inTenor,
        const std::vector<std::vector<Math::Vec>>& inDataForwardRate,
        std::size_t inDegZCB = 3 ) :
        ForwardRatesAbstract( inTerms, inTenor, inDataForwardRate, inDegZCB )
    {
    }
    template <class PayoffObject_>
    double calcExpectation( PayoffObject_ inPayoff ) const
    {
        double lResult           = 0.0;
        std::size_t lIndTenorPay = inPayoff.getIndexTenorPay();
        std::size_t lIndTermsPay = mTerms[lIndTenorPay];
        bool lIsTerminal         = ( lIndTenorPay == mTenor.size() );
        for ( std::size_t iPath = 0; iPath < mNPath; ++iPath )
        {
            double lValuePayoff = inPayoff( iPath );
            if ( !lIsTerminal && lValuePayoff != 0.0 )
            {
                Math::Vec lZCB =
                    mTenor.getTauVec() *
                        ( *msDataForwardRates )[iPath][lIndTermsPay] +
                    1.0;
                for ( std::size_t iZCB = lIndTenorPay; iZCB < mTenor.size();
                      ++iZCB )
                {
                    lValuePayoff *= lZCB( iZCB );
                }
            }
            lResult += lValuePayoff;
        }
        lResult *= ( *msZCB )( mTenor.term( mTenor.size() ) ) / mNPath;
        return lResult;
    }
};

template <typename T_, typename ElementState_>
concept C_State = requires( T_ inObj, std::size_t inIndPath,
                            std::size_t inIndTerms, ElementState_ inVal ) {
    inObj.setStateElement( inIndPath, inIndTerms, inVal );
    inObj.setValueElement( inIndPath, inIndTerms, inVal );
    {
        inObj.getStates()
    } -> std::same_as<const std::vector<std::vector<ElementState_>>&>;
    {
        inObj.getValues()
    } -> std::same_as<const std::vector<std::vector<ElementState_>>&>;
};
template <typename T_>
concept C_LogExpElement = requires( T_ inObj ) {
    log( inObj );
    exp( inObj );
};

template <typename ElementState_> class StatesPlain;
template <C_LogExpElement ElementState_> class StatesLog;

}  // namespace ModelData
}  // namespace Process

#include "process/model_data.tpp"

#endif