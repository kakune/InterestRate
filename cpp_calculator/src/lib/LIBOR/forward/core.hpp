/**
 * @file core.hpp
 * @brief This defines class for forward model.
 * @author kakune
 * @date 3/29/2024
 */

#ifndef MARKET_MODEL_FORWARD_CORE_HPP
#define MARKET_MODEL_FORWARD_CORE_HPP

#include <memory>
#include <vector>

#include "math/matrix.hpp"
#include "process/market_data.hpp"
#include "process/model_data.hpp"
#include "process/random.hpp"

namespace LIBOR
{
namespace Forward
{

/**
 * @brief This is the abstract class for market forward rate models.
 */
class ModelAbstract
{
protected:
    std::size_t mNFR;                    //! the number of forward rate
    std::size_t mNPath;                  //! the number of Path
    Process::MarketData::Terms mTerms;   //! terms
    std::vector<std::size_t> mIndTenor;  //! indices of tenor
    std::vector<double> mTauTenor;       //! the intervals of tenor
    Math::Vec mInitFRs;                  //! initial forward rate
    std::unique_ptr<Process::RandomVec::StdBrownAbstract>
        muStdBrown;   //! random vec
    Math::Mat mCorr;  //! model correlation
    Math::Mat mRho;   //! correlation of brownian motion

    /**
     * @brief The drift term in SDE of f[inIndPath][inIndTerm]
     * @param inIndPath the index of path
     * @param inIndTerm the index of term
     * @param inStates the set of states
     * @return Math::Vec the coefficient
     */
    virtual Math::Vec driftTerm(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<Math::Vec>>& inStates ) const = 0;

    /**
     * @brief The shock due to dW in SDE of f[inIndPath][inIndTerm]
     * @param inIndPath the index of path
     * @param inIndTerm the index of term
     * @param inStates the set of states
     * @param inRandomVec the vector of dW
     * @return Math::Vec the shock
     */
    virtual Math::Vec volTerm(
        std::size_t inIndPath, std::size_t inIndTerm,
        const std::vector<std::vector<Math::Vec>>& inStates,
        const Math::Vec& inRandomVec ) const = 0;

    /**
     * @brief transf function FR to State in SDE.
     * @param inFR forward rate
     * @param inIndTime index of time
     * @return Math::Vec state vector
     */
    virtual Math::Vec transfFRToState( const Math::Vec& inFR,
                                       std::size_t inIndTime ) const;

    /**
     * @brief transf function State to FR.
     * @param inState state vector in SDE
     * @param inIndTime index of time
     * @return Math::Vec forward rate
     */
    virtual Math::Vec transfStateToFR( const Math::Vec& inState,
                                       std::size_t inIndTime ) const;

    std::vector<double> calcTauTenor(
        const std::vector<std::size_t>& inIndTenor,
        const Process::MarketData::Terms& inTerms ) const;

public:
    ModelAbstract(
        std::size_t inNPath, const Process::MarketData::Terms& inTerms,
        const std::vector<std::size_t>& inIndTenor, const Math::Vec& inInitFRs,
        std::unique_ptr<Process::RandomVec::StdBrownAbstract> inuStdBrown,
        const Math::Mat& inCorr ) :
        mNFR( inIndTenor.size() - 1 ),
        mNPath( inNPath ),
        mTerms( inTerms ),
        mIndTenor( inIndTenor ),
        mTauTenor( calcTauTenor( inIndTenor, inTerms ) ),
        mInitFRs( inInitFRs ),
        muStdBrown( std::move( inuStdBrown ) ),
        mCorr( inCorr ),
        mRho( dot( inCorr.transpose(), inCorr ) )
    {
    }
    virtual ~ModelAbstract() = default;
    virtual Process::ModelData::ForwardRates createForwardRates() const;
};

/**
 * @brief This build the object of ModelAbstract.
 */
class ModelAbstractBuilder
{
protected:
    std::size_t mNPath;                                   //! the number of Path
    std::unique_ptr<Process::MarketData::Terms> muTerms;  //! term structure
    std::vector<std::size_t> mIndTenor;                   //! indices of tenor
    std::unique_ptr<Process::RandomVec::StdBrownAbstract>
        muStdBrown;  //! random path
    Math::Vec mInitFRs = Math::Vec( 0 );
    Math::Mat mCorr    = Math::Mat( 0, 0 );

public:
    ModelAbstractBuilder& setNPath( std::size_t inNPath )
    {
        mNPath = inNPath;
        return *this;
    }
    ModelAbstractBuilder& setTerms(
        std::shared_ptr<const std::vector<double>> insTerms )
    {
        muTerms = std::make_unique<Process::MarketData::Terms>( insTerms );
        return *this;
    }
    ModelAbstractBuilder& setTerms( const Process::MarketData::Terms& inTerms )
    {
        muTerms = std::make_unique<Process::MarketData::Terms>( inTerms );
        return *this;
    }
    ModelAbstractBuilder& setIndTenor(
        const std::vector<std::size_t>& inIndTenor )
    {
        mIndTenor = inIndTenor;
        return *this;
    }
    ModelAbstractBuilder& setInitFRs( const Math::Vec& inInitFRs )
    {
        mInitFRs = inInitFRs;
        return *this;
    }
    ModelAbstractBuilder& setRandom(
        std::unique_ptr<Process::RandomVec::StdBrownAbstract> inuStdBrown,
        const Math::Mat& inCorr )
    {
        muStdBrown = std::move( inuStdBrown );
        mCorr      = inCorr;
        return *this;
    }
    virtual ~ModelAbstractBuilder() = default;
};

}  // namespace Forward

}  // namespace LIBOR

#endif