/**
 * @file market.hpp
 * @brief This defines the struct that stores the market data
 * @author kakune
 * @date 3/8/2024
 */

#ifndef PROCESS_MARKET_HPP
#define PROCESS_MARKET_HPP

#include <limits>
#include <memory>
#include <vector>

#include "math/interpolate_1d.hpp"

namespace Process
{
namespace Market
{

struct Data
{
    std::shared_ptr<const std::vector<double> > msTerms;  //! term structure
    std::size_t mNMesh;  //! the number of mesh in making ZCB from forward rate
                         //! or vice versa.
    std::vector<double> mFineTerms;  //! term structure cut by mNMesh
    Math::Interpolate1d::NewtonSpline
        mInterpZCB;  //! interpolated function of ZCB
    Math::Interpolate1d::NewtonSpline
        mInterpInstantaneousForwardRate;  //! interpolated function of ZCB
    void setZCB( std::vector<double> inZCB );
    void setForwardRate( std::vector<double> inForwardRate );
    Data( std::shared_ptr<const std::vector<double> > insTerms,
          std::size_t inNMesh = 10 ) :
        msTerms( insTerms ),
        mNMesh( inNMesh ),
        mFineTerms( ( insTerms->size() - 1 ) * inNMesh + 1 ),
        mInterpZCB( 3 ),
        mInterpInstantaneousForwardRate( 3 )
    {
        for ( std::size_t iTerm = 1; iTerm < msTerms->size(); ++iTerm )
        {
            double lTmpDif =
                ( msTerms->at( iTerm ) - msTerms->at( iTerm - 1 ) ) /
                double( mNMesh );
            std::size_t lStartInd      = ( iTerm - 1 ) * mNMesh;
            mFineTerms.at( lStartInd ) = msTerms->at( iTerm - 1 );
            for ( std::size_t iMesh = 1; iMesh < mNMesh; ++iMesh )
            {
                mFineTerms.at( lStartInd + iMesh ) =
                    mFineTerms.at( lStartInd + iMesh - 1 ) + lTmpDif;
            }
        }
        mFineTerms.back() = msTerms->back();
    }
};

}  // namespace Market
}  // namespace Process

#endif