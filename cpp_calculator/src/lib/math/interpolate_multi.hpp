/**
 * @file interpolate_multi.hpp
 * @brief This defines interpolation function for multi dimensional dataset.
 * @author kakune
 * @date 3/13/2024
 */

#ifndef MATH_INTERPOLATE_MULTI_HPP
#define MATH_INTERPOLATE_MULTI_HPP

#include <omp.h>

#include <memory>
#include <vector>

#include "math/matrix.hpp"

namespace Math
{
namespace InterpolateMulti
{

class RBFAbstract
{
protected:
    bool mIsBuilt;
    std::size_t mNCoeff;
    std::size_t mNDim;
    std::vector<Math::Vec> mPoints;
    Math::Vec mCoeffs;
    double mFactorDistance;
    double mMinDistance;
    double mFactorDecay;
    virtual double distance( const Math::Vec& inX1,
                             const Math::Vec& inX2 ) const;
    virtual double radial( double inDist ) const = 0;

public:
    RBFAbstract( double inMinDistance = 0.2, double inFactorDecay = 0.001 ) :
        mFactorDistance( 1.0 ),
        mMinDistance( inMinDistance ),
        mFactorDecay( inFactorDecay ),
        mPoints( 0, Math::Vec( 0 ) ),
        mCoeffs( 0 )
    {
    }
    void build(
        std::vector<std::shared_ptr<const std::vector<double>>> insRefVars,
        std::shared_ptr<const std::vector<double>> insRefVals );

    double operator()( const std::vector<double>& inVar ) const;
};

class RBFGaussian : public RBFAbstract
{
private:
    double radial( double inDist ) const override;

public:
    RBFGaussian( double inMinDistance = 0.2, double inFactorDecay = 0.001 ) :
        RBFAbstract( inMinDistance, inFactorDecay )
    {
    }
};

}  // namespace InterpolateMulti
}  // namespace Math

#endif