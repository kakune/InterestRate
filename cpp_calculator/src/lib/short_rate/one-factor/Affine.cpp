/**
 * @file Affine.cpp
 * @brief This implements one-factor Affine short-rate models. (see 3.2.4 of
 * Brigo)
 * @author kakune
 * @date 3/14/2024
 */

#include "short_rate/one-factor/Affine.hpp"

#include <vector>

namespace ShortRate
{
namespace OneFactor
{

double ConstantAffine::driftCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots ) const
{
    return mLambda * inSpots[inIndPath][inIndTerm - 1] + mEta;
}
double ConstantAffine::volCoeff(
    std::size_t inIndPath, std::size_t inIndTerm,
    const std::vector<std::vector<double>>& inSpots ) const
{
    return std::sqrt( mGamma * inSpots[inIndPath][inIndTerm - 1] + mDelta );
}

}  // namespace OneFactor
}  // namespace ShortRate
