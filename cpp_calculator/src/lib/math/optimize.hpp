/**
 * @file optimize.hpp
 * @brief This defines optimizers.
 * @author kakune
 * @date 4/30/2024
 */

#ifndef MATH_OPTIMIZE_HPP
#define MATH_OPTIMIZE_HPP

#include <concepts>

#include "math/matrix.hpp"

namespace Math::Optimize
{

template <typename Func_>
concept C_ObjectiveFunction = requires( Func_ inFunc, const Math::Vec& inVec ) {
    { inFunc( inVec ) } -> std::convertible_to<double>;
};

/**
 * @brief This try to find Vec that minimizes inObjectiveFunc.
 * @param inObjectiveFunc objective function
 * @param inInitialValue vector of initial value of finding
 * @param inAlpha
 * @param inBeta1
 * @param inBeta2
 * @param inEpsilon
 * @param inMaxIter
 * @param inEpsGradNorm
 * @return Math::Vec found Vec minimizes inObjectiveFunc
 */
Math::Vec adam( C_ObjectiveFunction auto inObjectiveFunc,
                const Math::Vec& inInitialValue, double inAlpha = 0.01,
                double inBeta1 = 0.9, double inBeta2 = 0.999,
                double inEpsilon = 1e-8, std::size_t inMaxIter = 10000000,
                double inEpsGradNorm = 1e-12 );

}  // namespace Math::Optimize

#ifndef NINCLUDE_TPP
#include "math/optimize.tpp"
#endif

#endif