/**
 * @file integral_1d.hpp
 * @brief This defines the 1d-integral functions
 * @author kakune
 * @date 3/8/2024
 */

#ifndef MATH_INTEGRAL_1D_HPP
#define MATH_INTEGRAL_1D_HPP
#include <cstddef>

namespace Math::Integral
{
namespace FiniteInterval
{

template <typename TypeX_ = double, typename TypeY_ = double>
TypeY_ trapz( auto inFunc, TypeX_ inMin, TypeX_ inMax,
              std::size_t inNDivision );

template <typename TypeY_ = double>
TypeY_ DEFormula( auto inFunc, double inMin, double inMax,
                  double inTolRel = 1e-6 );

template <typename TypeY_ = double>
TypeY_ doublyAdaptiveNewtonCotes( auto inFunc, double inMin, double inMax,
                                  double inTolAbs = 0.0,
                                  double inTolRel = 1e-6 );

}  // namespace FiniteInterval

namespace InfiniteInterval
{

template <typename TypeY_ = double>
TypeY_ DEFormula( auto inFunc, double inTolRel = 1e-6 );

template <typename TypeY_ = double>
TypeY_ doublyAdaptiveNewtonCotes( auto inFunc, double inTolAbs = 0.0,
                                  double inTolRel = 1e-6 );

}  // namespace InfiniteInterval

namespace UpperInfiniteInterval
{

template <typename TypeY_ = double>
TypeY_ DEFormula( auto inFunc, double inMin, double inTolRel = 1e-6 );

template <typename TypeY_ = double>
TypeY_ DEFormulaForExp( auto inFunc, double inMin, double inTolRel = 1e-6 );

template <typename TypeY_ = double>
TypeY_ DEFormulaForGaussian( auto inFunc, double inMin,
                             double inTolRel = 1e-6 );

template <typename TypeY_ = double>
TypeY_ DEFormulaForVibration( auto inFunc, double inMin, double inTolRel = 1e-6,
                              double inFactorK = 6.0 );

template <typename TypeY_ = double>
TypeY_ doublyAdaptiveNewtonCotes( auto inFunc, double inMin,
                                  double inTolAbs = 0.0,
                                  double inTolRel = 1e-6 );

}  // namespace UpperInfiniteInterval

}  // namespace Math::Integral

#ifndef NINCLUDE_TPP
#include "math/integral_1d_adaptive.tpp"
#include "math/integral_1d_de.tpp"
#endif

#endif