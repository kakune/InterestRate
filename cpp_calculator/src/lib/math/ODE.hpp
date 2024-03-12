/**
 * @file ODE.hpp
 * @brief This defines 1-dim ode solver.
 * @author kakune
 * @date 3/11/2024
 */

#ifndef MATH_ODE_HPP
#define MATH_ODE_HPP

namespace Math
{
namespace ODE
{

/**
 * @brief Solves an ordinary differential equation (ODE).
 * @tparam Func_ y'(x,y) double(double, double) A callable object representing
 * the ODE.
 * @param inInitX The initial value of the independent variable.
 * @param inInitY The initial value of the dependent variable.
 * @param inEndX The end value of the independent variable.
 * @param inRelTol The relative tolerance for the solver (default is 1e-6).
 * @param inMaxDif The maximum allowed difference between successive
 * approximations (default is 1.0).
 * @return A pair of vectors, [independent, dependent].
 */
template <auto Func_>
std::pair<std::vector<double>, std::vector<double>> solveRungeKutta45(
    double inInitX, double inInitY, double inEndX, double inRelTol = 1e-6,
    double inMaxDif = 1.0 );

/**
 * @brief Solves a system of ordinary differential equations (ODEs).
 * @tparam Func_ y'(x,y) vector<double>(double, vector<double>) A callable
 * object representing the system of ODEs.
 * @param inInitX The initial value of the independent variable.
 * @param inInitY A vector containing the initial values of the dependent
 * variables.
 * @param inEndX The end value of the independent variable.
 * @param inRelTol The relative tolerance for the solver (default is 1e-6).
 * @param inMaxDif The maximum allowed difference between successive
 * approximations (default is 1.0).
 * @return A vector of vectors, [independent, dependent1, dependent2, ...].
 */
template <auto Func_>
std::vector<std::vector<double>> solveSIMLRungeKutta45(
    double inInitX, std::vector<double> inInitY, double inEndX,
    double inRelTol = 1e-6, double inMaxDif = 1.0 );

/**
 * @brief Solves a second-order ordinary differential equation (ODE).
 * @tparam Func_ y''(x,y,y') double(double, double, double) A callable object
 * representing the ODE.
 * @param inInitX The initial value of the independent variable.
 * @param inInitY The initial value of the dependent variable.
 * @param inInitDY The initial value of the derivative of the dependent
 * variable.
 * @param inEndX The end value of the independent variable.
 * @param inRelTol The relative tolerance for the solver (default is 1e-6).
 * @param inMaxDif The maximum allowed difference between successive
 * approximations (default is 1.0).
 * @return A tuple containing vectors of the independent variable, the first
 * dependent variable, and the second dependent variable at each step.
 */
template <auto Func_>
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
solveSecondOrderRungeKutta45( double inInitX, double inInitY, double inInitDY,
                              double inEndX, double inRelTol,
                              double inMaxDif = 1.0 );

/**
 * @brief Solves a second-order system of ordinary differential equations
 * (ODEs).
 * @tparam Func_ y''(x,y,y') vector<double>(double, vector<double>,
 * vector<double>) A callable object representing the system of ODEs.
 * @param inInitX The initial value of the independent variable.
 * @param inInitY A vector containing the initial values of the dependent
 * variables.
 * @param inInitDY A vector containing the initial derivatives of the dependent
 * variable.s
 * @param inEndX The end value of the independent variable.
 * @param inRelTol The relative tolerance for the solver (default is 1e-6).
 * @param inMaxDif The maximum allowed difference between successive
 * approximations (default is 1.0).
 * @return A vector of vectors, where each inner vector represents the values of
 * the dependent variables at a step.
 */
template <auto Func_>
std::vector<std::vector<double>> solveSIMLSecondOrderRungeKutta45(
    double inInitX, std::vector<double> inInitY, std::vector<double> inInitDY,
    double inEndX, double inRelTol, double inMaxDif = 1.0 );

}  // namespace ODE
}  // namespace Math

#include "math/ODE.tpp"

#endif