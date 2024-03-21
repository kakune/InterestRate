/**
 * @file interpolate_1d.hpp
 * @brief This defines interpolation function for 1d.
 * @author kakune
 * @date 1/29/2024
 */

#ifndef MATH_INTERPOLATE_1D_HPP
#define MATH_INTERPOLATE_1D_HPP

#include <memory>
#include <vector>

namespace Math
{
namespace Interpolate1D
{

/**
 * @brief This is the abstract class for 1d-interpolation classes.
 */
class GeneratorAbstract
{
protected:
    /**
     * @brief This calculates coefficients of interpolation function.
     * @param insRefXs vector of X ( must be > mNDeg )
     * @param insRefYs vector of Y corresponding to f(msRefXs)
     */
    virtual void build(
        std::shared_ptr<const std::vector<double> > insRefXs,
        std::shared_ptr<const std::vector<double> > insRefYs ) = 0;

public:
    /**
     * @brief This calculates the value f(inX) using interpolation.
     * @param inX
     * @return double f(inX)
     */
    virtual double operator()( double inX ) const = 0;
    virtual ~GeneratorAbstract()                  = default;
};

/**
 * @brief This class computes the Newton spline of arbitary degree.
 */
class NewtonSpline : public GeneratorAbstract
{
private:
    std::size_t mNDeg;
    std::size_t mSizeIndex;
    std::shared_ptr<const std::vector<double> > msRefXs;  //! vector of X
    std::vector<std::vector<double> > mCoeff;
    std::vector<std::vector<double> > mCoeffIntegral;
    std::vector<double> mSumIntegral;
    void build( std::shared_ptr<const std::vector<double> > insRefXs,
                std::shared_ptr<const std::vector<double> > insRefYs ) override;
    /**
     * @brief This prepares the coefficients for integral.
     */
    void buildIntegral();

public:
    /**
     * @brief This constructs a new NewtonSpline.
     * @param insRefXs vector of X
     * @param insRefYs vector of Y corresponding to f(msRefXs)
     * @param inNDeg degree of spline
     */
    NewtonSpline( std::shared_ptr<const std::vector<double> > insRefXs,
                  std::shared_ptr<const std::vector<double> > insRefYs,
                  std::size_t inNDeg = 1 ) :
        mNDeg( inNDeg ),
        mCoeff( std::vector<std::vector<double> >( inNDeg + 1 ) ),
        msRefXs( insRefXs )
    {
        build( insRefXs, insRefYs );
        buildIntegral();
    }

    /**
     * @brief This calculates the value f(inX) using mNDeg-degree Newton spline.
     * @param inX
     * @return double f(inX)
     */
    double operator()( double inX ) const override;
    /**
     * @brief This calculates the derivative of f(inX) using mNDeg-degree Newton
     * spline.
     * @param inX
     * @param inOrder the order of derivative
     * @return double f^(inOrder)(inX)
     */
    double deriv( double inX, std::size_t inOrder = 1 ) const;
    /**
     * @brief This calculates the integral in the area [x0, inX].
     * @param inX
     * @return double
     */
    double integral( double inX ) const;
    /**
     * @brief This calculates the integral in the area [inA, inB].
     * @param inX
     * @return double
     */
    double integral( double inA, double inB ) const;
};

}  // namespace Interpolate1D
}  // namespace Math

#endif