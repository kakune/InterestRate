/**
 * @file interpolate_1d.hpp
 * @brief This defines interpolation function for 1d.
 * @author kakune
 * @date 1/29/2024
 */

#include <memory>
#include <vector>

namespace Math
{
namespace Interpolate1d
{

/**
 * @brief This is the abstract class for 1d-interpolation classes.
 */
class GeneratorAbstract
{
protected:
    bool mIsBuilt = false;
    std::shared_ptr< const std::vector< double > > msRefXs;  //! vector of X
    std::shared_ptr< const std::vector< double > >
        msRefYs;  //! vector of Y corresponding to f(msRefXs)
public:
    /**
     * @brief This calculates coefficients of interpolation function.
     * @param insRefXs vector of X
     * @param insRefYs vector of Y corresponding to f(msRefXs)
     */
    virtual void build(
        std::shared_ptr< const std::vector< double > > insRefXs,
        std::shared_ptr< const std::vector< double > > insRefYs ) = 0;
    /**
     * @brief This calculates the value f(inX) using interpolation.
     * @param inX
     * @return double f(inX)
     */
    virtual double operator()( double inX ) = 0;
    virtual ~GeneratorAbstract()            = default;
};

/**
 * @brief This class computes the Newton spline of arbitary degree.
 */
class NewtonSpline : public GeneratorAbstract
{
private:
    std::size_t mNDeg;
    std::vector< std::vector< double > > mCoeff;

public:
    /**
     * @brief This constructs a new NewtonSpline.
     * @param inNDeg degree of spline
     */
    NewtonSpline( std::size_t inNDeg = 1 ) :
        mNDeg( inNDeg ),
        mCoeff( std::vector< std::vector< double > >( inNDeg ) )
    {
    }
    /**
     * @brief This calculates coefficients of interpolation function.
     * @param insRefXs vector of X ( must be > mNDeg )
     * @param insRefYs vector of Y corresponding to f(msRefXs)
     */
    void build(
        std::shared_ptr< const std::vector< double > > insRefXs,
        std::shared_ptr< const std::vector< double > > insRefYs ) override;
    /**
     * @brief This calculates the value f(inX) using mNDeg-degree Newton spline.
     * @param inX
     * @return double f(inX)
     */
    double operator()( double inX ) override;
};

}  // namespace Interpolate1d
}  // namespace Math