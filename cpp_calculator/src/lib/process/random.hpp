/**
 * @file random.hpp
 * @brief This defines random path classes.
 * @author kakune
 * @date 1/29/2024
 */

#ifndef PROCESS_RANDOM_HPP
#define PROCESS_RANDOM_HPP

#include <memory>
#include <random>
#include <vector>

#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <curand_kernel.h>
#endif

#include "market_data.hpp"
#include "math/matrix.hpp"

namespace Process
{
namespace Random
{

template <class T_>
concept C_StdBrown = requires( T_ inObj ) {
    inObj.initialize();
    {
        inObj()
    } -> std::convertible_to<double>;
};

/**
 * @brief This is the abstract class for standard brownian motion generator
 * classes.
 */
class StdBrownAbstract
{
public:
    /**
     * @brief This initializes the variance reduction method.
     */
    virtual void initialize() {}
    /**
     * @brief This generates the random value.
     * @return double random value
     */
    virtual double operator()() = 0;
    virtual ~StdBrownAbstract() = default;
};

/**
 * @brief This is standard brownian motion generator without using variance
 * reduction.
 */
class StdBrownPlain : public StdBrownAbstract
{
private:
    std::random_device mDevice;
    std::mt19937 mGenerator;
    std::normal_distribution<double> mDistribution;

public:
    StdBrownPlain() : mGenerator( mDevice() ), mDistribution( 0.0, 1.0 ) {}
    double operator()() override;
};

/**
 * @brief This is standard brownian motion generator with antithetic variates
 * method. path[i+1][t] = -path[i][t] (i:even)
 */
class StdBrownAntithetic : public StdBrownAbstract
{
private:
    bool mIsNextNew = true;
    double mPrevRandomValue;
    std::random_device mDevice;
    std::mt19937 mGenerator;
    std::normal_distribution<double> mDistribution;

public:
    StdBrownAntithetic() :
        mIsNextNew( true ), mGenerator( mDevice() ), mDistribution( 0.0, 1.0 )
    {
    }
    void initialize() override;
    double operator()() override;
};

}  // namespace Random

namespace RandomVec
{

template <class T_>
concept C_StdBrown = requires( T_ inObj ) {
    inObj.initialize();
    {
        inObj()
    } -> std::convertible_to<Math::Vec>;
};

/**
 * @brief This is the abstract class for standard brownian motion generator.
 */
class StdBrownAbstract
{
protected:
    std::size_t mDim;       //! the dimension of Vec
    std::size_t mIndStart;  //! components of random vec whose (i < mStartInd)
                            //! become zero.

public:
    /**
     * @brief This constructs a new StdBrownAbstract.
     * @param inDim the dimension of generated vectors
     */
    StdBrownAbstract( std::size_t inDim ) : mDim( inDim ), mIndStart( 0 ) {}
    /**
     * @brief This initializes the variance reduction method.
     */
    virtual void initialize( std::size_t inIndStart = 0 )
    {
        mIndStart = inIndStart;
    }
    /**
     * @brief This generates the random value.
     * @return double random value
     */
    virtual Math::Vec operator()() = 0;
    virtual ~StdBrownAbstract()    = default;
};

/**
 * @brief This is a standard brownian motion generator without using variance
 * reduction.
 */
class StdBrownPlain : public StdBrownAbstract
{
private:
    std::random_device mDevice;
    std::mt19937 mGenerator;
    std::normal_distribution<double> mDistribution;

public:
    StdBrownPlain( std::size_t inDim ) :
        StdBrownAbstract( inDim ),
        mGenerator( mDevice() ),
        mDistribution( 0.0, 1.0 )
    {
    }
    Math::Vec operator()() override;
};

/**
 * @brief This is a standard brownian motion generator with antithetic variates
 * method. path[i+1][t] = -path[i][t] (i:even)
 */
class StdBrownAntithetic : public StdBrownAbstract
{
private:
    bool mIsNextNew;
    Math::Vec mPrevRandomValue;
    std::random_device mDevice;
    std::mt19937 mGenerator;
    std::normal_distribution<double> mDistribution;

public:
    StdBrownAntithetic( std::size_t inDim ) :
        StdBrownAbstract( inDim ),
        mIsNextNew( true ),
        mGenerator( mDevice() ),
        mDistribution( 0.0, 1.0 ),
        mPrevRandomValue( inDim )
    {
    }
    void initialize( std::size_t inIndStart = 0 ) override;
    Math::Vec operator()() override;
};

}  // namespace RandomVec
}  // namespace Process

#endif