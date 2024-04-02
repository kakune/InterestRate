/**
 * @file random.hpp
 * @brief This defines random path classes on GPU.
 * @author kakune
 * @date 2/13/2024
 */

#ifndef PROCESS_CUDA_RANDOM_HPP
#define PROCESS_CUDA_RANDOM_HPP

#include <cuda.h>
#include <cuda_runtime_api.h>

#include <memory>
#include <random>
#include <vector>

namespace ProcessCUDA
{
namespace Random
{

/**
 * @brief This is the abstract class for random path classes.
 */
class StdBrownAbstract
{
protected:
    std::size_t mNPath, mNTerm;  //! the number of Path
    double* mpcTerms;            //! term structure
    double* mpcRandomValues;     //! generated random values

public:
    /**
     * @brief This constructs a new StdBrownAbstract.
     * @param inNPath the number of Path
     * @param insTerms shared pointer of the term structure
     */
    StdBrownAbstract( std::size_t inNPath,
                      std::shared_ptr<const std::vector<double> > insTerms ) :
        mNPath( inNPath ),
        mNTerm( insTerms->size() ),
        mpcTerms( nullptr ),
        mpcRandomValues( nullptr )
    {
        cudaMalloc( (void**)&mpcTerms, insTerms->size() * sizeof( double ) );
        cudaMemcpy( mpcTerms, insTerms->data(),
                    insTerms->size() * sizeof( double ),
                    cudaMemcpyHostToDevice );
    }

    StdBrownAbstract( std::size_t inNPath, std::size_t inNTerm,
                      double* inpcTerms ) :
        mNPath( inNPath ),
        mNTerm( inNTerm ),
        mpcTerms( inpcTerms ),
        mpcRandomValues( nullptr )
    {
    }
    /**
     * @brief This makes random path on GPU.
     */
    virtual void makeRandomVals( unsigned long long seed = 5556 ) = 0;
    double* getPtrTerms() { return mpcTerms; }
    double* getPtrRandomValues()
    {
        if ( mpcRandomValues == nullptr ) { makeRandomVals(); }
        return mpcRandomValues;
    }
    virtual ~StdBrownAbstract()
    {
        if ( mpcRandomValues != nullptr ) { cudaFree( mpcRandomValues ); }
    };
};

/**
 * @brief This is path of Brownian motion without using variance reduction.
 */
class PathBrownPlain : public StdBrownAbstract
{
public:
    using StdBrownAbstract::StdBrownAbstract;
    void makeRandomVals( unsigned long long seed = 5556 ) override;
};

class PathBrownAntithetic : public StdBrownAbstract
{
public:
    PathBrownAntithetic( std::size_t inNPath, std::size_t inNTerm,
                         double* inpcTerms ) :
        StdBrownAbstract( ( inNPath + 1 ) / 2 * 2, inNTerm, inpcTerms )
    {
    }
    PathBrownAntithetic(
        std::size_t inNPath,
        std::shared_ptr<const std::vector<double> > insTerms ) :
        StdBrownAbstract( ( inNPath + 1 ) / 2 * 2, insTerms )
    {
    }
    void makeRandomVals( unsigned long long seed = 5556 ) override;
};

// /**
//  * @brief This is path of Brownian motion with antithetic variates method.
//  * path[i+1][t] = -path[i][t] (i:even)
//  */
// class PathBrownAntithetic : public StdBrownAbstract
// {
// private:
//     double mTmpSqrtInterval;
//     bool mIsNextNew = true;
//     double mPrevRandomValue;
//     std::random_device mDevice;
//     std::mt19937 mGenerator;
//     std::normal_distribution< double > mDistribution;

// public:
//     PathBrownAntithetic(
//         std::size_t inNPath,
//         std::shared_ptr< const std::vector< double > > insTerms ) :
//         StdBrownAbstract( inNPath, insTerms ),
//         mIsNextNew( true ),
//         mGenerator( mDevice() ),
//         mDistribution( 0.0, 1.0 )
//     {
//     }
//     void setIndexTime( std::size_t inIndex ) override
//     {
//         StdBrownAbstract::setIndexTime( inIndex );
//         mTmpSqrtInterval =
//             std::sqrt( msTerms->at( inIndex ) - msTerms->at( inIndex - 1 ) );
//         mIsNextNew = true;
//     }
//     double generateRandomVal() override;
// };

}  // namespace Random
}  // namespace ProcessCUDA

#endif