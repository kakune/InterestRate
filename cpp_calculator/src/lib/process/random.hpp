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

namespace Process
{
namespace Random
{

/**
 * @brief This is the abstract class for random path classes.
 */
class PathAbstract
{
protected:
    std::size_t mNPath;  //! the number of Path
    std::shared_ptr< const std::vector< double > > msTerms;  //! term structure
    std::vector< std::vector< double > >
        mRandomValues;        //! generated random values
    std::size_t mIndTmpTime;  //! the index of temporary time
public:
    /**
     * @brief This constructs a new PathAbstract.
     * @param inNPath the number of Path
     * @param insTerms shared pointer of the term structure
     */
    PathAbstract( std::size_t inNPath,
                  std::shared_ptr< const std::vector< double > > insTerms ) :
        mNPath( inNPath ), msTerms( insTerms )
    {
    }
    /**
     * @brief This accesses the nth random path.
     * @param inIndPath index of path
     * @return const std::vector<double>& nth random path
     */
    const std::vector< double >& operator[]( std::size_t inIndPath ) const
    {
        return mRandomValues.at( inIndPath );
    }
    /**
     * @brief This accesses the nth random path.
     * @param inIndPath index of path
     * @return const std::vector<double>& nth random path
     */
    const std::vector< double >& at( std::size_t inIndPath ) const
    {
        return mRandomValues.at( inIndPath );
    }
    /**
     * @brief This initialize mRandomValues.
     */
    void initRandomValues();
    /**
     * @brief This sets the Time Index before generating random variable.
     * @param inIndex index of time
     */
    virtual void setIndexTime( std::size_t inIndex ) { mIndTmpTime = inIndex; }
    /**
     * @brief This generates the random value.
     * @return double random value
     */
    virtual double generateRandomVal() = 0;
    /**
     * @brief This makes random path.
     */
    virtual void makePath();
    /**
     * @brief This makes random variables, NOT path.
     */
    virtual void makeRandomVariables();
    virtual ~PathAbstract() = default;
};

/**
 * @brief This is path of Brownian motion without using variance reduction.
 */
class PathBrownPlain : public PathAbstract
{
private:
    double mTmpSqrtInterval;
    std::random_device mDevice;
    std::mt19937 mGenerator;
    std::normal_distribution< double > mDistribution;

public:
    PathBrownPlain( std::size_t inNPath,
                    std::shared_ptr< const std::vector< double > > insTerms ) :
        PathAbstract( inNPath, insTerms ),
        mGenerator( mDevice() ),
        mDistribution( 0.0, 1.0 )
    {
    }
    void setIndexTime( std::size_t inIndex ) override
    {
        PathAbstract::setIndexTime( inIndex );
        mTmpSqrtInterval =
            std::sqrt( msTerms->at( inIndex ) - msTerms->at( inIndex - 1 ) );
    }
    double generateRandomVal() override;
};

/**
 * @brief This is path of Brownian motion with antithetic variates method.
 * path[i+1][t] = -path[i][t] (i:even)
 */
class PathBrownAntithetic : public PathAbstract
{
private:
    double mTmpSqrtInterval;
    bool mIsNextNew = true;
    double mPrevRandomValue;
    std::random_device mDevice;
    std::mt19937 mGenerator;
    std::normal_distribution< double > mDistribution;

public:
    PathBrownAntithetic(
        std::size_t inNPath,
        std::shared_ptr< const std::vector< double > > insTerms ) :
        PathAbstract( inNPath, insTerms ),
        mIsNextNew( true ),
        mGenerator( mDevice() ),
        mDistribution( 0.0, 1.0 )
    {
    }
    void setIndexTime( std::size_t inIndex ) override
    {
        PathAbstract::setIndexTime( inIndex );
        mTmpSqrtInterval =
            std::sqrt( msTerms->at( inIndex ) - msTerms->at( inIndex - 1 ) );
        mIsNextNew = true;
    }
    double generateRandomVal() override;
};

}  // namespace Random
}  // namespace Process

#endif