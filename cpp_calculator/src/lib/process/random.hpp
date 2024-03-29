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

#include "market_data.hpp"
#include "math/matrix.hpp"

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
    std::vector<std::vector<double> >
        mRandomValues;         //! generated random values
    std::size_t mIndTmpTime;   //! the index of temporary time
    MarketData::Terms mTerms;  //! term structure

public:
    /**
     * @brief This constructs a new PathAbstract.
     * @param inNPath the number of Path
     * @param insTerms shared pointer of the term structure
     */
    PathAbstract( std::size_t inNPath,
                  std::shared_ptr<const std::vector<double> > insTerms ) :
        mNPath( inNPath ), mTerms( insTerms )
    {
    }
    PathAbstract( std::size_t inNPath, const MarketData::Terms& inTerms ) :
        mNPath( inNPath ), mTerms( inTerms )
    {
    }
    /**
     * @brief This accesses the nth random path.
     * @param inIndPath index of path
     * @return const std::vector<double>& nth random path
     */
    const std::vector<double>& operator[]( std::size_t inIndPath ) const
    {
        return mRandomValues.at( inIndPath );
    }
    /**
     * @brief This accesses the nth random path.
     * @param inIndPath index of path
     * @return const std::vector<double>& nth random path
     */
    const std::vector<double>& at( std::size_t inIndPath ) const
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
    std::normal_distribution<double> mDistribution;

public:
    PathBrownPlain( std::size_t inNPath,
                    std::shared_ptr<const std::vector<double> > insTerms ) :
        PathAbstract( inNPath, insTerms ),
        mGenerator( mDevice() ),
        mDistribution( 0.0, 1.0 )
    {
    }
    PathBrownPlain( std::size_t inNPath, const MarketData::Terms& inTerms ) :
        PathAbstract( inNPath, inTerms ),
        mGenerator( mDevice() ),
        mDistribution( 0.0, 1.0 )
    {
    }
    void setIndexTime( std::size_t inIndex ) override
    {
        PathAbstract::setIndexTime( inIndex );
        mTmpSqrtInterval = std::sqrt( mTerms[inIndex] - mTerms[inIndex - 1] );
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
    std::normal_distribution<double> mDistribution;

public:
    PathBrownAntithetic(
        std::size_t inNPath,
        std::shared_ptr<const std::vector<double> > insTerms ) :
        PathAbstract( inNPath, insTerms ),
        mIsNextNew( true ),
        mGenerator( mDevice() ),
        mDistribution( 0.0, 1.0 )
    {
    }
    PathBrownAntithetic( std::size_t inNPath,
                         const MarketData::Terms& inTerms ) :
        PathAbstract( inNPath, inTerms ),
        mIsNextNew( true ),
        mGenerator( mDevice() ),
        mDistribution( 0.0, 1.0 )
    {
    }
    void setIndexTime( std::size_t inIndex ) override
    {
        PathAbstract::setIndexTime( inIndex );
        mTmpSqrtInterval = std::sqrt( mTerms[inIndex] - mTerms[inIndex - 1] );
        mIsNextNew       = true;
    }
    double generateRandomVal() override;
};

}  // namespace Random

namespace RandomVec
{

/**
 * @brief This is the abstract class for random vector classes.
 */
class PathAbstract
{
protected:
    std::size_t mNPath;  //! the number of Path
    std::size_t mDim;    //! the dimension of Vec
    std::vector<std::vector<Math::Vec> >
        mRandomValues;         //! generated random values
    std::size_t mIndTmpTime;   //! the index of temporary time
    MarketData::Terms mTerms;  //! term structurez

public:
    /**
     * @brief This constructs a new PathAbstract.
     * @param inNPath the number of Path
     * @param insTerms shared pointer of the term structure
     * @param inDim the dimension of generated vectors
     */
    PathAbstract( std::size_t inNPath,
                  std::shared_ptr<const std::vector<double> > insTerms,
                  std::size_t inDim ) :
        mNPath( inNPath ), mTerms( insTerms ), mDim( inDim )
    {
    }
    PathAbstract( std::size_t inNPath, const MarketData::Terms& inTerms,
                  std::size_t inDim ) :
        mNPath( inNPath ), mTerms( inTerms ), mDim( inDim )
    {
    }
    /**
     * @brief This accesses the nth random path.
     * @param inIndPath index of path
     * @return const std::vector<double>& nth random path
     */
    const std::vector<Math::Vec>& operator[]( std::size_t inIndPath ) const
    {
        return mRandomValues.at( inIndPath );
    }
    /**
     * @brief This accesses the nth random path.
     * @param inIndPath index of path
     * @return const std::vector<double>& nth random path
     */
    const std::vector<Math::Vec>& at( std::size_t inIndPath ) const
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
    virtual Math::Vec generateRandomVal() = 0;
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
    std::normal_distribution<double> mDistribution;

public:
    PathBrownPlain( std::size_t inNPath,
                    std::shared_ptr<const std::vector<double> > insTerms,
                    std::size_t inDim ) :
        PathAbstract( inNPath, insTerms, inDim ),
        mGenerator( mDevice() ),
        mDistribution( 0.0, 1.0 )
    {
    }
    PathBrownPlain( std::size_t inNPath, const MarketData::Terms& inTerms,
                    std::size_t inDim ) :
        PathAbstract( inNPath, inTerms, inDim ),
        mGenerator( mDevice() ),
        mDistribution( 0.0, 1.0 )
    {
    }
    void setIndexTime( std::size_t inIndex ) override
    {
        PathAbstract::setIndexTime( inIndex );
        mTmpSqrtInterval = std::sqrt( mTerms[inIndex] - mTerms[inIndex - 1] );
    }
    Math::Vec generateRandomVal() override;
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
    Math::Vec mPrevRandomValue;
    std::random_device mDevice;
    std::mt19937 mGenerator;
    std::normal_distribution<double> mDistribution;

public:
    PathBrownAntithetic( std::size_t inNPath,
                         std::shared_ptr<const std::vector<double> > insTerms,
                         std::size_t inDim ) :
        PathAbstract( inNPath, insTerms, inDim ),
        mIsNextNew( true ),
        mGenerator( mDevice() ),
        mDistribution( 0.0, 1.0 ),
        mPrevRandomValue( inDim )
    {
    }
    PathBrownAntithetic( std::size_t inNPath, const MarketData::Terms& inTerms,
                         std::size_t inDim ) :
        PathAbstract( inNPath, inTerms, inDim ),
        mIsNextNew( true ),
        mGenerator( mDevice() ),
        mDistribution( 0.0, 1.0 ),
        mPrevRandomValue( inDim )
    {
    }
    void setIndexTime( std::size_t inIndex ) override
    {
        PathAbstract::setIndexTime( inIndex );
        mTmpSqrtInterval = std::sqrt( mTerms[inIndex] - mTerms[inIndex - 1] );
        mIsNextNew       = true;
    }
    Math::Vec generateRandomVal() override;
};

}  // namespace RandomVec
}  // namespace Process

#endif