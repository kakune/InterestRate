/**
 * @file random.hpp
 * @brief This defines random path classes.
 * @author kakune
 * @date 1/29/2024
 */

#ifndef PROCESS_RANDOM_HPP
#define PROCESS_RANDOM_HPP

#include <memory>
#include <vector>

namespace Process
{
namespace Random
{

/**
 * @brief This is the interface of random paths.
 */
class PathAbstract
{
protected:
    std::size_t mNPath;  //! the number of Path
    std::shared_ptr< const std::vector< double > > msTerms;  //! term structure
    std::vector< std::vector< double > >
        mRandomValues;  //! generated random values
public:
    /**
     * @brief This constructs a new PathAbstract.
     * @param inNPath the number of Path
     * @param insTerms shared pointer of the term structure
     */
    PathAbstract( std::size_t inNPath,
                  std::shared_ptr< const std::vector< double > > insTerms ) :
        mNPath( inNPath ),
        msTerms( insTerms ),
        mRandomValues( std::vector< std::vector< double > >(
            inNPath, std::vector< double >( insTerms->size(), 0 ) ) )
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
     * @brief This calculate the nth path value at arbital time.
     * @param inIndPath index of path
     * @param inTime time that must be within Terms.
     * @return double the nth path value at time
     */
    double operator()( std::size_t inIndPath, double inTime );
    /**
     * @brief This makes random path.
     */
    virtual void makePath() = 0;
};

/**
 * @brief This is path of Brownian motion without using variance reduction.
 */
class PathBrownPlain : public PathAbstract
{
public:
    PathBrownPlain( std::size_t inNPath,
                    std::shared_ptr< const std::vector< double > > insTerms ) :
        PathAbstract( inNPath, insTerms )
    {
    }
    void makePath() override;
};

/**
 * @brief This is path of Brownian motion with antithetic variates method.
 * path[i+1][t] = -path[i][t] (i:even)
 */
class PathBrownAntithetic : public PathAbstract
{
public:
    PathBrownAntithetic(
        std::size_t inNPath,
        std::shared_ptr< const std::vector< double > > insTerms ) :
        PathAbstract( inNPath, insTerms )
    {
    }
    void makePath() override;
};

}  // namespace Random
}  // namespace Process

#endif