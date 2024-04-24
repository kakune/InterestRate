/**
 * @file core.hpp
 * @brief This defines classes for short model calculation with PDE.
 * @author kakune
 * @date 3/21/2024
 */

#ifndef SHORT_RATE_PDE_CORE_HPP
#define SHORT_RATE_PDE_CORE_HPP

#include <algorithm>

namespace ShortRate
{
namespace PDE
{

/**
 * @brief This is the abstract class for short rate models.
 */
class ModelAbstract
{
protected:
    const double mStartTime, mEndTime;  //! the time region for PDE
    const double mStepTime;             //! the step of time in PDE

public:
    /**
     * @brief This constructs a new ModelAbstract.
     * @param inStartTime the start-time of PDE
     * @param inEndTime  the end-time of PDE
     * @param inStepTime the step of time in PDE
     */
    ModelAbstract( double inStartTime, double inEndTime, double inStepTime ) :
        mStartTime( std::min( { inStartTime, inEndTime } ) ),
        mEndTime( std::max( { inStartTime, inEndTime } ) ),
        mStepTime( inStepTime > 0.0
                       ? inStepTime
                       : std::abs( inEndTime - inStartTime ) / 10.0 )
    {
    }
    virtual ~ModelAbstract() = default;
};

/**
 * @brief This build the object of ModelAbstract.
 */
class ModelAbstractBuilder
{
protected:
    double mStartTime, mEndTime;
    double mStepTime = 0.0;

public:
    ModelAbstractBuilder& setRegTime( double inStartTime, double inEndTime )
    {
        mStartTime = inStartTime;
        mEndTime   = inEndTime;
        return *this;
    }
    ModelAbstractBuilder& setStepTime( double inStepTime )
    {
        mStepTime = inStepTime;
        return *this;
    }
    virtual ~ModelAbstractBuilder() = default;
};

}  // namespace PDE
}  // namespace ShortRate

#endif