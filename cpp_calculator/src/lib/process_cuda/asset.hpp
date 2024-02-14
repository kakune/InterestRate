/**
 * @file asset.hpp
 * @brief This defines calculator using CUDA of assets.
 * @author kakune
 * @date 2/14/2024
 */

#ifndef PROCESS_CUDA_ASSET_HPP
#define PROCESS_CUDA_ASSET_HPP

#include <cuda.h>
#include <cuda_runtime.h>

#include <cmath>
#include <memory>
#include <vector>

namespace ProcessCUDA
{
namespace Asset
{

/**
 * @brief This is the abstract class for SDE Model of forward price.
 */
class ModelForwardAbstract
{
protected:
    std::size_t mNPath;       //! the number of Path
    std::size_t mNTerm;       //! the number of Term
    double* mpcSqDts;         //! square roots of time interval
    double* mpcForwardPrice;  //! calcurated forward price
    double mInitPrice;        //! initial forward price
    double mTimeMaturity;     //! maturity time
    bool mIsMadeDt;  //! flag whether mpcSqDts was created by this class.

    /**
     * @brief This calcurate interest rates of each random path.
     */
    virtual void calcEachForwardPrice() = 0;

public:
    /**
     * @brief This constructs a ModelAbstract.
     * @param inNPath the number of Path
     * @param insTerms term structure
     * @param inInitPrice the initial price of forward price
     */
    ModelForwardAbstract( std::size_t inNPath,
                          std::shared_ptr<const std::vector<double> > insTerms,
                          double inInitPrice ) :
        mNPath( inNPath ),
        mNTerm( insTerms->size() ),
        mpcSqDts( nullptr ),
        mTimeMaturity( insTerms->back() - insTerms->front() ),
        mpcForwardPrice( nullptr ),
        mInitPrice( inInitPrice ),
        mIsMadeDt( true )
    {
        std::vector<double> lSqDts( insTerms->size(), 0.0 );
        for ( std::size_t i = 1; i < insTerms->size(); ++i )
        {
            lSqDts.at( i ) =
                std::sqrt( insTerms->at( i ) - insTerms->at( i - 1 ) );
        }
        cudaMalloc( (void**)&mpcSqDts, mNTerm * sizeof( double ) );
        cudaMemcpy( mpcSqDts, lSqDts.data(), mNTerm * sizeof( double ),
                    cudaMemcpyHostToDevice );
    }
    /**
     * @brief This constructs a ModelAbstract.
     * @param inNPath the number of Path
     * @param inNTerm the number of Term
     * @param inpcSqDts square roots of time interval
     * @param inTimeMaturity maturity time
     * @param inInitPrice the initial price of forward price
     */
    ModelForwardAbstract( std::size_t inNPath, std::size_t inNTerm,
                          double* inpcSqDts, double inTimeMaturity,
                          double inInitPrice ) :
        mNPath( inNPath ),
        mNTerm( inNTerm ),
        mpcSqDts( inpcSqDts ),
        mTimeMaturity( inTimeMaturity ),
        mpcForwardPrice( nullptr ),
        mInitPrice( inInitPrice ),
        mIsMadeDt( false )
    {
    }
    ~ModelForwardAbstract()
    {
        if ( mpcForwardPrice != nullptr ) { cudaFree( mpcForwardPrice ); }
        if ( mIsMadeDt && mpcSqDts != nullptr ) { cudaFree( mpcSqDts ); }
    };
    /**
     * @brief This computes the price of call option WITHOUT discount factor.
     * @param inStrike the price of strike
     * @return double price of call option
     */
    double priceCallOption( double inStrike );
    /**
     * @brief This computes the price of put option WITHOUT discount factor.
     * @param inStrike the price of strike
     * @return double price of call option
     */
    double pricePutOption( double inStrike );
    /**
     * @brief This computes the implied volatility.
     * @param inStrike the price of strike
     * @return double implied volatility
     */
    double impliedVolatility( double inStrike );
};

/**
 * @brief This build the object of ModelAbstract
 */
class ModelForwardAbstractBuilder
{
protected:
    std::size_t mNPath;          //! the number of Path
    std::size_t mNTerm;          //! the number of Term
    double* mpcSqDts = nullptr;  //! square roots of time interval
    double mTimeMaturity;        //! maturity time
    std::shared_ptr<const std::vector<double> > msTerms;  //! term structure
    double mInitPrice;  //! initial forward price

public:
    ModelForwardAbstractBuilder& setNPath( std::size_t inNPath )
    {
        mNPath = inNPath;
        return *this;
    }
    ModelForwardAbstractBuilder& setSqDts( std::size_t inNTerm,
                                           double* inpcSqDts,
                                           double inTimeMaturity )
    {
        mNTerm        = inNTerm;
        mpcSqDts      = inpcSqDts;
        mTimeMaturity = inTimeMaturity;
        return *this;
    }
    ModelForwardAbstractBuilder& setTerms(
        std::shared_ptr<const std::vector<double> > insTerms )
    {
        msTerms = insTerms;
        return *this;
    }
    ModelForwardAbstractBuilder& setInitPrice( double inInitPrice )
    {
        mInitPrice = inInitPrice;
        return *this;
    }
    virtual ~ModelForwardAbstractBuilder() = default;
};

/**
 * @brief This is SDE calculator for Black-Sholes model.
 */
class BlackSholesForward : public ModelForwardAbstract
{
private:
    double mVol;  //! constant volatility
    void calcEachForwardPrice() override;

public:
    BlackSholesForward( std::size_t inNPath,
                        std::shared_ptr<const std::vector<double> > insTerms,
                        double inInitPrice, double inVol ) :
        ModelForwardAbstract( inNPath, insTerms, inInitPrice ), mVol( inVol )
    {
        calcEachForwardPrice();
    }
    BlackSholesForward( std::size_t inNPath, std::size_t inNTerm,
                        double* inpcSqDts, double inTimeMaturity,
                        double inInitPrice, double inVol ) :
        ModelForwardAbstract( inNPath, inNTerm, inpcSqDts, inTimeMaturity,
                              inInitPrice ),
        mVol( inVol )
    {
        calcEachForwardPrice();
    }
};

/**
 * @brief This build the object of BlackSholesForward
 */
class BlackSholesForwardBuilder : public ModelForwardAbstractBuilder
{
private:
    double mVol;  //! constant volatility

public:
    BlackSholesForwardBuilder& setVol( double inVol )
    {
        mVol = inVol;
        return *this;
    }
    BlackSholesForward build()
    {
        if ( mpcSqDts == nullptr )
        {
            return BlackSholesForward( mNPath, msTerms, mInitPrice, mVol );
        }
        return BlackSholesForward( mNPath, mNTerm, mpcSqDts, mTimeMaturity,
                                   mInitPrice, mVol );
    }
};

/**
 * @brief This is calculator for SABR.
 */
class SABRForward : public ModelForwardAbstract
{
private:
    double mExponent;  //! exponent of forward price in SDE
    double mVolvol;    //! volvol
    double mInitVol;   //! initial volatility
    double mCorr;      //! correlation between brownians
    void calcEachForwardPrice() override;

public:
    SABRForward( std::size_t inNPath,
                 std::shared_ptr<const std::vector<double> > insTerms,
                 double inInitPrice, double inInitVol, double inCorr,
                 double inExponent, double inVolvol ) :
        ModelForwardAbstract( inNPath, insTerms, inInitPrice ),
        mInitVol( inInitVol ),
        mCorr( inCorr ),
        mExponent( inExponent ),
        mVolvol( inVolvol )
    {
        calcEachForwardPrice();
    }
    SABRForward( std::size_t inNPath, std::size_t inNTerm, double* inpcSqDts,
                 double inTimeMaturity, double inInitPrice, double inInitVol,
                 double inCorr, double inExponent, double inVolvol ) :
        ModelForwardAbstract( inNPath, inNTerm, inpcSqDts, inTimeMaturity,
                              inInitPrice ),
        mInitVol( inInitVol ),
        mCorr( inCorr ),
        mExponent( inExponent ),
        mVolvol( inVolvol )
    {
        calcEachForwardPrice();
    }
};

/**
 * @brief This build the object of SABRForward
 */
class SABRForwardBuilder : public ModelForwardAbstractBuilder
{
private:
    double mInitVol;   //! initial volatility
    double mCorr;      //! correlation between brownians
    double mExponent;  //! exponent of forward price in SDE
    double mVolvol;    //! volvol

public:
    SABRForwardBuilder& setInitVol( double inInitVol )
    {
        mInitVol = inInitVol;
        return *this;
    }
    SABRForwardBuilder& setCorr( double inCorr )
    {
        mCorr = inCorr;
        return *this;
    }
    SABRForwardBuilder& setExponent( double inExponent )
    {
        mExponent = inExponent;
        return *this;
    }
    SABRForwardBuilder& setVolvol( double inVolvol )
    {
        mVolvol = inVolvol;
        return *this;
    }
    SABRForward build()
    {
        if ( mpcSqDts == nullptr )
        {
            return SABRForward( mNPath, msTerms, mInitPrice, mInitVol, mCorr,
                                mExponent, mVolvol );
        }
        return SABRForward( mNPath, mNTerm, mpcSqDts, mTimeMaturity, mInitPrice,
                            mInitVol, mCorr, mExponent, mVolvol );
    }
};

__global__ void priceCallOptionKernel( double* inResult, std::size_t inNPath,
                                       double inStrike, double* inPrice );
__global__ void pricePutOptionKernel( double* inResult, std::size_t inNPath,
                                      double inStrike, double* inPrice );
__global__ void calcBlackEachForwardKernel( double* inResult,
                                            std::size_t inNPath,
                                            std::size_t inNTerm, double inVol,
                                            double inInitPrice, double* inSqDts,
                                            unsigned long long inSeed = 5556 );
__global__ void calcBlackEachForwardWithLogKernel(
    double* inResult, std::size_t inNPath, std::size_t inNTerm, double inVol,
    double inInitPrice, double* inSqDts, unsigned long long inSeed = 5556 );

__global__ void calcSABREachForwardKernel(
    double* inResult, std::size_t inNPath, std::size_t inNTerm,
    double inInitPrice, double inInitVol, double inCorr, double inExponent,
    double inVolvol, double* inSqDts, unsigned long long inSeed = 5556 );
__global__ void calcSABREachForwardWithLogKernel(
    double* inResult, std::size_t inNPath, std::size_t inNTerm,
    double inInitPrice, double inInitVol, double inCorr, double inExponent,
    double inVolvol, double* inSqDts, unsigned long long inSeed = 5556 );

}  // namespace Asset
}  // namespace ProcessCUDA

#endif