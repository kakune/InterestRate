/**
 * @file random.cu
 * @brief This implements random generators with CUDA.
 * @author kakune
 * @date 4/2/2024
 */

#include <stdio.h>

#include "short_rate/one-factor/core.hpp"

namespace ShortRate
{
namespace OneFactor
{
typedef void ( *createShortRates )( Process::ModelData::SpotRatesForDevice, void* );

struct DeviceSpotRateCreator
{
    void* mParam;
    Process::ModelData::SpotRatesForDevice mSpotRates;
    createShortRates mCreateFunc;
};

__host__ __device__ DeviceStdBrownGenerator::DeviceStdBrownGenerator(
    generateRandom inFunc, initializeRandom inInitFunc, void* inParam ) :
    mGenerateFunc( inFunc ), mInitializeFunc( inInitFunc ), mParam( inParam )
{
}
__device__ void DeviceStdBrownGenerator::initialize( std::size_t inSequence,
                                                     std::size_t inSeed )
{
    mInitializeFunc( inSequence, inSeed, &mRandState, mParam );
}
__device__ void DeviceStdBrownGenerator::operator()( double* inpResult )
{
    mGenerateFunc( inpResult, &mRandState, mParam );
}

__device__ void generateRandomPlain( double* inpResult,
                                     curandState* inRandState, void* inParam )
{
    *inpResult = curand_normal_double( inRandState );
}
__device__ void initializeRandomPlain( std::size_t inSequence,
                                       std::size_t inSeed,
                                       curandState* inRandState, void* inParam )
{
    curand_init( inSeed, inSequence, 0, inRandState );
}

__device__ void generateRandomAntithetic( double* inpResult,
                                          curandState* inRandState,
                                          void* inParam )
{
    DeviceStdBrownAntitheticParam* lParam =
        reinterpret_cast<DeviceStdBrownAntitheticParam*>( inParam );
    if ( lParam->mIsNextNew )
    {
        lParam->mIsNextNew       = false;
        lParam->mPrevRandomValue = curand_normal_double( inRandState );
        *inpResult               = lParam->mPrevRandomValue;
    }
    else
    {
        lParam->mIsNextNew = true;
        *inpResult         = -lParam->mPrevRandomValue;
    }
}
__device__ void initializeRandomAntithetic( std::size_t inSequence,
                                            std::size_t inSeed,
                                            curandState* inRandState,
                                            void* inParam )
{
    DeviceStdBrownAntitheticParam* lParam =
        reinterpret_cast<DeviceStdBrownAntitheticParam*>( inParam );
    lParam->mIsNextNew = true;
    curand_init( inSeed, inSequence, 0, inRandState );
}

__device__ DeviceStdBrownGenerator
createDeviceStdBrownGenerator( std::size_t inNumberRand, void* inParam )
{
    if ( inNumberRand == PROCESSCUDA_RANDOM_ANTITHETIC )
    {
        return DeviceStdBrownGenerator( generateRandomAntithetic,
                                        initializeRandomAntithetic, inParam );
    }
    return DeviceStdBrownGenerator( generateRandomPlain, initializeRandomPlain,
                                    inParam );
}

__device__ std::size_t sizeofDeviceParam( std::size_t inNumberRand )
{
    if ( inNumberRand == PROCESSCUDA_RANDOM_ANTITHETIC )
    {
        return sizeof( DeviceStdBrownAntitheticParam );
    }
    return 0;
}

}  // namespace OneFactor

}  // namespace ShortRate