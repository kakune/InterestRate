#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

#include "math/matrix.hpp"

double testDotVecVec( Math::Vec inLhs, Math::Vec inRhs )
{
    return Math::dot( inLhs, inRhs );
}
double testDotVecMat( Math::Vec inLhs, Math::Mat inRhs, std::size_t inInd )
{
    return Math::dot( inLhs, inRhs )[inInd];
}
double testDotMatVec( Math::Mat inLhs, Math::Vec inRhs, std::size_t inInd )
{
    return Math::dot( inLhs, inRhs )[inInd];
}
double testDotMatMat( Math::Mat inLhs, Math::Mat inRhs, std::size_t inIndLeft,
                      std::size_t inIndRight )
{
    Math::print( inLhs );
    std::cout << inLhs( 2, 1 ) << std::endl;
    return Math::dot( inLhs, inRhs )( inIndLeft, inIndRight );
}
double testCholeskyDecomposition( Math::Mat inMat, std::size_t inIndLeft,
                                  std::size_t inIndRight )
{
    inMat.choleskyDecompose();
    return inMat( inIndLeft, inIndRight );
}

double testSolveEqPositiveDefinite( Math::Mat inMat, Math::Vec inVec,
                                    std::size_t inInd )
{
    solveEqPositiveDefinite( inMat, inVec );
    return inVec[inInd];
}

TEST( MatrixTest, Dot )
{
    EXPECT_NEAR( 10.0, testDotVecVec( { 1.0, 2.0, 3.0 }, { 3.0, 2.0, 1.0 } ),
                 1e-6 );
    EXPECT_NEAR( 7.0,
                 testDotVecMat( { 1.0, 2.0, 3.0 },
                                {
                                    { 3.0, 2.0, 1.0, 5.0 },
                                    { 2.0, 2.0, 2.0, 3.0 },
                                    { 0.0, 1.0, 2.0, 0.0 },
                                },
                                0 ),
                 1e-6 );
    EXPECT_NEAR( 10.0,
                 testDotMatVec(
                     {
                         { 3.0, 2.0, 1.0 },
                         { 2.0, 2.0, 2.0 },
                         { 0.0, 1.0, 2.0 },
                         { 0.0, 1.0, 2.0 },
                     },
                     { 1.0, 2.0, 3.0 }, 0 ),
                 1e-6 );
    EXPECT_NEAR( 13.0,
                 testDotMatMat(
                     {
                         { 3.0, 2.0, 1.0 },
                         { 2.0, 2.0, 2.0 },
                         { 0.0, 1.0, 2.0 },
                         { 0.0, 1.0, 2.0 },
                     },
                     {
                         { 3.0, 2.0, 1.0, 5.0 },
                         { 2.0, 2.0, 2.0, 3.0 },
                         { 0.0, 1.0, 2.0, 0.0 },
                     },
                     0, 0 ),
                 1e-6 );
}

TEST( MatrixTest, CholeskyDecomposition )
{
    EXPECT_NEAR( 2.236,
                 testCholeskyDecomposition(
                     {
                         { 5.0, 2.0, 3.0 },
                         { 2.0, 7.0, -4.0 },
                         { 3.0, -4.0, 9.0 },
                     },
                     0, 0 ),
                 1e-3 );
}
TEST( MatrixTest, SolveEqTriangular )
{
    EXPECT_NEAR( 3.0,
                 testSolveEqPositiveDefinite(
                     {
                         { 5.0, 2.0, 3.0 },
                         { 2.0, 7.0, -4.0 },
                         { 3.0, -4.0, 9.0 },
                     },
                     { 18.0, 4.0, 22.0 }, 2 ),
                 1e-3 );
}
TEST( MatrixTest, Eigen )
{
    Math::Mat lMat{ { 1.0, 2.0, 3.0 }, { 2.0, 4.0, 1.0 }, { 3.0, 1.0, -2.0 } };
    Math::print( lMat );
    auto [lEigenVec, lEigenMat] = lMat.symLargeEigens( 2 );
    Math::print( lEigenVec );
    Math::print( lEigenMat );
}