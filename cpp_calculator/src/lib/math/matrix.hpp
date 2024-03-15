/**
 * @file matrix.hpp
 * @brief This defines matrix.
 * @author kakune
 * @date 3/13/2024
 */

#ifndef MATH_MATRIX_HPP
#define MATH_MATRIX_HPP

#include <iostream>
#include <numeric>
#include <valarray>

#ifndef NUSE_MKL
#include "mkl.h"
#else
#include "lapacke.h"
#endif

namespace Math
{

class Vec;
class Mat;

class Vec
{
public:
    std::valarray<double> mData;
    int mNSize;

    Vec( int inNSize );
    Vec( std::initializer_list<double> inVal );
    Vec( const std::vector<double>& inVal );

    double& operator()( int i );
    const double& operator()( int i ) const;
    Vec& operator+=( const Vec& inVec );
    Vec operator+( const Vec& inRhs ) const;
    Vec& operator-=( const Vec& inVec );
    Vec operator-( const Vec& inRhs ) const;
    Vec& operator*=( const Vec& inVec );
    Vec operator*( const Vec& inRhs ) const;
    Vec& operator/=( const Vec& inVec );
    Vec operator/( const Vec& inRhs ) const;
    const Vec& print() const;
    Vec& print();

    Vec& solveEqLCholesky( const Mat& inL );
};
class Mat
{
public:
    std::valarray<double> mData;
    int mNRow, mNCol;

    Mat( int inNRow, int inNCol );
    Mat( std::initializer_list<std::initializer_list<double>> inVal );
    Mat( const std::vector<std::vector<double>>& inVal );

    double& operator()( int i, int j );
    const double& operator()( int i, int j ) const;
    Mat& operator+=( const Mat& inMat );
    Mat operator+( const Mat& inRhs ) const;
    Mat& operator-=( const Mat& inMat );
    Mat operator-( const Mat& inRhs ) const;
    Mat& operator*=( const Mat& inMat );
    Mat operator*( const Mat& inRhs ) const;
    Mat& operator/=( const Mat& inMat );
    Mat operator/( const Mat& inRhs ) const;
    Mat transpose() const;
    const Mat& print() const;
    Mat& print();

    Mat& choleskyDecompose();
};

void solveEqPositiveDefinite( Mat& inMat, Vec& inVec );
double dot( const Vec& inLhs, const Vec& inRhs );
Vec dot( const Mat& inLhs, const Vec& inRhs );
Vec dot( const Vec& inLhs, const Mat& inRhs );
Mat dot( const Mat& inLhs, const Mat& inRhs );
Mat choleskyDecompose( const Mat& inMat );
Vec solveEqLowerTriangular( const Mat& inMat, const Vec& inVec );

}  // namespace Math

#endif