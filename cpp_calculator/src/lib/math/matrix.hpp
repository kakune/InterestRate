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

namespace Math
{

class Vec;
class Mat;

class Vec
{
private:
    std::valarray<double> mData;
    std::size_t mNSize;
    friend Mat;

public:
    Vec( std::size_t inNSize, double inVal = 0 );
    Vec( std::initializer_list<double> inVal );
    Vec( const std::vector<double>& inVal );
    Vec( const std::valarray<double>& inVal );

    double& operator()( std::size_t i );
    const double& operator()( std::size_t i ) const;

    const Vec& operator+() const&;
    Vec operator+() &&;
    Vec operator-() const&;
    Vec operator-() &&;

    Vec& operator+=( const Vec& inVec );
    Vec operator+( const Vec& inRhs ) const;
    Vec& operator-=( const Vec& inVec );
    Vec operator-( const Vec& inRhs ) const;
    Vec& operator*=( const Vec& inVec );
    Vec operator*( const Vec& inRhs ) const;
    Vec& operator/=( const Vec& inVec );
    Vec operator/( const Vec& inRhs ) const;

    Vec& operator+=( double inVal );
    Vec operator+( double inRhs ) const;
    Vec& operator-=( double inVal );
    Vec operator-( double inRhs ) const;
    Vec& operator*=( double inVal );
    Vec operator*( double inRhs ) const;
    Vec& operator/=( double inVal );
    Vec operator/( double inRhs ) const;

    friend Vec operator+( double inLhs, const Vec& inRhs );
    friend Vec operator-( double inLhs, const Vec& inRhs );
    friend Vec operator*( double inLhs, const Vec& inRhs );
    friend Vec operator/( double inLhs, const Vec& inRhs );

    std::size_t size() const;
    double sum() const;
    double sum( std::size_t inIndBegin, std::size_t inIndEnd ) const;
    double min() const;
    double max() const;
    const Vec& print() const;
    Vec& print();

    Vec& solveEqLCholesky( const Mat& inL );

    friend Vec sqrt( const Vec& inVec );
    friend Vec exp( const Vec& inVec );
    friend Vec log( const Vec& inVec );
    friend Vec abs( const Vec& inVec );
    friend Vec pow( const Vec& inVec, double inExponent );

    friend void solveEqPositiveDefinite( Mat& inMat, Vec& inVec );
    friend double dot( const Vec& inLhs, const Vec& inRhs );
    friend Vec dot( const Mat& inLhs, const Vec& inRhs );
    friend Vec dot( const Vec& inLhs, const Mat& inRhs );
    friend Vec solveEqLowerTriangular( const Mat& inMat, const Vec& inVec );
    friend Mat dotVecVecToMat( const Vec& inLhs, const Vec& inRhs );

    Vec& multiplyUpperMatFromLeft( const Mat& inLhs );
    Vec& multiplyLowerMatFromLeft( const Mat& inLhs );
    friend Vec dotUpperMatVec( const Mat& inLhs, Vec inRhs );
    friend Vec dotLowerMatVec( const Mat& inLhs, Vec inRhs );
};

class Mat
{
private:
    std::valarray<double> mData;
    std::size_t mNRow, mNCol;
    friend Vec;

public:
    Mat( std::size_t inNRow, std::size_t inNCol, double inVal = 0 );
    Mat( std::size_t inNRow, std::size_t inNCol,
         const std::valarray<double>& inVal );
    Mat( std::initializer_list<std::initializer_list<double>> inVal );
    Mat( const std::vector<std::vector<double>>& inVal );

    double& operator()( std::size_t i, std::size_t j );
    const double& operator()( std::size_t i, std::size_t j ) const;
    Vec getRow(std::size_t inIndRow ) const;
    Vec getCol(std::size_t inIndCol ) const;

    const Mat& operator+() const&;
    Mat operator+() &&;
    Mat operator-() const&;
    Mat operator-() &&;

    Mat& operator+=( const Mat& inMat );
    Mat operator+( const Mat& inRhs ) const;
    Mat& operator-=( const Mat& inMat );
    Mat operator-( const Mat& inRhs ) const;
    Mat& operator*=( const Mat& inMat );
    Mat operator*( const Mat& inRhs ) const;
    Mat& operator/=( const Mat& inMat );
    Mat operator/( const Mat& inRhs ) const;

    Mat& operator+=( double inVal );
    Mat operator+( double inRhs ) const;
    Mat& operator-=( double inVal );
    Mat operator-( double inRhs ) const;
    Mat& operator*=( double inVal );
    Mat operator*( double inRhs ) const;
    Mat& operator/=( double inVal );
    Mat operator/( double inRhs ) const;

    friend Mat operator+( double inLhs, const Mat& inRhs );
    friend Mat operator-( double inLhs, const Mat& inRhs );
    friend Mat operator*( double inLhs, const Mat& inRhs );
    friend Mat operator/( double inLhs, const Mat& inRhs );

    std::size_t sizeRow() const;
    std::size_t sizeCol() const;
    Mat transpose() const;
    const Mat& print() const;
    Mat& print();

    std::pair<Vec, Mat> symLargeEigens( std::size_t inNumEigen );

    friend Vec& Vec::solveEqLCholesky( const Mat& inL );
    Mat& choleskyDecompose();

    friend void solveEqPositiveDefinite( Mat& inMat, Vec& inVec );
    friend Vec dot( const Mat& inLhs, const Vec& inRhs );
    friend Vec dot( const Vec& inLhs, const Mat& inRhs );
    friend Mat dot( const Mat& inLhs, const Mat& inRhs );
    friend Mat choleskyDecompose( const Mat& inMat );
    friend Vec solveEqLowerTriangular( const Mat& inMat, const Vec& inVec );

    friend Vec& Vec::multiplyUpperMatFromLeft( const Mat& inLhs );
    friend Vec& Vec::multiplyLowerMatFromLeft( const Mat& inLhs );
    friend Vec dotUpperMatVec( const Mat& inLhs, Vec inRhs );
    friend Vec dotLowerMatVec( const Mat& inLhs, Vec inRhs );
};

Mat unitMat( std::size_t inSize, double inVal = 1.0 );
Mat diagMat( const Vec& inVec );

}  // namespace Math

#endif