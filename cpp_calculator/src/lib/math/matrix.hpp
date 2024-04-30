/**
 * @file matrix.hpp
 * @brief This defines matrix.
 * @author kakune
 * @date 3/13/2024
 */

#ifndef MATH_MATRIX_HPP
#define MATH_MATRIX_HPP

#include <valarray>
#include <vector>

namespace Math
{

using Vec = std::valarray<double>;
class Mat;

inline Vec makeVec( std::size_t inNSize ) { return Vec( inNSize ); }
inline Vec makeVec( std::size_t inNSize, double inVal )
{
    return std::valarray<double>( inVal, inNSize );
}
inline Vec makeVec( const std::vector<double>& inVec )
{
    return std::valarray<double>( inVec.data(), inVec.size() );
}
inline Vec makeVec( std::valarray<double> inVec ) { return inVec; }

class Mat
{
private:
    std::size_t mNRow, mNCol;
    Vec mData;
    double& operator[]( std::size_t i );
    const double& operator[]( std::size_t i ) const;

public:
    Mat( std::size_t inNRow, std::size_t inNCol, double inVal = 0 );
    Mat( std::size_t inNRow, std::size_t inNCol,
         const std::valarray<double>& inVal );
    Mat( std::initializer_list<std::initializer_list<double>> inVal );
    Mat( const std::vector<std::vector<double>>& inVal );

    double& operator()( std::size_t i, std::size_t j );
    const double& operator()( std::size_t i, std::size_t j ) const;

    const Mat& operator+() const&;
    Mat operator+() &&;
    Mat operator-() const&;
    Mat operator-() &&;

    Mat& operator+=( const Mat& inMat );
    Mat& operator-=( const Mat& inMat );
    Mat& operator*=( const Mat& inMat );
    Mat& operator/=( const Mat& inMat );

    // Mat operator+( const Mat& inRhs ) const;
    // Mat operator-( const Mat& inRhs ) const;
    // Mat operator*( const Mat& inRhs ) const;
    // Mat operator/( const Mat& inRhs ) const;

    Mat& operator+=( double inVal );
    Mat& operator-=( double inVal );
    Mat& operator*=( double inVal );
    Mat& operator/=( double inVal );

    // Mat operator+( double inRhs ) const;
    // Mat operator-( double inRhs ) const;
    // Mat operator*( double inRhs ) const;
    // Mat operator/( double inRhs ) const;

    // friend Mat operator+( double inLhs, const Mat& inRhs );
    // friend Mat operator-( double inLhs, const Mat& inRhs );
    // friend Mat operator*( double inLhs, const Mat& inRhs );
    // friend Mat operator/( double inLhs, const Mat& inRhs );

    std::size_t sizeRow() const;
    std::size_t sizeCol() const;
    Mat transpose() const;

    std::pair<Vec, Mat> symLargeEigens( std::size_t inNumEigen );

    friend Vec solveEqLCholesky( const Mat& inLowerMat, Vec inVec );
    Mat& choleskyDecompose();

    friend void solveEqPositiveDefinite( Mat& inMat, Vec& inVec );
    friend Vec dot( const Mat& inLhs, const Vec& inRhs );
    friend Vec dot( const Vec& inLhs, const Mat& inRhs );
    friend Mat dot( const Mat& inLhs, const Mat& inRhs );
    friend Mat choleskyDecompose( const Mat& inMat );
    friend Vec solveEqLowerTriangular( const Mat& inMat, Vec inVec );

    friend Vec dotUpperMatVec( const Mat& inLhs, Vec inRhs );
    friend Vec dotLowerMatVec( const Mat& inLhs, Vec inRhs );
};

Mat unitMat( std::size_t inSize, double inVal = 1.0 );
Mat diagMat( const Vec& inVec );

double dot( const Vec& inLhs, const Vec& inRhs );
Vec dot( const Mat& inLhs, const Vec& inRhs );
Vec dot( const Vec& inLhs, const Mat& inRhs );
Mat dot( const Mat& inLhs, const Mat& inRhs );

const Vec& print( const Vec& inVec );
const std::vector<double>& print( const std::vector<double>& inVec );
const Mat& print( const Mat& inMat );
const std::vector<std::vector<double>>& print(
    const std::vector<std::vector<double>>& inMat );

}  // namespace Math

#endif