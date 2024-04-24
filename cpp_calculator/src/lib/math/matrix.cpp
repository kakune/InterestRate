/**
 * @file matrix.hpp
 * @brief This implements matrix.
 * @author kakune
 * @date 3/14/2024
 */

#include "math/matrix.hpp"

#ifndef NUSE_MKL
#include "mkl.h"
#else
#include <cblas.h>
#include <lapacke.h>
#endif
#include <cassert>
#include <iostream>
#include <numeric>
namespace Math
{
Vec::Vec( std::size_t inNSize, double inVal ) :
    mNSize( inNSize ), mData( inVal, inNSize )
{
}
Vec::Vec( std::initializer_list<double> inVal ) :
    mNSize( inVal.size() ), mData( inVal.size() )
{
    std::size_t i = 0;
    for ( const auto& val : inVal ) { mData[i++] = val; }
}
Vec::Vec( const std::vector<double>& inVal ) :
    mNSize( inVal.size() ), mData( &inVal.data()[0], inVal.size() )
{
}
Vec::Vec( const std::valarray<double>& inVal ) :
    mNSize( inVal.size() ), mData( &inVal[0], inVal.size() )
{
}
double& Vec::operator()( std::size_t i ) { return mData[i]; }
const double& Vec::operator()( std::size_t i ) const { return mData[i]; }

const Vec& Vec::operator+() const& { return *this; }
Vec Vec::operator+() &&
{
    Vec lResult = std::move( *this );
    return lResult;
}
Vec Vec::operator-() const& { return Vec( -mData ); }
Vec Vec::operator-() &&
{
    Vec lResult   = std::move( *this );
    lResult.mData = -lResult.mData;
    return lResult;
}
Vec& Vec::operator+=( const Vec& inVec )
{
    assert( mNSize == inVec.mNSize );
    mData += inVec.mData;
    return *this;
}
Vec Vec::operator+( const Vec& inRhs ) const
{
    assert( mNSize == inRhs.mNSize );
    Vec lResult = *this;
    lResult += inRhs;
    return lResult;
}
Vec& Vec::operator-=( const Vec& inVec )
{
    assert( mNSize == inVec.mNSize );
    mData -= inVec.mData;
    return *this;
}
Vec Vec::operator-( const Vec& inRhs ) const
{
    assert( mNSize == inRhs.mNSize );
    Vec lResult = *this;
    lResult -= inRhs;
    return lResult;
}
Vec& Vec::operator*=( const Vec& inVec )
{
    assert( mNSize == inVec.mNSize );
    mData *= inVec.mData;
    return *this;
}
Vec Vec::operator*( const Vec& inRhs ) const
{
    assert( mNSize == inRhs.mNSize );
    Vec lResult = *this;
    lResult *= inRhs;
    return lResult;
}
Vec& Vec::operator/=( const Vec& inVec )
{
    assert( mNSize == inVec.mNSize );
    mData /= inVec.mData;
    return *this;
}
Vec Vec::operator/( const Vec& inRhs ) const
{
    assert( mNSize == inRhs.mNSize );
    Vec lResult = *this;
    lResult /= inRhs;
    return lResult;
}
Vec& Vec::operator+=( double inVal )
{
    mData += inVal;
    return *this;
}
Vec Vec::operator+( double inRhs ) const
{
    Vec lResult = *this;
    lResult += inRhs;
    return lResult;
}
Vec& Vec::operator-=( double inVal )
{
    mData -= inVal;
    return *this;
}
Vec Vec::operator-( double inRhs ) const
{
    Vec lResult = *this;
    lResult -= inRhs;
    return lResult;
}
Vec& Vec::operator*=( double inVal )
{
    mData *= inVal;
    return *this;
}
Vec Vec::operator*( double inRhs ) const
{
    Vec lResult = *this;
    lResult *= inRhs;
    return lResult;
}
Vec& Vec::operator/=( double inVal )
{
    mData /= inVal;
    return *this;
}
Vec Vec::operator/( double inRhs ) const
{
    Vec lResult = *this;
    lResult /= inRhs;
    return lResult;
}
Vec operator+( double inLhs, const Vec& inRhs )
{
    Vec lResult( inRhs );
    lResult.mData += inLhs;
    return lResult;
}
Vec operator-( double inLhs, const Vec& inRhs )
{
    Vec lResult( -inRhs.mData );
    lResult.mData += inLhs;
    return lResult;
}
Vec operator*( double inLhs, const Vec& inRhs )
{
    Vec lResult( inRhs );
    lResult.mData *= inLhs;
    return lResult;
}
Vec operator/( double inLhs, const Vec& inRhs )
{
    Vec lResult( inLhs / inRhs.mData );
    return lResult;
}
std::size_t Vec::size() const { return mNSize; }
double Vec::sum() const { return mData.sum(); }
double Vec::sum( std::size_t inIndBegin, std::size_t inIndEnd ) const
{
    std::slice lSlice( inIndBegin, inIndEnd - inIndBegin, 1 );
    return mData[lSlice].sum();
}
double Vec::min() const { return mData.min(); }
double Vec::max() const { return mData.max(); }
const Vec& Vec::print() const
{
    for ( std::size_t i = 0; i < mNSize; i++ ) { std::cout << mData[i] << " "; }
    std::cout << std::endl;
    return *this;
}
Vec& Vec::print()
{
    for ( std::size_t i = 0; i < mNSize; i++ ) { std::cout << mData[i] << " "; }
    std::cout << std::endl;
    return *this;
}
Vec sqrt( const Vec& inVec ) { return Vec( sqrt( inVec.mData ) ); }
Vec exp( const Vec& inVec ) { return Vec( exp( inVec.mData ) ); }
Vec log( const Vec& inVec ) { return Vec( log( inVec.mData ) ); }
Vec abs( const Vec& inVec ) { return Vec( abs( inVec.mData ) ); }
Vec pow( const Vec& inVec, double inExponent )
{
    return Vec( pow( inVec.mData, inExponent ) );
}

Mat::Mat( std::size_t inNRow, std::size_t inNCol, double inVal ) :
    mNRow( inNRow ), mNCol( inNCol ), mData( inVal, inNRow * inNCol )
{
}
Mat::Mat( std::size_t inNRow, std::size_t inNCol,
          const std::valarray<double>& inVal ) :
    mNRow( inNRow ), mNCol( inNCol ), mData( inVal )
{
}
Mat::Mat( std::initializer_list<std::initializer_list<double>> inVal ) :
    mNRow( inVal.size() ),
    mNCol( inVal.begin()->size() ),
    mData( inVal.size() * inVal.begin()->size() )
{
    std::size_t iRow = 0;
    for ( const auto& row : inVal )
    {
        assert( row.size() == mNCol );
        std::size_t iCol = 0;
        for ( const auto& val : row )
        {
            mData[iRow + ( iCol++ * mNRow )] = val;
        }
        ++iRow;
    }
}
Mat::Mat( const std::vector<std::vector<double>>& inVal ) :
    mNRow( inVal.size() ),
    mNCol( inVal.begin()->size() ),
    mData( inVal.size() * inVal.begin()->size() )
{
    std::size_t i = 0;
    for ( const auto& row : inVal )
    {
        assert( row.size() == mNCol );
        for ( const auto& val : row ) { mData[i++] = val; }
    }
}

double& Mat::operator()( std::size_t i, std::size_t j )
{
    return mData[i + j * mNRow];
}
const double& Mat::operator()( std::size_t i, std::size_t j ) const
{
    return mData[i + j * mNRow];
}
Vec Mat::getRow( std::size_t inIndRow ) const
{
    return Vec(
        std::valarray<double>( mData[std::slice( inIndRow, mNCol, mNRow )] ) );
}
Vec Mat::getCol( std::size_t inIndCol ) const
{
    return Vec( std::valarray<double>(
        mData[std::slice( inIndCol * mNRow, mNRow, 1 )] ) );
}
const Mat& Mat::operator+() const& { return *this; }
Mat Mat::operator+() &&
{
    Mat lResult = std::move( *this );
    return lResult;
}
Mat Mat::operator-() const&
{
    Mat lResult( *this );
    lResult.mData = -mData;
    return lResult;
}
Mat Mat::operator-() &&
{
    Mat lResult   = std::move( *this );
    lResult.mData = -lResult.mData;
    return lResult;
}
Mat& Mat::operator+=( const Mat& inMat )
{
    assert( mNRow == inMat.mNRow );
    assert( mNCol == inMat.mNCol );
    mData += inMat.mData;
    return *this;
}
Mat Mat::operator+( const Mat& inRhs ) const
{
    assert( mNRow == inRhs.mNRow );
    assert( mNCol == inRhs.mNCol );
    Mat lResult = *this;
    lResult += inRhs;
    return lResult;
}
Mat& Mat::operator-=( const Mat& inMat )
{
    assert( mNRow == inMat.mNRow );
    assert( mNCol == inMat.mNCol );
    mData -= inMat.mData;
    return *this;
}
Mat Mat::operator-( const Mat& inRhs ) const
{
    assert( mNRow == inRhs.mNRow );
    assert( mNCol == inRhs.mNCol );
    Mat lResult = *this;
    lResult -= inRhs;
    return lResult;
}
Mat& Mat::operator*=( const Mat& inMat )
{
    assert( mNRow == inMat.mNRow );
    assert( mNCol == inMat.mNCol );
    mData *= inMat.mData;
    return *this;
}
Mat Mat::operator*( const Mat& inRhs ) const
{
    assert( mNRow == inRhs.mNRow );
    assert( mNCol == inRhs.mNCol );
    Mat lResult = *this;
    lResult *= inRhs;
    return lResult;
}
Mat& Mat::operator/=( const Mat& inMat )
{
    assert( mNRow == inMat.mNRow );
    assert( mNCol == inMat.mNCol );
    mData /= inMat.mData;
    return *this;
}
Mat Mat::operator/( const Mat& inRhs ) const
{
    assert( mNRow == inRhs.mNRow );
    assert( mNCol == inRhs.mNCol );
    Mat lResult = *this;
    lResult /= inRhs;
    return lResult;
}
Mat& Mat::operator+=( double inVal )
{
    mData += inVal;
    return *this;
}
Mat Mat::operator+( double inRhs ) const
{
    Mat lResult( *this );
    lResult.mData += inRhs;
    return lResult;
}
Mat& Mat::operator-=( double inVal )
{
    mData -= inVal;
    return *this;
}
Mat Mat::operator-( double inRhs ) const
{
    Mat lResult( *this );
    lResult.mData -= inRhs;
    return lResult;
}
Mat& Mat::operator*=( double inVal )
{
    mData *= inVal;
    return *this;
}
Mat Mat::operator*( double inRhs ) const
{
    Mat lResult( *this );
    lResult.mData *= inRhs;
    return lResult;
}
Mat& Mat::operator/=( double inVal )
{
    mData /= inVal;
    return *this;
}
Mat Mat::operator/( double inRhs ) const
{
    Mat lResult( *this );
    lResult.mData /= inRhs;
    return lResult;
}
Mat operator+( double inLhs, const Mat& inRhs )
{
    Mat lResult( inRhs );
    lResult.mData += inLhs;
    return lResult;
}
Mat operator-( double inLhs, const Mat& inRhs )
{
    Mat lResult( inRhs.mNRow, inRhs.mNCol );
    lResult.mData = ( inLhs - inRhs.mData );
    return lResult;
}
Mat operator*( double inLhs, const Mat& inRhs )
{
    Mat lResult( inRhs );
    lResult.mData *= inLhs;
    return lResult;
}
Mat operator/( double inLhs, const Mat& inRhs )
{
    Mat lResult( inRhs.mNRow, inRhs.mNCol );
    lResult.mData = ( inLhs / inRhs.mData );
    return lResult;
}

std::size_t Mat::sizeRow() const { return mNRow; }
std::size_t Mat::sizeCol() const { return mNCol; }
Mat Mat::transpose() const
{
    Mat lResult( mNCol, mNRow );
    for ( std::size_t i = 0; i < mNRow; ++i )
    {
        for ( std::size_t j = 0; j < mNCol; ++j )
        {
            lResult( j, i ) = ( *this )( i, j );
        }
    }
    return lResult;
}

Mat& Mat::choleskyDecompose()
{
    assert( mNRow == mNCol );
    int info = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'L', mNCol, &mData[0], mNRow );
    if ( info != 0 )
    {
        std::cerr << "Error: Math::choleskyDecompose()" << std::endl
                  << "Error code : " << info << std::endl;
        exit( 1 );
    }
    return *this;
}

const Mat& Mat::print() const
{
    for ( std::size_t i = 0; i < mNRow; i++ )
    {
        for ( std::size_t j = 0; j < mNCol; j++ )
        {
            std::cout << mData[i + j * mNRow] << " ";
        }
        std::cout << std::endl;
    }
    return *this;
}
Mat& Mat::print()
{
    for ( std::size_t i = 0; i < mNRow; i++ )
    {
        for ( std::size_t j = 0; j < mNCol; j++ )
        {
            std::cout << mData[i + j * mNRow] << " ";
        }
        std::cout << std::endl;
    }
    return *this;
}

std::pair<Vec, Mat> Mat::symLargeEigens( std::size_t inNumEigen )
{
    Vec lEigenValues( inNumEigen );
    Mat lEigenVectors( mNRow, inNumEigen );
    int lNFound;
    std::vector<int> lIndFailed( inNumEigen );

    int lInfo = LAPACKE_dsyevx(
        LAPACK_COL_MAJOR, 'V', 'I', 'U', mNRow, &mData[0], mNRow, 0.0, 0.0,
        mNRow - inNumEigen + 1, mNRow, 0.0, &lNFound, &lEigenValues.mData[0],
        &lEigenVectors.mData[0], lEigenVectors.mNRow, &lIndFailed[0] );
    if ( lInfo != 0 )
    {
        std::cerr << "Error: Math::Vec::symLargeEigens()" << std::endl
                  << "Error code : " << lInfo << std::endl;
        exit( 1 );
    }

    return std::make_pair( lEigenValues, lEigenVectors );
}

Vec& Vec::solveEqLCholesky( const Mat& inL )
{
    assert( inL.mNCol == inL.mNRow && inL.mNRow == mNSize );
    int info = LAPACKE_dpotrs( LAPACK_COL_MAJOR, 'L', inL.mNRow, 1,
                               &inL.mData[0], inL.mNRow, &mData[0], mNSize );
    if ( info != 0 )
    {
        std::cerr << "Error: Math::Vec::solveEqLCholesky()" << std::endl
                  << "Error code : " << info << std::endl;
        exit( 1 );
    }
    return *this;
}

void solveEqPositiveDefinite( Mat& inMat, Vec& inVec )
{
    assert( inMat.mNCol == inMat.mNRow && inMat.mNRow == inVec.mNSize );
    inMat.choleskyDecompose();
    inVec.solveEqLCholesky( inMat );
}
double dot( const Vec& inLhs, const Vec& inRhs )
{
    return std::inner_product( &inLhs.mData[0], &inLhs.mData[0] + inLhs.mNSize,
                               &inRhs.mData[0], 0.0 );
}
Vec dot( const Mat& inLhs, const Vec& inRhs )
{
    assert( inLhs.mNCol == inRhs.mNSize );
    Vec lResult( inLhs.mNRow );
    cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, inLhs.mNRow, 1,
                 inLhs.mNCol, 1.0, &inLhs.mData[0], inLhs.mNRow,
                 &inRhs.mData[0], inRhs.mNSize, 0.0, &lResult.mData[0],
                 lResult.mNSize );
    return lResult;
}
Vec dot( const Vec& inLhs, const Mat& inRhs )
{
    assert( inLhs.mNSize == inRhs.mNRow );
    Vec lResult( inRhs.mNCol );
    cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, 1, inRhs.mNCol,
                 inLhs.mNSize, 1.0, &inLhs.mData[0], 1, &inRhs.mData[0],
                 inRhs.mNRow, 0.0, &lResult.mData[0], 1 );
    return lResult;
}
Mat dot( const Mat& inLhs, const Mat& inRhs )
{
    assert( inLhs.mNCol == inRhs.mNRow );
    Mat lResult( inLhs.mNRow, inRhs.mNCol );
    cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, inLhs.mNRow,
                 inRhs.mNCol, inLhs.mNCol, 1.0, &inLhs.mData[0], inLhs.mNRow,
                 &inRhs.mData[0], inRhs.mNRow, 0.0, &lResult.mData[0],
                 lResult.mNRow );
    return lResult;
}
Mat dotVecVecToMat( const Vec& inLhs, const Vec& inRhs )
{
    std::valarray<double> lResult( 0.0, inLhs.size() * inRhs.size() );
    cblas_dger( CblasColMajor, inLhs.size(), inRhs.size(), 1.0, &inLhs.mData[0],
                1, &inRhs.mData[0], 1, &lResult[0], inLhs.size() );
    return Mat( inLhs.size(), inRhs.size(), lResult );
}

Mat choleskyDecompose( const Mat& inMat )
{
    assert( inMat.mNRow == inMat.mNCol );
    Mat lResult = inMat;
    return lResult.choleskyDecompose();
}

Vec solveEqLowerTriangular( const Mat& inMat, const Vec& inVec )
{
    Vec lResult = inVec;
    return lResult.solveEqLCholesky( inMat );
}
Vec& Vec::multiplyUpperMatFromLeft( const Mat& inLhs )
{
    cblas_dtrmv( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                 inLhs.mNRow, &inLhs.mData[0], inLhs.mNRow, &mData[0], 1 );
    return *this;
}
Vec& Vec::multiplyLowerMatFromLeft( const Mat& inLhs )
{
    cblas_dtrmv( CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                 inLhs.mNRow, &inLhs.mData[0], inLhs.mNRow, &mData[0], 1 );
    return *this;
}
Vec dotUpperMatVec( const Mat& inLhs, Vec inRhs )
{
    return inRhs.multiplyUpperMatFromLeft( inLhs );
}
Vec dotLowerMatVec( const Mat& inLhs, Vec inRhs )
{
    return inRhs.multiplyLowerMatFromLeft( inLhs );
}

Mat unitMat( std::size_t inSize, double inVal )
{
    Mat lResult( inSize, inSize, 0.0 );
    for ( std::size_t i = 0; i < inSize; ++i ) { lResult( i, i ) = inVal; }
    return lResult;
}
Mat diagMat( const Vec& inVec )
{
    Mat lResult( inVec.size(), inVec.size(), 0.0 );
    for ( std::size_t i = 0; i < inVec.size(); ++i )
    {
        lResult( i, i ) = inVec( i );
    }
    return lResult;
}

}  // namespace Math