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
namespace Math
{
Vec::Vec( int inNSize ) : mNSize( inNSize ), mData( inNSize ) {}
Vec::Vec( std::initializer_list<double> inVal ) :
    mNSize( inVal.size() ), mData( inVal.size() )
{
    int i = 0;
    for ( const auto& val : inVal ) { mData[i++] = val; }
}
Vec::Vec( const std::vector<double>& inVal ) :
    mNSize( inVal.size() ), mData( inVal.size() )
{
    int i = 0;
    for ( const auto& val : inVal ) { mData[i++] = val; }
}
double& Vec::operator()( int i ) { return mData[i]; }
const double& Vec::operator()( int i ) const { return mData[i]; }
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
const Vec& Vec::print() const
{
    for ( int i = 0; i < mNSize; i++ ) { std::cout << mData[i] << " "; }
    std::cout << std::endl;
    return *this;
}
Vec& Vec::print()
{
    for ( int i = 0; i < mNSize; i++ ) { std::cout << mData[i] << " "; }
    std::cout << std::endl;
    return *this;
}

Mat::Mat( int inNRow, int inNCol ) :
    mNRow( inNRow ), mNCol( inNCol ), mData( inNRow * inNCol )
{
}
Mat::Mat( std::initializer_list<std::initializer_list<double>> inVal ) :
    mNRow( inVal.size() ),
    mNCol( inVal.begin()->size() ),
    mData( inVal.size() * inVal.begin()->size() )
{
    int iRow = 0;
    for ( const auto& row : inVal )
    {
        assert( row.size() == mNCol );
        int iCol = 0;
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
    int i = 0;
    for ( const auto& row : inVal )
    {
        assert( row.size() == mNCol );
        for ( const auto& val : row ) { mData[i++] = val; }
    }
}
double& Mat::operator()( int i, int j ) { return mData[i + j * mNRow]; }
const double& Mat::operator()( int i, int j ) const
{
    return mData[i + j * mNRow];
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

Mat Mat::transpose() const
{
    Mat lResult( mNCol, mNRow );
    for ( int i = 0; i < mNRow; ++i )
    {
        for ( int j = 0; j < mNCol; ++j )
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
    for ( int i = 0; i < mNRow; i++ )
    {
        for ( int j = 0; j < mNCol; j++ )
        {
            std::cout << mData[i + j * mNRow] << " ";
        }
        std::cout << std::endl;
    }
    return *this;
}
Mat& Mat::print()
{
    for ( int i = 0; i < mNRow; i++ )
    {
        for ( int j = 0; j < mNCol; j++ )
        {
            std::cout << mData[i + j * mNRow] << " ";
        }
        std::cout << std::endl;
    }
    return *this;
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

}  // namespace Math