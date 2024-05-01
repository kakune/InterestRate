/**
 * @file matrix.hpp
 * @brief This implements matrix.
 * @author kakune
 * @date 3/14/2024
 */

#include "math/matrix.hpp"

#ifndef NUSE_MKL
#include "mkl_cblas.h"
#include "mkl_lapacke.h"
#else
#include <cblas.h>
#include <lapacke.h>
#endif
#include <cassert>
#include <iostream>
namespace Math
{

double& Mat::operator[]( std::size_t i ) { return mData[i]; }
const double& Mat::operator[]( std::size_t i ) const { return mData[i]; }

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
Mat& Mat::operator-=( const Mat& inMat )
{
    assert( mNRow == inMat.mNRow );
    assert( mNCol == inMat.mNCol );
    mData -= inMat.mData;
    return *this;
}
Mat& Mat::operator*=( const Mat& inMat )
{
    assert( mNRow == inMat.mNRow );
    assert( mNCol == inMat.mNCol );
    mData *= inMat.mData;
    return *this;
}
Mat& Mat::operator/=( const Mat& inMat )
{
    assert( mNRow == inMat.mNRow );
    assert( mNCol == inMat.mNCol );
    mData /= inMat.mData;
    return *this;
}

// Mat Mat::operator+( const Mat& inRhs ) const
// {
//     assert( mNRow == inRhs.mNRow );
//     assert( mNCol == inRhs.mNCol );
//     Mat lResult = *this;
//     lResult += inRhs;
//     return lResult;
// }
// Mat Mat::operator-( const Mat& inRhs ) const
// {
//     assert( mNRow == inRhs.mNRow );
//     assert( mNCol == inRhs.mNCol );
//     Mat lResult = *this;
//     lResult -= inRhs;
//     return lResult;
// }
// Mat Mat::operator*( const Mat& inRhs ) const
// {
//     assert( mNRow == inRhs.mNRow );
//     assert( mNCol == inRhs.mNCol );
//     Mat lResult = *this;
//     lResult *= inRhs;
//     return lResult;
// }
// Mat Mat::operator/( const Mat& inRhs ) const
// {
//     assert( mNRow == inRhs.mNRow );
//     assert( mNCol == inRhs.mNCol );
//     Mat lResult = *this;
//     lResult /= inRhs;
//     return lResult;
// }

Mat& Mat::operator+=( double inVal )
{
    mData += inVal;
    return *this;
}
Mat& Mat::operator-=( double inVal )
{
    mData -= inVal;
    return *this;
}
Mat& Mat::operator*=( double inVal )
{
    mData *= inVal;
    return *this;
}
Mat& Mat::operator/=( double inVal )
{
    mData /= inVal;
    return *this;
}

// Mat Mat::operator+( double inRhs ) const
// {
//     Mat lResult( *this );
//     lResult.mData += inRhs;
//     return lResult;
// }
// Mat Mat::operator-( double inRhs ) const
// {
//     Mat lResult( *this );
//     lResult.mData -= inRhs;
//     return lResult;
// }
// Mat Mat::operator*( double inRhs ) const
// {
//     Mat lResult( *this );
//     lResult.mData *= inRhs;
//     return lResult;
// }
// Mat Mat::operator/( double inRhs ) const
// {
//     Mat lResult( *this );
//     lResult.mData /= inRhs;
//     return lResult;
// }

// Mat operator+( double inLhs, const Mat& inRhs )
// {
//     Mat lResult( inRhs );
//     lResult.mData += inLhs;
//     return lResult;
// }
// Mat operator-( double inLhs, const Mat& inRhs )
// {
//     Mat lResult( inRhs.mNRow, inRhs.mNCol );
//     lResult.mData = ( inLhs - inRhs.mData );
//     return lResult;
// }
// Mat operator*( double inLhs, const Mat& inRhs )
// {
//     Mat lResult( inRhs );
//     lResult.mData *= inLhs;
//     return lResult;
// }
// Mat operator/( double inLhs, const Mat& inRhs )
// {
//     Mat lResult( inRhs.mNRow, inRhs.mNCol );
//     lResult.mData = ( inLhs / inRhs.mData );
//     return lResult;
// }

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

std::pair<Vec, Mat> Mat::symLargeEigens( std::size_t inNumEigen )
{
    Vec lEigenValues( inNumEigen );
    Mat lEigenVectors( mNRow, inNumEigen );
    int lNFound;
    std::vector<int> lIndFailed( inNumEigen );

    int lInfo = LAPACKE_dsyevx(
        LAPACK_COL_MAJOR, 'V', 'I', 'U', mNRow, &mData[0], mNRow, 0.0, 0.0,
        mNRow - inNumEigen + 1, mNRow, 0.0, &lNFound, &lEigenValues[0],
        &lEigenVectors[0], lEigenVectors.mNRow, &lIndFailed[0] );
    if ( lInfo != 0 )
    {
        std::cerr << "Error: Math::Vec::symLargeEigens()" << std::endl
                  << "Error code : " << lInfo << std::endl;
        exit( 1 );
    }

    return std::make_pair( lEigenValues, lEigenVectors );
}

Vec solveEqLCholesky( const Mat& inLowerMat, Vec inVec )
{
    assert( inLowerMat.mNCol == inLowerMat.mNRow &&
            inLowerMat.mNRow == inVec.size() );
    int info = LAPACKE_dpotrs( LAPACK_COL_MAJOR, 'L', inLowerMat.mNRow, 1,
                               &inLowerMat[0], inLowerMat.mNRow, &inVec[0],
                               inVec.size() );
    if ( info != 0 )
    {
        std::cerr << "Error: Math::Vec::solveEqLCholesky()" << std::endl
                  << "Error code : " << info << std::endl;
        exit( 1 );
    }
    return inVec;
}

void solveEqPositiveDefinite( Mat& inMat, Vec& inVec )
{
    assert( inMat.mNCol == inMat.mNRow && inMat.mNRow == inVec.size() );
    inMat.choleskyDecompose();
    inVec = solveEqLCholesky( inMat, inVec );
}
double dot( const Vec& inLhs, const Vec& inRhs )
{
    return ( inLhs * inRhs ).sum();
}
Vec dot( const Mat& inLhs, const Vec& inRhs )
{
    assert( inLhs.mNCol == inRhs.size() );
    Vec lResult( inLhs.mNRow );
    cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, inLhs.mNRow, 1,
                 inLhs.mNCol, 1.0, &inLhs[0], inLhs.mNRow, &inRhs[0],
                 inRhs.size(), 0.0, &lResult[0], lResult.size() );
    return lResult;
}
Vec dot( const Vec& inLhs, const Mat& inRhs )
{
    assert( inLhs.size() == inRhs.mNRow );
    Vec lResult( inRhs.mNCol );
    cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, 1, inRhs.mNCol,
                 inLhs.size(), 1.0, &inLhs[0], 1, &inRhs[0], inRhs.mNRow, 0.0,
                 &lResult[0], 1 );
    return lResult;
}
Mat dot( const Mat& inLhs, const Mat& inRhs )
{
    assert( inLhs.mNCol == inRhs.mNRow );
    Mat lResult( inLhs.mNRow, inRhs.mNCol );
    cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, inLhs.mNRow,
                 inRhs.mNCol, inLhs.mNCol, 1.0, &inLhs[0], inLhs.mNRow,
                 &inRhs[0], inRhs.mNRow, 0.0, &lResult[0], lResult.mNRow );
    return lResult;
}
Mat dotVecVecToMat( const Vec& inLhs, const Vec& inRhs )
{
    std::valarray<double> lResult( 0.0, inLhs.size() * inRhs.size() );
    cblas_dger( CblasColMajor, inLhs.size(), inRhs.size(), 1.0, &inLhs[0], 1,
                &inRhs[0], 1, &lResult[0], inLhs.size() );
    return Mat( inLhs.size(), inRhs.size(), lResult );
}

Mat choleskyDecompose( const Mat& inMat )
{
    assert( inMat.mNRow == inMat.mNCol );
    Mat lResult = inMat;
    return lResult.choleskyDecompose();
}

Vec solveEqLowerTriangular( const Mat& inMat, Vec inVec )
{
    return solveEqLCholesky( inMat, inVec );
}

Vec dotUpperMatVec( const Mat& inLhs, Vec inRhs )
{
    cblas_dtrmv( CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                 inLhs.mNRow, &inLhs[0], inLhs.mNRow, &inRhs[0], 1 );
    return inRhs;
}
Vec dotLowerMatVec( const Mat& inLhs, Vec inRhs )
{
    cblas_dtrmv( CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
                 inLhs.mNRow, &inLhs[0], inLhs.mNRow, &inRhs[0], 1 );
    return inRhs;
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
        lResult( i, i ) = inVec[i];
    }
    return lResult;
}
const Vec& print( const Vec& inVec )
{
    for ( std::size_t i = 0; i < inVec.size(); i++ )
    {
        std::cout << inVec[i] << " ";
    }
    std::cout << std::endl;
    return inVec;
}
const std::vector<double>& print( const std::vector<double>& inVec )
{
    for ( std::size_t i = 0; i < inVec.size(); i++ )
    {
        std::cout << inVec[i] << " ";
    }
    std::cout << std::endl;
    return inVec;
}
const Mat& print( const Mat& inMat )
{
    for ( std::size_t i = 0; i < inMat.sizeRow(); i++ )
    {
        for ( std::size_t j = 0; j < inMat.sizeCol(); j++ )
        {
            std::cout << inMat( i, j ) << " ";
        }
        std::cout << std::endl;
    }
    return inMat;
}
const std::vector<std::vector<double>>& print(
    const std::vector<std::vector<double>>& inMat )
{
    for ( const auto& lVec : inMat )
    {
        for ( const auto& lVal : lVec ) { std::cout << lVal << " "; }
        std::cout << std::endl;
    }
    return inMat;
}

}  // namespace Math