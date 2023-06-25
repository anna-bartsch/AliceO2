// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file  BandMatrixSolverNew.h
/// \brief Definition of BandMatrixSolverNew class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_BANDMATRIXSOLVERNEW_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_BANDMATRIXSOLVERNEW_H

#include "GPUCommonDef.h"
#include "GPUCommonRtypes.h"
#include <vector>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <iostream>
#include <math.h>

namespace GPUCA_NAMESPACE
{
namespace gpu
{

/// Linear Equation Solver for a symmetric positive-definite band matrix A[n x n].
///
/// The matrix has a pattern of BandWidthT adjacent non-zero entries right to the diagonal in each row
/// Here is an example with n==10, BandWidthT==4.  (*) means non-zero element, (+) means symmetric element):
///     (****      )
///     (+****     )
///     (++****    )
///     (+++****   )
/// A = ( +++****  )
///     (  +++**** )
///     (   +++****)
///     (    +++***)
///     (     +++**)
///     (      +++*)
///
/// The non-zero matrix elements are stored in [n x BandWidthT] array mA
///
/// The equation to sove is A[n][n] x X[n][Bdim] = B[n][Bdim].
/// During calculations, the initial values of mA and mB get lost, so one can call solve() only once.
/// The solution X is stored in mB.
///
class BandMatrixSolverNew
{
 public:

  /// Consructor
  BandMatrixSolverNew(int N, int Bdim, int BW) : mN(N), mBdim(Bdim), mBandWidth(BW)
  {
    assert(N > 0 && Bdim > 0);
    mA.resize(mN * mBandWidth, 0.);
    mB.resize(mN * mBdim, 0.);
  }

  /// debug tool: init arrays with NaN's
  void initWithNaN()
  {
    // Assign NaN's to ensure that uninitialized elements (for the matrix type 1) are not used in calculations.
    mA.assign(mA.size(), std::numeric_limits<double>::signaling_NaN());
    mB.assign(mB.size(), std::numeric_limits<double>::signaling_NaN());
  }

  /// access to A elements
  double& A(int i, int j)
  {
    auto ij = std::minmax(i, j);
    assert(ij.first >= 0 && ij.second < mN);
    int k = ij.second - ij.first;
    assert(k < mBandWidth);
    return mA[ij.first * mBandWidth + k];
  }

  /// access to B elements
  double& B(int i, int j)
  {
    assert(i >= 0 && i < mN && j >= 0 && j < mBdim);
    return mB[i * mBdim + j];
  }

  /// solve the equation
  void solve();


  bool isFailed = false;

 private:
  
  void triangulateBlock(double AA[], double bb[], int nRows);

  void dioganalizeBlock(double A[], double b[], int nCols);

 private:
  int mN = 0;
  int mBdim = 0;
  int mBandWidth = 0;
  std::vector<double> mA;
  std::vector<double> mB;

#ifndef GPUCA_ALIROOT_LIB
//  ClassDefNV(BandMatrixSolverNew, 0);
#endif
};


inline void BandMatrixSolverNew::triangulateBlock(double AA[], double bb[], int nRows)
{
  {
    int m = mBandWidth;
    double* A = AA;
    for (int rows = 0; rows < nRows; rows++) {
      double c = 1. / A[0];
      //std::cout<<A[0]<<std::endl;
      if (!std::isfinite(c) || fabs(A[0])<1.e-15 ) {
        isFailed = true;
        //std::cout << "\n\nBandMatrixSover: Can not construct spline! \n\n" << std::endl;
        c = 0.;
      }
      A[0] = c; // store 1/a[0][0]
      double* rowi = A + mBandWidth - 1;
      for (int i = 1; i < m; i++) { // row 0+i
        double ai = c * A[i];       // A[0][i] die neuen Werte 
        for (int j = i; j < m; j++) {
          rowi[j] -= ai * A[j]; // A[i][j] -= A[0][j]/A[0][0]*A[i][0]
        }
        A[i] = ai; // A[0][i] /= A[0][0]
        rowi += mBandWidth - 1;
      }
      m--;
      A += mBandWidth;
    }
  }

  for (int k = 0; k < mBdim; k++) {
    int m = mBandWidth;
    double* A = AA;
    double* b = bb;
    for (int rows = 0; rows < nRows; rows++) {
      double bk = b[k];
      for (int i = 1; i < m; i++) {
        b[mBdim * i + k] -= A[i] * bk;
      }
      b[k] *= A[0];
      m--;
      A += mBandWidth;
      b += mBdim;
    }
  }
}

inline void BandMatrixSolverNew::dioganalizeBlock(double AA[], double bb[], int nCols)
{
  for (int k = 0; k < mBdim; k++) {
    int rows = mBandWidth;
    double* A = AA;
    double* b = bb;
    for (int col = 0; col < nCols; col++) {
      double bk = b[k];
      for (int i = 1; i < rows; i++) {
        b[-i * mBdim + k] -= A[mBandWidth * (-i) + i] * bk;
      }
      A -= mBandWidth;
      b -= mBdim;
      rows--;
    }
  }
}

inline void BandMatrixSolverNew::solve()
{
  /// Solution slover

  isFailed = false;

  const int stepA = mBandWidth;
  const int stepB = mBdim;
  // Upper Triangulization
  {
    int k = 0;
    double* Ak = &mA[0];
    double* bk = &mB[0];
    for (; k < mN - mBandWidth; k += 1, Ak += stepA, bk += stepB) { // for each row k
      triangulateBlock(Ak, bk, 1);
    }
    // last m rows
    triangulateBlock(Ak, bk, mBandWidth);
  }

  // Diagonalization
  {
    int k = mN - 1;
    double* Ak = &mA[mBandWidth * k];
    double* bk = &mB[mBdim * k];
    for (; k > mBandWidth - 1; k -= 1, Ak -= stepA, bk -= stepB) { // for each row k
      dioganalizeBlock(Ak, bk, 1);
    }
    // first m rows
    dioganalizeBlock(Ak, bk, mBandWidth);
  }
}


} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
