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

/// \file  BandMatrix2dSolver.h
/// \brief Definition of BandMatrix2dSolver class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#ifndef ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_BandMatrix2dSolver_H
#define ALICEO2_GPUCOMMON_TPCFASTTRANSFORMATION_BandMatrix2dSolver_H

#include "GPUCommonDef.h"
#include "GPUCommonRtypes.h"
#include <vector>
#include <cassert>
#include <algorithm>

namespace GPUCA_NAMESPACE
{
namespace gpu
{

/// Linear Equation Solver for a symmetric positive-definite matrix A[n x n].
///
/// A[n x n] * X [n x m] = B[n x m]
///
/// A elements are stored in the upper triangle of A.
/// Thus A(i,j) and A(j,i) access the same element.
///
class BandMatrix2dSolver
{
 public:
  BandMatrix2dSolver(int N, int M, int BS) : mN(N), mM(M), mBandShift(BS), mShift(mN + mM)
  {
    assert(mN > 0 && mM > 0 && mBandShift>0);
    mA.resize(mN * mShift, 0.);
  }

  /// access to A elements
  double& A(int i, int j)
  {
    auto ij = std::minmax(i, j);
    assert(ij.first >= 0 && ij.second < mN);
    return mA[ij.first * mShift + ij.second];
  }

  /// access to B elements
  double& B(int i, int j)
  {
    assert(i >= 0 && i < mN && j >= 0 && j < mM);
    return mA[i * mShift + mN + j];
  }

  ///
  void solve();

  ///
  void print();

  /// Test the class functionality. Returns 1 when ok, 0 when not ok
  static int test(bool prn = 0);

 private:
 private:
  int mN = 0;
  int mM = 0;
  int mBandShift = 0; // a gap - the distance between the bands
  int mShift = 0;
  std::vector<double> mA;

#ifndef GPUCA_ALIROOT_LIB
  ClassDefNV(BandMatrix2dSolver, 0);
#endif
};

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif