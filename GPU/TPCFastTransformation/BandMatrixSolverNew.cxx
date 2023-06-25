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

/// \file  BandMatrixSolverNew.cxx
/// \brief Implementation of BandMatrixSolverNew class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#include "BandMatrixSolverNew.h"
#include "GPUCommonLogger.h"
#include <random>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace GPUCA_NAMESPACE::gpu;

//templateClassImp(GPUCA_NAMESPACE::gpu::BandMatrixSolverNew);

