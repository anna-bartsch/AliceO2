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

/// \file  BandMatrix2dSolver.cxx
/// \brief Implementation of BandMatrix2dSolver class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>

#include "BandMatrix2dSolver.h"
#include "GPUCommonLogger.h"

#include <iostream>
#include <random>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace GPUCA_NAMESPACE::gpu;

ClassImp(GPUCA_NAMESPACE::gpu::BandMatrix2dSolver);

/*
//Annas Code

void BandMatrix2dSolver::solve(){

  // Upper Triangulization
  for (int i = 0; i < mN; i++) { //wir gehen die Spalten durch
    double* rowI = &mA[i * mShift]; //pointer to the first A element in row i
    double* rowIb = &mA[i * mShift + mN]; //pointer to the first B element in row i
    double c = 1. / rowI[i]; // 1/das Diagonalelement
    double* rowJ = rowI + mShift; //row beneath row i
    //Band 1
    for (int j = i + 1; j<i+8+mBandShift+12 && j<mN; j++, rowJ += mShift) { // NEU: keine Bänder mehr (j<i+8) (statt nur j<mN, nKnots muss noch berücksichtigt werden..) i+1 and rowJ goes one row down, die Reihen werden durchgegangen
     // if (rowI[j] != 0.) {
      double aij = c * rowI[j]; // A[i][j] / A[i][i] , der neue Wert von rowI[j]
      for (int k = j; k<0.5*mShift+8+mBandShift+12 && k<mShift; k++) { //NEU: kein k<mShift sondern auch nur bis zum Ende des letzten Bandes -> wir wissen auch hier, wo die Nullen sind
        rowJ[k] -= aij * rowI[k]; // A[j][k] -= A[i][k]/A[i][i]*A[j][i] //muss das nicht rowJ[i] sein? oder egal weil symmetrisch? von den Elementen unter row i werden die gewichteten Werte von i abgezogen, um vorne eine 0 zu bekommen 
      }
      rowI[j] = aij; // A[i][j] /= A[i][i] //der neue Wert für die Zahl in Zeile i, bei der ihr symmetrischer Wert gerade null gemacht wurde
     // }
    }
      //Band 2
    //*rowJ = rowI + (mBandShift-4)*mShift;
    //*for (int j = i + mBandShift-4 ; j<i+mBandShift+12 && j<mN; j++, rowJ+=mShift) { //NEU: keine Lücke mehr, nicht j = i+mBandshift-4 (statt nur j<mN j<4*7, nKnots muss noch berücksichtigt werden..) i+1 and rowJ goes one row down, die Reihen werden durchgegangen      
      // rowJ = rowI + (j-1)*mShift; // da ist was falsch -----> geht nicht zu Ende
     // if (rowI[j] != 0.) {
      //if weggenommen , wenn das Element rechts neben i nicht null ist (symm Matrix) -> nicht mehr so oft auf nullen püfen, wir wissen wann es nicht null ist
     //* double aij = c * rowI[j]; // A[i][j] / A[i][i] , der neue Wert von rowI[j]
     //* for (int k = j; k<mShift; k++) { //-> wir wissen auch hier, wo die Nullen sind
     //*   rowJ[k] -= aij * rowI[k]; // A[j][k] -= A[i][k]/A[i][i]*A[j][i] //muss das nicht rowJ[i] sein? oder egal weil symmetrisch? von den Elementen unter row i werden die gewichteten Werte von i abgezogen, um vorne eine 0 zu bekommen 
     //* }
     //* rowI[j] = aij; // A[i][j] /= A[i][i] //der neue Wert für die Zahl in Zeile i, bei der ihr symmetrischer Wert gerade null gemacht wurde
     // }
    //*}   

    for (int k = 0; k < mM; k++) { //die neuen b Werte mittels Teilen durch Diagonalwert
      rowIb[k] *= c;
    }
  } 

  // Diagonalization
  for (int i = mN - 1; i >= 0; i--) {
    double* rowIb = &mA[i * mShift + mN]; //pointer to the first B element in row i
    double* rowJb = rowIb - mShift; //pointer to the first B element one row above row i
    //Band 1
    for (int j = i - 1; j >= 0 && j >i-8-mBandShift-12 ; j--, rowJb -= mShift) { // NEU: nicht mehr j>=i-8  row j
      double aji = mA[j * mShift + i]; //die Schleife geht ab dem Diagonalelement und die Spalte darüber ab
      //if (aji != 0.) {
        for (int k = 0; k < mM; k++) {
          rowJb[k] -= aji * rowIb[k];
        }
      //}
    }
    //* //Band 2
    //*rowJb = rowIb - (mBandShift-4)*mShift;
    //*for (int j = i - mBandShift +4; j >= 0 && j > i-mBandShift-12 ; j--, rowJb -= mShift) { // row j
   //*   double aji = mA[j * mShift + i]; //die Schleife geht ab dem Diagonalelement und die Spalte darüber ab
      //if (aji != 0.) {
     //*   for (int k = 0; k < mM; k++) {
     //*     rowJb[k] -= aji * rowIb[k];
    //*    }
      //}
   //* }

  }
}

*/


   




/*
// original solver in easy language
void BandMatrix2dSolver::solve()
{
  // Upper Triangulization
  for (int i = 0; i < mN; i++) {
    double c = 1. / A(i,i);
    for (int j = i + 1; j < mN; j++) { // row j
      double aij = c * A(i,j); // A[i][j] / A[i][i]
      for (int k = j; k < mShift; k++) {
        A(j,k) -= aij * A(i,k); // A[j][k] -= A[i][k]/A[i][i]*A[j][i]
      }
      A(i,j) = aij; // A[i][j] /= A[i][i]
    }    
    for (int k = 0; k < mM; k++) {
      B(i,k) *= c;
    }    
  }
  
  // Diagonalization
  for (int i = mN - 1; i >= 0; i--) {
    for (int j = i - 1; j >= 0; j--) { // row j
      double aji = A(j,i);
      for (int k = 0; k < mM; k++) {
        B(j,k) -= aji * B(i,k);
      }
    }
  }
}
*/


  
 // original solver
void BandMatrix2dSolver::solve()
{
  //printing the matrix with 0 and touched entries

  double Matrix[mN][mN];
  for (int i = 0; i < mN; i++) {
    for (int j = 0; j < mN; j++) {
      Matrix[i][j]=0;
    }
    Matrix[i][i]= 2;
  }

  // Upper Triangulization
  for (int i = 0; i < mN; i++) {
    double* rowI = &mA[i * mShift];
    double* rowIb = &mA[i * mShift + mN];
    double c = 1. / rowI[i];
    double* rowJ = rowI + mShift;
    for (int j = i + 1; j < mN; j++, rowJ += mShift) { // row j
      //if (rowI[j] != 0.) {
        if (A(i,j)!=0 && Matrix[i][j]==0) Matrix[i][j] = j; //Test für Ausgabe
        double aij = c * rowI[j]; // A[i][j] / A[i][i]
        for (int k = j; k < mShift; k++) {
          rowJ[k] -= aij * rowI[k]; // A[j][k] -= A[i][k]/A[i][i]*A[j][i]
        }
        rowI[j] = aij; // A[i][j] /= A[i][i]
      //}
    }    
    for (int k = 0; k < mM; k++) {
      rowIb[k] *= c;
    }    
  }
  
  // Diagonalization
  for (int i = mN - 1; i >= 0; i--) {
    double* rowIb = &mA[i * mShift + mN];
    double* rowJb = rowIb - mShift;
    for (int j = i - 1; j >= 0; j--, rowJb -= mShift) { // row j
    if (A(i,j)!=0 && Matrix[i][j]==0) Matrix[i][j] = j; //Test für Ausgabe
      double aji = mA[j * mShift + i];
      //if (aji != 0.) {
        for (int k = 0; k < mM; k++) {
          rowJb[k] -= aji * rowIb[k];
        }
      //}
    }
  }
  for (int i = 0; i < mN; i++) {
    for (int j = 0; j < mN; j++) {
      cout << Matrix[i][j] << " ";
    }
    cout << endl;
  }
  
}
 


void BandMatrix2dSolver::print()
{
  for (int i = 0; i < mN; i++) {
    LOG(info) << "";
    for (int j = 0; j < mN; j++) {
      LOG(info) << std::fixed << std::setw(5) << std::setprecision(2) << A(i, j) << " ";
    }
    LOG(info) << " | ";
    for (int j = 0; j < mM; j++) {
      LOG(info) << std::fixed << std::setw(5) << std::setprecision(2) << B(i, j) << " ";
    }
  }
  LOG(info) << std::setprecision(-1);
}

int BandMatrix2dSolver::test(bool prn)
{
  constexpr int n = 40; //vorher 30
  constexpr int d = 3;
  constexpr int bandShift = 8;

  // std::random_device rd;  // Will be used to obtain a seed for the random
  std::mt19937 gen(1); // Standard mersenne_twister_engine seeded with 1
  std::uniform_real_distribution<> uniform(-.999, .999);

  double maxDiff = 0.;
  int nTries = 10000;

  auto tmpTime = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(tmpTime - tmpTime);
  auto durationMult = duration;

  for (int iter = 0; iter < nTries; iter++) {

    double x[n][d];
    double A[n][n];
    {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < d; j++) {
          x[i][j] = 1. * uniform(gen);
        }
      }
      for (int i = 0; i < n; i++) {
        A[i][i] = fabs(2. + uniform(gen));
      }
      for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
          A[i][j] = A[i][i] * A[j][j] * uniform(gen);
          if (j >= i + 8 && j < i+8+bandShift) A[i][j] = 0; //null setzen ---------+8 bis bandShift und bandShift +12-----------------------------------------
          if (j >= i+bandShift + 12 +8) A[i][j] = 0;
          A[j][i] = A[i][j];
        }
      }

      for (int i = 0; i < n; i++) {
        A[i][i] = A[i][i] * A[i][i];
      }
      if (prn && iter == nTries - 1) {
        for (int i = 0; i < n; i++) {
          LOG(info) << "";
          for (int j = 0; j < n; j++) {
            std::cout << std::fixed << std::setw(5) << std::setprecision(2) << A[i][j] << " ";
          }
        }
        LOG(info) << "";
      }
    }
    double b[n][d];
    auto startMult = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < d; k++) {
        b[i][k] = 0.;
      }
      for (int j = 0; j < n; j++) {
        for (int k = 0; k < d; k++) {
          b[i][k] += x[j][k] * A[i][j];
        }
      }
    }
    auto stopMult = std::chrono::high_resolution_clock::now();
    durationMult += std::chrono::duration_cast<std::chrono::nanoseconds>(stopMult - startMult);

    BandMatrix2dSolver sym(n, d, bandShift);

    for (int i = 0; i < n; i++) {
      for (int k = 0; k < d; k++) {
        sym.B(i, k) = b[i][k];
      }
      for (int j = i; j < n; j++) {
        sym.A(i, j) = A[i][j];
      }
    }

    auto start = std::chrono::high_resolution_clock::now();
    sym.solve();
    auto stop = std::chrono::high_resolution_clock::now();
    duration += std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    double diff = 0.;
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < d; k++) {
        double t = fabs(x[i][k] - sym.B(i, k));
        if (diff < t) {
          diff = t;
        }
      }
    }
    if (maxDiff < diff) {
      maxDiff = diff;
    }
    // LOG(info) << std::defaultfloat ;
    // LOG(info) << "\n\n max diff " <<diff << "\n";
  }

  int ok = (maxDiff < 1.e-7);

  if (prn || !ok) {
    LOG(info) << std::defaultfloat;
    LOG(info) << "\n\n Overall max diff " << maxDiff << "\n";
    LOG(info) << " time " << duration.count() / nTries;
    LOG(info) << " time multiplication " << durationMult.count() / nTries;
  }
  return ok;
}
