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
#include "BandMatrixSolver.h"
#include "BandMatrixSolverNew.h"
#include "GPUCommonLogger.h"

#include <iostream>
#include <random>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace GPUCA_NAMESPACE::gpu;

ClassImp(GPUCA_NAMESPACE::gpu::BandMatrix2dSolver);




//Annas Code ohne Division

/*
void BandMatrix2dSolver::solve(){

  int fullBandSize = 8+mBandShift+12; // size of the band, including the diagonal element
  // Upper Triangulization
  for (int i = 0; i < mN; i++) { //wir gehen die Spalten durch
    double* rowI = &mA[i * mShift]; //pointer to the first A element in row i
    double* rowIb = &mA[i * mShift + mN]; //pointer to the first B element in row i
    double aii = rowI[i]; // diagonal element
    std::cout << "i aii " << i << " " << aii << std::endl;
    double* rowJ = rowI + mShift; //row beneath row i

    for (int j = i+1; j<mN; j++, rowJ += mShift) { // wir lassen die Zeile, die Null werden soll, aus NEU: keine Bänder mehr (j<i+8) (statt nur j<mN, nKnots muss noch berücksichtigt werden..) i+1 and rowJ goes one row down, die Reihen werden durchgegangen
      double aji = rowJ[i]; // A[i][j] / A[i][i] , der neue Wert von rowI[j]
      for (int k = 0; k<mShift; k++) { //NEU: kein k<mShift sondern auch nur bis zum Ende des letzten Bandes -> wir wissen auch hier, wo die Nullen sind        
        rowJ[k] = rowJ[k] - (aji/aii) * rowI[k];
        //rowJ[k] = rowJ[k]  - aij *rowI[k] / aii; // A[j][k] -= A[i][k]/A[i][i]*A[j][i] //muss das nicht rowJ[i] sein? oder egal weil symmetrisch? von den Elementen unter row i werden die gewichteten Werte von i abgezogen, um vorne eine 0 zu bekommen 
        //rowJ[k]*=aii;
      }      
    }  
    for (int k = i; k<mShift; k++) {
      //rowI[k] /=aii; 
    }
  }

  // Diagonalization
  for (int i = mN - 1; i >= 0; i--) { 

    double* rowIb = &mA[i * mShift + mN]; //pointer to the first B element in row i
    double* rowJb = rowIb - mShift; //pointer to the first B element one row above row i

    double aii = mA[i * mShift + i]; 
    for (int k = 0; k < mM; k++) {
      rowIb[k] = rowIb[k]/aii;
    } 
    mA[i * mShift + i] = 1.; 

    //Band 1
    for (int j = i - 1; j >= 0; j--, rowJb -= mShift) { // NEU: nicht mehr j>=i-8  row j
      double aji = mA[j * mShift + i]; //die Schleife geht ab dem Diagonalelement und die Spalte darüber ab
      for (int k = 0; k < mM; k++) {
        rowJb[k] = rowJb[k] - aji * rowIb[k];
      }
    }

  }
}
*/



///*
//Annas Code ohne Dreieck


void BandMatrix2dSolver::solve(){
  isFailed = false;
  int fullBandSize = 8+mBandShift+12; // size of the band, including the diagonal elementfabs
  // Upper Triangulization
  for (int i = 0; i < mN; i++) { //wir gehen die Spalten durch
    double* rowI = &mA[i * mShift]; //pointer to the first A element in row i
    double* rowIb = &mA[i * mShift + mN]; //pointer to the first B element in row i

    double c = 1. / rowI[i]; // 1 / diagonal element
    //std::cout<<A[0]<<std::endl;
    if (!std::isfinite(c) || !std::isfinite(rowI[i]) || fabs(rowI[i]) < 1.e-15 ) {
      isFailed = true;
      //std::cout << "\n\nBandMatrixSover: Can not construct spline! \n\n" << std::endl;
      c = 0.;
    }

    double* rowJ = rowI + mShift; //row beneath row i

    for (int j = i + 1; j<i+fullBandSize && j<mN; j++, rowJ += mShift) { // NEU: keine Bänder mehr (j<i+8) (statt nur j<mN, nKnots muss noch berücksichtigt werden..) i+1 and rowJ goes one row down, die Reihen werden durchgegangen
     // if (rowI[j] != 0.) {
      double aij = c * rowI[j]; // A[i][j] / A[i][i] , der neue Wert von rowI[j]
      for (int k = j; k<j+fullBandSize && k<mN; k++) { //NEU: kein k<mShift sondern auch nur bis zum Ende des letzten Bandes -> wir wissen auch hier, wo die Nullen sind
        rowJ[k] -= aij * rowI[k]; // A[j][k] -= A[i][k]/A[i][i]*A[j][i] //muss das nicht rowJ[i] sein? oder egal weil symmetrisch? von den Elementen unter row i werden die gewichteten Werte von i abgezogen, um vorne eine 0 zu bekommen 
      }
      for (int k = mN; k<mShift; k++) { //NEU: kein k<mShift sondern auch nur bis zum Ende des letzten Bandes -> wir wissen auch hier, wo die Nullen sind
        rowJ[k] -= aij * rowI[k]; // A[j][k] -= A[i][k]/A[i][i]*A[j][i] //muss das nicht rowJ[i] sein? oder egal weil symmetrisch? von den Elementen unter row i werden die gewichteten Werte von i abgezogen, um vorne eine 0 zu bekommen 
      }
      rowI[j] = aij; // A[i][j] /= A[i][i] //der neue Wert für die Zahl in Zeile i, bei der ihr symmetrischer Wert gerade null gemacht wurde
     // }
    }
    for (int k = 0; k < mM; k++) { //die neuen b Werte mittels Teilen durch Diagonalwert
      rowIb[k] *= c;
    }
  }

  // Diagonalization
  for (int i = mN - 1; i >= 0; i--) {
    double* rowIb = &mA[i * mShift + mN]; //pointer to the first B element in row i
    double* rowJb = rowIb - mShift; //pointer to the first B element one row above row i
    //Band 1
    for (int j = i - 1; j >= 0 && j >i-fullBandSize ; j--, rowJb -= mShift) { // NEU: nicht mehr j>=i-8  row j
      double aji = mA[j * mShift + i]; //die Schleife geht ab dem Diagonalelement und die Spalte darüber ab
      //if (aji != 0.) {
        for (int k = 0; k < mM; k++) {
          rowJb[k] -= aji * rowIb[k]; 
        }
      //}
    }
  }
}

//*/


//Annas Code mit Dreieck
/*

void BandMatrix2dSolver::solve(){

  // Upper Triangulization
  for (int i = 0; i < mN; i++) { //wir gehen die Spalten durch
    double* rowI = &mA[i * mShift]; //pointer to the first A element in row i
    double* rowIb = &mA[i * mShift + mN]; //pointer to the first B element in row i
    double c = 1. / rowI[i]; // 1 / diagonal element
    double* rowJ = rowI + mShift; //row beneath row i

    //Band 1
    int triangle = 0; //Lücke zwischen Bändern wird kleiner bis 0
    if (mBandShift > i){
      triangle = mBandShift - i; //could also do steps of 4 
    }
    

    for (int j = i + 1; j<i+8 && j<mN; j++, rowJ += mShift) { // NEU: keine Bänder mehr (j<i+8) (statt nur j<mN, nKnots muss noch berücksichtigt werden..) i+1 and rowJ goes one row down, die Reihen werden durchgegangen
     // if (rowI[j] != 0.) {
      double aij = c * rowI[j]; // A[i][j] / A[i][i] , der neue Wert von rowI[j]
      for (int k = j; k<0.5*mShift+8+mBandShift+12 && k<mShift; k++) { //NEU: kein k<mShift sondern auch nur bis zum Ende des letzten Bandes -> wir wissen auch hier, wo die Nullen sind
        rowJ[k] -= aij * rowI[k]; // A[j][k] -= A[i][k]/A[i][i]*A[j][i] //muss das nicht rowJ[i] sein? oder egal weil symmetrisch? von den Elementen unter row i werden die gewichteten Werte von i abgezogen, um vorne eine 0 zu bekommen 
      }
      rowI[j] = aij; // A[i][j] /= A[i][i] //der neue Wert für die Zahl in Zeile i, bei der ihr symmetrischer Wert gerade null gemacht wurde
     // }
    }

    // Band 2
    if (i+8>=mN) { //?where is the mistake
    //rowJ = rowI + (triangle+1)*mShift; //darunter i + 8 + triangle?
    for (int j = i + 8 ; j<i+8+mBandShift+12 && j<mN; j++, rowJ+=mShift) { //NEU: keine Lücke mehr, nicht j = i+mBandshift-4 (statt nur j<mN j<4*7, nKnots muss noch berücksichtigt werden..) i+1 and rowJ goes one row down, die Reihen werden durchgegangen      
      // rowJ = rowI + (j-1)*mShift; // da ist was falsch -----> geht nicht zu Ende
     // if (rowI[j] != 0.) {   //if weggenommen , wenn das Element rechts neben i nicht null ist (symm Matrix) -> nicht mehr so oft auf nullen püfen, wir wissen wann es nicht null ist
      double aij = c * rowI[j]; // A[i][j] / A[i][i] , der neue Wert von rowI[j]
      for (int k = j; k<0.5*mShift+8+mBandShift+12 && k<mShift; k++) { //NEU: kein k<mShift sondern auch nur bis zum Ende des letzten Bandes -> wir wissen auch hier, wo die Nullen sind
        rowJ[k] -= aij * rowI[k]; // A[j][k] -= A[i][k]/A[i][i]*A[j][i] //muss das nicht rowJ[i] sein? oder egal weil symmetrisch? von den Elementen unter row i werden die gewichteten Werte von i abgezogen, um vorne eine 0 zu bekommen 
      }
      rowI[j] = aij; // A[i][j] /= A[i][i] //der neue Wert für die Zahl in Zeile i, bei der ihr symmetrischer Wert gerade null gemacht wurde
     // }
    }   

    } // von if 105 
 
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
// original solver in easy language without division
void BandMatrix2dSolver::solve()
{
  // Upper Triangulization
  for (int i = 0; i < mN; i++) {
    double aii =  A(i,i);
    for (int j = i + 1; j < mN; j++) { // row j
      double aij =  A(i,j); // A[i][j] / A[i][i]
      for (int k = j; k < mShift; k++) {
       // A(j,k) = A(j,k)*aii - aij * A(i,k); // *= aii - aij * A(i,k); // A[j][k] -= A[i][k]/A[i][i]*A[j][i] 
      double ajk = A(j,k);
      double aik = A(i,k);
        A(j,k) = ajk - aik * aij /aii ;
        //A(j,k)*= aii;
      }
      //A(i,j) = aij; // A[i][j] /= A[i][i]
    }    
    //for (int k = 0; k < mM; k++) {
      //B(i,k) *= c;
   // }    
  }
  
  // Diagonalization
  for (int i = mN - 1; i >= 0; i--) {
    double aii = A(i,i);
    
    
    for (int j = i - 1; j >= 0; j--) { // row j
      double aji = A(j,i);
      for (int k = 0; k < mM; k++) {
        B(j,k) = B(j,k)*aii - aji* B(i,k); //B(j,k) -= aji * B(i,k) / aii ?
      }
    }
  }
  for (int i = mN-1; i>= 0; i--){
    double aii = A(i,i);
    for (int k = 0; k < mM; k++) {
      B(i,k) /= aii;
    }  
  }
}
*/


  /*
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
        if (A(i,j)!=0 && Matrix[i][j]==0) Matrix[i][j] = 1; //Test für Ausgabe
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
    if (A(i,j)!=0 && Matrix[i][j]==0) Matrix[i][j] = 1; //Test für Ausgabe
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
 */


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
  constexpr int n = 160; //vorher 30 bzw 160 40
  constexpr int d = 3;
  constexpr int bandShift = 28; //vorher 8 bzw 28 8 

  // std::random_device rd;  // Will be used to obtain a seed for the random
  std::mt19937 gen(1); // Standard mersenne_twister_engine seeded with 1
  std::uniform_real_distribution<> uniform(-.999, .999);

  double maxDiff = -1.;
  double maxDiff1d = 0.;
  double maxDiff1dNew = 0.;
  int nTries = 10000;

  auto tmpTime = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(tmpTime - tmpTime); //to get ns
  auto durationMult = duration;
  auto duration1d = duration;
  auto duration1dNew = duration;

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
          A[j][i] = A[i][j];  //symm
        }
      }

      for (int i = 0; i < n; i++) {
        A[i][i] = A[i][i] * A[i][i];
      }
      if (prn && iter == nTries - 1) {  //beim letzten Durchgang die Matrix ausgeben
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
    BandMatrixSolver<8+bandShift+12> sym1d(n,d);
    BandMatrixSolverNew sym1dNew(n,d,8+bandShift+12);

    for (int i = 0; i < n; i++) {
      for (int k = 0; k < d; k++) {
        sym.B(i, k) = b[i][k];
        sym1dNew.B(i,k) = b[i][k];
        sym1d.B(i,k) = b[i][k];
      }
      for (int j = i; j < n; j++) {
        sym.A(i, j) = A[i][j];
        sym.A(j, i) = A[i][j];
        sym1d.A(i, j) = A[i][j];
        if(j<i+8+bandShift+12){  //-------wieso nochmal? ist doch schon ne Bandmatrix
          sym1dNew.A(i,j)= A[i][j];
        }
      }
    }

    auto start = std::chrono::high_resolution_clock::now(); //Zeitstoppen 
    sym.solve();
    auto stop = std::chrono::high_resolution_clock::now();
    duration += std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

    auto start1d = std::chrono::high_resolution_clock::now();
    sym1d.solve();
    auto stop1d = std::chrono::high_resolution_clock::now();
    duration1d += std::chrono::duration_cast<std::chrono::nanoseconds>(stop1d - start1d);

    auto start1dNew = std::chrono::high_resolution_clock::now();
    sym1dNew.solve();
    auto stop1dNew = std::chrono::high_resolution_clock::now();
    duration1dNew += std::chrono::duration_cast<std::chrono::nanoseconds>(stop1dNew - start1dNew);


    double diff = -1.;
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < d; k++) {
        double t = fabs(x[i][k] - sym.B(i, k));
        if (diff < t || !finite(t)) {
          diff = t;
        }
      }
    }
    if (maxDiff < diff || !finite(diff)) {
      maxDiff = diff;
    }

    double diff1d = -1.;
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < d; k++) {
        double t = fabs(x[i][k] - sym1d.B(i, k));
        if (diff1d < t || !finite(t)) {
          diff1d = t;
        }
      }
    }
    if (maxDiff1d < diff1d || !finite(diff1d)) {
      maxDiff1d = diff1d;
    }


    double diff1dNew = 0.;
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < d; k++) {
        double t = fabs(x[i][k] - sym1dNew.B(i, k));
        if (diff1dNew < t || !finite(t)) {
          diff1dNew = t;
        }
      }
    }

    if (maxDiff1dNew < diff1dNew || !finite(diff1dNew)) {
      maxDiff1dNew = diff1dNew;
    }

    // LOG(info) << std::defaultfloat ;
    // LOG(info) << "\n\n max diff " <<diff << "\n";
  }

  int ok = (maxDiff < 1.e-7);
  int ok1d = (maxDiff1d < 1.e-7);

  if (prn || !ok || !ok1d) {
    LOG(info) << std::defaultfloat;
    LOG(info) << "\n\n Overall max diff " << maxDiff << "\n";
    LOG(info) << "\n\n Overall max diff 1d " << maxDiff1d << "\n";
    LOG(info) << "\n\n Overall max diff 1d New " << maxDiff1dNew << "\n";
    LOG(info) << " time " << duration.count() / nTries;
   LOG(info) << " time 1d " << duration1d.count() / nTries;
   LOG(info) << " time 1d New " << duration1dNew.count() / nTries;
    LOG(info) << " time multiplication " << durationMult.count() / nTries;
  }
  return ok;
}
