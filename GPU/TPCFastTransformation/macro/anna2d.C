/*
   root -l -q anna.C
 */

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "TFile.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TH1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TLine.h"
#include "GPU/Spline2D.h"
#include "Math/Functor.h"
#include "Math/ChebyshevApprox.h"
#include "GPU/Spline2DHelper.h"
#endif

#include <thread>
#include <mutex>

constexpr int nThreads = 1;

constexpr int nKnotsU = 4;
constexpr int nKnotsV = 4;
constexpr int nKnots = nKnotsU * nKnotsV;

constexpr int nCellsU = nKnotsU - 1;
constexpr int nCellsV = nKnotsV - 1;
constexpr int nCells = nCellsU * nCellsV;

constexpr int nPointsMax = 16;
constexpr int modulo = nPointsMax + 1;

double F(double u)
{
  return 0;
}

int printMode = 1;

bool ask()
{
  std::cout << "type: 'q'-exit";
  std::cout << ", 'a'- draw all cases";
  std::cout << ", 'c'- contraversal cases";

  std::cout << std::endl;

  std::string str;
  std::getline(std::cin, str);
  if (str == "a") {
    printMode = 2;
  } else if (str == "c") {
    printMode = 1;
  } else if (str == "q" || str == ".q") {
    printMode = 0;
  }
  return (printMode != 0);
}

int anna2d()
{

  using namespace o2::gpu;

  std::cout << "Test interpolation.." << std::endl;

  long double combTotal = 1.;
  for (int i = 0; i < nCells; i++) {
    combTotal *= modulo;
  }

  int cellCut[nCellsU][nCellsV];

  for (int i = 0, iexp = 1; i < nCellsU; i++) {
    for (int j = 0; j < nCellsV; j++, iexp *= modulo) {
      cellCut[i][j] = iexp;
    }
  }

  long unsigned int combStart = (long unsigned int)0.25 * combTotal;
  std::mutex mtx;

  auto myThread = [&](int iThread) {

    vector<double> pu, pv, pF(10000, 0.);

    // spline with many data points
    o2::gpu::Spline2D<float, 1> spline(nKnotsU, nKnotsV);
    o2::gpu::Spline2DHelper<float> helper;

    // loop over the combination index comb

    for (long unsigned int comb = combStart + iThread; printMode != 0; comb += nThreads) {

      if (iThread == 0 && (comb - combStart) % 1000000 * nThreads == 0) {
        std::cout << "processed " << (100. * comb) / combTotal << " \% of " << combTotal << std::endl;
      }

      // construct the function
      //gRandom->SetSeed(seed);
      //std::cout << "Random seed: " << seed << " " << gRandom->GetSeed() << std::endl;

      // create data points

      int nPointsTotal = 0;

      for (int i = 0; i < nCellsU; i++) {
        for (int j = 0; j < nCellsV; j++) {
          nPointsTotal += (comb / cellCut[i][j]) % modulo;
        }
      }

      if (comb > 0 && nPointsTotal == 0) {
        break;
      }

      if (nPointsTotal < nKnots * 4) {
        continue;
      }

      pu.clear();
      pv.clear();

      for (int iCell = 0; iCell < nCellsU; iCell++) {
        for (int jCell = 0; jCell < nCellsV; jCell++) {

          int nPoints = (comb / cellCut[iCell][jCell]) % modulo;
          int nPoints1D = (int)ceil(sqrt(nPoints));

          assert(nPoints <= nPoints1D * nPoints1D);

          for (int ip = 0, k = 0; ip < nPoints1D && k < nPoints; ip++) {
            for (int jp = 0; jp < nPoints1D && k < nPoints; jp++, k++) {
              double u = iCell + (ip + 1.) / (nPoints1D + 1.);
              double v = jCell + (jp + 1.) / (nPoints1D + 1.);
              pu.push_back(u);
              pv.push_back(v);
            }
          }
        }
      }

      assert(nPointsTotal == (int)pu.size());

      bool constructionOk = helper.approximateDataPoints(spline, 0., nKnotsU - 1, 0., nKnotsV - 1, &pu[0], &pv[0], &pF[0], pu.size());

      //----------------------------------------------anna test ---------------------------------------------

      int testAnna = 2; // 0 : geht nicht, 1: geht, 2: unbekannt

      {
        int nDataPerCell[nCellsU][nCellsV] = {0}; //bei drei plus 1 anstatt plus 2? Testprogramm schreiben, verschiedene Szenarien durchrechnen
        int nDataPerKnot[nKnotsU][nKnotsV] = {0};             // <---

        for (int i = 0; i < pu.size(); i++) {                      //wir gehen die Punkte durch
          int iu = spline.getGridX1().getLeftKnotIndexForU(pu[i]); //is the data point in this Cell? like this, a border point is only counted for the first Cell
          int iv = spline.getGridX2().getLeftKnotIndexForU(pv[i]); //is the data point in this Cell? like this, a border point is only counted for the first Cell
          nDataPerCell[iu][iv]++;
        }

        int cornerCells[4][2] = {{0, 0}, {nCellsU - 1, 0}, {0, nCellsV - 1}, {nCellsU - 1, nCellsV - 1}};
        int cornerKnoten[4][2] = {{0, 0}, {nKnotsU - 1, 0}, {0, nKnotsV - 1}, {nKnotsU - 1, nKnotsV - 1}};

        for (int i = 0; i < 4; i++) {
          if (nDataPerCell[cornerCells[i][0]][cornerCells[i][1]] < 4) {
            testAnna = 0;
            cout << "the corners do not have enough knots" << endl; //Problem with threads?
            break;
          }
          nDataPerKnot[cornerKnoten[i][0]][cornerKnoten[i][1]] = 4;
          nDataPerCell[cornerCells[i][0]][cornerCells[i][1]] -= 4;
        }
        //damit keine wertvollen data points verschwendet werden, wenn nebenan eine Cell >16 data points hat
        for (int i=0; i<nCellsU; i++){ // -> aber die Corner haben schon 4 weniger
          for (int j=0; j<nCellsV; j++){
            if (nDataPerCell[i][j] >= 16) { //umliegende Knoten füllen
              nDataPerKnot[i][j] = 4;
              nDataPerKnot[i][j+1] = 4;
              nDataPerKnot[i+1][j] = 4;
              nDataPerKnot[i+1][j+1] = 4;
              }                                //nDataPerCell ändern? 
          }
        }
        for (int i = 0; i < nCellsU; i++){
          for (int j = 0; j < nCellsV; j++){
            if (nDataPerCell<4){
              n
            }
            
          }
          
        }
        



        /*
        for (int i = 0; i < 2 * nKnotsU + 2 * nKnotsV; i++) { //den Rand abgehen
          int s0 = 4 - nDataPerKnot[i];
          int n0 = std::min(s0, nDataPerCell[i]);
          int n1 = std::min(4, nDataPerCell[i] - n0);
          nDataPerKnot[i] += n0;
          nDataPerKnot[i + 1] = n1;
        }

        for (int i = 0; i < nKnots; i++) {
          if (nDataPerKnot[i] < 2) {
            testAnna = 0;
            cout << "Knot " << i << " has only " << nDataPerKnot[i] << " data points" << endl;
            //vector mit tuple Problemstelle und Anzahl Punkte?
          }
        }
        */
      }

      //-----------------------------------------------end of anna's test -------------------------------------------------

      if (printMode != 1 || testAnna == constructionOk || testAnna == 2) {
        continue;
      }

      // print

      mtx.lock();

      if (printMode == 0) {
        mtx.unlock();
        break;
      }

      std::cout << "thread " << iThread << " combination " << comb << " n points " << nPointsTotal << std::endl;

      for (int i = 0, i6 = 1; i < nCells; i++, i6 *= modulo) {
        int iCell = i / nCellsU;
        int jCell = i % nCellsU;
        int nPoints = (comb / i6) % modulo;
        std::cout << nPoints << " ";
        if (jCell == nCellsU - 1) {
          std::cout << std::endl;
        }
      }
      std::cout << endl;

      helper.approximateDataPoints(spline, 0, nKnotsU - 1, 0, nKnotsV - 1, &pu[0], &pv[0], &pF[0], pu.size(), true);

      std::cout << " matrix solver: " << constructionOk << ", anna's test: " << testAnna << std::endl;

      ask();

      mtx.unlock();

    } // loop over combinations
  };  // end of myThread()

  std::thread threads[nThreads];

  for (int i = 0; i < nThreads; i++) {
    threads[i] = std::thread(myThread, i);
  }

  for (int i = 0; i < nThreads; i++) {
    threads[i].join();
  }

  return 0;
}
