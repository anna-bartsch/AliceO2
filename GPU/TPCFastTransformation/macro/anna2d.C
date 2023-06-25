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

constexpr int nThreads = 60;

constexpr int nKnotsU = 4;
constexpr int nKnotsV = 4;
constexpr int nKnots = nKnotsU * nKnotsV;

constexpr int nIntervalsU = nKnotsU - 1;
constexpr int nIntervalsV = nKnotsV - 1;
constexpr int nIntervals = nIntervalsU * nIntervalsV;

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
  for (int i = 0; i < nIntervals; i++) {
    combTotal *= modulo;
  }

  int cellCut[nIntervalsU][nIntervalsV];

  for (int i = 0, iexp = 1; i < nIntervalsU; i++) {
    for (int j = 0; j < nIntervalsV; j++, iexp *= modulo) {
      cellCut[i][j] = iexp;
    }
  }

  long unsigned int combStart = (long unsigned int) 0.25 * combTotal;
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

      for (int i = 0; i < nIntervalsU; i++) {
        for (int j = 0; j < nIntervalsV; j++) {
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

      for (int iCell = 0; iCell < nIntervalsU; iCell++) {
        for (int jCell = 0; jCell < nIntervalsV; jCell++) {

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

      assert(nPointsTotal == (int) pu.size());

      bool constructionOk = helper.approximateDataPoints(spline, 0., nKnotsU - 1, 0., nKnotsV - 1, &pu[0], &pv[0], &pF[0], pu.size());

      //----------------------------------------------anna test ---------------------------------------------
      bool testAnna = 0;
      {
        constructionOk = helper.approximateDataPoints(spline, 0, nKnots - 1, &vu[0], &vy[0], vu.size());

      std::vector<int> nDataPerInterval(nIntervals, 0);            //bei drei plus 1 anstatt plus 2? Testprogramm schreiben, verschiedene Szenarien durchrechnen
      std::vector<int> nDataPerKnot(spline.getNumberOfKnots(), 0); // <---

      for (int i = 0; i < vu.size(); i++) { //wir gehen die Punkte durch
        int interval = spline.getLeftKnotIndexForU(vu[i]); //is the data point in this interval? like this, a border point is only counted for the first interval
        nDataPerInterval[interval]++;
      }
      int corner[] = [0, nKnotsU-1, nKnots-nKnotsU-1, nKnots-1];
      for (int i = 0; i < 4; i++) {
      if (nDataPerInterval[corner[i]] < 4)
      {
        testAnna = 0;
        cout << "the corners do not have enough knots" << endl; //Problem with threads?
      }
      }

      for (int i = 0; i < nIntervals; i++) {
        int s0 = 2 - nDataPerKnot[i];
        int n0 = std::min(s0, nDataPerInterval[i]);
        int n1 = std::min(2, nDataPerInterval[i] - n0);
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


      }

      //-----------------------------------------------end of anna's test -------------------------------------------------

      if (printMode != 1 || testAnna == constructionOk) {
        continue;
      }

      // print

      mtx.lock();

      if (printMode == 0) {
        mtx.unlock();
        break;
      }

      std::cout << "thread " << iThread << " combination " << comb << " n points " << nPointsTotal << std::endl;

      for (int i = 0, i6 = 1; i < nIntervals; i++, i6 *= modulo) {
        int iCell = i / nIntervalsU;
        int jCell = i % nIntervalsU;
        int nPoints = (comb / i6) % modulo;
        std::cout << nPoints << " ";
        if (jCell == nIntervalsU - 1) {
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
