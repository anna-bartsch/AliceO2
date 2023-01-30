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
#include "GPU/Spline1D.h"
#include "GPU/IrregularSpline1D.h"
#include "GPU/Spline1DHelperOld.h"
#include "Math/Functor.h"
#include "Math/ChebyshevApprox.h"

#endif

const int Fdegree = 5;
int nKnots = 6;

static double Fcoeff[2 * (Fdegree + 1)];

double F(double u)
{
  double uu = u * TMath::Pi() / (nKnots - 1);
  double f = 0; // Fcoeff[0]/2;
  for (int i = 1; i <= Fdegree; i++) {
    f += Fcoeff[2 * i] * TMath::Cos(i * uu) + Fcoeff[2 * i + 1] * TMath::Sin(i * uu);
  }
  return f;
}

double Flocal(double u)
{
  double f = F(u * (nKnots - 1));
  return f;
}

TCanvas* canv = new TCanvas("cQA", "Spline Demo", 1600, 800);

int drawOption = 0;

bool ask()
{
  canv->Update();
  std::cout << "type: 'q'-exit";
  std::cout << ", 'a'- draw all cases";
  std::cout << ", 'c'- contraversal cases";

  std::cout << std::endl;

  std::string str;
  std::getline(std::cin, str);
  if (str == "a") {
    drawOption = 0;
  } else if (str == "c") {
    drawOption = 1;
  }

  return (str != "q" && str != ".q");
}

int anna()
{

  using namespace o2::gpu;

  std::cout << "Test interpolation.." << std::endl;

  // TCanvas* canv = new TCanvas("cQA", "Spline1D  QA", 2000, 1000);

  gRandom->SetSeed(0);
  int seed = 12;
  long int comb = 0; // combination index

  for (;; comb++, seed++) {

    // construct the function
    gRandom->SetSeed(seed);
    std::cout << "Random seed: " << seed << " " << gRandom->GetSeed() << std::endl;

    for (int i = 0; i < 2 * (Fdegree + 1); i++) {
      Fcoeff[i] = gRandom->Uniform(-1, 1);
    }

    // create data points

    vector<double> vu, vy;
    int nPointsTotal = 0;

    for (int i = 0, i5 = 1; i < nKnots - 1; i++, i5 *= 5) {
      int nPoints = (comb / i5) % 5; // i5 is 5^i
      nPointsTotal += nPoints;
      std::cout << nPoints << " ";
      for (int ip = 0; ip < nPoints; ip++) {
        double u = i + (ip + 1.) / (nPoints + 1.);
        vu.push_back(u);
        vy.push_back(F(u));
      }
    }
    std::cout << endl;

    if (comb > 0 && nPointsTotal == 0)
      break;

    TNtuple* drawPoints = new TNtuple("drawPoints", "drawPoints", "type:u:f");
    for (int i = 0; i < vu.size(); i++) {
      drawPoints->Fill(1, vu[i], vy[i]);
    }

    // spline with many data points
    o2::gpu::Spline1D<float, 1> spline(nKnots);
    o2::gpu::Spline1DHelper<float> helper;

    bool testAnna = 1;
    bool constructionOk = 1;

    /*{     //---------------------------------------------test anna 1--------------------------

      constructionOk = helper.approximateDataPoints(spline, 0, nKnots - 1, &vu[0], &vy[0], vu.size());


      int nIntervals = spline.getNumberOfKnots() - 1;

      std::vector<int> nDataPerInterval(nIntervals, 0); //bei drei plus 1 anstatt plus 2? Testprogramm schreiben, verschiedene Szenarien durchrechnen

      for (int i = 0; i < vu.size(); i++) {
        int interval = spline.getLeftKnotIndexForU(vu[i]); //is the data point in this interval? like this, a border point is only counted for the first interval
        nDataPerInterval[interval]++;
      }

      bool isChanged = 1;
      while (isChanged) {
        isChanged = 0;
        for (int j = 0; j < nIntervals; j++) {
          if (nDataPerInterval[j] >= 4) {
            if (j > 0) { //if it exists
              nDataPerInterval[j - 1] += 2;
            }
            if (j < nIntervals - 1) { //if it exists
              nDataPerInterval[j + 1] += 2;
            }
            nDataPerInterval[j] = -10; //because we don't want it to increase the neighbours again; a value can be incremented twice per run
            isChanged = 1;
          }
        }
        for (int j = 0; j < nIntervals - 1; j++) {
          if (nDataPerInterval[j] == 3 && nDataPerInterval[j + 1] == 3) {
            nDataPerInterval[j] = 4;
            nDataPerInterval[j + 1] = 4;
            isChanged = 1;
          }
        }
      }

      for (int i = 0; i < nIntervals; i++) {
        if (nDataPerInterval[i] < 4 && nDataPerInterval[i] >= 0) {
          testAnna = 0;
          cout << "Interval " << i << " has only " << nDataPerInterval[i] << " data points" << endl;
          //vector mit tuple Problemstelle und Anzahl Punkte?
        }
      }

      //-------------------------------------end of test 1-------------------------------------------
    } */

    //----------------------------------------------anna test 2---------------------------------------------

    {
      constructionOk = helper.approximateDataPoints(spline, 0, nKnots - 1, &vu[0], &vy[0], vu.size());

      int nKnots = spline.getNumberOfKnots();
      int nIntervals = nKnots - 1;

      std::vector<int> nDataPerInterval(nIntervals, 0);            //bei drei plus 1 anstatt plus 2? Testprogramm schreiben, verschiedene Szenarien durchrechnen
      std::vector<int> nDataPerKnot(spline.getNumberOfKnots(), 0); // <---

      for (int i = 0; i < vu.size(); i++) {
        int interval = spline.getLeftKnotIndexForU(vu[i]); //is the data point in this interval? like this, a border point is only counted for the first interval
        nDataPerInterval[interval]++;
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

    //-----------------------------------------------end of test 2 -------------------------------------------------

    std::cout << " matrix solver: " << constructionOk << ", anna's test: " << testAnna << std::endl;

    if (drawOption == 1 && testAnna == constructionOk) {
      continue;
    }

    // spline.print();

    canv->Draw();

    for (int i = 0; i < nKnots; i++) {
      double u = spline.getKnot(i).u;
      double fs = spline.interpolate(spline.convUtoX(u));
      drawPoints->Fill(2, u, fs);
    }

    TNtuple* nt = new TNtuple("nt", "nt", "u:f0:fSpline");

    float stepS = 1.e-4;
    int nSteps = (int)(1. / stepS + 1);

    double drawMax = -1.e20;
    double drawMin = 1.e20;

    for (float s = 0; s < 1. + stepS; s += stepS) {
      double u = s * (nKnots - 1);
      double f0 = F(u);
      double fSpline = spline.interpolate(spline.convUtoX(u));
      nt->Fill(u, f0, fSpline);
      drawMax = std::max(drawMax, std::max(f0, fSpline));
      drawMin = std::min(drawMin, std::min(f0, fSpline));
    }

    /*
      canv->cd(1);
      qaX->Draw();
      canv->cd(2);
    */

    // nt->SetMarkerColor(kBlack);
    // nt->Draw("f0:u","","");

    {
      TNtuple* ntRange = new TNtuple("ntRange", "nt", "u:f");
      double L = drawMax - drawMin;
      drawMin -= 0.0 * L;
      drawMax += 0.1 * L;

      ntRange->Fill(-0.001, drawMin);
      ntRange->Fill(-0.001, drawMax);
      ntRange->Fill(nKnots - 1 - 0.005, drawMin);
      ntRange->Fill(nKnots - 1 - 0.005, drawMax);
      ntRange->SetMarkerColor(kWhite);
      ntRange->SetMarkerSize(0.1);
      ntRange->Draw("f:u", "", "");
      delete ntRange;
    }

    //auto legend = new TLegend(0.1, 0.72, 0.4, 0.95);
    // legend->SetHeader("Splines of the same size:","C"); // option "C" allows to center the header

    nt->SetMarkerColor(kGray);
    nt->SetMarkerStyle(8);
    nt->SetMarkerSize(2.);
    nt->Draw("f0:u", "", "P,same");

    TH1* htemp = (TH1*)gPad->GetPrimitive("htemp");
    htemp->SetTitle("Splines of the same size");
    htemp->GetXaxis()->SetTitle("x");
    htemp->GetYaxis()->SetTitle("R");

    TLine* l0 = new TLine();
    l0->SetLineWidth(7);
    l0->SetLineColor(kGray);
    // legend->AddEntry(l0, "Input function", "L");
    //legend->AddEntry(l0, "Function to approximate", "L");
    //legend->Draw();

    drawPoints->SetMarkerStyle(8);
    drawPoints->SetMarkerSize(1.5);

    nt->SetMarkerSize(1.);
    nt->SetMarkerColor(kGreen + 2);
    nt->Draw("fSpline:u", "", "P,same");

    drawPoints->SetMarkerStyle(21);
    drawPoints->SetMarkerColor(kGreen + 2);
    drawPoints->SetMarkerSize(2.5);             // 5.
    drawPoints->Draw("f:u", "type==2", "same"); // Spline
    TNtuple* l1 = new TNtuple();
    l1->SetMarkerStyle(21);
    l1->SetMarkerColor(kGreen + 2);
    l1->SetMarkerSize(1.5); // 3.5
    l1->SetLineColor(kGreen + 2);
    l1->SetLineWidth(2.); // 5.
    // legend->AddEntry(l1, Form("Interpolation spline (%d drawPoints + %d slopes)", nKnots, nKnots), "LP");
    //legend->AddEntry(l1, "Best-fit spline", "LP");
    //legend->Draw();

    {
      drawPoints->SetMarkerColor(kBlack);
      drawPoints->SetMarkerSize(1.5);
      drawPoints->SetMarkerStyle(8);

      drawPoints->Draw("f:u", "type==1", "same"); // best-fit data points
      // drawPoints->Draw("f:u", "type==5", "same"); // chebyshev, data points
      TMarker* l5 = new TMarker;
      l5->SetMarkerStyle(8);
      l5->SetMarkerSize(1.5);
      l5->SetMarkerColor(kBlack);
      //legend->AddEntry(l5, "Construction points", "P");
      //legend->Draw();
    }

    if (!ask()) {
      break;
    }

    //delete legend;
  }

  return 0;
}
