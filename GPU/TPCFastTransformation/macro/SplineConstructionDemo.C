/*
   root -l -q SplineDemo.C
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
int nKnots = 4;

static double Fcoeff[2 * (Fdegree + 1)];

void F(double u, double f[])
{
  double uu = u * TMath::Pi() / (nKnots - 1);
  f[0] = 0; // Fcoeff[0]/2;
  for (int i = 1; i <= Fdegree; i++) {
    f[0] += Fcoeff[2 * i] * TMath::Cos(i * uu) + Fcoeff[2 * i + 1] * TMath::Sin(i * uu);
  }
}

/*
void F(float u, float f[])
{
  double uu = u * 2 / (nKnots - 1)-1.;
  double t0=1;
  double t1 = uu;
  f[0] = 0;
  for (int i = 1; i <= Fdegree*2; i++) {
     double t = t = 2*uu*t1-t0;
     f[0] += Fcoeff[i]*t;
     t0 = t1;
     t1 = t;
  }
}
*/

double F1D(double u)
{
  double f = 0;
  F(u, &f);
  return f;
}

double Flocal(double u)
{
  double f = 0;
  F(u * (nKnots - 1), &f);
  return f;
}

TCanvas* canv = new TCanvas("cQA", "Spline Demo", 1600, 800);

bool doAskSteps = 1;
bool drawConstruction = 1;

bool ask()
{
  canv->Update();
  std::cout << "type: 'q'-exit";
  std::cout << ", 's'-individual steps";
  std::cout << ", 'p'-construction points";

  std::cout << std::endl;

  std::string str;
  std::getline(std::cin, str);
  if (str == "s") {
    doAskSteps = !doAskSteps;
  } else if (str == "p") {
    drawConstruction = !drawConstruction;
  }

  return (str != "q" && str != ".q");
}

bool askStep()
{
  return (!doAskSteps) ? 1 : ask();
}

int SplineConstructionDemo()
{

  const int nAxiliaryPoints = 10;

  using namespace o2::gpu;

  std::cout << "Test interpolation.." << std::endl;

  // TCanvas* canv = new TCanvas("cQA", "Spline1D  QA", 2000, 1000);

  gRandom->SetSeed(0);

  TH1F* histDfSpline = new TH1F("histDfSpline", "Df Spline", 100, -1., 1.);
  TH1F* histDfSpline1 = new TH1F("histDfSpline1", "Df Spline1", 100, -1., 1.);
 
  for (int seed = 12;; seed++) {

    // seed = gRandom->Integer(100000); // 605

    gRandom->SetSeed(seed);
    std::cout << "Random seed: " << seed << " " << gRandom->GetSeed() << std::endl;

    for (int i = 0; i < 2 * (Fdegree + 1); i++) {
      Fcoeff[i] = gRandom->Uniform(-1, 1);
    }

    TNtuple* drawPoints = new TNtuple("drawPoints", "drawPoints", "type:u:f");

    o2::gpu::Spline1D<float, 1> spline(nKnots);
    spline.approximateFunction(0, nKnots - 1, F, nAxiliaryPoints);

    o2::gpu::Spline1D<float, 1> spline1(nKnots); // spline with missing data points          ---         with how many knots does it crash? ----
    o2::gpu::Spline1DHelper<float> helper;
    {
      vector<double> vu, vy;
      for (int i = 0; i < nKnots; i += 1) {
        double u = spline.getKnot(i).u;
        if (i >= 1) {
          vu.push_back(u);
          vy.push_back(F1D(u));

         /* vu.push_back(u + 0.6);
          vy.push_back(F1D(u + 0.6)); */
        }
        if (i >= 0) {

          vu.push_back(u + 0.5);
          vy.push_back(F1D(u + 0.5));

          /*
          vu.push_back(u + 0.6);
          vy.push_back(F1D(u + 0.6));
          */
        }
      }
      helper.approximateDataPoints(spline1, 0, nKnots - 1, &vu[0], &vy[0], vu.size());
      for(int i=0; i<vu.size(); i++){
        drawPoints->Fill(1, vu[i], vy[i]); 
      }
    }


    spline.print();
    spline1.print();

    canv->Draw();

    for (int i = 0; i < nKnots; i++) {
      double u = spline1.getKnot(i).u;
      double fs = spline1.interpolate(spline1.convUtoX(u));
      drawPoints->Fill(2, u, fs);
    }

   for (int i = 0; i < nKnots; i++) {
      double u = spline.getKnot(i).u;
      double fs = spline.interpolate(spline.convUtoX(u));
      drawPoints->Fill(4, u, fs);
    }

    o2::gpu::Spline1DHelperOld<float> helperOld;
    helperOld.setSpline(spline, 1, nAxiliaryPoints);
    for (int j = 0; j < helperOld.getNumberOfDataPoints(); j++) {
      const typename Spline1DHelperOld<float>::DataPoint& p = helperOld.getDataPoint(j);
      double f0;
      F(p.u, &f0);
      double fs = spline.interpolate(spline.convUtoX(p.u));    
      drawPoints->Fill(3, p.u, f0); // spline knos
    }

    TNtuple* nt = new TNtuple("nt", "nt", "u:f0:fSpline:fSpline1");

    float stepS = 1.e-4;
    int nSteps = (int)(1. / stepS + 1);

    double statDfSpline = 0;
    double statDfSpline1 = 0;
  
    double statMinMaxSpline = 0;
    double statMinMaxSpline1 = 0;
   
    double drawMax = -1.e20;
    double drawMin = 1.e20;
    int statN = 0;
    for (float s = 0; s < 1. + stepS; s += stepS) {
      double u = s * (nKnots - 1);
      double f0;
      F(u, &f0);
      double fSpline = spline.interpolate(spline.convUtoX(u));
      double fSpline1 = spline1.interpolate(spline1.convUtoX(u));
       nt->Fill(u, f0, fSpline, fSpline1);
      drawMax = std::max(drawMax, (double)f0);
      drawMin = std::min(drawMin, (double)f0);
      drawMax = std::max(drawMax, std::max(fSpline, fSpline1));
      drawMin = std::min(drawMin, std::min(fSpline, fSpline1));
      statDfSpline += (fSpline - f0) * (fSpline - f0);
      statDfSpline1 += (fSpline1 - f0) * (fSpline1 - f0);
      statN++;
      histDfSpline->Fill(fSpline - f0);
      histDfSpline1->Fill(fSpline1 - f0);

      statMinMaxSpline = std::max(statMinMaxSpline, fabs(fSpline - f0));
      statMinMaxSpline1 = std::max(statMinMaxSpline1, fabs(fSpline1 - f0));
     }

    //histMinMaxSpline->Fill(statMinMaxSpline);
    //histMinMaxSpline1->Fill(statMinMaxSpline1);

    std::cout << "\n"
              << std::endl;
    std::cout << "\nRandom seed: " << seed << " " << gRandom->GetSeed() << std::endl;
    std::cout << "Spline : std dev " << sqrt(statDfSpline / statN) << " minmax " << statMinMaxSpline << std::endl;
    std::cout << "Spline1     : std dev " << sqrt(statDfSpline1 / statN) << " minmax " << statMinMaxSpline1 << std::endl;

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

    auto legend = new TLegend(0.1, 0.72, 0.4, 0.95);
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
    legend->AddEntry(l0, "Function to approximate", "L");
    legend->Draw();

    drawPoints->SetMarkerStyle(8);
    drawPoints->SetMarkerSize(1.5);

    if (!askStep()) {
      break;
    }

    nt->SetMarkerSize(1.);
    nt->SetMarkerColor(kGreen + 2);
    nt->Draw("fSpline:u", "", "P,same");

    drawPoints->SetMarkerStyle(21);
    drawPoints->SetMarkerColor(kGreen + 2);
    drawPoints->SetMarkerSize(2.5);             // 5.
    drawPoints->Draw("f:u", "type==4", "same"); // Spline
    TNtuple* l1 = new TNtuple();
    l1->SetMarkerStyle(21);
    l1->SetMarkerColor(kGreen + 2);
    l1->SetMarkerSize(1.5); // 3.5
    l1->SetLineColor(kGreen + 2);
    l1->SetLineWidth(2.); // 5.
    // legend->AddEntry(l1, Form("Interpolation spline (%d drawPoints + %d slopes)", nKnots, nKnots), "LP");
    legend->AddEntry(l1, "Best-fit spline", "LP");
    legend->Draw();

    if (!askStep()) {
      break;
    }

    nt->SetMarkerColor(kRed);
    nt->Draw("fSpline1:u", "", "P,same");

    drawPoints->SetMarkerStyle(20);
    drawPoints->SetMarkerColor(kRed);
    drawPoints->SetMarkerSize(2.5);             // 5.
    drawPoints->Draw("f:u", "type==2", "same"); // best-fit splines

    // TMarker * l3 = new TMarker();
    TNtuple* l3 = new TNtuple();
    l3->SetMarkerStyle(20);
    l3->SetMarkerColor(kRed);
    l3->SetMarkerSize(2.5); // 3.5
    l3->SetLineColor(kRed);
    l3->SetLineWidth(5.);
    // legend->AddEntry(l3, Form("Best-fit spline (%d drawPoints + %d slopes)", nKnots, nKnots), "PL");
    legend->AddEntry(l3, "Best-fit spline 1", "PL");
    legend->Draw();

    drawPoints->SetMarkerStyle(8);

    if (!askStep()) {
      break;
    }

    if (drawConstruction) {
      drawPoints->SetMarkerColor(kYellow);
      drawPoints->SetMarkerSize(1.5);
      drawPoints->SetMarkerStyle(8);

      drawPoints->Draw("f:u", "type==3", "same"); // best-fit data points
      // drawPoints->Draw("f:u", "type==5", "same"); // chebyshev, data points
      TMarker* l4 = new TMarker;
      l4->SetMarkerStyle(8);
      l4->SetMarkerSize(1.5);
      l4->SetMarkerColor(kYellow);
      legend->AddEntry(l4, "Construction points", "P");
      legend->Draw();

      drawPoints->SetMarkerColor(kBlack);
      drawPoints->SetMarkerSize(1.5);
      drawPoints->SetMarkerStyle(8);

      drawPoints->Draw("f:u", "type==1", "same"); // best-fit data points
      // drawPoints->Draw("f:u", "type==5", "same"); // chebyshev, data points
      TMarker* l5 = new TMarker;
      l4->SetMarkerStyle(8);
      l4->SetMarkerSize(1.5);
      l4->SetMarkerColor(kBlack);
      legend->AddEntry(l5, "Construction points", "P");
      legend->Draw();
 
      if (!askStep()) {
        break;
      }
    }

    if (!doAskSteps && !ask()) {
      break;
    }

    delete legend;
  }

  return 0;
}
