#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TList.h>
#include "../include/StyleProposal.h"
#include "../include/PaperColors.h"
using namespace std;

//--------Important info---------------------//
//=====In data:
// 1. K*0, K*+- is averaged
// 2. Pi,K,p is averaged
//=====In EPOS:
// 1. Only K*0 is averaged.
//===== In Pythia (Local):
// 1. Nothing is averaged. All as summed.

// Important: Do not scale any graph without cloning it first, otherwise the original graph will be scaled and it will affect other plots.
TGraphErrors *SmoothGraph(const TGraphErrors *gr);
TGraphErrors *MakeRatio(const TGraphErrors *numerator, const TGraphErrors *denominator, bool isModel = false, float CorrFactorDen = 1);
TGraphErrors **MakeRatioUncorr(TGraphErrors **numerator, TGraphErrors **denominator, float CorrFactorDen = 1, int nGraphs = 3);
TGraphErrors *GetGraph(TFile *f, const string &name);
TFile *OpenFile(const string &path);

void ScaleGraph(TGraph *gr, double scale);
TGraphErrors **DivideByMult(TGraphErrors **gr, double WhichMultPoint, double tolerance = 0.5, int noOfGraphs = 3);
TGraphErrors *DivideByMultModel(TGraphErrors *gr, double WhichMultPoint, double tolerance = 0.5);
void RestrictModelXaxis(TGraphErrors *gr, double xMin, double xMax);

struct SuppressionSignificanceResult
{
  double highMultiplicity;
  double lowMultiplicity;
  double systematicSignificance;
  double uncorrelatedSignificance;
  double uncorrelatedPlusStatSignificance;
  double highToLowRatio;
};

// Graph ordering follows the existing macro convention: point 0 is the
// highest-multiplicity class and the last point is the lowest-multiplicity class.
SuppressionSignificanceResult CalculateSuppressionSignificance(
    TGraphErrors **ratioGraphs,
    const string &ratioLabel,
    bool includeStatInSystematic = false,
    bool reportUncorrelatedPlusStat = false,
    bool reportHighToLowRatio = false)
{
  const double invalid = std::numeric_limits<double>::quiet_NaN();
  SuppressionSignificanceResult result{invalid, invalid, invalid, invalid, invalid, invalid};
  if (!ratioGraphs || !ratioGraphs[0] || !ratioGraphs[1] || !ratioGraphs[2])
  {
    cout << "Cannot calculate " << ratioLabel << " suppression significance: missing ratio/error graph." << endl;
    return result;
  }

  const int nPoints = ratioGraphs[0]->GetN();
  if (nPoints < 1 || ratioGraphs[1]->GetN() != nPoints || ratioGraphs[2]->GetN() != nPoints)
  {
    cout << "Cannot calculate " << ratioLabel << " suppression significance: empty or inconsistent ratio graphs." << endl;
    return result;
  }

  double ignoredX;
  ratioGraphs[1]->GetPoint(0, ignoredX, result.highMultiplicity);
  ratioGraphs[1]->GetPoint(nPoints - 1, ignoredX, result.lowMultiplicity);
  const double sysHigh = ratioGraphs[1]->GetErrorY(0);
  const double sysLow = ratioGraphs[1]->GetErrorY(nPoints - 1);
  const double statHigh = ratioGraphs[0]->GetErrorY(0);
  const double statLow = ratioGraphs[0]->GetErrorY(nPoints - 1);
  const double uncorrHigh = ratioGraphs[2]->GetErrorY(0);
  const double uncorrLow = ratioGraphs[2]->GetErrorY(nPoints - 1);
  const double difference = std::fabs(result.highMultiplicity - result.lowMultiplicity);

  double systematicVariance = sysHigh * sysHigh + sysLow * sysLow;
  if (includeStatInSystematic)
    systematicVariance += statHigh * statHigh + statLow * statLow;
  const double uncorrelatedVariance = uncorrHigh * uncorrHigh + uncorrLow * uncorrLow;
  const double uncorrelatedStatVariance = uncorrelatedVariance + statHigh * statHigh + statLow * statLow;

  result.systematicSignificance = systematicVariance > 0.0 ? difference / std::sqrt(systematicVariance) : invalid;
  result.uncorrelatedSignificance = uncorrelatedVariance > 0.0 ? difference / std::sqrt(uncorrelatedVariance) : invalid;
  result.uncorrelatedPlusStatSignificance = uncorrelatedStatVariance > 0.0 ? difference / std::sqrt(uncorrelatedStatVariance) : invalid;
  result.highToLowRatio = result.lowMultiplicity != 0.0 ? result.highMultiplicity / result.lowMultiplicity : invalid;

  cout << ratioLabel << " HM-LM significance ("
       << (includeStatInSystematic ? "systematic + statistical" : "systematic")
       << "): " << result.systematicSignificance << " sigma" << endl;
  cout << ratioLabel << " HM-LM significance (uncorrelated): "
       << result.uncorrelatedSignificance << " sigma" << endl;
  if (reportUncorrelatedPlusStat)
    cout << ratioLabel << " HM-LM significance (uncorrelated + statistical): "
         << result.uncorrelatedPlusStatSignificance << " sigma" << endl;
  if (reportHighToLowRatio)
    cout << ratioLabel << " HM/LM ratio: " << result.highToLowRatio << endl;
  return result;
}

void ParticleRatios()
{
  gROOT->SetBatch(kTRUE); // Avoid opening the intermediate canvases while the source plots are prepared.
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(0);
  // gStyle->SetLineScalePS(3.0);

  Int_t Run2Color = PaperColors::GetROOTColor(PaperColors::kRun2);
  Int_t Run3Color = PaperColors::GetROOTColor(PaperColors::kRun3);

  //============================================================
  // USER CONFIGURATION: edit paths and plot output settings here
  //============================================================
  bool isSavePlots = true;
  const string plotOutputDirectory = "../figures";

  // This work
  const string inputKstarRun3 = "../inputs/Kstar0_136/Results.root";
  const string inputPhiRun3 = "../inputs/PhiResults/Sawan/ResultsPhi.root";

  const string inputPionRun3 = "../inputs/PiKp/Pi_results.root";
  const string inputProtonRun3 = "../inputs/PiKp/Pr_results.root";
  const string inputKaonRun3 = "../inputs/PiKp/Ka_results.root";
  const string inputPiKpCombined = "../inputs/PreviousALICEMeasurements/ppRun3_PiKpHEP.root";
  const string inputKstarPhiRun2 = "../inputs/PreviousALICEMeasurements/pp13TeV_KstarPhi.root";
  const string inputChargedKstarRun2 = "../inputs/PreviousALICEMeasurements/pp13TeV_ChKstar.root";
  const string inputMultStrangeRun2 = "../inputs/PreviousALICEMeasurements/pp13TeV_MultStrange.root";
  const string inputXiStarSigmaRun2 = "../inputs/PreviousALICEMeasurements/pp13TeV_XiStarSigma.root";
  const string inputLambda1520Run2 = "../inputs/PreviousALICEMeasurements/pp13TeV_LambdaStar.root";
  const string inputLambda1520Run3 = "../inputs/LambdaRun3/Sawan/ResultsLambda1520.root";
  const string inputKshortRun3 = "../inputs/K0s_Run3_Results/Sawan/ResultsK0s.root";
  const string inputALICE7TeV = "../inputs/PreviousALICEMeasurements/pp7TeVALICE.root";

  // Model inputs
  const string inputPythiaModels = "../inputs/Pythia/ModelResultsXi.root";
  const string inputEPOS = "../inputs/EPOS/FinalEPOSResults_BongFix4.root";
  const string inputThermalModel = "../inputs/ThermalModel/thermal_results_SCE_new2.root";

  // DATA INPUTS

  // Run 3
  TFile *fKstar = OpenFile(inputKstarRun3);
  TFile *fPhiRun3 = OpenFile(inputPhiRun3);

  TFile *fPion = OpenFile(inputPionRun3);
  TFile *fProton = OpenFile(inputProtonRun3);
  TFile *fKaon = OpenFile(inputKaonRun3);
  TFile *fLambda1520 = OpenFile(inputLambda1520Run3);
  TFile *fKshortRun3 = OpenFile(inputKshortRun3);
  TFile *fXiRun3 = OpenFile("../inputs/XiOmegaResults/Sawan/ResultsXi.root");

  TFile *fPiKp = OpenFile(inputPiKpCombined);

  TFile *fKstarPhiRun2 = OpenFile(inputKstarPhiRun2);
  TFile *fChKstar = OpenFile(inputChargedKstarRun2);
  TFile *fMultStrange = OpenFile(inputMultStrangeRun2);
  TFile *fXiStarSigma = OpenFile(inputXiStarSigmaRun2);
  TFile *fLambda1520Run2 = OpenFile(inputLambda1520Run2);

  // Available in run3

  // ALICE Run2 results
  TFile *fpp7TeV = OpenFile(inputALICE7TeV);

  // Run 3 data graphs: mean pT
  TGraphErrors *gMPtKstar[3];
  TGraphErrors *gMPtPion[3];
  TGraphErrors *gMPtProton[3];
  TGraphErrors *gMPtKaon[3];
  TGraphErrors *gMPtKshortRun3[3];
  TGraphErrors *gMPtPhiRun3[3];

  // Run 3 data graphs: yields
  TGraphErrors *gMYieldKstar[3];
  TGraphErrors *gMYieldPion[3];
  TGraphErrors *gMYieldProton[3];
  TGraphErrors *gMYieldKaon[3];
  TGraphErrors *gMYieldKshortRun3[3];
  TGraphErrors *gMYieldPhiRun3[3];
  TGraphErrors *gMYieldXiRun3[3];

  // Run 2 data graphs: mean pT
  TGraphErrors *gMPtKstarRun2[3];
  TGraphErrors *gMPtPhiRun2[3];
  TGraphErrors *gMPtChKstarRun2[3];
  TGraphErrors *gMPtKshortRun2[3];
  TGraphErrors *gMPtLambdaRun2[3];
  TGraphErrors *gMPtXiRun2[3];
  TGraphErrors *gMPtOmegaRun2[3];
  TGraphErrors *gMPtXiStarRun2[3];
  TGraphErrors *gMPtSigmaRun2[3];
  TGraphErrors *gMPtLambda1520[3];
  TGraphErrors *gMPtRho[3];

  // Run 2 data graphs: yields
  TGraphErrors *gMYieldKstarRun2[3];
  TGraphErrors *gMYieldPhiRun2[3];
  TGraphErrors *gMYieldChKstarRun2[3];
  TGraphErrors *gMYieldKshortRun2[3];
  TGraphErrors *gMYieldLambdaRun2[3];
  TGraphErrors *gMYieldXiRun2[3];
  TGraphErrors *gMYieldOmegaRun2[3];
  TGraphErrors *gMYieldXiStarRun2[3];
  TGraphErrors *gMYieldSigmaRun2[3];
  TGraphErrors *gMYieldLambda1520[3];
  TGraphErrors *gMYieldRho[3];

  for (int i = 0; i < 3; i++)
  {
    string suffix;
    if (i == 0)
      suffix = "_stat";
    else if (i == 1)
      suffix = "_sys";
    else
      suffix = "_sysuncorr";

    // Run 3 mean pT graphs
    gMPtKstar[i] = GetGraph(fKstar, Form("gMeanpTRun3%s", suffix.c_str()));
    gMPtPhiRun3[i] = GetGraph(fPhiRun3, Form("gMeanpTRun3%s", suffix.c_str()));
    gMPtPion[i] = GetGraph(fPion, Form("gMeanpTRun3%s", suffix.c_str()));
    gMPtKaon[i] = GetGraph(fKaon, Form("gMeanpTRun3%s", suffix.c_str()));
    gMPtProton[i] = GetGraph(fProton, Form("gMeanpTRun3%s", suffix.c_str()));
    gMPtKshortRun3[i] = GetGraph(fKshortRun3, Form("gMeanpTRun3%s", suffix.c_str()));

    // Run 3 yield graphs
    gMYieldKstar[i] = GetGraph(fKstar, Form("gMeanYieldRun3%s", suffix.c_str()));
    gMYieldPhiRun3[i] = GetGraph(fPhiRun3, Form("gMeanYieldRun3%s", suffix.c_str()));
    gMYieldPion[i] = GetGraph(fPion, Form("gMeanYieldRun3%s", suffix.c_str()));
    gMYieldKaon[i] = GetGraph(fKaon, Form("gMeanYieldRun3%s", suffix.c_str()));
    gMYieldProton[i] = GetGraph(fProton, Form("gMeanYieldRun3%s", suffix.c_str()));
    gMYieldKshortRun3[i] = GetGraph(fKshortRun3, Form("gMeanYieldRun3%s", suffix.c_str()));
    gMYieldXiRun3[i] = GetGraph(fXiRun3, Form("gMeanYieldRun3%s", suffix.c_str()));

    // Run 2 mean pT graphs
    gMPtKstarRun2[i] = GetGraph(fKstarPhiRun2, Form("gKstar_MeanpT%s", suffix.c_str()));
    gMPtPhiRun2[i] = GetGraph(fKstarPhiRun2, Form("gPhi_MeanpT%s", suffix.c_str()));
    gMPtChKstarRun2[i] = GetGraph(fChKstar, Form("gChKstar_MeanpT%s", suffix.c_str()));
    gMPtKshortRun2[i] = GetGraph(fMultStrange, Form("gKshort_MeanpT%s", suffix.c_str()));
    gMPtLambdaRun2[i] = GetGraph(fMultStrange, Form("gLambda_MeanpT%s", suffix.c_str()));
    gMPtXiRun2[i] = GetGraph(fMultStrange, Form("gXi_MeanpT%s", suffix.c_str()));
    gMPtOmegaRun2[i] = GetGraph(fMultStrange, Form("gOmega_MeanpT%s", suffix.c_str()));
    gMPtXiStarRun2[i] = GetGraph(fXiStarSigma, Form("gXi1530_MeanpT%s", suffix.c_str()));
    gMPtSigmaRun2[i] = GetGraph(fXiStarSigma, Form("gSigma1385_MeanpT%s", suffix.c_str()));
    gMPtLambda1520[i] = GetGraph(fLambda1520Run2, Form("gLambda_MeanpT%s", suffix.c_str()));

    // Run 2 yield graphs
    gMYieldKstarRun2[i] = GetGraph(fKstarPhiRun2, Form("gKstar_MeanYield%s", suffix.c_str()));
    gMYieldPhiRun2[i] = GetGraph(fKstarPhiRun2, Form("gPhi_MeanYield%s", suffix.c_str()));
    gMYieldChKstarRun2[i] = GetGraph(fChKstar, Form("gChKstar_MeanYield%s", suffix.c_str()));
    gMYieldKshortRun2[i] = GetGraph(fMultStrange, Form("gKshort_MeanYield%s", suffix.c_str()));
    gMYieldLambdaRun2[i] = GetGraph(fMultStrange, Form("gLambda_MeanYield%s", suffix.c_str()));
    gMYieldXiRun2[i] = GetGraph(fMultStrange, Form("gXi_MeanYield%s", suffix.c_str()));
    gMYieldOmegaRun2[i] = GetGraph(fMultStrange, Form("gOmega_MeanYield%s", suffix.c_str()));
    gMYieldXiStarRun2[i] = GetGraph(fXiStarSigma, Form("gXi1530_MeanYield%s", suffix.c_str()));
    gMYieldSigmaRun2[i] = GetGraph(fXiStarSigma, Form("gSigma1385_MeanYield%s", suffix.c_str()));
    gMYieldLambda1520[i] = GetGraph(fLambda1520Run2, Form("gLambda_MeanYield%s", suffix.c_str()));
  }

  // Published ratio graphs: 13 TeV and 7 TeV
  // 13 TeV ratios are stored in the K* / Phi Run 2 results file.
  TGraphErrors *gKstarKshortRatio_13TeV[2];
  TGraphErrors *gPhiKshortRatio_13TeV[2];
  TGraphErrors *gPhiPiRatio_13TeV[2];
  TGraphErrors *gKstarKaRatio_13TeV[2];
  TGraphErrors *gKstarPiRatio_13TeV[2];
  TGraphErrors *gPhiKaRatio_13TeV[2];
  TGraphErrors *gPhiPrRatio_13TeV[2];
  TGraphErrors *gPhiLambdaRatio_13TeV[2];

  // 7 TeV ratios are stored in the ALICE 7 TeV results file.
  TGraphErrors *gKstarPiRatio_7TeV[2];
  TGraphErrors *gPhiPiRatio_7TeV[2];

  for (int i = 0; i < 2; i++)
  {
    const char *uncertaintySuffix = (i == 0) ? "_stat" : "_sys";

    // 13 TeV ratios
    gKstarKshortRatio_13TeV[i] = GetGraph(fKstarPhiRun2, Form("gKstar_KshortRatio%s", uncertaintySuffix));
    gPhiKshortRatio_13TeV[i] = GetGraph(fKstarPhiRun2, Form("gPhi_KshortRatio%s", uncertaintySuffix));
    gPhiPiRatio_13TeV[i] = GetGraph(fKstarPhiRun2, Form("gPhi_PiRatio%s", uncertaintySuffix));
    gKstarKaRatio_13TeV[i] = GetGraph(fKstarPhiRun2, Form("gKstar_KaRatio%s", uncertaintySuffix));
    gKstarPiRatio_13TeV[i] = GetGraph(fKstarPhiRun2, Form("gKstar_PiRatio%s", uncertaintySuffix));
    gPhiKaRatio_13TeV[i] = GetGraph(fKstarPhiRun2, Form("gPhi_KaRatio%s", uncertaintySuffix));
    gPhiPrRatio_13TeV[i] = GetGraph(fKstarPhiRun2, Form("gPhi_PrRatio%s", uncertaintySuffix));
    gPhiLambdaRatio_13TeV[i] = GetGraph(fKstarPhiRun2, Form("gPhi_LambdaRatio%s", uncertaintySuffix));

    // 7 TeV ratios
    gKstarPiRatio_7TeV[i] = GetGraph(fpp7TeV, Form("gKstar_PiRatio%s", uncertaintySuffix));
    gPhiPiRatio_7TeV[i] = GetGraph(fpp7TeV, Form("gPhi_PiRatio%s", uncertaintySuffix));
  }

  // MODEL INPUTS AND LOADING: Pythia Hyperloop
  //================================================
  //==============Models from hyperloop outputs======
  //=================================================

  enum Models
  {
    kPythiaMonash,
    kPythiaRopes,
    kPythiaShoving,
    kEPOSWithUrQMD,
    kEPOSWithoutUrQMD,
    kThermalModel,
    kNModels
  };
  enum Particle
  {
    kKstar,
    kPhi,
    kPion,
    kKaon,
    kProton,
    kPionMinus,
    kKaonMinus,
    kAntiProton,
    kXi1530,
    kKshort,
    kKstarPM,
    kXi,
    kXiMinus,
    kNParticles
  };

  struct Model
  {
    Int_t color;
    int style;
    int width;
    const char *label;
  };

  Model models[kNModels] = {
      {PaperColors::GetROOTColor(PaperColors::kMonash), 1, 2, "Pythia Monash"},    // Pythia Monash
      {PaperColors::GetROOTColor(PaperColors::kRopes), 2, 2, "Pythia Ropes"},      // Pythia Ropes
      {PaperColors::GetROOTColor(PaperColors::kShoving), 3, 2, "Pythia Shoving"},  // Pythia Shoving
      {PaperColors::GetROOTColor(PaperColors::kEPOS), 7, 2, "EPOS4"},              // EPOS4
      {PaperColors::GetROOTColor(PaperColors::kEPOSUrQMD), 9, 2, "EPOS4 + UrQMD"}, // EPOS4 with UrQMD
      {PaperColors::GetROOTColor(PaperColors::kThermal), 5, 2, "Thermal model"},   // Thermal model
  };

  TFile *ModelsHyperloop = new TFile(inputPythiaModels.c_str(), "read");
  if (ModelsHyperloop->IsZombie())
  {
    cout << "Error: Hyperloop model file not found" << endl;
    return;
  }

  constexpr int nPythiaModels = kEPOSWithUrQMD;

  string hyperloopModels[nPythiaModels] = {"Pythia_Monash", "Pythia_Ropes", "Pythia_Shoving"};
  string particles[kNParticles] = {"Kstar", "Phi", "Pion", "Kaon", "Proton", "PionMinus", "KaonMinus", "AntiProton", "Xi1530", "Kshort", "KstarPM", "XiMinus", "XiPlus"};
  vector<int> modelsToPlot = {kPythiaMonash, kPythiaShoving, kPythiaRopes};

  TGraphErrors *gMYield[nPythiaModels][kNParticles];
  TGraphErrors *gMeanPt[nPythiaModels][kNParticles];
  TH1D *hSpectra[nPythiaModels][kNParticles];

  for (int iModel = 0; iModel < nPythiaModels; iModel++)
  {
    for (int iParticle = 0; iParticle < kNParticles; iParticle++)
    {
      gMYield[iModel][iParticle] = GetGraph(ModelsHyperloop, hyperloopModels[iModel] + "/" + particles[iParticle] + "/gMeanYield_" + particles[iParticle]);
      gMeanPt[iModel][iParticle] = GetGraph(ModelsHyperloop, hyperloopModels[iModel] + "/" + particles[iParticle] + "/gMeanpT_" + particles[iParticle]);
      hSpectra[iModel][iParticle] = (TH1D *)ModelsHyperloop->Get((hyperloopModels[iModel] + "/" + particles[iParticle] + "/hPt_" + particles[iParticle] + "_MinBias").c_str());

      if (gMYield[iModel][iParticle] == nullptr || gMeanPt[iModel][iParticle] == nullptr || hSpectra[iModel][iParticle] == nullptr)
      {
        cout << "Error: Graphs or histogram for model " << hyperloopModels[iModel] << " and particle " << particles[iParticle] << " not found" << endl;
        return;
      }

      SetGraphStyle(gMYield[iModel][iParticle]);
      SetGraphStyle(gMeanPt[iModel][iParticle]);
      SetHistStyle(hSpectra[iModel][iParticle]);

      // Apply model style once
      gMYield[iModel][iParticle]->SetLineColor(models[iModel].color);
      gMYield[iModel][iParticle]->SetLineStyle(models[iModel].style);
      gMYield[iModel][iParticle]->SetLineWidth(models[iModel].width);

      gMeanPt[iModel][iParticle]->SetLineColor(models[iModel].color);
      gMeanPt[iModel][iParticle]->SetLineStyle(models[iModel].style);
      gMeanPt[iModel][iParticle]->SetLineWidth(models[iModel].width);

      int nPoints = gMYield[iModel][iParticle]->GetN();
    }

    // Add pion with pionMinus, kaon with kaonMinus, proton with anti-proton
    for (int i = 0; i < gMYield[iModel][kPion]->GetN(); i++)
    {
      double x, yPion, yPionMinus;
      gMYield[iModel][kPion]->GetPoint(i, x, yPion);
      gMYield[iModel][kPionMinus]->GetPoint(i, x, yPionMinus);
      double yPionTotal = yPion + yPionMinus;
      gMYield[iModel][kPion]->SetPoint(i, x, yPionTotal);

      // Kaon
      double yKaon, yKaonMinus;
      gMYield[iModel][kKaon]->GetPoint(i, x, yKaon);
      gMYield[iModel][kKaonMinus]->GetPoint(i, x, yKaonMinus);
      double yKaonTotal = yKaon + yKaonMinus;
      gMYield[iModel][kKaon]->SetPoint(i, x, yKaonTotal);

      // Proton
      double yProton, yAntiProton;
      gMYield[iModel][kProton]->GetPoint(i, x, yProton);
      gMYield[iModel][kAntiProton]->GetPoint(i, x, yAntiProton);
      double yProtonTotal = yProton + yAntiProton;
      gMYield[iModel][kProton]->SetPoint(i, x, yProtonTotal);

      // Xi
      double yXi, yXiMinus;
      gMYield[iModel][kXi]->GetPoint(i, x, yXi);
      gMYield[iModel][kXiMinus]->GetPoint(i, x, yXiMinus);
      double yXiTotal = yXi + yXiMinus;
      gMYield[iModel][kXi]->SetPoint(i, x, yXiTotal);
    }
  }

  // MODEL INPUTS AND LOADING: EPOS
  //================================================
  //==============EPOS Model========================
  //================================================

  TFile *fEPOS = OpenFile(inputEPOS);

  TGraphErrors *gEPOS_KstardNdy_UrQMDOff = GetGraph(fEPOS, "g_dNdy_opt0");
  TGraphErrors *gEPOS_KstardNdy_UrQMDOn = GetGraph(fEPOS, "g_dNdy_opt1");
  TGraphErrors *gEPOS_KstarMeanpT_UrQMDOff = GetGraph(fEPOS, "g_pT_opt0");
  TGraphErrors *gEPOS_KstarMeanpT_UrQMDOn = GetGraph(fEPOS, "g_pT_opt1");

  TGraphErrors *gEPOS_PhidNdy_UrQMDOff = GetGraph(fEPOS, "g_phi_dNdy_opt0");
  TGraphErrors *gEPOS_PhidNdy_UrQMDOn = GetGraph(fEPOS, "g_phi_dNdy_opt1");
  TGraphErrors *gEPOS_PhiMeanpT_UrQMDOff = GetGraph(fEPOS, "g_phi_pT_opt0");
  TGraphErrors *gEPOS_PhiMeanpT_UrQMDOn = GetGraph(fEPOS, "g_phi_pT_opt1");

  for (int i = gEPOS_PhidNdy_UrQMDOff->GetN() - 1; i >= 0; --i)
  {
    gEPOS_PhidNdy_UrQMDOff->SetPointError(i, 0, 0);   // Remove error
    gEPOS_PhidNdy_UrQMDOn->SetPointError(i, 0, 0);    // Remove error
    gEPOS_PhiMeanpT_UrQMDOff->SetPointError(i, 0, 0); // Remove error
    gEPOS_PhiMeanpT_UrQMDOn->SetPointError(i, 0, 0);  // Remove error
  }

  for (auto model : modelsToPlot)
  {

    for (auto particle : {kKstar, kPhi, kPion, kKaon, kProton, kKshort})
    {
      gMYield[model][particle]->SetLineColor(models[model].color);
      gMYield[model][particle]->SetLineStyle(models[model].style);
      gMYield[model][particle]->SetLineWidth(models[model].width);

      gMeanPt[model][particle]->SetLineColor(models[model].color);
      gMeanPt[model][particle]->SetLineStyle(models[model].style);
      gMeanPt[model][particle]->SetLineWidth(models[model].width);
    }
  }

  // MODEL INPUTS AND LOADING: thermal model
  //================================================
  // ============Thermal model========================
  //================================================
  TFile *fThermal = OpenFile(inputThermalModel);
  TGraphErrors *gRatioKstarPhi_Thermal = GetGraph(fThermal, "ratioKstarPhi");
  TGraphErrors *gRatioKstarKa_Thermal = GetGraph(fThermal, "ratioKstarKaon");
  TGraphErrors *gRatioKstarPi_Thermal = GetGraph(fThermal, "ratioKstarPion");
  TGraphErrors *gRatioKstarK0s_Thermal = GetGraph(fThermal, "ratioKstarK0s");
  TGraphErrors *gRatioPhiKa_Thermal = GetGraph(fThermal, "ratioPhiKaon");
  TGraphErrors *gRatioPhiPi_Thermal = GetGraph(fThermal, "ratioPhiPion");
  TGraphErrors *gRatioPhiK0s_Thermal = GetGraph(fThermal, "ratioPhiK0s");
  TGraphErrors *gRatioPhiXi_Thermal = GetGraph(fThermal, "ratioPhiXi");

  // Scale the graphs for averaging
  ScaleGraph(gRatioKstarPhi_Thermal, 0.5);
  ScaleGraph(gRatioKstarK0s_Thermal, 0.5);
  ScaleGraph(gRatioPhiKa_Thermal, 2.0);
  ScaleGraph(gRatioPhiPi_Thermal, 2.0);
  ScaleGraph(gRatioPhiXi_Thermal, 1.0);

  // Calculate the ratios for EPOS

  TGraphErrors *gdNdy_Kstar_On = GetGraph(fEPOS, "g_dNdy_opt1");
  TGraphErrors *gdNdy_Kstar_Off = GetGraph(fEPOS, "g_dNdy_opt0");
  TGraphErrors *gdNdy_Phi_On = GetGraph(fEPOS, "g_phi_dNdy_opt1");
  TGraphErrors *gdNdy_Phi_Off = GetGraph(fEPOS, "g_phi_dNdy_opt0");
  TGraphErrors *gdNdy_K0s_On = GetGraph(fEPOS, "g_stable_dNdy_opt1_K0s");
  TGraphErrors *gdNdy_K0s_Off = GetGraph(fEPOS, "g_stable_dNdy_opt0_K0s");
  TGraphErrors *gdNdy_Kaon_On = GetGraph(fEPOS, "g_stable_dNdy_opt1_kaon");
  TGraphErrors *gdNdy_Kaon_Off = GetGraph(fEPOS, "g_stable_dNdy_opt0_kaon");
  TGraphErrors *gdNdy_Pion_On = GetGraph(fEPOS, "g_stable_dNdy_opt1_pi");
  TGraphErrors *gdNdy_Pion_Off = GetGraph(fEPOS, "g_stable_dNdy_opt0_pi");
  TGraphErrors *gdNdy_Xi_On = GetGraph(fEPOS, "g_stable_dNdy_opt1_Xi");
  TGraphErrors *gdNdy_Xi_Off = GetGraph(fEPOS, "g_stable_dNdy_opt0_Xi");

  TGraphErrors *gKstarKaRatio_EPOSNew_UrQMDOff = new TGraphErrors();
  TGraphErrors *gKstarKaRatio_EPOSNew_UrQMDOn = new TGraphErrors();
  TGraphErrors *gKstarPiRatio_EPOSNew_UrQMDOff = new TGraphErrors();
  TGraphErrors *gKstarPiRatio_EPOSNew_UrQMDOn = new TGraphErrors();
  TGraphErrors *gKstarK0sRatio_EPOSNew_UrQMDOff = new TGraphErrors();
  TGraphErrors *gKstarK0sRatio_EPOSNew_UrQMDOn = new TGraphErrors();

  TGraphErrors *gPhiKaRatio_EPOSNew_UrQMDOff = new TGraphErrors();
  TGraphErrors *gPhiKaRatio_EPOSNew_UrQMDOn = new TGraphErrors();
  TGraphErrors *gPhiPiRatio_EPOSNew_UrQMDOff = new TGraphErrors();
  TGraphErrors *gPhiPiRatio_EPOSNew_UrQMDOn = new TGraphErrors();
  TGraphErrors *gPhiK0sRatio_EPOSNew_UrQMDOff = new TGraphErrors();
  TGraphErrors *gPhiK0sRatio_EPOSNew_UrQMDOn = new TGraphErrors();
  TGraphErrors *gPhiXiRatio_EPOSNew_UrQMDOff = new TGraphErrors();
  TGraphErrors *gPhiXiRatio_EPOSNew_UrQMDOn = new TGraphErrors();

  TGraphErrors *gKstarPhiRatio_EPOSNew_UrQMDOff = new TGraphErrors();
  TGraphErrors *gKstarPhiRatio_EPOSNew_UrQMDOn = new TGraphErrors();

  int graphPointsEPOSNew = gdNdy_Kaon_On->GetN();
  for (int i = 0; i < graphPointsEPOSNew; i++)
  {
    double xK0s, yK0s, yKaon, yPion, yKstarOff, yKstarOn;
    double yK0sON, yKaonON, yPionON, yPhiOff, yPhiOn;

    gdNdy_Phi_Off->GetPoint(i, xK0s, yPhiOff);
    gdNdy_Phi_On->GetPoint(i, xK0s, yPhiOn);

    gdNdy_Kstar_Off->GetPoint(i, xK0s, yKstarOff);
    gdNdy_Kstar_On->GetPoint(i, xK0s, yKstarOn);

    gdNdy_K0s_Off->GetPoint(i, xK0s, yK0s);
    gdNdy_Kaon_Off->GetPoint(i, xK0s, yKaon);
    gdNdy_Pion_Off->GetPoint(i, xK0s, yPion);

    gdNdy_K0s_On->GetPoint(i, xK0s, yK0sON);
    gdNdy_Kaon_On->GetPoint(i, xK0s, yKaonON);
    gdNdy_Pion_On->GetPoint(i, xK0s, yPionON);

    gKstarKaRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yKaon);
    gKstarKaRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yKaon);
    gKstarPiRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yPion);
    gKstarPiRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yPion);
    gKstarK0sRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yK0s);
    gKstarK0sRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yK0s);

    gKstarPhiRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yKstarOff / yPhiOff);
    gKstarPhiRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yKstarOn / yPhiOn);
  }

  int graphPointsEPOSNewPhi = gdNdy_Phi_On->GetN();
  for (int i = 0; i < graphPointsEPOSNewPhi; i++)
  {
    double xK0s, yK0s, yKaon, yPion, yPhiOff, yPhiOn;
    double yK0sON, yKaonON, yPionON, yXiOff, yXiOn;

    gdNdy_Phi_Off->GetPoint(i, xK0s, yPhiOff);
    gdNdy_Phi_On->GetPoint(i, xK0s, yPhiOn);

    gdNdy_K0s_Off->GetPoint(i, xK0s, yK0s);
    gdNdy_Kaon_Off->GetPoint(i, xK0s, yKaon);
    gdNdy_Pion_Off->GetPoint(i, xK0s, yPion);

    gdNdy_K0s_On->GetPoint(i, xK0s, yK0sON);
    gdNdy_Kaon_On->GetPoint(i, xK0s, yKaonON);
    gdNdy_Pion_On->GetPoint(i, xK0s, yPionON);
    gdNdy_Xi_Off->GetPoint(i, xK0s, yXiOff);
    gdNdy_Xi_On->GetPoint(i, xK0s, yXiOn);

    gPhiKaRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yPhiOff / yKaon);
    gPhiKaRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yPhiOn / yKaon);
    gPhiPiRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yPhiOff / yPion);
    gPhiPiRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yPhiOn / yPion);
    gPhiK0sRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yPhiOff / yK0s);
    gPhiK0sRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yPhiOn / yK0s);
    gPhiXiRatio_EPOSNew_UrQMDOff->SetPoint(i, xK0s, yPhiOff / yXiOff);
    gPhiXiRatio_EPOSNew_UrQMDOn->SetPoint(i, xK0s, yPhiOn / yXiOff);
  }

  TGraphErrors *gKstarKaRatio_UrQMDOffSmooth = SmoothGraph(gKstarKaRatio_EPOSNew_UrQMDOff);
  TGraphErrors *gKstarKaRatio_UrQMDOnSmooth = SmoothGraph(gKstarKaRatio_EPOSNew_UrQMDOn);
  TGraphErrors *gKstarPiRatio_UrQMDOffSmooth = SmoothGraph(gKstarPiRatio_EPOSNew_UrQMDOff);
  TGraphErrors *gKstarPiRatio_UrQMDOnSmooth = SmoothGraph(gKstarPiRatio_EPOSNew_UrQMDOn);
  TGraphErrors *gKstarK0sRatio_UrQMDOffSmooth = SmoothGraph(gKstarK0sRatio_EPOSNew_UrQMDOff);
  TGraphErrors *gKstarK0sRatio_UrQMDOnSmooth = SmoothGraph(gKstarK0sRatio_EPOSNew_UrQMDOn);

  TGraphErrors *gPhiKaRatio_UrQMDOffSmooth = SmoothGraph(gPhiKaRatio_EPOSNew_UrQMDOff);
  TGraphErrors *gPhiKaRatio_UrQMDOnSmooth = SmoothGraph(gPhiKaRatio_EPOSNew_UrQMDOn);
  TGraphErrors *gPhiPiRatio_UrQMDOffSmooth = SmoothGraph(gPhiPiRatio_EPOSNew_UrQMDOff);
  TGraphErrors *gPhiPiRatio_UrQMDOnSmooth = SmoothGraph(gPhiPiRatio_EPOSNew_UrQMDOn);
  TGraphErrors *gPhiK0sRatio_UrQMDOffSmooth = SmoothGraph(gPhiK0sRatio_EPOSNew_UrQMDOff);
  TGraphErrors *gPhiK0sRatio_UrQMDOnSmooth = SmoothGraph(gPhiK0sRatio_EPOSNew_UrQMDOn);
  TGraphErrors *gPhiXiRatio_UrQMDOffSmooth = SmoothGraph(gPhiXiRatio_EPOSNew_UrQMDOff);
  TGraphErrors *gPhiXiRatio_UrQMDOnSmooth = SmoothGraph(gPhiXiRatio_EPOSNew_UrQMDOn);

  TGraphErrors *gKstarPhiRatio_UrQMDOffSmooth = SmoothGraph(gKstarPhiRatio_EPOSNew_UrQMDOff);
  TGraphErrors *gKstarPhiRatio_UrQMDOnSmooth = SmoothGraph(gKstarPhiRatio_EPOSNew_UrQMDOn);

  // Low to high multiplicity suppression in EPOS UrQMD Off and On
  double EPOSHM_Off, EPOSLM_Off, EPOSHM_On, EPOSLM_On, xLM2, xHM2;
  gKstarK0sRatio_UrQMDOffSmooth->GetPoint(0, xLM2, EPOSLM_Off);
  gKstarK0sRatio_UrQMDOffSmooth->GetPoint(gKstarK0sRatio_UrQMDOffSmooth->GetN() - 1, xHM2, EPOSHM_Off);
  gKstarK0sRatio_UrQMDOnSmooth->GetPoint(0, xLM2, EPOSLM_On);
  gKstarK0sRatio_UrQMDOnSmooth->GetPoint(gKstarK0sRatio_UrQMDOnSmooth->GetN() - 1, xHM2, EPOSHM_On);
  cout << "EPOS UrQMD Off K*0/K0s ratio suppression: " << EPOSHM_Off / EPOSLM_Off << endl;
  cout << "EPOS UrQMD On K*0/K0s ratio suppression: " << EPOSHM_On / EPOSLM_On << endl;

  TGraphErrors **gRatioKstarKaon = MakeRatioUncorr(gMYieldKstar, gMYieldKaon, 1.0, 3);

  CalculateSuppressionSignificance(gRatioKstarKaon, "K*0/Kaon", true, true, false);

  TGraphErrors **gRatioKstarPhi = MakeRatioUncorr(gMYieldKstar, gMYieldPhiRun3, 1.0, 3);

  CalculateSuppressionSignificance(gRatioKstarPhi, "K*0/Phi", false, true, false);

  TGraphErrors **gRatioPhiK = MakeRatioUncorr(gMYieldPhiRun3, gMYieldKaon, 1.0, 3);

  CalculateSuppressionSignificance(gRatioPhiK, "Phi/K", false, true, true);

  TGraphErrors **gRatioKstarPion = MakeRatioUncorr(gMYieldKstar, gMYieldPion, 1.0, 3);

  CalculateSuppressionSignificance(gRatioKstarPion, "K*0/Pion", false, true, true);

  TGraphErrors **gRatioPhiPion = MakeRatioUncorr(gMYieldPhiRun3, gMYieldPion, 1.0, 3);

  CalculateSuppressionSignificance(gRatioPhiPion, "Phi/Pion", false, true, true);

  TGraphErrors **gRatioPhiXi = MakeRatioUncorr(gMYieldPhiRun3, gMYieldXiRun3, 1.0, 3);
  ScaleGraph(gPhiXiRatio_EPOSNew_UrQMDOff, 0.5); // To average the Xi yield
  ScaleGraph(gPhiXiRatio_EPOSNew_UrQMDOn, 0.5);

  TGraphErrors **gRatioKstarKshort = MakeRatioUncorr(gMYieldKstar, gMYieldKshortRun3, 1.0, 3);

  CalculateSuppressionSignificance(gRatioKstarKshort, "K*0/Kshort", true, false, false);

  TGraphErrors **gRatioPhiK0s = MakeRatioUncorr(gMYieldPhiRun3, gMYieldKshortRun3, 1.0, 3);

  CalculateSuppressionSignificance(gRatioPhiK0s, "Phi/K0s", true, false, true);

  //======================================================================
  //  Double Yield Ratio (Mult/LM): Kstar/K, Kstar/Pi
  //======================================================================

  TGraphErrors **gYieldKstarKaLMRatio = DivideByMult(gRatioKstarKaon, -1, 0.5, 3);
  TGraphErrors **gYieldKstarPiLMRatio = DivideByMult(gRatioKstarPion, -1, 0.5, 3);

  for (auto model : modelsToPlot)
  {
    ScaleGraph(gMYield[model][kKstar], 0.5);
    gMYield[model][kKstar]->Draw("l same");
  }

  // // ======================================================
  // //       FOUR-PANEL PARTICLE-RATIO COMPARISON
  // //
  // //       (a) K*0 / K       (b) Phi / K
  // //       (c) K*0 / pi      (d) Phi / pi
  // // ======================================================

  // TCanvas *cParticleRatios = new TCanvas("cParticleRatios", "Particle ratios", 900, 900);
  // SetPad(cParticleRatios);
  // cParticleRatios->Divide(2, 2, 0.0, 0.0);

  // // Pad 1: Kstar / K
  // cParticleRatios->cd(1);

  // gPad->SetLeftMargin(0.16);
  // gPad->SetRightMargin(0.005);
  // gPad->SetTopMargin(0.08);
  // gPad->SetBottomMargin(0.0);
  // gPad->SetTicks(1, 1);

  // // Pad 2: Phi / K

  // cParticleRatios->cd(2);

  // gPad->SetLeftMargin(0.005);
  // gPad->SetRightMargin(0.16);
  // gPad->SetTopMargin(0.08);
  // gPad->SetBottomMargin(0.0);
  // gPad->SetTicks(1, 1);

  // // Pad 3: Kstar / pi

  // cParticleRatios->cd(3);

  // gPad->SetLeftMargin(0.16);
  // gPad->SetRightMargin(0.005);
  // gPad->SetTopMargin(0.00);
  // gPad->SetBottomMargin(0.15);
  // gPad->SetTicks(1, 1);

  // // Pad 4: Phi / pi

  // cParticleRatios->cd(4);

  // gPad->SetLeftMargin(0.005);
  // gPad->SetRightMargin(0.16);
  // gPad->SetTopMargin(0.0);
  // gPad->SetBottomMargin(0.15);
  // gPad->SetTicks(1, 1);

  // // ------------------------------------------------------
  // // Helper to draw one ratio panel
  // // ------------------------------------------------------

  // auto drawRatioPanel = [&](int pad,
  //                           const char *frameName,
  //                           const char *panelLabel,
  //                           TGraphErrors **data,
  //                           TGraphErrors *eposOff,
  //                           TGraphErrors *eposOn,
  //                           int numerator,
  //                           int denominator,
  //                           double modelScale,
  //                           double yMin,
  //                           double yMax,
  //                           bool rightAxis,
  //                           bool showXaxis,
  //                           bool showYaxis)
  // {
  //   cParticleRatios->cd(pad);

  //   // Frame
  //   TH1D *hFrame = new TH1D(frameName, "", 100, 1.8, 23.5);
  //   SetHistStyle(hFrame);

  //   if (!showXaxis)
  //   {
  //     hFrame->GetXaxis()->SetTitleSize(0.0);
  //   }
  //   if (!showYaxis)
  //   {
  //     hFrame->GetYaxis()->SetTitleSize(0.0);
  //   }

  //   hFrame->GetXaxis()->SetTitle("#LT d#it{N_{#it{ch}}/#it{d}#eta} #GT_{|#it{#eta}|<0.5}");
  //   hFrame->GetYaxis()->SetTitle("d#it{N}/d#it{y} ratio");
  //   hFrame->GetXaxis()->SetLimits(1.8, 23.5);
  //   hFrame->GetYaxis()->SetRangeUser(yMin, yMax);
  //   hFrame->GetYaxis()->SetMaxDigits(3);
  //   hFrame->GetYaxis()->SetTitleOffset(rightAxis ? 1.25 : 1.35);

  //   // Draw the right-hand y axis for the right column
  //   if (rightAxis)
  //     hFrame->Draw("Y+");
  //   else
  //     hFrame->Draw();

  //   // --------------------------------------------------
  //   // Run 3 data
  //   // --------------------------------------------------

  //   SetGraphStyle(data[0], colors[5], 1.0, 20, 1.0, 1, 1);
  //   data[0]->Draw("PE SAME");

  //   // Total systematic uncertainty
  //   SetGraphStyle(data[1], colors[5], 1.0, 20, 1.0, 1, 1);
  //   data[1]->SetFillStyle(0);
  //   data[1]->Draw("5 SAME");

  //   // Uncorrelated systematic uncertainty
  //   SetGraphStyle(data[2], colors[5], 1.0, 20, 1.0, 1, 1);
  //   data[2]->SetFillStyle(1001);
  //   data[2]->SetFillColorAlpha(colors[5], 0.35);
  //   data[2]->Draw("E2 SAME");

  //   // --------------------------------------------------
  //   // EPOS4 models
  //   // --------------------------------------------------

  //   SetGraphStyle(eposOff, models[kEPOSWithoutUrQMD].color, 1.0, 1, 1.0, models[kEPOSWithoutUrQMD].width, models[kEPOSWithoutUrQMD].style);
  //   eposOff->Draw("LX SAME");

  //   SetGraphStyle(eposOn, models[kEPOSWithUrQMD].color, 1.0, 1, 1.0, models[kEPOSWithUrQMD].width, models[kEPOSWithUrQMD].style);
  //   eposOn->Draw("LX SAME");

  //   // --------------------------------------------------
  //   // PYTHIA models
  //   // --------------------------------------------------

  //   std::vector<TGraphErrors *> modelGraphs;

  //   for (auto model : modelsToPlot)
  //   {
  //     TGraphErrors *gRatioModel = MakeRatio(
  //         gMYield[model][numerator],
  //         gMYield[model][denominator],
  //         true);

  //     ScaleGraph(gRatioModel, modelScale);

  //     SetGraphStyle(gRatioModel, models[model].color, 1.0, 1, 1.0, models[model].width, models[model].style);
  //     gRatioModel->Draw("L SAME");

  //     modelGraphs.push_back(gRatioModel);
  //   }

  //   return modelGraphs;
  // };

  // // ------------------------------------------------------
  // // Draw all four panels
  // // ------------------------------------------------------

  // // (a) K*0 / K -- top left
  // std::vector<TGraphErrors *> modelGraphsKstarKaon =
  //     drawRatioPanel(
  //         1,
  //         "hFrameFourRatioKstarKaon",
  //         "K*^{0} / K",
  //         gRatioKstarKaon,
  //         gKstarKaRatio_UrQMDOffSmooth,
  //         gKstarKaRatio_UrQMDOnSmooth,
  //         kKstar,
  //         kKaon,
  //         2.0,
  //         0.235, 0.44,
  //         false,
  //         false,
  //         true);

  // // (b) Phi / K -- top right
  // std::vector<TGraphErrors *> modelGraphsPhiKaon =
  //     drawRatioPanel(
  //         2,
  //         "hFrameFourRatioPhiKaon",
  //         "#phi / K",
  //         gRatioPhiK,
  //         gPhiKaRatio_UrQMDOffSmooth,
  //         gPhiKaRatio_UrQMDOnSmooth,
  //         kPhi,
  //         kKaon,
  //         2.0,
  //         0.07, 0.17,
  //         true,
  //         false,
  //         false);

  // // (c) K*0 / pi -- bottom left
  // std::vector<TGraphErrors *> modelGraphsKstarPion =
  //     drawRatioPanel(
  //         3,
  //         "hFrameFourRatioKstarPion",
  //         "K*^{0} / #pi",
  //         gRatioKstarPion,
  //         gKstarPiRatio_UrQMDOffSmooth,
  //         gKstarPiRatio_UrQMDOnSmooth,
  //         kKstar,
  //         kPion,
  //         2.0,
  //         0.031, 0.062,
  //         false,
  //         false,
  //         false);

  // // (d) Phi / pi -- bottom right
  // std::vector<TGraphErrors *> modelGraphsPhiPion =
  //     drawRatioPanel(
  //         4,
  //         "hFrameFourRatioPhiPion",
  //         "#phi / #pi",
  //         gRatioPhiPion,
  //         gPhiPiRatio_UrQMDOffSmooth,
  //         gPhiPiRatio_UrQMDOnSmooth,
  //         kPhi,
  //         kPion,
  //         2.0,
  //         0.007, 0.023,
  //         true,
  //         true,
  //         false);

  // // --------------------------------------------------
  // // Panel labels
  // // --------------------------------------------------

  // cParticleRatios->cd(1);

  // TLatex panelLatex;
  // panelLatex.SetNDC();
  // panelLatex.SetTextFont(42);
  // panelLatex.SetTextSize(0.045);
  // panelLatex.DrawLatex(0.92, 0.85, "(a)");

  // panelLatex.SetTextSize(0.05);
  // panelLatex.DrawLatex(0.78, 0.845, "K*^{0} / K");

  // panelLatex.SetTextSize(0.065);
  // panelLatex.DrawLatex(0.22, 0.82, "#bf{ALICE}");

  // TLegend *legendRun3 = new TLegend(0.2, 0.74, 0.54, 0.8);
  // SetLegend(legendRun3);
  // legendRun3->AddEntry(gRatioKstarKaon[0], "pp, #sqrt{s} = 13.6 TeV, |y| < 0.5", "p");
  // legendRun3->Draw();

  // // ------------------------------------------------------
  // // Common model legend
  // // ------------------------------------------------------

  // cParticleRatios->cd(2);

  // TLegend *legendFourRatios = new TLegend(0.05, 0.53, 0.54, 0.89);

  // SetLegend(legendFourRatios);
  // legendFourRatios->SetNColumns(1);

  // legendFourRatios->AddEntry(gKstarKaRatio_UrQMDOffSmooth, models[kEPOSWithoutUrQMD].label, "l");
  // legendFourRatios->AddEntry(gKstarKaRatio_UrQMDOnSmooth, models[kEPOSWithUrQMD].label, "l");

  // for (int i = 0; i < (int)modelGraphsKstarKaon.size(); i++)
  // {
  //   auto model = modelsToPlot[i];
  //   legendFourRatios->AddEntry(modelGraphsKstarKaon[i], models[model].label, "l");
  // }

  // legendFourRatios->Draw();

  // panelLatex.SetTextSize(0.045);
  // panelLatex.DrawLatex(0.76, 0.85, "(b)");

  // panelLatex.SetTextSize(0.05);
  // panelLatex.DrawLatex(0.65, 0.85, "#phi / K");

  // cParticleRatios->cd(3);
  // panelLatex.SetTextSize(0.045);
  // panelLatex.DrawLatex(0.92, 0.92, "(c)");

  // panelLatex.SetTextSize(0.05);
  // panelLatex.DrawLatex(0.78, 0.92, "K*^{0} / #pi");

  // cParticleRatios->cd(4);
  // panelLatex.SetTextSize(0.045);
  // panelLatex.DrawLatex(0.76, 0.92, "(d)");

  // panelLatex.SetTextSize(0.05);
  // panelLatex.DrawLatex(0.65, 0.92, "#phi / #pi");

  // if (isSavePlots)
  // {
  //   cParticleRatios->SaveAs(Form("%s/ParticleRatios_Combined_Run3.pdf", plotOutputDirectory.c_str()));
  // }

  //======================================================
  // COMBINED K*0 AND PHI OBSERVABLE PLOTS
  //======================================================
  // These panels reuse the measured, EPOS, and Pythia graphs already loaded
  // and drawn above. No additional physics scaling is applied here.

  auto makeHorizontalPad = [](const char *name, const char *title,
                              double xMin, double xMax, bool isLeft)
  {
    TPad *pad = new TPad(name, title, xMin, 0.0, xMax, 1.0);
    pad->SetLeftMargin(isLeft ? 0.18 : 0.005);
    pad->SetRightMargin(isLeft ? 0.005 : 0.18);
    pad->SetBottomMargin(0.15);
    pad->SetTopMargin(0.05);
    pad->SetTickx(1);
    pad->SetTicky(1);
    pad->Draw();
    return pad;
  };

  struct PanelAnnotations
  {
    double particleLabelX;
    double particleLabelY;
    double panelLabelX;
    double panelLabelY;
    int legendKind; // 0: none, 1: data, 2: models
    double legendX1;
    double legendY1;
    double legendX2;
    double legendY2;
  };

  // // Annotation coordinates are specified per pad in NDC, so titles, panel
  // // letters, and legends can be positioned independently in each layout.
  // auto drawObservablePanel = [&](TVirtualPad *pad,
  //                                const char *frameName,
  //                                TGraphErrors **run3Data,
  //                                TGraphErrors **run2Data,
  //                                TGraphErrors *eposOff,
  //                                TGraphErrors *eposOn,
  //                                int modelParticle,
  //                                bool isMeanPt,
  //                                double yMin,
  //                                double yMax,
  //                                bool rightAxis,
  //                                bool showXTitle,
  //                                const char *panelLabel,
  //                                const char *particleLabel,
  //                                const PanelAnnotations &annotations)
  // {
  //   pad->cd();

  //   TH1D *frame = new TH1D(frameName, "", 100, 2.0, 28.0);
  //   SetHistStyle(frame);
  //   frame->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
  //   frame->GetYaxis()->SetTitleSize(0.065);
  //   frame->GetYaxis()->SetLabelSize(0.055);
  //   frame->GetYaxis()->SetTitle(isMeanPt? "#LT #it{p}_{T} #GT (GeV/#it{c})": "dN/dy");
  //   frame->GetXaxis()->SetTitleSize(0.065);
  //   frame->GetXaxis()->SetLabelSize(0.055);
  //   frame->GetXaxis()->SetLimits(2.0, 28.0);
  //   frame->GetYaxis()->SetRangeUser(yMin, yMax);

  //   frame->GetYaxis()->SetMaxDigits(3);
  //   frame->GetYaxis()->SetTitleOffset(rightAxis ? 1.25 : 1.35);
  //   if (!showXTitle)
  //     frame->GetXaxis()->SetTitleSize(0.0);
  //   if (rightAxis)
  //   {
  //     frame->Draw("Y+");
  //   }
  //   else
  //     frame->Draw();

  //   // Clone the shared input graphs so this panel's styling does not
  //   // change the graphs used by the earlier plots in the macro.
  //   TGraphErrors *run3Stat = (TGraphErrors *)run3Data[0]->Clone(Form("%s_Run3Stat", frameName));
  //   TGraphErrors *run3Sys = (TGraphErrors *)run3Data[1]->Clone(Form("%s_Run3Sys", frameName));
  //   TGraphErrors *run3Uncorr = (TGraphErrors *)run3Data[2]->Clone(Form("%s_Run3Uncorr", frameName));
  //   TGraphErrors *run2Stat = (TGraphErrors *)run2Data[0]->Clone(Form("%s_Run2Stat", frameName));
  //   TGraphErrors *run2Sys = (TGraphErrors *)run2Data[1]->Clone(Form("%s_Run2Sys", frameName));
  //   TGraphErrors *run2Uncorr = (TGraphErrors *)run2Data[2]->Clone(Form("%s_Run2Uncorr", frameName));

  //   // Run 3 data: statistical, total systematic, and uncorrelated systematic.
  //   SetGraphStyle(run3Stat, Run3Color, 1.0, 20, 0.8, 1);
  //   run3Stat->Draw("PE SAME");
  //   SetGraphStyle(run3Sys, Run3Color, 1.0, 20, 0.8, 1);
  //   run3Sys->SetFillStyle(0);
  //   run3Sys->Draw("5 SAME");
  //   SetGraphStyle(run3Uncorr, Run3Color, 1.0, 20, 0.8, 1);
  //   run3Uncorr->SetFillStyle(1001);
  //   run3Uncorr->SetFillColorAlpha(Run3Color, 0.30);
  //   run3Uncorr->Draw("E2 SAME");

  //   // Run 2 data.
  //   SetGraphStyle(run2Stat, Run2Color, 1.0, 21, 0.8, 1);
  //   run2Stat->Draw("PE SAME");
  //   SetGraphStyle(run2Sys, Run2Color, 1.0, 21, 0.8, 1);
  //   run2Sys->SetFillStyle(0);
  //   run2Sys->Draw("5 SAME");
  //   SetGraphStyle(run2Uncorr, Run2Color, 1.0, 21, 0.8, 1);
  //   run2Uncorr->SetFillStyle(1001);
  //   run2Uncorr->SetFillColorAlpha(Run2Color, 0.30);
  //   run2Uncorr->Draw("E2 SAME");

  //   // EPOS with and without UrQMD.
  //   TGraphErrors *eposOffCopy = (TGraphErrors *)eposOff->Clone(Form("%s_EPOSOff", frameName));
  //   TGraphErrors *eposOnCopy = (TGraphErrors *)eposOn->Clone(Form("%s_EPOSOn", frameName));
  //   SetGraphStyle(eposOffCopy, models[kEPOSWithoutUrQMD].color, 1.0, 1, 1.0, models[kEPOSWithoutUrQMD].width, models[kEPOSWithoutUrQMD].style);
  //   eposOffCopy->Draw("L SAME");
  //   SetGraphStyle(eposOnCopy, models[kEPOSWithUrQMD].color, 1.0, 1, 1.0, models[kEPOSWithUrQMD].width, models[kEPOSWithUrQMD].style);
  //   eposOnCopy->Draw("L SAME");

  //   // Pythia model curves. The K* yield graphs have already been scaled
  //   // for particle/antiparticle averaging by the original yield plot.
  //   std::vector<TGraphErrors *> drawnModelGraphs;
  //   for (auto model : modelsToPlot)
  //   {
  //     TGraphErrors *modelSource = isMeanPt
  //                                     ? gMeanPt[model][modelParticle]
  //                                     : gMYield[model][modelParticle];
  //     TGraphErrors *modelGraph = (TGraphErrors *)modelSource->Clone(Form("%s_Model%d", frameName, (int)model));
  //     SetGraphStyle(modelGraph, models[model].color, 1.0, 1, 0.0, models[model].width, models[model].style);
  //     modelGraph->Draw("LX SAME");
  //     drawnModelGraphs.push_back(modelGraph);
  //   }

  //   TLatex panelText;
  //   panelText.SetNDC();
  //   panelText.SetTextFont(42);
  //   panelText.SetTextSize(0.045);
  //   panelText.DrawLatex(annotations.panelLabelX, annotations.panelLabelY, panelLabel);
  //   panelText.SetTextSize(0.055);
  //   panelText.DrawLatex(annotations.particleLabelX, annotations.particleLabelY, particleLabel);

  //   if (annotations.legendKind == 1)
  //   {
  //     TLegend *dataLegend = new TLegend(annotations.legendX1, annotations.legendY1,
  //                                       annotations.legendX2, annotations.legendY2);
  //     SetLegend(dataLegend);
  //     dataLegend->SetTextAlign(12);
  //     dataLegend->AddEntry(run3Stat, "pp, #sqrt{s} = 13.6 TeV", "p");
  //     dataLegend->AddEntry(run2Stat, "pp, #sqrt{s} = 13 TeV", "p");
  //     dataLegend->Draw();
  //   }
  //   else if (annotations.legendKind == 2)
  //   {
  //     TLegend *modelLegend = new TLegend(annotations.legendX1, annotations.legendY1,
  //                                        annotations.legendX2, annotations.legendY2);
  //     SetLegend(modelLegend);
  //     modelLegend->AddEntry(eposOffCopy, models[kEPOSWithoutUrQMD].label, "l");
  //     modelLegend->AddEntry(eposOnCopy, models[kEPOSWithUrQMD].label, "l");
  //     for (size_t iModel = 0; iModel < modelsToPlot.size(); ++iModel)
  //     {
  //       const auto model = modelsToPlot[iModel];
  //       modelLegend->AddEntry(drawnModelGraphs[iModel], models[model].label, "l");
  //     }
  //     modelLegend->Draw();
  //   }
  // };

  // TCanvas *cMeanPtKstarPhi = new TCanvas("cMeanPtKstarPhi", "Mean pT: K*0 and phi", 850, 450);
  // SetPad(cMeanPtKstarPhi);
  // cMeanPtKstarPhi->cd();
  // TPad *padMeanPtKstar = makeHorizontalPad("padMeanPtKstar", "K*0", 0.0, 0.5, true);
  // TPad *padMeanPtPhi = makeHorizontalPad("padMeanPtPhi", "phi", 0.5, 1.0, false);

  // const PanelAnnotations meanPtKstarAnnotations{
  //     0.27, // left x for particle label
  //     0.87, // left y for particle label
  //     0.2,  // left x for panel label
  //     0.88, // left y for panel label
  //     1,    // legend kind: 1 = data, 2 = models
  //     0.2,  // x1 for legend
  //     0.7,  // y1 for legend
  //     0.5,  // x2 for legend
  //     0.85  // y2 for legend
  // };

  // const PanelAnnotations meanPtPhiAnnotations{0.11, 0.88, 0.04, 0.88, 2, 0.4, 0.18, 0.92, 0.48};

  // drawObservablePanel(padMeanPtKstar, "hFrameMeanPtKstarPhi_Kstar",
  //                     gMPtKstar, gMPtKstarRun2,
  //                     gEPOS_KstarMeanpT_UrQMDOff, gEPOS_KstarMeanpT_UrQMDOn,
  //                     kKstar, true, 0.58, 1.76, false, true, "(a)", "K*^{0}", meanPtKstarAnnotations);
  // drawObservablePanel(padMeanPtPhi, "hFrameMeanPtKstarPhi_Phi",
  //                     gMPtPhiRun3, gMPtPhiRun2,
  //                     gEPOS_PhiMeanpT_UrQMDOff, gEPOS_PhiMeanpT_UrQMDOn,
  //                     kPhi, true, 0.58, 1.76, true, true, "(b)", "#phi", meanPtPhiAnnotations);

  // // K*0 and phi dN/dy side-by-side.
  // TCanvas *cYieldKstarPhi = new TCanvas("cYieldKstarPhi", "dN/dy: K*0 and phi", 850, 450);
  // SetPad(cYieldKstarPhi);
  // cYieldKstarPhi->cd();

  // if (isSavePlots)
  // {
  //   cMeanPtKstarPhi->SaveAs(Form("%s/MeanPt_Kstar_Phi_Combined_Run3.pdf", plotOutputDirectory.c_str()));
  // }

  // TPad *padYieldKstar = makeHorizontalPad("padYieldKstar", "K*0", 0.0, 0.5, true);
  // TPad *padYieldPhi = makeHorizontalPad("padYieldPhi", "phi", 0.5, 1.0, false);
  // const PanelAnnotations yieldKstarAnnotations{
  //     0.27, // left x for particle label
  //     0.87, // left y for particle label
  //     0.2,  // left x for panel label
  //     0.88, // left y for panel label
  //     1,    // legend kind: 1 = data, 2 = models
  //     0.2,  // x1 for legend
  //     0.7,  // y1 for legend
  //     0.5,  // x2 for legend
  //     0.85  // y2 for legend
  // };
  // const PanelAnnotations yieldPhiAnnotations{0.11, 0.88, 0.04, 0.88, 2, 0.04, 0.55, 0.5, 0.85};
  // drawObservablePanel(padYieldKstar, "hFrameYieldKstarPhi_Kstar",
  //                     gMYieldKstar, gMYieldKstarRun2,
  //                     gEPOS_KstardNdy_UrQMDOff, gEPOS_KstardNdy_UrQMDOn,
  //                     kKstar, false, 0.0, 0.71, false, true, "(a)", "K*^{0}", yieldKstarAnnotations);
  // drawObservablePanel(padYieldPhi, "hFrameYieldKstarPhi_Phi",
  //                     gMYieldPhiRun3, gMYieldPhiRun2,
  //                     gEPOS_PhidNdy_UrQMDOff, gEPOS_PhidNdy_UrQMDOn,
  //                     kPhi, false, 0.0, 0.31, true, true, "(b)", "#phi", yieldPhiAnnotations);

  // ======================================================
  // Four-panel mean pT and dN/dy: K*0 and phi
  //       (a) top left: K*0 mean pT
  //       (b) top right: phi mean pT
  //       (c) bottom left: K*0 dN/dy
  //       (d) bottom right: phi dN/dy
  // ======================================================

  // Annotation coordinates are specified per pad in NDC, so titles, panel
  // letters, and legends can be positioned independently in each layout.
  auto drawObservablePanel = [&](TVirtualPad *pad,
                                 const char *frameName,
                                 TGraphErrors **run3Data,
                                 TGraphErrors **run2Data,
                                 TGraphErrors *eposOff,
                                 TGraphErrors *eposOn,
                                 int modelParticle,
                                 bool isMeanPt,
                                 double yMin,
                                 double yMax,
                                 bool rightAxis,
                                 bool showXTitle,
                                 bool left,
                                 bool top)
  {
    pad->cd();

    TH1D *frame = new TH1D(frameName, "", 100, 2.0, 28.0);
    SetHistStyle(frame);
    frame->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    frame->GetYaxis()->SetTitleSize(top ? 0.073 : 0.068);
    frame->GetYaxis()->SetLabelSize(top ? 0.065 : 0.055);
    frame->GetYaxis()->SetTitle(isMeanPt ? "#LT #it{p}_{T} #GT (GeV/#it{c})" : "dN/dy");
    frame->GetXaxis()->SetTitleSize(0.065);
    frame->GetXaxis()->SetLabelSize(0.055);
    frame->GetXaxis()->SetLimits(2.0, 28.0);
    frame->GetYaxis()->SetRangeUser(yMin, yMax);

    frame->GetYaxis()->SetMaxDigits(3);
    if (top)
      frame->GetYaxis()->SetTitleOffset(rightAxis ? 1.22 : 1.32);
    else
      frame->GetYaxis()->SetTitleOffset(rightAxis ? 1.3 : 1.4);
    if (!showXTitle)
      frame->GetXaxis()->SetTitleSize(0.0);
    if (rightAxis)
    {
      frame->Draw("Y+");
    }
    else
      frame->Draw();

    // Clone the shared input graphs so this panel's styling does not
    // change the graphs used by the earlier plots in the macro.
    TGraphErrors *run3Stat = (TGraphErrors *)run3Data[0]->Clone(Form("%s_Run3Stat", frameName));
    TGraphErrors *run3Sys = (TGraphErrors *)run3Data[1]->Clone(Form("%s_Run3Sys", frameName));
    TGraphErrors *run3Uncorr = (TGraphErrors *)run3Data[2]->Clone(Form("%s_Run3Uncorr", frameName));
    TGraphErrors *run2Stat = (TGraphErrors *)run2Data[0]->Clone(Form("%s_Run2Stat", frameName));
    TGraphErrors *run2Sys = (TGraphErrors *)run2Data[1]->Clone(Form("%s_Run2Sys", frameName));
    TGraphErrors *run2Uncorr = (TGraphErrors *)run2Data[2]->Clone(Form("%s_Run2Uncorr", frameName));

    // Run 3 data: statistical, total systematic, and uncorrelated systematic.
    SetGraphStyle(run3Stat, Run3Color, 1.0, 20, 0.8, 1);
    run3Stat->Draw("PE SAME");
    SetGraphStyle(run3Sys, Run3Color, 1.0, 20, 0.8, 1);
    run3Sys->SetFillStyle(0);
    run3Sys->Draw("5 SAME");
    SetGraphStyle(run3Uncorr, Run3Color, 1.0, 20, 0.8, 1);
    run3Uncorr->SetFillStyle(1001);
    run3Uncorr->SetFillColorAlpha(Run3Color, 0.30);
    run3Uncorr->Draw("E2 SAME");

    // Run 2 data.
    SetGraphStyle(run2Stat, Run2Color, 1.0, 21, 0.8, 1);
    run2Stat->Draw("PE SAME");
    SetGraphStyle(run2Sys, Run2Color, 1.0, 21, 0.8, 1);
    run2Sys->SetFillStyle(0);
    run2Sys->Draw("5 SAME");
    SetGraphStyle(run2Uncorr, Run2Color, 1.0, 21, 0.8, 1);
    run2Uncorr->SetFillStyle(1001);
    run2Uncorr->SetFillColorAlpha(Run2Color, 0.30);
    run2Uncorr->Draw("E2 SAME");

    // EPOS with and without UrQMD.
    TGraphErrors *eposOffCopy = (TGraphErrors *)eposOff->Clone(Form("%s_EPOSOff", frameName));
    TGraphErrors *eposOnCopy = (TGraphErrors *)eposOn->Clone(Form("%s_EPOSOn", frameName));
    SetGraphStyle(eposOffCopy, models[kEPOSWithoutUrQMD].color, 1.0, 1, 1.0, models[kEPOSWithoutUrQMD].width, models[kEPOSWithoutUrQMD].style);
    eposOffCopy->Draw("L SAME");
    SetGraphStyle(eposOnCopy, models[kEPOSWithUrQMD].color, 1.0, 1, 1.0, models[kEPOSWithUrQMD].width, models[kEPOSWithUrQMD].style);
    eposOnCopy->Draw("L SAME");

    // Pythia model curves. The K* yield graphs have already been scaled
    // for particle/antiparticle averaging by the original yield plot.
    std::vector<TGraphErrors *> drawnModelGraphs;
    for (auto model : modelsToPlot)
    {
      TGraphErrors *modelSource = isMeanPt
                                      ? gMeanPt[model][modelParticle]
                                      : gMYield[model][modelParticle];
      TGraphErrors *modelGraph = (TGraphErrors *)modelSource->Clone(Form("%s_Model%d", frameName, (int)model));
      SetGraphStyle(modelGraph, models[model].color, 1.0, 1, 0.0, models[model].width, models[model].style);
      modelGraph->Draw("LX SAME");
      drawnModelGraphs.push_back(modelGraph);
    }

    if (left && top)
    {
      TLegend *dataLegend = new TLegend(0.23, 0.62, 0.5, 0.78);
      SetLegend(dataLegend);
      dataLegend->SetTextSize(0.05);
      dataLegend->AddEntry(run3Stat, "pp, #sqrt{s} = 13.6 TeV", "p");
      dataLegend->AddEntry(run2Stat, "pp, #sqrt{s} = 13 TeV", "p");
      dataLegend->Draw();

      TLegend *modelLegend = new TLegend(0.55, 0.05, 0.95, 0.2);
      SetLegend(modelLegend);
      modelLegend->SetTextSize(0.05);
      modelLegend->AddEntry(eposOffCopy, models[kEPOSWithoutUrQMD].label, "l");
      modelLegend->AddEntry(eposOnCopy, models[kEPOSWithUrQMD].label, "l");
      modelLegend->Draw();
    }
    if (!left && top)
    {
      TLegend *modelLegend = new TLegend(0.35, 0.05, 0.9, 0.3);
      SetLegend(modelLegend);
      modelLegend->SetTextSize(0.05);
      // modelLegend->AddEntry(eposOffCopy, models[kEPOSWithoutUrQMD].label, "l");
      // modelLegend->AddEntry(eposOnCopy, models[kEPOSWithUrQMD].label, "l");
      for (size_t iModel = 0; iModel < modelsToPlot.size(); ++iModel)
      {
        const auto model = modelsToPlot[iModel];
        modelLegend->AddEntry(drawnModelGraphs[iModel], models[model].label, "l");
      }
      modelLegend->Draw();
    }
  };

  TCanvas *cMeanPtYieldKstarPhi = new TCanvas("cMeanPtYieldKstarPhi", "Mean pT and dN/dy: K*0 and phi", 900, 860);
  SetPad(cMeanPtYieldKstarPhi);
  TPad *pMeanPtKstar = new TPad("pMeanPtKstar", "K*0 mean pT", 0.0, 0.55, 0.5, 1.0);
  pMeanPtKstar->SetLeftMargin(0.2);
  pMeanPtKstar->SetRightMargin(0.005);
  pMeanPtKstar->SetBottomMargin(0.0);
  pMeanPtKstar->SetTopMargin(0.07);
  pMeanPtKstar->SetTickx(1);
  pMeanPtKstar->SetTicky(1);
  pMeanPtKstar->Draw();

  TPad *pMeanPtPhi = new TPad("pMeanPtPhi", "phi mean pT", 0.5, 0.55, 1.0, 1.0);
  pMeanPtPhi->SetLeftMargin(0.005);
  pMeanPtPhi->SetRightMargin(0.2);
  pMeanPtPhi->SetBottomMargin(0.0);
  pMeanPtPhi->SetTopMargin(0.07);
  pMeanPtPhi->SetTickx(1);
  pMeanPtPhi->SetTicky(1);
  pMeanPtPhi->Draw();

  TPad *pYieldKstar = new TPad("pYieldKstar", "K*0 dN/dy", 0.0, 0.0, 0.5, 0.55);
  pYieldKstar->SetLeftMargin(0.2);
  pYieldKstar->SetRightMargin(0.005);
  pYieldKstar->SetBottomMargin(0.22);
  pYieldKstar->SetTopMargin(0.0);
  pYieldKstar->SetTickx(1);
  pYieldKstar->SetTicky(1);
  pYieldKstar->Draw();

  TPad *pYieldPhi = new TPad("pYieldPhi", "phi dN/dy", 0.5, 0.0, 1.0, 0.55);
  pYieldPhi->SetLeftMargin(0.005);
  pYieldPhi->SetRightMargin(0.2);
  pYieldPhi->SetBottomMargin(0.22);
  pYieldPhi->SetTopMargin(0.0);
  pYieldPhi->SetTickx(1);
  pYieldPhi->SetTicky(1);
  pYieldPhi->Draw();

  pMeanPtKstar->cd();
  drawObservablePanel(gPad, "hFrameMeanPtYieldKstarPhi_KstarPt",
                      gMPtKstar, gMPtKstarRun2,
                      gEPOS_KstarMeanpT_UrQMDOff, gEPOS_KstarMeanpT_UrQMDOn,
                      kKstar, true, 0.58, 1.76, false, false, true, true);

  TLatex panelA;
  panelA.SetNDC();
  panelA.SetTextFont(42);
  panelA.SetTextSize(0.08);
  panelA.DrawLatex(0.75, 0.81, "#bf{ALICE}");
  panelA.SetTextSize(0.065);
  panelA.DrawLatex(0.25, 0.82, "(a) K*^{0} ");

  pMeanPtPhi->cd();
  drawObservablePanel(gPad, "hFrameMeanPtYieldKstarPhi_PhiPt",
                      gMPtPhiRun3, gMPtPhiRun2,
                      gEPOS_PhiMeanpT_UrQMDOff, gEPOS_PhiMeanpT_UrQMDOn,
                      kPhi, true, 0.58, 1.76, true, false, false, true);

  TLatex panelB;
  panelB.SetNDC();
  panelB.SetTextFont(42);
  panelB.SetTextSize(0.065);
  panelB.DrawLatex(0.05, 0.82, "(b) #phi");

  pYieldKstar->cd();
  drawObservablePanel(gPad, "hFrameMeanPtYieldKstarPhi_KstarYield",
                      gMYieldKstar, gMYieldKstarRun2,
                      gEPOS_KstardNdy_UrQMDOff, gEPOS_KstardNdy_UrQMDOn,
                      kKstar, false, 0.0, 0.71, false, true, true, false);

  TLatex panelC;
  panelC.SetNDC();
  panelC.SetTextFont(42);
  panelC.SetTextSize(0.055);
  panelC.DrawLatex(0.25, 0.9, "(c) K*^{0} ");

  pYieldPhi->cd();
  drawObservablePanel(gPad, "hFrameMeanPtYieldKstarPhi_PhiYield",
                      gMYieldPhiRun3, gMYieldPhiRun2,
                      gEPOS_PhidNdy_UrQMDOff, gEPOS_PhidNdy_UrQMDOn,
                      kPhi, false, 0.0, 0.26, true, true, false, false);
  TLatex panelD;
  panelD.SetNDC();
  panelD.SetTextFont(42);
  panelD.SetTextSize(0.055);
  panelD.DrawLatex(0.05, 0.9, "(d) #phi");

  if (isSavePlots)
  {
    cMeanPtYieldKstarPhi->SaveAs(Form("%s/ParticleRatios_Observables4_Run3.pdf", plotOutputDirectory.c_str()));
  }

  return;

  //======================================================
  // 1. TWO-PANEL RATIO PLOT (2,1): K*0 / Phi and K*0 / K0s
  //======================================================
  TCanvas *cRatioKstarPhiK0s = new TCanvas("cRatioKstarPhiK0s", "Ratios: K*0/Phi and K*0/K0s", 850, 450);
  SetPad(cRatioKstarPhiK0s);
  cRatioKstarPhiK0s->cd();

  TPad *padRatioKstarPhi = makeHorizontalPad("padRatioKstarPhi", "K*0/Phi", 0.0, 0.5, true);
  TPad *padRatioKstarK0s = makeHorizontalPad("padRatioKstarK0s", "K*0/K0s", 0.5, 1.0, false);

  padRatioKstarPhi->SetLeftMargin(0.2);
  padRatioKstarPhi->SetTopMargin(0.03);
  padRatioKstarPhi->SetBottomMargin(0.17);
  padRatioKstarK0s->SetRightMargin(0.2);
  padRatioKstarK0s->SetTopMargin(0.03);
  padRatioKstarK0s->SetBottomMargin(0.17);

  auto drawTwoPanelRatio = [&](TVirtualPad *pad,
                               const char *frameName,
                               TGraphErrors **data,
                               TGraphErrors *eposOff,
                               TGraphErrors *eposOn,
                               TGraphErrors *thermal,
                               int numParticle,
                               int denParticle,
                               double modelScale,
                               double yMin,
                               double yMax,
                               bool rightAxis,
                               const char *panelLabel,
                               const char *ratioLabel,
                               bool showDataLegend,
                               bool showModelLegend)
  {
    pad->cd();

    TH1D *frame = new TH1D(frameName, "", 100, 1.8, 23.5);
    SetHistStyle(frame);
    frame->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    frame->GetXaxis()->SetTitleSize(0.065);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetTitle(ratioLabel);
    frame->GetYaxis()->SetTitleSize(0.065);
    frame->GetYaxis()->SetLabelSize(0.05);

    frame->GetXaxis()->SetLimits(1.8, 23.5);
    frame->GetYaxis()->SetRangeUser(yMin, yMax);
    frame->GetYaxis()->SetMaxDigits(3);
    frame->GetYaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetTitleOffset(rightAxis ? 1.35 : 1.25);

    if (rightAxis)
      frame->Draw("Y+");
    else
      frame->Draw();

    // Run 3 Data (cloned to keep original styling intact)
    TGraphErrors *dStat = (TGraphErrors *)data[0]->Clone(Form("%s_stat", frameName));
    TGraphErrors *dSys = (TGraphErrors *)data[1]->Clone(Form("%s_sys", frameName));
    TGraphErrors *dUncorr = (TGraphErrors *)data[2]->Clone(Form("%s_uncorr", frameName));

    SetGraphStyle(dStat, kRed, 1.0, 20, 0.8, 1);
    dStat->Draw("PE SAME");
    SetGraphStyle(dSys, kRed, 1.0, 20, 0.8, 1);
    dSys->SetFillStyle(0);
    dSys->Draw("5 SAME");
    SetGraphStyle(dUncorr, kRed, 1.0, 20, 0.8, 1);
    dUncorr->SetFillStyle(1001);
    dUncorr->SetFillColorAlpha(kRed, 0.35);
    dUncorr->Draw("E2 SAME");
    // dStat->Draw("PE SAME"); // Draw again to ensure points are on top

    // EPOS4 & Thermal-FIST
    TGraphErrors *eposOffCopy = (TGraphErrors *)eposOff->Clone(Form("%s_eposOff", frameName));
    TGraphErrors *eposOnCopy = (TGraphErrors *)eposOn->Clone(Form("%s_eposOn", frameName));
    TGraphErrors *thermalCopy = (TGraphErrors *)thermal->Clone(Form("%s_thermal", frameName));

    SetGraphStyle(eposOffCopy, models[kEPOSWithoutUrQMD].color, 1.0, 1, 1.0, models[kEPOSWithoutUrQMD].width, models[kEPOSWithoutUrQMD].style);
    eposOffCopy->Draw("LX SAME");
    SetGraphStyle(eposOnCopy, models[kEPOSWithUrQMD].color, 1.0, 1, 1.0, models[kEPOSWithUrQMD].width, models[kEPOSWithUrQMD].style);
    eposOnCopy->Draw("LX SAME");

    thermalCopy->SetLineColor(models[kThermalModel].color);
    thermalCopy->SetLineStyle(models[kThermalModel].style);
    thermalCopy->SetLineWidth(models[kThermalModel].width);
    thermalCopy->Draw("LX SAME");

    // Pythia models
    std::vector<TGraphErrors *> pythiaGraphs;
    std::vector<int> pythiaModelIndices;
    for (auto model : modelsToPlot)
    {
      TGraphErrors *gMod = MakeRatio(gMYield[model][numParticle], gMYield[model][denParticle], true);
      ScaleGraph(gMod, modelScale);
      SetGraphStyle(gMod, models[model].color, 1.0, 1, 1.0, models[model].width, models[model].style);
      gMod->Draw("L SAME");
      pythiaGraphs.push_back(gMod);
      pythiaModelIndices.push_back(model);
    }

    TLatex txt;
    txt.SetNDC();
    txt.SetTextFont(42);
    txt.SetTextSize(0.045);
    txt.DrawLatex(rightAxis ? 0.05 : 0.22, 0.87, panelLabel);
    txt.SetTextSize(0.052);

    if (showDataLegend)
    {
      txt.SetTextSize(0.065);
      txt.DrawLatex(0.25, 0.88, "#bf{ALICE}");
      txt.SetTextSize(0.05);
      txt.DrawLatex(0.25, 0.81, "pp, #sqrt{s} = 13.6 TeV, |y| < 0.5");
    }

    if (showModelLegend)
    {
      // Create a legend for the model graphs
      TLegend *legMod = new TLegend(0.3, 0.59, 0.90, 0.94);
      SetLegend(legMod);
      legMod->SetTextSize(0.04);
      legMod->AddEntry(eposOffCopy, models[kEPOSWithoutUrQMD].label, "l");
      legMod->AddEntry(eposOnCopy, models[kEPOSWithUrQMD].label, "l");
      legMod->AddEntry(thermalCopy, models[kThermalModel].label, "l");
      for (size_t i = 0; i < pythiaGraphs.size(); ++i)
      {
        legMod->AddEntry(pythiaGraphs[i], models[pythiaModelIndices[i]].label, "l");
      }
      legMod->Draw();
    }
  };

  // Panel (a): K*0 / Phi (already averaged in both numerator & denominator -> scale = 1.0)
  drawTwoPanelRatio(padRatioKstarPhi, "hFrame2Panel_KstarPhi",
                    gRatioKstarPhi,
                    gKstarPhiRatio_EPOSNew_UrQMDOff, gKstarPhiRatio_EPOSNew_UrQMDOn, gRatioKstarPhi_Thermal,
                    kKstar, kPhi, 1.0,
                    2.1, 5.15, false, "", "K*^{0} / #phi", true, false);

  // Panel (b): K*0 / K0s (K*0 averaged, K0s not -> scale = 1.0)
  drawTwoPanelRatio(padRatioKstarK0s, "hFrame2Panel_KstarK0s",
                    gRatioKstarKshort,
                    gKstarK0sRatio_EPOSNew_UrQMDOff, gKstarK0sRatio_EPOSNew_UrQMDOn, gRatioKstarK0s_Thermal,
                    kKstar, kKshort, 1.0,
                    0.245, 0.54, true, "", "K*^{0} / K_{S}^{0}", false, true);

  if (isSavePlots)
  {
    cRatioKstarPhiK0s->SaveAs(Form("%s/ParticleRatios_Kstar2_Run3.pdf", plotOutputDirectory.c_str()));
  }

  //======================================================
  // 2. THREE-PANEL RATIO PLOT (3,1) WITH SEPARATE Y-AXES:
  //    Phi / Pi, Phi / K0s, Phi / Xi
  //======================================================

  gStyle->SetLineScalePS(1.7);

  TCanvas *cRatioPhiThreePanel = new TCanvas("cRatioPhiThreePanel", "Ratios: Phi/Pi, Phi/K0s, Phi/Xi", 1350, 500);
  SetPad(cRatioPhiThreePanel);
  cRatioPhiThreePanel->cd();

  // Three equal-width pads with comfortable left margins so every y-axis label is visible
  TPad *padPhiPi = new TPad("padPhiPi", "Phi/Pi", 0.00, 0.0, 0.33, 1.0);
  TPad *padPhiK0s = new TPad("padPhiK0s", "Phi/K0s", 0.33, 0.0, 0.667, 1.0);
  TPad *padPhiXi = new TPad("padPhiXi", "Phi/Xi", 0.667, 0.0, 1.00, 1.0);

  for (TPad *p : {padPhiPi, padPhiK0s, padPhiXi})
  {
    if (p == padPhiPi)
    {
      p->SetLeftMargin(0.16);
      p->SetRightMargin(0.02);
    }
    else if (p == padPhiK0s)
    {
      p->SetLeftMargin(0.19);
      p->SetRightMargin(0.01);
    }
    else if (p == padPhiXi)
    {
      p->SetLeftMargin(0.14);
      p->SetRightMargin(0.05);
    }
    p->SetTopMargin(0.05);
    p->SetBottomMargin(0.18);

    p->SetTickx(1);
    p->SetTicky(1);
    p->Draw();
  }

  auto drawSpacedRatioPanel = [&](TVirtualPad *pad,
                                  const char *frameName,
                                  TGraphErrors **data,
                                  TGraphErrors *eposOff,
                                  TGraphErrors *eposOn,
                                  TGraphErrors *thermal,
                                  int numParticle,
                                  int denParticle,
                                  double modelScale,
                                  double yMin,
                                  double yMax,
                                  const char *yTitle,
                                  const char *panelLabel,
                                  const char *ratioLabel,
                                  bool showInfo,
                                  bool showInfo2,
                                  bool showModelLegend,
                                  double yTitleOffset)
  {
    pad->cd();

    TH1D *frame = new TH1D(frameName, "", 100, 1.8, 23.5);
    SetHistStyle(frame);
    frame->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    frame->GetXaxis()->SetTitleSize(0.065);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetTitle(yTitle);
    frame->GetYaxis()->SetTitleSize(0.065);
    frame->GetYaxis()->SetTitleOffset(yTitleOffset);
    frame->GetYaxis()->SetLabelSize(0.05);
    frame->GetXaxis()->SetLimits(1.8, 23.5);
    frame->GetYaxis()->SetRangeUser(yMin, yMax);
    frame->GetYaxis()->SetMaxDigits(3);
    frame->GetYaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetTitleOffset(yTitleOffset);
    frame->Draw();

    if (data != nullptr)
    {
      TGraphErrors *dStat = (TGraphErrors *)data[0]->Clone(Form("%s_stat", frameName));
      TGraphErrors *dSys = (TGraphErrors *)data[1]->Clone(Form("%s_sys", frameName));
      TGraphErrors *dUncorr = (TGraphErrors *)data[2]->Clone(Form("%s_uncorr", frameName));

      SetGraphStyle(dStat, kRed, 1.0, 20, 0.8, 1);
      dStat->Draw("PE SAME");
      SetGraphStyle(dSys, kRed, 1.0, 20, 0.8, 1);
      dSys->SetFillStyle(0);
      dSys->Draw("5 SAME");
      SetGraphStyle(dUncorr, kRed, 1.0, 20, 0.8, 1);
      dUncorr->SetFillStyle(1001);
      dUncorr->SetFillColorAlpha(kRed, 0.35);
      dUncorr->Draw("E2 SAME");
      // dStat->Draw("PE SAME"); // Draw again to ensure points are on top
    }

    TGraphErrors *eposOffCopy = (TGraphErrors *)eposOff->Clone(Form("%s_eposOff", frameName));
    TGraphErrors *eposOnCopy = (TGraphErrors *)eposOn->Clone(Form("%s_eposOn", frameName));
    TGraphErrors *thermalCopy = (TGraphErrors *)thermal->Clone(Form("%s_thermal", frameName));

    SetGraphStyle(eposOffCopy, models[kEPOSWithoutUrQMD].color, 1.0, 1, 1.0, models[kEPOSWithoutUrQMD].width, models[kEPOSWithoutUrQMD].style);
    eposOffCopy->Draw("LX SAME");
    SetGraphStyle(eposOnCopy, models[kEPOSWithUrQMD].color, 1.0, 1, 1.0, models[kEPOSWithUrQMD].width, models[kEPOSWithUrQMD].style);
    eposOnCopy->Draw("LX SAME");

    thermalCopy->SetLineColor(models[kThermalModel].color);
    thermalCopy->SetLineStyle(models[kThermalModel].style);
    thermalCopy->SetLineWidth(models[kThermalModel].width);
    thermalCopy->Draw("LX SAME");

    std::vector<TGraphErrors *> pythiaGraphs;
    std::vector<int> pythiaModelIndices;
    for (auto model : modelsToPlot)
    {
      TGraphErrors *gMod = MakeRatio(gMYield[model][numParticle], gMYield[model][denParticle], true);
      ScaleGraph(gMod, modelScale);
      SetGraphStyle(gMod, models[model].color, 1.0, 1, 1.0, models[model].width, models[model].style);
      gMod->Draw("L SAME");
      pythiaGraphs.push_back(gMod);
      pythiaModelIndices.push_back(model);
    }

    TLatex txt;
    txt.SetNDC();
    txt.SetTextFont(42);
    txt.SetTextSize(0.048);
    // txt.DrawLatex(0.22, 0.87, panelLabel);
    txt.SetTextSize(0.055);
    txt.DrawLatex(0.27, 0.87, ratioLabel);

    if (showInfo)
    {
      txt.SetTextSize(0.065);
      txt.DrawLatex(0.72, 0.86, "#bf{ALICE}");
    }

    if (showInfo2)
    {

      TLegend *legData = new TLegend(0.48, 0.79, 0.88, 0.92);
      SetLegend(legData);
      legData->SetTextSize(0.05);
      legData->AddEntry((TObject *)0, "pp, #sqrt{s} = 13.6 TeV", "");
      legData->AddEntry((TObject *)0, "|y| < 0.5", "");
      legData->Draw();
    }

    if (showModelLegend)
    {
      TLegend *legMod = new TLegend(0.48, 0.57, 1., 0.92);
      SetLegend(legMod);
      legMod->SetTextSize(0.04);
      legMod->AddEntry(eposOffCopy, models[kEPOSWithoutUrQMD].label, "l");
      legMod->AddEntry(eposOnCopy, models[kEPOSWithUrQMD].label, "l");
      legMod->AddEntry(thermalCopy, models[kThermalModel].label, "l");
      for (size_t i = 0; i < pythiaGraphs.size(); ++i)
      {
        legMod->AddEntry(pythiaGraphs[i], models[pythiaModelIndices[i]].label, "l");
      }
      legMod->Draw();
    }
  };

  // Panel (a): Phi / Pi (pion is summed pi+ + pi- in Pythia -> scale = 2.0)
  drawSpacedRatioPanel(padPhiPi, "hFrame3Panel_PhiPi",
                       gRatioPhiPion,
                       gPhiPiRatio_EPOSNew_UrQMDOff, gPhiPiRatio_EPOSNew_UrQMDOn, gRatioPhiPi_Thermal,
                       kPhi, kPion, 2.0,
                       0.007, 0.028, "#phi / #pi", "(a)", "", true, false, false, 1.0);

  padPhiPi->cd();
  TLatex a;
  a.SetNDC();
  a.SetTextFont(42);
  a.SetTextSize(0.048);
  a.DrawLatex(0.2, 0.87, "(a)");

  // Panel (b): Phi / K0s
  drawSpacedRatioPanel(padPhiK0s, "hFrame3Panel_PhiK0s",
                       gRatioPhiK0s,
                       gPhiK0sRatio_EPOSNew_UrQMDOff, gPhiK0sRatio_EPOSNew_UrQMDOn, gRatioPhiK0s_Thermal,
                       kPhi, kKshort, 1.0,
                       0.07, 0.16, "#phi / K_{S}^{0}", "(b)", "", false, true, false, 1.40);

  padPhiK0s->cd();
  TLatex b;
  b.SetNDC();
  b.SetTextFont(42);
  b.SetTextSize(0.048);
  b.DrawLatex(0.23, 0.87, "(b)");

  // Panel (c): Phi / Xi
  // Note: gPhiXiRatio_EPOSNew_UrQMDOff/On were already scaled by 0.5 earlier in the macro.
  // Pass nullptr for data while fXiRun3 / gMYieldXiRun3 are commented out;
  // once uncommented, pass MakeRatioUncorr(gMYieldPhiRun3, gMYieldXiRun3, 1.0, 3) instead.
  drawSpacedRatioPanel(padPhiXi, "hFrame3Panel_PhiXi",
                       gRatioPhiXi,
                       gPhiXiRatio_EPOSNew_UrQMDOff, gPhiXiRatio_EPOSNew_UrQMDOn, gRatioPhiXi_Thermal,
                       kPhi, kXi, 1.0,
                       1.1, 4.3, "#phi / #Xi", "(c)", "", false, false, true, 0.95);

  padPhiXi->cd();
  TLatex c;
  c.SetNDC();
  c.SetTextFont(42);
  c.SetTextSize(0.048);
  c.DrawLatex(0.18, 0.87, "(c)");

  if (isSavePlots)
  {
    cRatioPhiThreePanel->SaveAs(Form("%s/ParticleRatios_Phi3_Run3.pdf", plotOutputDirectory.c_str()));
  }

  //======================================================
  // 3. ADDITIONAL TWO-PANEL RATIO PLOT: K*0 / Phi and K*0 / K
  //    The original K*0 / K0s canvas above is kept unchanged.
  //======================================================
  // TCanvas *cRatioKstarPhiKaon = new TCanvas("cRatioKstarPhiKaon", "Ratios: K*0/Phi and K*0/K", 850, 450);
  // SetPad(cRatioKstarPhiKaon);
  // cRatioKstarPhiKaon->cd();

  // TPad *padRatioKstarPhiKaon = makeHorizontalPad("padRatioKstarPhiKaon", "K*0/Phi", 0.0, 0.5, true);
  // TPad *padRatioKstarKaon = makeHorizontalPad("padRatioKstarKaon", "K*0/K", 0.5, 1.0, false);

  // drawTwoPanelRatio(padRatioKstarPhiKaon, "hFrame2Panel_KstarPhi_ChargedKaon",
  //                   gRatioKstarPhi,
  //                   gKstarPhiRatio_EPOSNew_UrQMDOff, gKstarPhiRatio_EPOSNew_UrQMDOn, gRatioKstarPhi_Thermal,
  //                   kKstar, kPhi, 1.0,
  //                   2.1, 5.15, false, "", "K*^{0} / #phi", true, false);

  // drawTwoPanelRatio(padRatioKstarKaon, "hFrame2Panel_KstarKaon",
  //                   gRatioKstarKaon,
  //                   gKstarKaRatio_UrQMDOffSmooth, gKstarKaRatio_UrQMDOnSmooth, gRatioKstarKa_Thermal,
  //                   kKstar, kKaon, 2.0,
  //                   0.252, 0.48, true, "", "K*^{0} / K", false, true);

  //======================================================
  // 4. ADDITIONAL THREE-PANEL RATIO PLOT:
  //    Phi / Pi, Phi / K, Phi / Xi
  //    The original Phi / K0s canvas above is kept unchanged.
  //======================================================
  // TCanvas *cRatioPhiKaonThreePanel = new TCanvas("cRatioPhiKaonThreePanel", "Ratios: Phi/Pi, Phi/K, Phi/Xi", 1350, 460);
  // SetPad(cRatioPhiKaonThreePanel);
  // cRatioPhiKaonThreePanel->cd();

  // TPad *padPhiPiKaon = new TPad("padPhiPiKaon", "Phi/Pi", 0.00, 0.0, 0.334, 1.0);
  // TPad *padPhiKaon = new TPad("padPhiKaon", "Phi/K", 0.334, 0.0, 0.667, 1.0);
  // TPad *padPhiXiKaon = new TPad("padPhiXiKaon", "Phi/Xi", 0.667, 0.0, 1.00, 1.0);

  // for (TPad *p : {padPhiPiKaon, padPhiKaon, padPhiXiKaon})
  // {
  //   if (p == padPhiPiKaon)
  //     p->SetLeftMargin(0.14);
  //   else if (p == padPhiKaon)
  //     p->SetLeftMargin(0.18);
  //   else if (p == padPhiXiKaon)
  //     p->SetLeftMargin(0.13);

  //   p->SetRightMargin(0.03);
  //   p->SetTopMargin(0.07);
  //   p->SetBottomMargin(0.15);
  //   p->SetTickx(1);
  //   p->SetTicky(1);
  //   p->Draw();
  // }

  // drawSpacedRatioPanel(padPhiPiKaon, "hFrame3Panel_PhiPi_ChargedKaon",
  //                      gRatioPhiPion,
  //                      gPhiPiRatio_EPOSNew_UrQMDOff, gPhiPiRatio_EPOSNew_UrQMDOn, gRatioPhiPi_Thermal,
  //                      kPhi, kPion, 2.0,
  //                      0.007, 0.028, "#phi / #pi", "(a)", "", true, false, false, 1.25);

  // drawSpacedRatioPanel(padPhiKaon, "hFrame3Panel_PhiKaon",
  //                      gRatioPhiK,
  //                      gPhiKaRatio_UrQMDOffSmooth, gPhiKaRatio_UrQMDOnSmooth, gRatioPhiKa_Thermal,
  //                      kPhi, kKaon, 2.0,
  //                      0.07, 0.16, "#phi / K", "(b)", "", false, true, false, 1.60);

  // drawSpacedRatioPanel(padPhiXiKaon, "hFrame3Panel_PhiXi_ChargedKaon",
  //                      gRatioPhiXi,
  //                      gPhiXiRatio_EPOSNew_UrQMDOff, gPhiXiRatio_EPOSNew_UrQMDOn, gRatioPhiXi_Thermal,
  //                      kPhi, kXi, 1.0,
  //                      1.1, 4.3, "#phi / #Xi", "(c)", "", false, false, true, 0.95);

  // Export only the three requested combined canvases.
  // if (isSavePlots)
  // {
  //   cMeanPtYieldKstarPhi->SaveAs((plotOutputDirectory + "/ParticleRatios_Observables4_Run3.pdf").c_str());
  //   cRatioKstarPhiK0s->SaveAs((plotOutputDirectory + "/ParticleRatios_Kstar2_Run3.pdf").c_str());
  //   cRatioPhiKaonThreePanel->SaveAs((plotOutputDirectory + "/ParticleRatios_Phi3_Run3.pdf").c_str());
  //   cParticleRatios->SaveAs((plotOutputDirectory + "/ParticleRatios_Combined_Run3.pdf").c_str());
  // }
}

// GRAPH HELPERS: ratios, input access, scaling, and multiplicity normalization
TGraphErrors *MakeRatio(const TGraphErrors *numerator, const TGraphErrors *denominator, bool isModel, float CorrFactorDen)
{
  if (numerator == nullptr || denominator == nullptr)
  {
    cout << "Error: null input graph while building ratio" << endl;
    return nullptr;
  }

  if (numerator->GetN() != denominator->GetN())
  {
    cout << "Error: ratio inputs have different number of points" << endl;
    return nullptr;
  }

  TGraphErrors *ratio = new TGraphErrors(numerator->GetN());
  for (int i = 0; i < numerator->GetN(); ++i)
  {
    double xNumerator = 0;
    double yNumerator = 0;
    double xDenominator = 0;
    double yDenominator = 0;
    numerator->GetPoint(i, xNumerator, yNumerator);
    denominator->GetPoint(i, xDenominator, yDenominator);

    yDenominator *= CorrFactorDen; // Apply correction factor to denominator if needed

    double yRatio = (yDenominator != 0) ? yNumerator / yDenominator : 0;
    ratio->SetPoint(i, xNumerator, yRatio);

    double xErrorNumerator = numerator->GetErrorX(i);
    double numeratorError = numerator->GetErrorY(i);
    double denominatorError = denominator->GetErrorY(i);
    denominatorError *= CorrFactorDen; // Apply correction factor to denominator error if needed
    double yRatioError = 0;
    if (yDenominator != 0)
    {
      yRatioError = sqrt(pow(numeratorError / yDenominator, 2) + pow(yNumerator * denominatorError / (yDenominator * yDenominator), 2));
    }
    if (isModel)
      ratio->SetPointError(i, 0, 0);
    else
      ratio->SetPointError(i, xErrorNumerator, yRatioError);
  }

  SetGraphStyle(ratio);

  return ratio;
}

TGraphErrors **MakeRatioUncorr(TGraphErrors **numerator, TGraphErrors **denominator, float CorrFactorDen, int nGraphs)
{
  if (!numerator || !denominator)
  {
    cout << "Error: null input graph array" << endl;
    return nullptr;
  }

  if (nGraphs < 3)
  {
    cout << "Error: nGraphs must be at least 3" << endl;
    return nullptr;
  }

  TGraphErrors **ratio = new TGraphErrors *[3];

  for (int ig = 0; ig < 3; ig++)
  {
    if (!numerator[ig] || !denominator[ig])
    {
      cout << "Error: null graph at index " << ig << endl;
      delete[] ratio;
      return nullptr;
    }

    if (numerator[ig]->GetN() != denominator[ig]->GetN())
    {
      cout << "Error: different number of points in numerator " << "and denominator graphs" << endl;

      cout << "ig = " << ig << ", numerator = " << numerator[ig]->GetN() << ", denominator = " << denominator[ig]->GetN() << endl;

      delete[] ratio;
      return nullptr;
    }

    ratio[ig] = new TGraphErrors(numerator[ig]->GetN());

    for (int i = 0; i < numerator[ig]->GetN(); i++)
    {
      double xNum, yNum;
      double xDen, yDen;

      numerator[ig]->GetPoint(i, xNum, yNum);
      denominator[ig]->GetPoint(i, xDen, yDen);

      // Apply correction factor to denominator
      yDen *= CorrFactorDen;

      double yRatio = (yDen != 0.) ? yNum / yDen : 0.;

      ratio[ig]->SetPoint(i, xNum, yRatio);

      // ---------------------------------------------------------
      // Select the errors according to the requested definition
      // ---------------------------------------------------------

      double numErr = 0.;
      double denErr = 0.;

      if (ig == 0)
      {
        // Ratio 0:
        // stat(num) / stat(den)
        numErr = numerator[0]->GetErrorY(i);
        denErr = denominator[0]->GetErrorY(i);
      }
      else if (ig == 1)
      {
        // Ratio 1:
        // uncorr(num) / total(den) (wrong)
        // total (num) / total (den) (correct)
        numErr = numerator[1]->GetErrorY(i);
        denErr = denominator[1]->GetErrorY(i);
      }
      else if (ig == 2)
      {
        // Ratio 2:
        // uncorr(num) / uncorr(den)
        numErr = numerator[2]->GetErrorY(i);
        denErr = denominator[2]->GetErrorY(i);
      }

      // Correction factor also scales the denominator uncertainty
      denErr *= CorrFactorDen;

      double yRatioErr = 0.;

      if (yDen != 0.)
      {
        yRatioErr = sqrt(
            pow(numErr / yDen, 2) +
            pow(yNum * denErr / (yDen * yDen), 2));
      }

      double xErr = numerator[ig]->GetErrorX(i);

      ratio[ig]->SetPointError(i, xErr, yRatioErr);
    }

    SetGraphStyle(ratio[ig]);
  }

  return ratio;
}

TFile *OpenFile(const string &path)
{
  TFile *f = new TFile(path.c_str(), "read");
  if (f->IsZombie())
  {
    cout << "Error: File not found: " << path << endl;
    return nullptr;
  }
  return f;
}

TGraphErrors *GetGraph(TFile *f, const string &name)
{
  if (!f)
  {
    cout << "Error: null ROOT file while loading graph " << name << endl;
    return nullptr;
  }

  TObject *object = f->Get(name.c_str());
  TGraphErrors *graph = dynamic_cast<TGraphErrors *>(object);
  if (!graph)
  {
    // Run 3 phi graphs are stored as TGraphAsymmErrors. Convert their
    // asymmetric errors to the mean lower/upper uncertainty so the source
    // macro's existing symmetric ratio propagation remains applicable.
    TGraphAsymmErrors *asymm = dynamic_cast<TGraphAsymmErrors *>(object);
    if (asymm)
    {
      graph = new TGraphErrors(asymm->GetN());
      graph->SetName((name + "_symmetric").c_str());
      graph->SetTitle(asymm->GetTitle());
      for (int i = 0; i < asymm->GetN(); ++i)
      {
        double x, y;
        asymm->GetPoint(i, x, y);
        const double ex = 0.5 * (asymm->GetErrorXlow(i) + asymm->GetErrorXhigh(i));
        const double ey = 0.5 * (asymm->GetErrorYlow(i) + asymm->GetErrorYhigh(i));
        graph->SetPoint(i, x, y);
        graph->SetPointError(i, ex, ey);
      }
    }
  }

  if (!graph)
  {
    cout << "Error: graph " << name
         << " not found or has an unsupported class in file " << f->GetName() << endl;
    return nullptr;
  }

  SetGraphStyle(graph);
  graph->SetTitle("");
  return graph;
}

void ScaleGraph(TGraph *gr, double scale)
{
  if (!gr)
    return;

  for (int i = 0; i < gr->GetN(); ++i)
  {
    double x, y;
    gr->GetPoint(i, x, y);
    gr->SetPoint(i, x, y * scale);
  }

  if (auto *ge = dynamic_cast<TGraphErrors *>(gr))
  {
    for (int i = 0; i < ge->GetN(); ++i)
    {
      ge->SetPointError(i, ge->GetErrorX(i), ge->GetErrorY(i) * scale);
    }
  }
}

TGraphErrors **DivideByMult(TGraphErrors **gr, double WhichMultPoint, double tolerance, int noOfGraphs)
{
  if (!gr || !gr[0] || !gr[2])
  {
    cout << "Error: null graph array or required graphs" << endl;
    return nullptr;
  }

  if (noOfGraphs < 3)
  {
    cout << "Error: noOfGraphs must be at least 3" << endl;
    return nullptr;
  }

  // Output:
  // grCopy[0] = ratio with statistical uncertainty
  // grCopy[1] = ratio with uncorrelated uncertainty
  TGraphErrors **grCopy = new TGraphErrors *[2];

  // ------------------------------------------------------------
  // Find the reference multiplicity point
  // ------------------------------------------------------------

  int refIndex = -1;

  if (WhichMultPoint >= 0)
  {
    for (int i = 0; i < gr[0]->GetN(); i++)
    {
      double x, y;
      gr[0]->GetPoint(i, x, y);

      if (fabs(x - WhichMultPoint) < tolerance)
      {
        refIndex = i;
        break;
      }
    }

    if (refIndex < 0)
    {
      cout << "Error: Could not find multiplicity point "
           << WhichMultPoint << endl;

      delete[] grCopy;
      return nullptr;
    }
  }
  else
  {
    // Use the lowest multiplicity point
    double xFirst, yFirst;
    double xLast, yLast;

    gr[0]->GetPoint(0, xFirst, yFirst);
    gr[0]->GetPoint(gr[0]->GetN() - 1, xLast, yLast);

    refIndex = (xFirst < xLast) ? 0 : gr[0]->GetN() - 1;
  }

  // ------------------------------------------------------------
  // Reference value
  // ------------------------------------------------------------

  double xRef, yRef;
  gr[0]->GetPoint(refIndex, xRef, yRef);

  if (yRef == 0.)
  {
    cout << "Error: Reference value is zero" << endl;
    delete[] grCopy;
    return nullptr;
  }

  // Reference uncertainties
  double yRefStatErr = gr[0]->GetErrorY(refIndex);
  double yRefUncorrErr = gr[2]->GetErrorY(refIndex);

  // ------------------------------------------------------------
  // Create the two output graphs
  // ------------------------------------------------------------

  grCopy[0] = (TGraphErrors *)gr[0]->Clone();
  grCopy[1] = (TGraphErrors *)gr[2]->Clone();

  // ------------------------------------------------------------
  // Calculate ratios
  // ------------------------------------------------------------

  for (int i = 0; i < gr[0]->GetN(); i++)
  {
    double x, y;
    gr[0]->GetPoint(i, x, y);

    double xerr = gr[0]->GetErrorX(i);

    // ========================================================
    // 1. STATISTICAL RATIO
    // ========================================================

    double yStatErr = gr[0]->GetErrorY(i);

    double ratio = y / yRef;

    double ratioStatErr = sqrt(
        pow(yStatErr / yRef, 2) +
        pow(y * yRefStatErr / (yRef * yRef), 2));

    grCopy[0]->SetPoint(i, x, ratio);
    grCopy[0]->SetPointError(i, xerr, ratioStatErr);

    // ========================================================
    // 2. UNCORRELATED SYSTEMATIC RATIO
    // ========================================================

    double yUncorrErr = gr[2]->GetErrorY(i);

    double ratioUncorrErr = sqrt(pow(yUncorrErr / yRef, 2) + pow(y * yRefUncorrErr / (yRef * yRef), 2));

    grCopy[1]->SetPoint(i, x, ratio);
    grCopy[1]->SetPointError(i, xerr, ratioUncorrErr);
  }

  SetGraphStyle(grCopy[0]);
  SetGraphStyle(grCopy[1]);

  return grCopy;
}

void RestrictModelXaxis(TGraphErrors *gr, double xMin, double xMax)
{
  for (int i = gr->GetN() - 1; i >= 0; --i)
  {
    double x, y;
    gr->GetPoint(i, x, y);
    if (x < xMin || x > xMax)
    {
      gr->RemovePoint(i);
    }
  }
}

TGraphErrors *DivideByMultModel(TGraphErrors *gr, double WhichMultPoint, double tolerance)
{
  TGraphErrors *grCopy = (TGraphErrors *)gr->Clone();

  double yGivenMult = 1.0;
  double yGivenMultErr = 0.0;

  for (int i = 0; i < gr->GetN(); i++)
  {
    double x, y;
    gr->GetPoint(i, x, y);

    if (fabs(x - WhichMultPoint) < tolerance)
    {
      yGivenMult = y;
      yGivenMultErr = gr->GetErrorY(i);
      break;
    }
  }

  for (int i = 0; i < gr->GetN(); i++)
  {
    double x, y;
    gr->GetPoint(i, x, y);

    double yerr = gr->GetErrorY(i);
    double xerr = gr->GetErrorX(i);

    // Normalization point
    if (fabs(x - WhichMultPoint) < tolerance)
    {
      grCopy->SetPoint(i, x, 1.0);
      grCopy->SetPointError(i, xerr, 0.0);
      continue;
    }

    double ratio = y / yGivenMult;

    double ratioErr = sqrt(pow(yerr / yGivenMult, 2) + pow(yGivenMultErr * y / (yGivenMult * yGivenMult), 2));

    grCopy->SetPoint(i, x, ratio);
    grCopy->SetPointError(i, xerr, ratioErr);
  }

  return grCopy;
}

TGraphErrors *SmoothGraph(const TGraphErrors *gr)
{
  int n = gr->GetN();

  TGraphErrors *grSmooth = new TGraphErrors(n);

  for (int i = 0; i < n; i++)
  {

    double x, y;
    gr->GetPoint(i, x, y);

    double ySmooth = 0.0;
    double errSmooth = 0.0;
    int nPoints = 0;

    // Average neighbouring points
    for (int j = i - 1; j <= i + 1; j++)
    {

      if (j < 0 || j >= n)
        continue;

      double xj, yj;
      gr->GetPoint(j, xj, yj);

      double eyj = gr->GetErrorY(j);

      ySmooth += yj;
      errSmooth += eyj * eyj;

      nPoints++;
    }

    ySmooth /= nPoints;

    // Error propagation for the average
    errSmooth = std::sqrt(errSmooth) / nPoints;

    grSmooth->SetPoint(i, x, ySmooth);
    grSmooth->SetPointError(
        i,
        gr->GetErrorX(i),
        errSmooth);
  }

  return grSmooth;
}
