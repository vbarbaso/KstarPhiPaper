#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TString.h>
#include "StyleProposal.h"
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

void ParticleRatio_simple()
{
    gStyle->SetEndErrorSize(0);
    //============================================================
    // USER CONFIGURATION: edit paths and plot output settings here
    //============================================================
    bool isSavePlots = false;
    const string plotOutputDirectory = "../figures";

    const Color_t colorEPOSWithoutUrQMD = kGreen + 2;
    const int lineStyleEPOSWithoutUrQMD = 2;
    const Color_t colorEPOSWithUrQMD = kBlue - 6;
    const int lineStyleEPOSWithUrQMD = 9;

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
    // TFile *fRho = OpenFile("../inputs/Rho_Run3_Results/Sawan/ResultsRho.root");
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

        // Previous combined Pi/K/p input alternatives
        // gMPtPion[i] = GetGraph(fPiKp, Form("gPion_MeanpT%s", suffix.c_str()));
        // gMPtKaon[i] = GetGraph(fPiKp, Form("gKaon_MeanpT%s", suffix.c_str()));
        // gMPtProton[i] = GetGraph(fPiKp, Form("gProton_MeanpT%s", suffix.c_str()));
        // gMYieldPion[i] = GetGraph(fPiKp, Form("gPion_MeanYield%s", suffix.c_str()));
        // gMYieldKaon[i] = GetGraph(fPiKp, Form("gKaon_MeanYield%s", suffix.c_str()));
        // gMYieldProton[i] = GetGraph(fPiKp, Form("gProton_MeanYield%s", suffix.c_str()));

        // Rho Run 3 and Xi Run 3 input alternatives (kept as in the original)
        // gMPtRho[i] = GetGraph(fRho, Form("gMeanpTRun3%s", suffix.c_str()));
        // gMYieldRho[i] = GetGraph(fRho, Form("gMeanYieldRun3%s", suffix.c_str()));
        gMYieldXiRun3[i] = GetGraph(fXiRun3, Form("gMeanYieldRun3%s", suffix.c_str()));
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

    enum Model
    {
        kPythiaMonash,
        kPythiaRopes,
        kPythiaShoving,
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

    struct ModelStyle
    {
        Color_t color;
        int style;
    };

    ModelStyle modelStyle[kNModels] = {
        {kCyan + 1, 1},   // Pythia Monash
        {kYellow + 1, 1}, // Pythia Ropes
        {kMagenta + 1, 1} // Pythia Shoving
    };

    const char *modelLabel[kNModels] = {
        "Pythia Monash",
        "Pythia Ropes",
        "Pythia Shoving"};

    TFile *ModelsHyperloop = new TFile(inputPythiaModels.c_str(), "read");
    if (ModelsHyperloop->IsZombie())
    {
        cout << "Error: Hyperloop model file not found" << endl;
        return;
    }

    string hyperloopModels[kNModels] = {"Pythia_Monash", "Pythia_Ropes", "Pythia_Shoving"};
    string particles[kNParticles] = {"Kstar", "Phi", "Pion", "Kaon", "Proton", "PionMinus", "KaonMinus", "AntiProton", "Xi1530", "Kshort", "KstarPM", "XiMinus", "XiPlus"};
    vector<int> modelsToPlot = {kPythiaMonash, kPythiaShoving, kPythiaRopes};

    TGraphErrors *gMYield[kNModels][kNParticles];
    TGraphErrors *gMeanPt[kNModels][kNParticles];
    TH1D *hSpectra[kNModels][kNParticles];

    for (int iModel = 0; iModel < kNModels; iModel++)
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
            gMYield[iModel][iParticle]->SetLineColor(modelStyle[iModel].color);
            gMYield[iModel][iParticle]->SetLineStyle(modelStyle[iModel].style);

            gMeanPt[iModel][iParticle]->SetLineColor(modelStyle[iModel].color);
            gMeanPt[iModel][iParticle]->SetLineStyle(modelStyle[iModel].style);

            int nPoints = gMYield[iModel][iParticle]->GetN();

            // RestrictModelXaxis(gMYield[iModel][iParticle], 3.3, 22.5);
            // RestrictModelXaxis(gMeanPt[iModel][iParticle], 3.3, 22.5);
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
    // RestrictModelXaxis(gEPOS_KstardNdy_UrQMDOff, 3.3, 23.5);
    // RestrictModelXaxis(gEPOS_KstardNdy_UrQMDOn, 3.3, 24.5);
    // RestrictModelXaxis(gEPOS_KstarMeanpT_UrQMDOff, 3.3, 23.5);
    // RestrictModelXaxis(gEPOS_KstarMeanpT_UrQMDOn, 3.3, 24.5);

    TGraphErrors *gEPOS_PhidNdy_UrQMDOff = GetGraph(fEPOS, "g_phi_dNdy_opt0");
    TGraphErrors *gEPOS_PhidNdy_UrQMDOn = GetGraph(fEPOS, "g_phi_dNdy_opt1");
    TGraphErrors *gEPOS_PhiMeanpT_UrQMDOff = GetGraph(fEPOS, "g_phi_pT_opt0");
    TGraphErrors *gEPOS_PhiMeanpT_UrQMDOn = GetGraph(fEPOS, "g_phi_pT_opt1");
    // RestrictModelXaxis(gEPOS_PhidNdy_UrQMDOff, 3.3, 23.5);
    // RestrictModelXaxis(gEPOS_PhidNdy_UrQMDOn, 3.3, 24.5);
    // RestrictModelXaxis(gEPOS_PhiMeanpT_UrQMDOff, 3.3, 23.5);
    // RestrictModelXaxis(gEPOS_PhiMeanpT_UrQMDOn, 3.3, 24.5);

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
            gMYield[model][particle]->SetLineColor(modelStyle[model].color);
            gMYield[model][particle]->SetLineStyle(modelStyle[model].style);
            gMYield[model][particle]->SetLineWidth(2);

            gMeanPt[model][particle]->SetLineColor(modelStyle[model].color);
            gMeanPt[model][particle]->SetLineStyle(modelStyle[model].style);
            gMeanPt[model][particle]->SetLineWidth(2);
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

    // PLOTTING: yields, mean pT, and particle ratios
    ////======================================================
    //    ================Plots=====================
    //======================================================

    TCanvas *cdNdyKstar = new TCanvas("cdNdyKstar", "cdNdyKstar", 720, 720);
    SetPad(cdNdyKstar);
    gMYieldKstar[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMYieldKstar[0]->GetYaxis()->SetTitle("dN/dy");
    gMYieldKstar[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldKstar[0]->SetLineWidth(2);
    gMYieldKstar[0]->GetYaxis()->SetRangeUser(0, 0.71);
    gMYieldKstar[0]->SetMarkerColor(kRed + 1);
    gMYieldKstar[0]->SetLineColor(kRed + 1);
    gMYieldKstar[0]->Draw("APE");
    gMYieldKstar[1]->SetFillStyle(0);
    gMYieldKstar[1]->SetLineColor(kRed + 1);
    gMYieldKstar[1]->SetLineWidth(2);
    gMYieldKstar[1]->Draw("5 same");
    gMYieldKstar[2]->SetFillStyle(1001);
    gMYieldKstar[2]->SetFillColorAlpha(kRed + 1, 0.3);
    gMYieldKstar[2]->Draw("E2 same");
    gMYieldKstarRun2[0]->SetMarkerStyle(24);
    gMYieldKstarRun2[0]->SetMarkerColor(kBlue + 1);
    gMYieldKstarRun2[0]->SetLineColor(kBlue + 1);
    gMYieldKstarRun2[0]->SetLineWidth(2);
    gMYieldKstarRun2[0]->Draw("P same");
    gMYieldKstarRun2[1]->SetLineColor(kBlue + 1);
    gMYieldKstarRun2[1]->SetFillStyle(0);
    gMYieldKstarRun2[1]->SetLineWidth(2);
    gMYieldKstarRun2[1]->Draw("5 same");
    gMYieldKstarRun2[2]->SetFillStyle(1001);
    gMYieldKstarRun2[2]->SetFillColorAlpha(kBlue + 1, 0.3);
    gMYieldKstarRun2[2]->Draw("E2 same");

    int TotalPointsRun3 = gMYieldKstar[0]->GetN();
    int TotalPointsRun2 = gMYieldKstarRun2[0]->GetN();

    for (auto model : modelsToPlot)
    {
        ScaleGraph(gMYield[model][kKstar], 0.5);
        gMYield[model][kKstar]->Draw("l same");
    }

    int totalPointsEPOS = gEPOS_KstardNdy_UrQMDOff->GetN();
    for (int i = 0; i < totalPointsEPOS; i++)
    {
        double x, y;
        gEPOS_KstardNdy_UrQMDOff->SetPointError(i, 0, 0); // Remove error
        gEPOS_KstardNdy_UrQMDOn->SetPointError(i, 0, 0);  // Remove error
    }

    SetGraphStyle(gEPOS_KstardNdy_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gEPOS_KstardNdy_UrQMDOff->Draw("l same");
    SetGraphStyle(gEPOS_KstardNdy_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gEPOS_KstardNdy_UrQMDOn->Draw("l same");

    TLegend *legend = new TLegend(0.4, 0.82, 0.8, 0.91);
    SetLegend(legend);
    legend->SetTextSize(0.035);
    legend->AddEntry(gMYieldKstar[0], "pp, #sqrt{s} = 13.6 TeV", "P");
    legend->AddEntry(gMYieldKstarRun2[0], "pp, #sqrt{s} = 13 TeV", "P");
    legend->Draw();

    TLegend *legend2 = new TLegend(0.2, 0.66, 0.5, 0.92);
    SetLegend(legend2);
    legend2->SetTextSize(0.035);
    legend2->SetHeader("Models: pp #sqrt{s} = 13.6 TeV");

    legend2->AddEntry(gEPOS_KstardNdy_UrQMDOff, "EPOS4", "L");
    legend2->AddEntry(gEPOS_KstardNdy_UrQMDOn, "EPOS4 + UrQMD", "L");

    for (auto model : modelsToPlot)
    {
        legend2->AddEntry(gMYield[model][kKstar], modelLabel[model], "L");
    }
    legend->Draw();
    // legend2->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.22, 0.9, "ALICE");
    latex.DrawLatex(0.22, 0.78, "|y| < 0.5");
    latex.SetTextFont(22);
    latex.DrawLatex(0.22, 0.84, "K* (892)^{0}");
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.47, 0.3, "Uncertainties:");
    latex.SetTextFont(42);
    latex.DrawLatex(0.47, 0.25, "stat. (bars), total sys. (open box)");
    latex.DrawLatex(0.47, 0.2, "uncorr. sys. (shaded box)");

    if (isSavePlots)
    {
        cdNdyKstar->SaveAs((plotOutputDirectory + "/MeanYield_Kstar.pdf").c_str());
    }

    //===================================================
    //  ================dN/dy Phi======================
    //===================================================

    TCanvas *cdNdyPhi = new TCanvas("cdNdyPhi", "cdNdyPhi", 720, 720);
    SetPad(cdNdyPhi);
    gMYieldPhiRun3[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMYieldPhiRun3[0]->GetYaxis()->SetTitle("dN/dy");
    gMYieldPhiRun3[0]->GetXaxis()->SetLimits(0, 27);
    gMYieldPhiRun3[0]->SetLineWidth(2);
    gMYieldPhiRun3[0]->GetYaxis()->SetRangeUser(0, 0.26);
    gMYieldPhiRun3[0]->SetMarkerColor(kRed + 1);
    gMYieldPhiRun3[0]->SetLineColor(kRed + 1);
    gMYieldPhiRun3[0]->Draw("APE");
    gMYieldPhiRun3[1]->SetFillStyle(0);
    gMYieldPhiRun3[1]->SetLineColor(kRed + 1);
    gMYieldPhiRun3[1]->SetLineWidth(2);
    gMYieldPhiRun3[1]->Draw("5 same");
    gMYieldPhiRun3[2]->SetFillStyle(1001);
    gMYieldPhiRun3[2]->SetFillColorAlpha(kRed + 1, 0.3);
    gMYieldPhiRun3[2]->Draw("E2 same");
    gMYieldPhiRun2[0]->SetMarkerStyle(24);
    gMYieldPhiRun2[0]->SetMarkerColor(kBlue + 1);
    gMYieldPhiRun2[0]->SetLineColor(kBlue + 1);
    gMYieldPhiRun2[0]->SetLineWidth(2);
    gMYieldPhiRun2[0]->Draw("P same");
    gMYieldPhiRun2[1]->SetLineColor(kBlue + 1);
    gMYieldPhiRun2[1]->SetFillStyle(0);
    gMYieldPhiRun2[1]->SetLineWidth(2);
    gMYieldPhiRun2[1]->Draw("5 same");
    gMYieldPhiRun2[2]->SetFillStyle(1001);
    gMYieldPhiRun2[2]->SetFillColorAlpha(kBlue + 1, 0.3);
    gMYieldPhiRun2[2]->Draw("E2 same");

    int totalPointsEPOS_Phi = gEPOS_PhidNdy_UrQMDOff->GetN();
    for (int i = 0; i < totalPointsEPOS_Phi; i++)
    {
        gEPOS_PhidNdy_UrQMDOff->SetPointError(i, 0, 0); // Remove error
        gEPOS_PhidNdy_UrQMDOn->SetPointError(i, 0, 0);  // Remove error
    }

    SetGraphStyle(gEPOS_PhidNdy_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gEPOS_PhidNdy_UrQMDOff->Draw("l same");
    SetGraphStyle(gEPOS_PhidNdy_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gEPOS_PhidNdy_UrQMDOn->Draw("l same");

    for (auto model : modelsToPlot)
    {
        gMYield[model][kPhi]->Draw("l same");
    }

    legend2->Draw();
    latex.SetTextFont(22);
    latex.DrawLatex(0.76, 0.88, "#Phi (1020)");

    if (isSavePlots)
    {
        cdNdyPhi->SaveAs((plotOutputDirectory + "/MeanYield_Phi.pdf").c_str());
    }

    //===================================================
    //  ================<pT> K*0======================
    //===================================================
    TCanvas *cMeanPtKstar = new TCanvas("cMeanPtKstar", "cMeanPtKstar", 720, 720);
    SetPad(cMeanPtKstar);
    SetGraphStyle(gMPtKstar[0]);
    gMPtKstar[0]->SetTitle(0);
    gMPtKstar[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMPtKstar[0]->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/#it{c})");
    gMPtKstar[0]->GetXaxis()->SetLimits(0, 27);
    gMPtKstar[0]->GetYaxis()->SetRangeUser(0.58, 1.76);
    gMPtKstar[0]->SetMarkerColor(kRed);
    gMPtKstar[0]->SetLineColor(kRed);
    gMPtKstar[0]->SetLineWidth(2);
    gMPtKstar[0]->GetYaxis()->SetTitleOffset(1.3);
    gMPtKstar[0]->Draw("APE");
    gMPtKstar[1]->SetFillStyle(0);
    gMPtKstar[1]->SetLineColor(kRed);
    gMPtKstar[1]->SetLineWidth(2);
    gMPtKstar[1]->Draw("5 same");
    gMPtKstar[2]->SetFillStyle(1001);
    gMPtKstar[2]->SetFillColorAlpha(kRed, 0.3);
    gMPtKstar[2]->Draw("E2 same");

    gMPtKstarRun2[0]->SetMarkerStyle(24);
    gMPtKstarRun2[0]->SetMarkerColor(kBlue);
    gMPtKstarRun2[0]->SetLineColor(kBlue);
    gMPtKstarRun2[0]->SetLineWidth(2);
    gMPtKstarRun2[0]->Draw("P same");
    gMPtKstarRun2[1]->SetLineColor(kBlue);
    gMPtKstarRun2[1]->SetFillStyle(0);
    gMPtKstarRun2[1]->SetLineWidth(2);
    gMPtKstarRun2[1]->Draw("5 same");
    gMPtKstarRun2[2]->SetFillStyle(1001);
    gMPtKstarRun2[2]->SetFillColorAlpha(kBlue, 0.3);
    gMPtKstarRun2[2]->Draw("E2 same");

    int totalPointsEPOS_MeanPt = gEPOS_KstarMeanpT_UrQMDOff->GetN();
    for (int i = 0; i < totalPointsEPOS_MeanPt; i++)
    {
        gEPOS_KstarMeanpT_UrQMDOff->SetPointError(i, 0, 0); // Remove error
        gEPOS_KstarMeanpT_UrQMDOn->SetPointError(i, 0, 0);  // Remove error
    }

    SetGraphStyle(gEPOS_KstarMeanpT_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gEPOS_KstarMeanpT_UrQMDOff->Draw("l same");
    SetGraphStyle(gEPOS_KstarMeanpT_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gEPOS_KstarMeanpT_UrQMDOn->Draw("l same");

    for (auto model : modelsToPlot)
    {
        gMeanPt[model][kKstar]->Draw("l same");
    }

    legend->Draw();
    // legend2->Draw();
    latex.DrawLatex(0.22, 0.9, "ALICE");
    latex.DrawLatex(0.22, 0.78, "|y| < 0.5");
    latex.SetTextFont(22);
    latex.DrawLatex(0.22, 0.84, "K* (892)^{0}");
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.47, 0.3, "Uncertainties:");
    latex.SetTextFont(42);
    latex.DrawLatex(0.47, 0.25, "stat. (bars), total sys. (open box)");
    latex.DrawLatex(0.47, 0.2, "uncorr. sys. (shaded box)");
    if (isSavePlots)
    {
        cMeanPtKstar->SaveAs((plotOutputDirectory + "/MeanPt_Kstar.pdf").c_str());
    }

    //===================================================
    //  ================<pT> Phi======================
    //===================================================
    TCanvas *cMeanPtPhi = new TCanvas("cMeanPtPhi", "cMeanPtPhi", 720, 720);
    SetPad(cMeanPtPhi);
    SetGraphStyle(gMPtPhiRun3[0]);
    gMPtPhiRun3[0]->SetTitle(0);
    gMPtPhiRun3[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    gMPtPhiRun3[0]->GetYaxis()->SetTitle("#LT #it{p}_{T} #GT (GeV/#it{c})");
    gMPtPhiRun3[0]->GetXaxis()->SetLimits(0, 27);
    gMPtPhiRun3[0]->GetYaxis()->SetRangeUser(0.58, 1.76);
    gMPtPhiRun3[0]->SetMarkerColor(kRed);
    gMPtPhiRun3[0]->SetLineColor(kRed);
    gMPtPhiRun3[0]->SetLineWidth(2);
    gMPtPhiRun3[0]->GetYaxis()->SetTitleOffset(1.3);
    gMPtPhiRun3[0]->Draw("APE");
    gMPtPhiRun3[1]->SetFillStyle(0);
    gMPtPhiRun3[1]->SetLineColor(kRed);
    gMPtPhiRun3[1]->SetLineWidth(2);
    gMPtPhiRun3[1]->Draw("5 same");
    gMPtPhiRun3[2]->SetFillStyle(1001);
    gMPtPhiRun3[2]->SetFillColorAlpha(kRed, 0.3);
    gMPtPhiRun3[2]->Draw("E2 same");
    gMPtPhiRun2[0]->SetMarkerStyle(24);
    gMPtPhiRun2[0]->SetMarkerColor(kBlue);
    gMPtPhiRun2[0]->SetLineColor(kBlue);
    gMPtPhiRun2[0]->SetLineWidth(2);
    gMPtPhiRun2[0]->Draw("P same");
    gMPtPhiRun2[1]->SetLineColor(kBlue);
    gMPtPhiRun2[1]->SetFillStyle(0);
    gMPtPhiRun2[1]->SetLineWidth(2);
    gMPtPhiRun2[1]->Draw("5 same");
    gMPtPhiRun2[2]->SetFillStyle(1001);
    gMPtPhiRun2[2]->SetFillColorAlpha(kBlue, 0.3);
    gMPtPhiRun2[2]->Draw("E2 same");

    for (int i = 0; i < totalPointsEPOS_Phi; i++)
    {
        gEPOS_PhiMeanpT_UrQMDOff->SetPointError(i, 0, 0); // Remove error
        gEPOS_PhiMeanpT_UrQMDOn->SetPointError(i, 0, 0);  // Remove error
    }

    SetGraphStyle(gEPOS_PhiMeanpT_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gEPOS_PhiMeanpT_UrQMDOff->Draw("l same");
    SetGraphStyle(gEPOS_PhiMeanpT_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gEPOS_PhiMeanpT_UrQMDOn->Draw("l same");

    for (auto model : modelsToPlot)
    {
        gMeanPt[model][kKstar]->Draw("l same");
    }

    // legend->Draw();
    legend2->Draw();
    latex.SetTextSize(0.03);
    latex.SetTextFont(22);
    latex.DrawLatex(0.8, 0.90, "#Phi (1020)");
    if (isSavePlots)
        cMeanPtPhi->SaveAs((plotOutputDirectory + "/MeanPt_Phi.pdf").c_str());

    //====================================================
    //   ================Kstar/K Ratio ===================
    //====================================================

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
        // cout << "xPhi: " << xK0s << ", yPhiOff: " << yPhiOff << ", yPhiOn: " << yPhiOn << endl;

        gdNdy_Kstar_Off->GetPoint(i, xK0s, yKstarOff);
        gdNdy_Kstar_On->GetPoint(i, xK0s, yKstarOn);
        // cout << "xKstar: " << xK0s << ", yKstarOff: " << yKstarOff << ", yKstarOn: " << yKstarOn << endl;

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

    // // Low to high multiplicity suppression in EPOS UrQMD Off and On
    // double EPOSHM_Off, EPOSLM_Off, EPOSHM_On, EPOSLM_On, xLM2, xHM2;
    // gKstarK0sRatio_EPOSNew_UrQMDOff->GetPoint(0, xLM2, EPOSLM_Off);
    // gKstarK0sRatio_EPOSNew_UrQMDOff->GetPoint(gKstarK0sRatio_EPOSNew_UrQMDOff->GetN() - 1, xHM2, EPOSHM_Off);
    // gKstarK0sRatio_EPOSNew_UrQMDOn->GetPoint(0, xLM2, EPOSLM_On);
    // gKstarK0sRatio_EPOSNew_UrQMDOn->GetPoint(gKstarK0sRatio_EPOSNew_UrQMDOn->GetN() - 1, xHM2, EPOSHM_On);
    // cout << "EPOS UrQMD Off K*0/K0s ratio suppression: " << EPOSHM_Off / EPOSLM_Off << endl;
    // cout << "EPOS UrQMD On K*0/K0s ratio suppression: " << EPOSHM_On / EPOSLM_On << endl;

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

    TCanvas *cRatioKstarKaon = new TCanvas("cRatioKstarKaon", "cRatioKstarKaon", 720, 720);
    SetPad(cRatioKstarKaon);

    // for (int i3 = 0; i3 < TotalPointsRun3; i3++)
    // {
    //     double xkstar, ykstar, ysyskstar, yuncorrkstar;
    //     double xkaon, ykaon, ysyskaon, yuncorrkaon;
    //     gMYieldKstar[0]->GetPoint(i3, xkstar, ykstar);
    //     ysyskstar = gMYieldKstar[1]->GetErrorY(i3);
    //     yuncorrkstar = gMYieldKstar[2]->GetErrorY(i3);

    //     gMYieldKaon[0]->GetPoint(i3, xkaon, ykaon);
    //     ysyskaon = gMYieldKaon[1]->GetErrorY(i3);
    //     yuncorrkaon = gMYieldKaon[2]->GetErrorY(i3);

    //     cout << "Point " << i3 << ": K*0 yield = " << ykstar << " +/- " << ysyskstar << " (sys) +/- " << yuncorrkstar << " (uncorr)" << endl;
    //     cout << "Point " << i3 << ": Kaon yield = " << ykaon << " +/- " << ysyskaon << " (sys) +/- " << yuncorrkaon << " (uncorr)" << endl;
    //     double ratio = ykstar / ykaon;
    //     double ratioErrorsys = ratio * sqrt(pow(yuncorrkstar / ykstar, 2) + pow(ysyskaon / ykaon, 2));
    //     double ratioErroruncorr = ratio * sqrt(pow(yuncorrkstar / ykstar, 2) + pow(yuncorrkaon / ykaon, 2));
    //     cout << "Point " << i3 << ": K*0/Kaon ratio = " << ratio << " +/- " << ratioErrorsys << " (sys) +/- " << ratioErroruncorr << " (uncorr)" << endl;
    //     cout << endl;
    // }

    TGraphErrors **gRatioKstarKaon = MakeRatioUncorr(gMYieldKstar, gMYieldKaon, 1.0, 3);
    gRatioKstarKaon[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioKstarKaon[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphStyle(gRatioKstarKaon[0]);
    gRatioKstarKaon[0]->GetYaxis()->SetRangeUser(0.245, 0.44);
    gRatioKstarKaon[0]->GetXaxis()->SetLimits(1.8, 23.5);
    gRatioKstarKaon[0]->SetMarkerColor(kRed);
    gRatioKstarKaon[0]->SetLineColor(kRed);
    gRatioKstarKaon[0]->GetYaxis()->SetNdivisions(505);
    gRatioKstarKaon[0]->Draw("APE");
    gRatioKstarKaon[1]->SetFillStyle(0);
    gRatioKstarKaon[1]->SetLineColor(kRed);
    gRatioKstarKaon[1]->Draw("5 same");
    gRatioKstarKaon[2]->SetFillStyle(1001);
    gRatioKstarKaon[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioKstarKaon[2]->Draw("E2 same");

    // gKstarKaRatio_13TeV[0]->SetMarkerStyle(21);
    // gKstarKaRatio_13TeV[0]->SetMarkerColor(kBlue);
    // gKstarKaRatio_13TeV[0]->SetLineColor(kBlue);
    // gKstarKaRatio_13TeV[0]->Draw("P same");
    // gKstarKaRatio_13TeV[1]->SetLineColor(kBlue);
    // gKstarKaRatio_13TeV[1]->SetFillStyle(0);
    // gKstarKaRatio_13TeV[1]->Draw("5 same");

    SetGraphStyle(gKstarKaRatio_UrQMDOffSmooth, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gKstarKaRatio_UrQMDOffSmooth->Draw("LX same");
    SetGraphStyle(gKstarKaRatio_UrQMDOnSmooth, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gKstarKaRatio_UrQMDOnSmooth->Draw("LX same");

    gRatioKstarKa_Thermal->SetLineColor(kMagenta + 2);
    gRatioKstarKa_Thermal->SetLineStyle(5);
    gRatioKstarKa_Thermal->SetLineWidth(2);
    gRatioKstarKa_Thermal->Draw("LX same");

    ////Calculation of suppression significance between HM and LM for K*0/Kaon ratio
    int totalPointsKstarKaonRatio = gRatioKstarKaon[0]->GetN();
    double xLM, xHM, yHM, yLM, yHMError, yLMError, yHMStatError, yLMStatError, yLMUncorr, yHMUncorr;
    gRatioKstarKaon[1]->GetPoint(totalPointsKstarKaonRatio - 1, xLM, yLM);
    gRatioKstarKaon[1]->GetPoint(0, xHM, yHM);
    yLMError = gRatioKstarKaon[1]->GetErrorY(totalPointsKstarKaonRatio - 1);
    yHMError = gRatioKstarKaon[1]->GetErrorY(0);
    yLMStatError = gRatioKstarKaon[0]->GetErrorY(totalPointsKstarKaonRatio - 1);
    yHMStatError = gRatioKstarKaon[0]->GetErrorY(0);
    yLMUncorr = gRatioKstarKaon[2]->GetErrorY(totalPointsKstarKaonRatio - 1);
    yHMUncorr = gRatioKstarKaon[2]->GetErrorY(0);

    double sigmaDiffHM_LM = abs(yHM - yLM) / sqrt(pow(yHMError, 2) + pow(yLMError, 2) + pow(yHMStatError, 2) + pow(yLMStatError, 2));
    double sigmaDiffHM_LM_uncorr = abs(yHM - yLM) / sqrt(pow(yHMUncorr, 2) + pow(yLMUncorr, 2));
    double sigmaDiffHM_LM_uncorrStat = abs(yHM - yLM) / sqrt(pow(yHMStatError, 2) + pow(yLMStatError, 2) + pow(yHMUncorr, 2) + pow(yLMUncorr, 2));
    cout << "K*0/Kaon ratio difference between HM and LM: " << sigmaDiffHM_LM << " sigma" << endl;
    cout << "K*0/Kaon ratio difference between HM and LM (uncorrelated): " << sigmaDiffHM_LM_uncorr << " sigma" << endl;
    cout << "K*0/Kaon ratio difference between HM and LM (uncorrelated, statistical): " << sigmaDiffHM_LM_uncorrStat << " sigma" << endl;

    TLegend *legendRatio = new TLegend(0.65, 0.85, 0.9, 0.92);
    SetLegend(legendRatio);
    legendRatio->SetTextSize(0.04);
    legendRatio->AddEntry(gRatioKstarKaon[0], "K*^{0} / K", "p");

    TLegend *legendRatio2 = new TLegend(0.22, 0.66, 0.5, 0.92);
    SetLegend(legendRatio2);
    legendRatio2->SetTextSize(0.03);
    legendRatio2->AddEntry(gKstarKaRatio_UrQMDOffSmooth, "EPOS4", "l");
    legendRatio2->AddEntry(gKstarKaRatio_UrQMDOnSmooth, "EPOS4 + UrQMD", "l");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarKaModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kKaon], true);
        ScaleGraph(gRatioKstarKaModel, 2.0); // because Kstar yield was already average during dN/dy plotting, but kaon was remaining
        SetGraphStyle(gRatioKstarKaModel, modelStyle[model].color, 1.0, 1, 1.0, 3, modelStyle[model].style);
        gRatioKstarKaModel->Draw("l same");
        legendRatio2->AddEntry(gRatioKstarKaModel, modelLabel[model], "l");
    }
    legendRatio2->AddEntry(gRatioKstarKa_Thermal, "Thermal-FIST", "l");

    legendRatio->Draw();
    // legendRatio2->Draw();
    latex.DrawLatex(0.31, 0.9, "ALICE");
    latex.DrawLatex(0.31, 0.85, "pp, #sqrt{s} = 13.6 TeV");
    latex.DrawLatex(0.31, 0.80, "|y| < 0.5");
    if (isSavePlots)
    {
        cRatioKstarKaon->SaveAs((plotOutputDirectory + "/Ratio_KstarKaon_Run3.pdf").c_str());
    }

    //====================================================
    //   ================Kstar/Phi Ratio ===================
    //====================================================
    TCanvas *cRatioKstarPhi = new TCanvas("cRatioKstarPhi", "cRatioKstarPhi", 720, 720);
    SetPad(cRatioKstarPhi);
    TGraphErrors **gRatioKstarPhi = MakeRatioUncorr(gMYieldKstar, gMYieldPhiRun3, 1.0, 3);
    // TGraphErrors **gRatioKstarPhiRun2 = MakeRatioUncorr(gMYieldKstarRun2, gMYieldPhiRun2, 1.0, 3);
    gRatioKstarPhi[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioKstarPhi[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphStyle(gRatioKstarPhi[0]);
    gRatioKstarPhi[0]->GetYaxis()->SetRangeUser(2.1, 5.15);
    gRatioKstarPhi[0]->GetXaxis()->SetLimits(1.8, 23.5);
    gRatioKstarPhi[0]->SetMarkerColor(kRed);
    gRatioKstarPhi[0]->SetLineColor(kRed);
    gRatioKstarPhi[0]->GetYaxis()->SetNdivisions(505);
    gRatioKstarPhi[0]->Draw("APE");
    gRatioKstarPhi[1]->SetFillStyle(0);
    gRatioKstarPhi[1]->SetLineColor(kRed);
    gRatioKstarPhi[1]->Draw("5 same");
    gRatioKstarPhi[2]->SetFillStyle(1001);
    gRatioKstarPhi[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioKstarPhi[2]->Draw("E2 same");
    // gRatioKstarPhiRun2[0]->SetMarkerStyle(21);
    // gRatioKstarPhiRun2[0]->SetMarkerColor(kBlue);
    // gRatioKstarPhiRun2[0]->SetLineColor(kBlue);
    // gRatioKstarPhiRun2[0]->SetLineWidth(2);
    // gRatioKstarPhiRun2[0]->Draw("PE same");
    // gRatioKstarPhiRun2[1]->SetLineColor(kBlue);
    // gRatioKstarPhiRun2[1]->SetFillStyle(0);
    // gRatioKstarPhiRun2[1]->SetLineWidth(2);
    // gRatioKstarPhiRun2[1]->Draw("5 same");

    // calculate the significance of the difference between HM and LM for K*0/Phi ratio
    int totalPointsKstarPhiRatio = gRatioKstarPhi[0]->GetN();
    double xLM_KstarPhi, xHM_KstarPhi, yHM_KstarPhi, yLM_KstarPhi, yHMError_KstarPhi, yLMError_KstarPhi, yHMStatError_KstarPhi, yLMStatError_KstarPhi, yLMUncorr_KstarPhi, yHMUncorr_KstarPhi;
    gRatioKstarPhi[1]->GetPoint(totalPointsKstarPhiRatio - 1, xLM_KstarPhi, yLM_KstarPhi);
    gRatioKstarPhi[1]->GetPoint(0, xHM_KstarPhi, yHM_KstarPhi);
    yLMError_KstarPhi = gRatioKstarPhi[1]->GetErrorY(totalPointsKstarPhiRatio - 1);
    yHMError_KstarPhi = gRatioKstarPhi[1]->GetErrorY(0);
    yLMStatError_KstarPhi = gRatioKstarPhi[0]->GetErrorY(totalPointsKstarPhiRatio - 1);
    yHMStatError_KstarPhi = gRatioKstarPhi[0]->GetErrorY(0);
    yLMUncorr_KstarPhi = gRatioKstarPhi[2]->GetErrorY(totalPointsKstarPhiRatio - 1);
    yHMUncorr_KstarPhi = gRatioKstarPhi[2]->GetErrorY(0);
    double significance_KstarPhi = fabs(yHM_KstarPhi - yLM_KstarPhi) / sqrt(pow(yHMError_KstarPhi, 2) + pow(yLMError_KstarPhi, 2));
    cout << "Significance of K*0/Phi ratio difference: " << significance_KstarPhi << endl;
    double significance_KstarPhi_uncorr = fabs(yHM_KstarPhi - yLM_KstarPhi) / sqrt(pow(yHMUncorr_KstarPhi, 2) + pow(yLMUncorr_KstarPhi, 2));
    cout << "Significance of K*0/Phi ratio difference (uncorrelated): " << significance_KstarPhi_uncorr << endl;
    double significance_KstarPhi_uncorrStat = fabs(yHM_KstarPhi - yLM_KstarPhi) / sqrt(pow(yHMStatError_KstarPhi, 2) + pow(yLMStatError_KstarPhi, 2) + pow(yHMUncorr_KstarPhi, 2) + pow(yLMUncorr_KstarPhi, 2));
    cout << "Significance of K*0/Phi ratio difference (uncorrelated, statistical): " << significance_KstarPhi_uncorrStat << endl;

    // gKstarPhiRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    // gKstarPhiRatio_UrQMDOffSmooth->SetLineStyle(2);
    // gKstarPhiRatio_UrQMDOffSmooth->SetLineWidth(2);
    // gKstarPhiRatio_UrQMDOffSmooth->Draw("LX same");
    // gKstarPhiRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    // gKstarPhiRatio_UrQMDOnSmooth->SetLineStyle(9);
    // gKstarPhiRatio_UrQMDOnSmooth->SetLineWidth(2);
    // gKstarPhiRatio_UrQMDOnSmooth->Draw("LX same");

    SetGraphStyle(gKstarPhiRatio_EPOSNew_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gKstarPhiRatio_EPOSNew_UrQMDOff->Draw("LX same");
    SetGraphStyle(gKstarPhiRatio_EPOSNew_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gKstarPhiRatio_EPOSNew_UrQMDOn->Draw("LX same");

    gRatioKstarPhi_Thermal->SetLineColor(kMagenta + 2);
    gRatioKstarPhi_Thermal->SetLineStyle(5);
    gRatioKstarPhi_Thermal->SetLineWidth(2);
    gRatioKstarPhi_Thermal->Draw("LX same");

    TLegend *legendRatioKstarPhi = new TLegend(0.52, 0.73, 0.9, 0.94);
    SetLegend(legendRatioKstarPhi);
    legendRatioKstarPhi->SetTextSize(0.03);
    legendRatioKstarPhi->AddEntry(gKstarPhiRatio_EPOSNew_UrQMDOff, "EPOS4", "l");
    legendRatioKstarPhi->AddEntry(gKstarPhiRatio_EPOSNew_UrQMDOn, "EPOS4 + UrQMD", "l");
    legendRatioKstarPhi->AddEntry(gRatioKstarPhi_Thermal, "Thermal-FIST", "l");
    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarPhiModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kPhi], true);
        SetGraphStyle(gRatioKstarPhiModel, modelStyle[model].color, 1.0, 1, 1.0, 3, modelStyle[model].style);
        gRatioKstarPhiModel->Draw("l same");
        legendRatioKstarPhi->AddEntry(gRatioKstarPhiModel, modelLabel[model], "l");
    }

    // legendRatio->Draw();
    // legendRatio2->Draw();
    latex.DrawLatex(0.22, 0.9, "ALICE");
    latex.DrawLatex(0.22, 0.85, "pp, #sqrt{s} = 13.6 TeV");
    latex.DrawLatex(0.22, 0.80, "|y| < 0.5");
    latex.DrawLatex(0.22, 0.75, "K*^{0} / #Phi");

    legendRatioKstarPhi->Draw();

    if (isSavePlots)
    {
        cRatioKstarPhi->SaveAs((plotOutputDirectory + "/Ratio_KstarPhi_Run3.pdf").c_str());
    }

    //====================================================
    //   ================Phi/K Ratio ===================
    //====================================================
    TCanvas *cRatioPhiK = new TCanvas("cRatioPhiK", "cRatioPhiK", 720, 720);
    SetPad(cRatioPhiK);
    TGraphErrors **gRatioPhiK = MakeRatioUncorr(gMYieldPhiRun3, gMYieldKaon, 1.0, 3);
    gRatioPhiK[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioPhiK[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphStyle(gRatioPhiK[0]);
    gRatioPhiK[0]->GetYaxis()->SetRangeUser(0.07, 0.17);
    gRatioPhiK[0]->GetXaxis()->SetLimits(1.8, 23.5);
    // gRatioPhiK[0]->GetXaxis()->SetLimits(1.8, 26.5);
    gRatioPhiK[0]->SetMarkerColor(kRed);
    gRatioPhiK[0]->SetLineColor(kRed);
    gRatioPhiK[0]->GetYaxis()->SetNdivisions(505);
    gRatioPhiK[0]->Draw("APE");
    gRatioPhiK[1]->SetFillStyle(0);
    gRatioPhiK[1]->SetLineColor(kRed);
    gRatioPhiK[1]->Draw("5 same");
    gRatioPhiK[2]->SetFillStyle(1001);
    gRatioPhiK[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioPhiK[2]->Draw("E2 same");

    // Significance of Phi/K ratio difference between HM and LM
    int totalPointsPhiKaonRatio = gRatioPhiK[0]->GetN();
    double xLM_PhiKa, xHM_PhiKa, yHM_PhiKa, yLM_PhiKa, yHMError_PhiKa, yLMError_PhiKa, yHMStatError_PhiKa, yLMStatError_PhiKa, yLMUncorr_PhiKa, yHMUncorr_PhiKa;
    gRatioPhiK[1]->GetPoint(totalPointsPhiKaonRatio - 1, xLM_PhiKa, yLM_PhiKa);
    gRatioPhiK[1]->GetPoint(0, xHM_PhiKa, yHM_PhiKa);
    yLMError_PhiKa = gRatioPhiK[1]->GetErrorY(totalPointsPhiKaonRatio - 1);
    yHMError_PhiKa = gRatioPhiK[1]->GetErrorY(0);
    yLMStatError_PhiKa = gRatioPhiK[0]->GetErrorY(totalPointsPhiKaonRatio - 1);
    yHMStatError_PhiKa = gRatioPhiK[0]->GetErrorY(0);
    yLMUncorr_PhiKa = gRatioPhiK[2]->GetErrorY(totalPointsPhiKaonRatio - 1);
    yHMUncorr_PhiKa = gRatioPhiK[2]->GetErrorY(0);
    double significance_PhiKa = fabs(yHM_PhiKa - yLM_PhiKa) / sqrt(pow(yHMError_PhiKa, 2) + pow(yLMError_PhiKa, 2));
    cout << "Significance of Phi/K ratio difference: " << significance_PhiKa << endl;
    double significance_PhiKa_uncorr = fabs(yHM_PhiKa - yLM_PhiKa) / sqrt(pow(yHMUncorr_PhiKa, 2) + pow(yLMUncorr_PhiKa, 2));
    cout << "Significance of Phi/K ratio difference (uncorrelated): " << significance_PhiKa_uncorr << endl;
    double significance_PhiKa_uncorrStat = fabs(yHM_PhiKa - yLM_PhiKa) / sqrt(pow(yHMStatError_PhiKa, 2) + pow(yLMStatError_PhiKa, 2) + pow(yHMUncorr_PhiKa, 2) + pow(yLMUncorr_PhiKa, 2));
    cout << "Significance of Phi/K ratio difference (uncorrelated, statistical): " << significance_PhiKa_uncorrStat << endl;
    cout << "Phi/K ratio between HM and LM: " << yHM_PhiKa / yLM_PhiKa << endl;

    // gPhiKaRatio_13TeV[0]->SetMarkerStyle(21);
    // gPhiKaRatio_13TeV[0]->SetMarkerColor(kBlue);
    // gPhiKaRatio_13TeV[0]->SetLineColor(kBlue);
    // gPhiKaRatio_13TeV[0]->Draw("P same");
    // gPhiKaRatio_13TeV[1]->SetLineColor(kBlue);
    // gPhiKaRatio_13TeV[1]->SetFillStyle(0);
    // gPhiKaRatio_13TeV[1]->Draw("5 same");

    // gPhiKaRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    // gPhiKaRatio_UrQMDOffSmooth->SetLineStyle(2);
    // gPhiKaRatio_UrQMDOffSmooth->SetLineWidth(2);
    // gPhiKaRatio_UrQMDOffSmooth->Draw("LX same");
    // gPhiKaRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    // gPhiKaRatio_UrQMDOnSmooth->SetLineStyle(9);
    // gPhiKaRatio_UrQMDOnSmooth->SetLineWidth(2);
    // gPhiKaRatio_UrQMDOnSmooth->Draw("LX same");

    SetGraphStyle(gPhiKaRatio_EPOSNew_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gPhiKaRatio_EPOSNew_UrQMDOff->Draw("LX same");
    SetGraphStyle(gPhiKaRatio_EPOSNew_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gPhiKaRatio_EPOSNew_UrQMDOn->Draw("LX same");

    gRatioPhiKa_Thermal->SetLineColor(kMagenta + 2);
    gRatioPhiKa_Thermal->SetLineStyle(5);
    gRatioPhiKa_Thermal->SetLineWidth(2);
    gRatioPhiKa_Thermal->Draw("LX same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioPhiKaModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kKaon], true);
        SetGraphStyle(gRatioPhiKaModel, modelStyle[model].color, 1.0, 1, 1.0, 3, modelStyle[model].style);
        ScaleGraph(gRatioPhiKaModel, 2.0); // because Phi yield was already average during dN/dy plotting, but kaon was remaining
        gRatioPhiKaModel->Draw("l same");
    }

    legendRatio2->Draw();
    latex.SetTextSize(0.035);
    latex.SetTextFont(42);
    latex.DrawLatex(0.68, 0.88, "ALICE");
    latex.DrawLatex(0.68, 0.83, "pp, #sqrt{s} = 13.6 TeV");
    latex.DrawLatex(0.68, 0.78, "|y| < 0.5");
    latex.DrawLatex(0.68, 0.73, "#Phi / K");

    if (isSavePlots)
    {
        cRatioPhiK->SaveAs((plotOutputDirectory + "/Ratio_PhiKaon_Run3.pdf").c_str());
    }

    //==================================================
    // ===================Kstar / Pion Ratio ==================
    //==================================================
    TCanvas *cRatioKstarPion = new TCanvas("cRatioKstarPion", "cRatioKstarPion", 720, 720);
    SetPad(cRatioKstarPion);
    TGraphErrors **gRatioKstarPion = MakeRatioUncorr(gMYieldKstar, gMYieldPion, 1.0, 3);
    gRatioKstarPion[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioKstarPion[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphStyle(gRatioKstarPion[0]);
    // gRatioKstarPion[0]->GetYaxis()->SetRangeUser(0.029, 0.158);
    gRatioKstarPion[0]->GetYaxis()->SetRangeUser(0.031, 0.062);
    gRatioKstarPion[0]->GetXaxis()->SetLimits(1.8, 23.5);
    // gRatioKstarPion[0]->GetXaxis()->SetLimits(1.8, 26.5);
    gRatioKstarPion[0]->SetMarkerColor(kRed);
    gRatioKstarPion[0]->SetLineColor(kRed);
    gRatioKstarPion[0]->GetYaxis()->SetNdivisions(508);
    gRatioKstarPion[0]->GetYaxis()->SetMaxDigits(2);
    gRatioKstarPion[0]->Draw("APE");
    gRatioKstarPion[1]->SetFillStyle(0);
    gRatioKstarPion[1]->SetLineColor(kRed);
    gRatioKstarPion[1]->Draw("5 same");
    gRatioKstarPion[2]->SetFillStyle(1001);
    gRatioKstarPion[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioKstarPion[2]->Draw("E2 same");

    // significance of K*0/Pion ratio difference between HM and LM
    int totalPointsKstarPionRatio = gRatioKstarPion[0]->GetN();
    double xLM_KstarPi, xHM_KstarPi, yHM_KstarPi, yLM_KstarPi, yHMError_KstarPi, yLMError_KstarPi, yHMStatError_KstarPi, yLMStatError_KstarPi, yLMUncorr_KstarPi, yHMUncorr_KstarPi;
    gRatioKstarPion[1]->GetPoint(totalPointsKstarPionRatio - 1, xLM_KstarPi, yLM_KstarPi);
    gRatioKstarPion[1]->GetPoint(0, xHM_KstarPi, yHM_KstarPi);
    yLMError_KstarPi = gRatioKstarPion[1]->GetErrorY(totalPointsKstarPionRatio - 1);
    yHMError_KstarPi = gRatioKstarPion[1]->GetErrorY(0);
    yLMStatError_KstarPi = gRatioKstarPion[0]->GetErrorY(totalPointsKstarPionRatio - 1);
    yHMStatError_KstarPi = gRatioKstarPion[0]->GetErrorY(0);
    yLMUncorr_KstarPi = gRatioKstarPion[2]->GetErrorY(totalPointsKstarPionRatio - 1);
    yHMUncorr_KstarPi = gRatioKstarPion[2]->GetErrorY(0);
    double significance_KstarPi = fabs(yHM_KstarPi - yLM_KstarPi) / sqrt(pow(yHMError_KstarPi, 2) + pow(yLMError_KstarPi, 2));
    cout << "Significance of K*0/Pion ratio difference: " << significance_KstarPi << endl;
    double significance_KstarPi_uncorr = fabs(yHM_KstarPi - yLM_KstarPi) / sqrt(pow(yHMUncorr_KstarPi, 2) + pow(yLMUncorr_KstarPi, 2));
    cout << "Significance of K*0/Pion ratio difference (uncorrelated): " << significance_KstarPi_uncorr << endl;
    double significance_KstarPi_uncorrStat = fabs(yHM_KstarPi - yLM_KstarPi) / sqrt(pow(yHMStatError_KstarPi, 2) + pow(yLMStatError_KstarPi, 2) + pow(yHMUncorr_KstarPi, 2) + pow(yLMUncorr_KstarPi, 2));
    cout << "Significance of K*0/Pion ratio difference (uncorrelated, statistical): " << significance_KstarPi_uncorrStat << endl;
    cout << "K*0/Pion ratio between HM and LM: " << yHM_KstarPi / yLM_KstarPi << endl;

    // gKstarPiRatio_13TeV[0]->SetMarkerStyle(21);
    // gKstarPiRatio_13TeV[0]->SetMarkerColor(kBlue);
    // gKstarPiRatio_13TeV[0]->SetLineColor(kBlue);
    // gKstarPiRatio_13TeV[0]->Draw("P same");
    // gKstarPiRatio_13TeV[1]->SetLineColor(kBlue);
    // gKstarPiRatio_13TeV[1]->SetFillStyle(0);
    // gKstarPiRatio_13TeV[1]->Draw("5 same");

    // gKstarPiRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    // gKstarPiRatio_UrQMDOffSmooth->SetLineStyle(2);
    // gKstarPiRatio_UrQMDOffSmooth->SetLineWidth(2);
    // gKstarPiRatio_UrQMDOffSmooth->Draw("LX same");
    // gKstarPiRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    // gKstarPiRatio_UrQMDOnSmooth->SetLineStyle(9);
    // gKstarPiRatio_UrQMDOnSmooth->SetLineWidth(2);
    // gKstarPiRatio_UrQMDOnSmooth->Draw("LX same");

    SetGraphStyle(gKstarPiRatio_EPOSNew_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gKstarPiRatio_EPOSNew_UrQMDOff->Draw("LX same");
    SetGraphStyle(gKstarPiRatio_EPOSNew_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gKstarPiRatio_EPOSNew_UrQMDOn->Draw("LX same");

    gRatioKstarPi_Thermal->SetLineColor(kMagenta + 2);
    gRatioKstarPi_Thermal->SetLineStyle(5);
    gRatioKstarPi_Thermal->SetLineWidth(2);
    gRatioKstarPi_Thermal->Draw("LX same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarPiModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kPion], true);
        ScaleGraph(gRatioKstarPiModel, 2.0); // because Kstar yield was already average during dN/dy plotting, but pion was remaining
        SetGraphStyle(gRatioKstarPiModel, modelStyle[model].color, 1.0, 1, 1.0, 3, modelStyle[model].style);
        gRatioKstarPiModel->Draw("l same");
    }

    // latex.DrawLatex(0.28, 0.88, "#frac{K^{*0}}{#pi}");
    TLegend *legendRatio3 = new TLegend(0.65, 0.83, 0.75, 0.9);
    SetLegend(legendRatio3);
    legendRatio3->SetTextSize(0.04);
    legendRatio3->Clear();
    legendRatio3->AddEntry(gRatioKstarPion[0], "K*^{0} / #pi", "p");
    legendRatio3->Draw();
    // legendRatio2->Draw();
    latex.SetTextSize(0.035);
    latex.SetTextFont(22);
    latex.DrawLatex(0.22, 0.88, "ALICE");
    latex.DrawLatex(0.22, 0.83, "pp, #sqrt{s} = 13.6 TeV");
    latex.DrawLatex(0.22, 0.78, "|y| < 0.5");

    if (isSavePlots)
    {
        cRatioKstarPion->SaveAs((plotOutputDirectory + "/Ratio_KstarPion_Run3.pdf").c_str());
    }

    //==================================================
    // ===================Phi / Pion Ratio =============
    //==================================================
    TCanvas *cRatioPhiPion = new TCanvas("cRatioPhiPion", "cRatioPhiPion", 720, 720);
    SetPad(cRatioPhiPion);
    TGraphErrors **gRatioPhiPion = MakeRatioUncorr(gMYieldPhiRun3, gMYieldPion, 1.0, 3);
    gRatioPhiPion[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioPhiPion[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphStyle(gRatioPhiPion[0]);
    gRatioPhiPion[0]->GetYaxis()->SetRangeUser(7e-3, 28e-3);
    gRatioPhiPion[0]->GetXaxis()->SetLimits(1.8, 23.5);
    // gRatioPhiPion[0]->GetXaxis()->SetLimits(1.8, 26.5);
    gRatioPhiPion[0]->SetMarkerColor(kRed);
    gRatioPhiPion[0]->SetLineColor(kRed);
    gRatioPhiPion[0]->GetYaxis()->SetNdivisions(508);
    gRatioPhiPion[0]->GetYaxis()->SetMaxDigits(2);
    gRatioPhiPion[0]->Draw("APE");
    gRatioPhiPion[1]->SetFillStyle(0);
    gRatioPhiPion[1]->SetLineColor(kRed);
    gRatioPhiPion[1]->Draw("5 same");
    gRatioPhiPion[2]->SetFillStyle(1001);
    gRatioPhiPion[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioPhiPion[2]->Draw("E2 same");

    // Calculate the significance of the difference between HM and LM for Phi/Pion ratio
    double xLM_PhiPion, xHM_PhiPion, yHM_PhiPion, yLM_PhiPion, yHMError_PhiPion, yLMError_PhiPion, yHMStatError_PhiPion, yLMStatError_PhiPion, yLMUncorr_PhiPion, yHMUncorr_PhiPion;
    gRatioPhiPion[1]->GetPoint(totalPointsKstarPhiRatio - 1, xLM_PhiPion, yLM_PhiPion);
    gRatioPhiPion[1]->GetPoint(0, xHM_PhiPion, yHM_PhiPion);
    yLMError_PhiPion = gRatioPhiPion[1]->GetErrorY(totalPointsKstarPhiRatio - 1);
    yHMError_PhiPion = gRatioPhiPion[1]->GetErrorY(0);
    yLMStatError_PhiPion = gRatioPhiPion[0]->GetErrorY(totalPointsKstarPhiRatio - 1);
    yHMStatError_PhiPion = gRatioPhiPion[0]->GetErrorY(0);
    yLMUncorr_PhiPion = gRatioPhiPion[2]->GetErrorY(totalPointsKstarPhiRatio - 1);
    yHMUncorr_PhiPion = gRatioPhiPion[2]->GetErrorY(0);
    double significance_PhiPion = fabs(yHM_PhiPion - yLM_PhiPion) / sqrt(pow(yHMError_PhiPion, 2) + pow(yLMError_PhiPion, 2));
    cout << "Significance of Phi/Pion ratio difference: " << significance_PhiPion << endl;
    double significance_PhiPion_uncorr = fabs(yHM_PhiPion - yLM_PhiPion) / sqrt(pow(yHMUncorr_PhiPion, 2) + pow(yLMUncorr_PhiPion, 2));
    cout << "Significance of Phi/Pion ratio difference (uncorrelated): " << significance_PhiPion_uncorr << endl;
    double significance_PhiPion_uncorrStat = fabs(yHM_PhiPion - yLM_PhiPion) / sqrt(pow(yHMStatError_PhiPion, 2) + pow(yLMStatError_PhiPion, 2) + pow(yHMUncorr_PhiPion, 2) + pow(yLMUncorr_PhiPion, 2));
    cout << "Significance of Phi/Pion ratio difference (uncorrelated, statistical): " << significance_PhiPion_uncorrStat << endl;
    cout << "Phi/Pion ratio between HM and LM: " << yHM_PhiPion / yLM_PhiPion << endl;

    // ScaleGraph(gPhiPiRatio_13TeV[0], 2.0); // To average the pion yield
    // ScaleGraph(gPhiPiRatio_13TeV[1], 2.0);
    // gPhiPiRatio_13TeV[0]->SetMarkerStyle(21);
    // gPhiPiRatio_13TeV[0]->SetMarkerColor(kBlue);
    // gPhiPiRatio_13TeV[0]->SetLineColor(kBlue);
    // gPhiPiRatio_13TeV[0]->Draw("P same");
    // gPhiPiRatio_13TeV[1]->SetLineColor(kBlue);
    // gPhiPiRatio_13TeV[1]->SetFillStyle(0);
    // gPhiPiRatio_13TeV[1]->Draw("5 same");

    // gPhiPiRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    // gPhiPiRatio_UrQMDOffSmooth->SetLineStyle(2);
    // gPhiPiRatio_UrQMDOffSmooth->SetLineWidth(2);
    // gPhiPiRatio_UrQMDOffSmooth->Draw("LX same");
    // gPhiPiRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    // gPhiPiRatio_UrQMDOnSmooth->SetLineStyle(9);
    // gPhiPiRatio_UrQMDOnSmooth->SetLineWidth(2);
    // gPhiPiRatio_UrQMDOnSmooth->Draw("LX same");

    SetGraphStyle(gPhiPiRatio_EPOSNew_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gPhiPiRatio_EPOSNew_UrQMDOff->Draw("LX same");
    SetGraphStyle(gPhiPiRatio_EPOSNew_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gPhiPiRatio_EPOSNew_UrQMDOn->Draw("LX same");

    gRatioPhiPi_Thermal->SetLineColor(kMagenta + 2);
    gRatioPhiPi_Thermal->SetLineStyle(5);
    gRatioPhiPi_Thermal->SetLineWidth(2);
    gRatioPhiPi_Thermal->Draw("LX same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioPhiPiModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kPion], true);
        ScaleGraph(gRatioPhiPiModel, 2.0); // To average the pion yield
        SetGraphStyle(gRatioPhiPiModel, modelStyle[model].color, 1.0, 1, 1.0, 3, modelStyle[model].style);
        gRatioPhiPiModel->Draw("l same");
    }

    legendRatio2->Draw();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.82, 0.85, "#Phi/#pi");

    if (isSavePlots)
        cRatioPhiPion->SaveAs((plotOutputDirectory + "/Ratio_PhiPion_Run3.pdf").c_str());

    //==================================================
    // ===================Phi / Xi Ratio =============
    //==================================================
    TCanvas *cRatioPhiXi = new TCanvas("cRatioPhiXi", "cRatioPhiXi", 720, 720);
    SetPad(cRatioPhiXi);
    TGraphErrors **gRatioPhiXi = MakeRatioUncorr(gMYieldPhiRun3, gMYieldXiRun3, 1.0, 3);
    gRatioPhiXi[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioPhiXi[0]->GetYaxis()->SetTitle("#Phi / #Xi");
    SetGraphStyle(gRatioPhiXi[0]);
    // gRatioPhiXi[0]->GetYaxis()->SetRangeUser(1.5, 2.6);
    gRatioPhiXi[0]->GetYaxis()->SetRangeUser(1.1, 4.3);
    gRatioPhiXi[0]->GetXaxis()->SetLimits(1.8, 23.5);
    gRatioPhiXi[0]->SetMarkerColor(kRed);
    gRatioPhiXi[0]->SetLineColor(kRed);
    gRatioPhiXi[0]->GetYaxis()->SetNdivisions(505);
    gRatioPhiXi[0]->Draw("APE");
    gRatioPhiXi[1]->SetFillStyle(0);
    gRatioPhiXi[1]->SetLineColor(kRed);
    gRatioPhiXi[1]->Draw("5 same");
    gRatioPhiXi[2]->SetFillStyle(1001);
    gRatioPhiXi[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioPhiXi[2]->Draw("E2 same");

    // ScaleGraph(gPhiXiRatio_UrQMDOffSmooth, 0.5); // To average the Xi yield
    // ScaleGraph(gPhiXiRatio_UrQMDOnSmooth, 0.5);
    // gPhiXiRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    // gPhiXiRatio_UrQMDOffSmooth->SetLineStyle(2);
    // gPhiXiRatio_UrQMDOffSmooth->SetLineWidth(2);
    // gPhiXiRatio_UrQMDOffSmooth->Draw("LX same");
    // gPhiXiRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    // gPhiXiRatio_UrQMDOnSmooth->SetLineStyle(9);
    // gPhiXiRatio_UrQMDOnSmooth->SetLineWidth(2);
    // gPhiXiRatio_UrQMDOnSmooth->Draw("LX same");

    ScaleGraph(gPhiXiRatio_EPOSNew_UrQMDOff, 0.5); // To average the Xi yield
    ScaleGraph(gPhiXiRatio_EPOSNew_UrQMDOn, 0.5);
    SetGraphStyle(gPhiXiRatio_EPOSNew_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gPhiXiRatio_EPOSNew_UrQMDOff->Draw("LX same");
    SetGraphStyle(gPhiXiRatio_EPOSNew_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gPhiXiRatio_EPOSNew_UrQMDOn->Draw("LX same");

    gRatioPhiXi_Thermal->SetLineColor(kMagenta + 2);
    gRatioPhiXi_Thermal->SetLineStyle(5);
    gRatioPhiXi_Thermal->SetLineWidth(2);
    gRatioPhiXi_Thermal->Draw("LX same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioPhiXiModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kXi], true);
        ScaleGraph(gRatioPhiXiModel, 1.0); // To average the Xi yield
        SetGraphStyle(gRatioPhiXiModel, modelStyle[model].color, 1.0, 1, 1.0, 3, modelStyle[model].style);
        gRatioPhiXiModel->Draw("l same");
    }

    TLegend *legendRatioPhiXi = new TLegend(0.65, 0.75, 0.88, 0.9);
    SetLegend(legendRatioPhiXi);
    legendRatioPhiXi->SetTextSize(0.035);
    // legendRatioPhiXi->AddEntry(gRatioPhiXi[0], "#Phi / #Xi", "p");
    legendRatioPhiXi->AddEntry(gPhiXiRatio_EPOSNew_UrQMDOff, "EPOS4", "l");
    legendRatioPhiXi->AddEntry(gPhiXiRatio_EPOSNew_UrQMDOn, "EPOS4 + UrQMD", "l");
    legendRatioPhiXi->Draw();

    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.28, 0.88, "ALICE");
    latex.DrawLatex(0.28, 0.83, "pp, #sqrt{s} = 13.6 TeV");
    latex.DrawLatex(0.28, 0.78, "|y| < 0.5");

    if (isSavePlots)
        cRatioPhiXi->SaveAs((plotOutputDirectory + "/Ratio_PhiXi_Run3.pdf").c_str());

    //======================================================
    //================Kstar/Kshort Ratio====================
    //======================================================
    TCanvas *cRatioKstarKshort = new TCanvas("cRatioKstarKshort", "cRatioKstarKshort", 720, 720);
    SetPad(cRatioKstarKshort);
    TGraphErrors **gRatioKstarKshort = MakeRatioUncorr(gMYieldKstar, gMYieldKshortRun3, 1.0, 3);
    gRatioKstarKshort[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioKstarKshort[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphStyle(gRatioKstarKshort[0]);
    gRatioKstarKshort[0]->GetYaxis()->SetRangeUser(0.245, 0.54);
    gRatioKstarKshort[0]->GetXaxis()->SetLimits(1.8, 23.5);
    gRatioKstarKshort[0]->SetMarkerColor(kRed);
    gRatioKstarKshort[0]->SetLineColor(kRed);
    gRatioKstarKshort[0]->GetYaxis()->SetNdivisions(505);
    gRatioKstarKshort[0]->GetYaxis()->SetMaxDigits(3);
    gRatioKstarKshort[0]->Draw("APE");
    gRatioKstarKshort[1]->SetFillStyle(0);
    gRatioKstarKshort[1]->SetLineColor(kRed);
    gRatioKstarKshort[1]->Draw("5 same");
    gRatioKstarKshort[2]->SetFillStyle(1001);
    gRatioKstarKshort[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioKstarKshort[2]->Draw("E2 same");

    // cout << "\n Information for K*0 and Kshort yields and their errors " << endl;
    // for (int i3 = 0; i3 < TotalPointsRun3; i3++)
    // {
    //     double xkstar, ykstar, ysyskstar, yuncorrkstar;
    //     double xkaon, ykaon, ysyskaon, yuncorrkaon;
    //     gMYieldKstar[0]->GetPoint(i3, xkstar, ykstar);
    //     ysyskstar = gMYieldKstar[1]->GetErrorY(i3);
    //     yuncorrkstar = gMYieldKstar[2]->GetErrorY(i3);

    //     gMYieldKshortRun3[0]->GetPoint(i3, xkaon, ykaon);
    //     ysyskaon = gMYieldKshortRun3[1]->GetErrorY(i3);
    //     yuncorrkaon = gMYieldKshortRun3[2]->GetErrorY(i3);

    //     cout << "Point " << i3 << ": K*0 yield = " << ykstar << " +/- " << ysyskstar << " (sys) +/- " << yuncorrkstar << " (uncorr)" << endl;
    //     cout << "Point " << i3 << ": Kshort yield = " << ykaon << " +/- " << ysyskaon << " (sys) +/- " << yuncorrkaon << " (uncorr)" << endl;
    //     double ratio = ykstar / ykaon;
    //     double ratioErrorsys = ratio * sqrt(pow(ysyskstar / ykstar, 2) + pow(ysyskaon / ykaon, 2));
    //     double ratioErroruncorr = ratio * sqrt(pow(yuncorrkstar / ykstar, 2) + pow(yuncorrkaon / ykaon, 2));
    //     cout << "Point " << i3 << ": K*0/Kshort ratio = " << ratio << " +/- " << ratioErrorsys << " (sys) +/- " << ratioErroruncorr << " (uncorr)" << endl;
    //     cout << endl;
    // }

    ////Calculation of suppression significance between HM and LM for K*0/Kshort ratio
    int totalPointsKstarKshortRatio = gRatioKstarKshort[0]->GetN();
    double xLM_KstarKshort, xHM_KstarKshort, yHM_KstarKshort, yLM_KstarKshort, yHMError_KstarKshort, yLMError_KstarKshort, yHMStatError_KstarKshort, yLMStatError_KstarKshort, yLMUncorr_KstarKshort, yHMUncorr_KstarKshort;
    gRatioKstarKshort[1]->GetPoint(totalPointsKstarKshortRatio - 1, xLM_KstarKshort, yLM_KstarKshort);
    gRatioKstarKshort[1]->GetPoint(0, xHM_KstarKshort, yHM_KstarKshort);
    yLMError_KstarKshort = gRatioKstarKshort[1]->GetErrorY(totalPointsKstarKshortRatio - 1);
    yHMError_KstarKshort = gRatioKstarKshort[1]->GetErrorY(0);
    yLMStatError_KstarKshort = gRatioKstarKshort[0]->GetErrorY(totalPointsKstarKshortRatio - 1);
    yHMStatError_KstarKshort = gRatioKstarKshort[0]->GetErrorY(0);
    yLMUncorr_KstarKshort = gRatioKstarKshort[2]->GetErrorY(totalPointsKstarKshortRatio - 1);
    yHMUncorr_KstarKshort = gRatioKstarKshort[2]->GetErrorY(0);
    double sigmaDiffHM_LM_KstarKshort = abs(yHM_KstarKshort - yLM_KstarKshort) / sqrt(pow(yHMError_KstarKshort, 2) + pow(yLMError_KstarKshort, 2) + pow(yHMStatError_KstarKshort, 2) + pow(yLMStatError_KstarKshort, 2));
    double sigmaDiffHM_LM_uncorr_KstarKshort = abs(yHM_KstarKshort - yLM_KstarKshort) / sqrt(pow(yHMUncorr_KstarKshort, 2) + pow(yLMUncorr_KstarKshort, 2));
    cout << "K*0/Kshort ratio difference between HM and LM: " << sigmaDiffHM_LM_KstarKshort << " sigma" << endl;
    cout << "K*0/Kshort ratio difference between HM and LM (uncorrelated): " << sigmaDiffHM_LM_uncorr_KstarKshort << " sigma" << endl;

    // gKstarK0sRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    // gKstarK0sRatio_UrQMDOffSmooth->SetLineStyle(2);
    // gKstarK0sRatio_UrQMDOffSmooth->SetLineWidth(2);
    // gKstarK0sRatio_UrQMDOffSmooth->Draw("LX same");
    // gKstarK0sRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    // gKstarK0sRatio_UrQMDOnSmooth->SetLineStyle(9);
    // gKstarK0sRatio_UrQMDOnSmooth->SetLineWidth(2);
    // gKstarK0sRatio_UrQMDOnSmooth->Draw("LX same");

    SetGraphStyle(gKstarK0sRatio_EPOSNew_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gKstarK0sRatio_EPOSNew_UrQMDOff->Draw("LX same");
    SetGraphStyle(gKstarK0sRatio_EPOSNew_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gKstarK0sRatio_EPOSNew_UrQMDOn->Draw("LX same");

    gRatioKstarK0s_Thermal->SetLineColor(kMagenta + 2);
    gRatioKstarK0s_Thermal->SetLineStyle(5);
    gRatioKstarK0s_Thermal->SetLineWidth(2);
    gRatioKstarK0s_Thermal->Draw("LX same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioKstarKshortModel = MakeRatio(gMYield[model][kKstar], gMYield[model][kKshort], true);
        SetGraphStyle(gRatioKstarKshortModel, modelStyle[model].color, 1.0, 1, 1.0, 3, modelStyle[model].style);
        gRatioKstarKshortModel->Draw("l same");
    }

    TLegend *legendRatio4 = new TLegend(0.75, 0.83, 0.9, 0.9);
    SetLegend(legendRatio4);
    legendRatio4->SetTextSize(0.04);
    legendRatio4->AddEntry(gRatioKstarKshort[0], "K*^{0} / K_{S}^{0}", "p");
    legendRatio4->Draw();
    legendRatio2->Draw();

    if (isSavePlots)
    {
        cRatioKstarKshort->SaveAs((plotOutputDirectory + "/Ratio_KstarKshort_Run3.pdf").c_str());
    }

    //=====================================================================
    // ===================Phi/K0s ratio ================================
    //=====================================================================
    TCanvas *cRatioPhiK0s = new TCanvas("cRatioPhiK0s", "cRatioPhiK0s", 720, 720);
    SetPad(cRatioPhiK0s);
    TGraphErrors **gRatioPhiK0s = MakeRatioUncorr(gMYieldPhiRun3, gMYieldKshortRun3, 1.0, 3);
    gRatioPhiK0s[0]->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5} ");
    gRatioPhiK0s[0]->GetYaxis()->SetTitle("dN / dy ratio");
    SetGraphStyle(gRatioPhiK0s[0]);
    gRatioPhiK0s[0]->GetYaxis()->SetRangeUser(0.07, 0.18);
    gRatioPhiK0s[0]->GetXaxis()->SetLimits(1.8, 23.5);
    gRatioPhiK0s[0]->SetMarkerColor(kRed);
    gRatioPhiK0s[0]->SetLineColor(kRed);
    gRatioPhiK0s[0]->GetYaxis()->SetNdivisions(505);
    gRatioPhiK0s[0]->GetYaxis()->SetMaxDigits(3);
    gRatioPhiK0s[0]->Draw("APE");
    gRatioPhiK0s[1]->SetFillStyle(0);
    gRatioPhiK0s[1]->SetLineColor(kRed);
    gRatioPhiK0s[1]->Draw("5 same");
    gRatioPhiK0s[2]->SetFillStyle(1001);
    gRatioPhiK0s[2]->SetFillColorAlpha(kRed, 0.35);
    gRatioPhiK0s[2]->Draw("E2 same");

    // significance calculation for Phi/K0s ratio between HM and LM
    int totalPointsPhiK0sRatio = gRatioPhiK0s[0]->GetN();
    double xLM_PhiK0s, xHM_PhiK0s, yHM_PhiK0s, yLM_PhiK0s, yHMError_PhiK0s, yLMError_PhiK0s, yHMStatError_PhiK0s, yLMStatError_PhiK0s, yLMUncorr_PhiK0s, yHMUncorr_PhiK0s;
    gRatioPhiK0s[1]->GetPoint(totalPointsPhiK0sRatio - 1, xLM_PhiK0s, yLM_PhiK0s);
    gRatioPhiK0s[1]->GetPoint(0, xHM_PhiK0s, yHM_PhiK0s);
    yLMError_PhiK0s = gRatioPhiK0s[1]->GetErrorY(totalPointsPhiK0sRatio - 1);
    yHMError_PhiK0s = gRatioPhiK0s[1]->GetErrorY(0);
    yLMStatError_PhiK0s = gRatioPhiK0s[0]->GetErrorY(totalPointsPhiK0sRatio - 1);
    yHMStatError_PhiK0s = gRatioPhiK0s[0]->GetErrorY(0);
    yLMUncorr_PhiK0s = gRatioPhiK0s[2]->GetErrorY(totalPointsPhiK0sRatio - 1);
    yHMUncorr_PhiK0s = gRatioPhiK0s[2]->GetErrorY(0);
    double sigmaDiffHM_LM_PhiK0s = abs(yHM_PhiK0s - yLM_PhiK0s) / sqrt(pow(yHMError_PhiK0s, 2) + pow(yLMError_PhiK0s, 2) + pow(yHMStatError_PhiK0s, 2) + pow(yLMStatError_PhiK0s, 2));
    double sigmaDiffHM_LM_uncorr_PhiK0s = abs(yHM_PhiK0s - yLM_PhiK0s) / sqrt(pow(yHMUncorr_PhiK0s, 2) + pow(yLMUncorr_PhiK0s, 2));
    cout << "Phi/K0s ratio difference between HM and LM: " << sigmaDiffHM_LM_PhiK0s << " sigma" << endl;
    cout << "Phi/K0s ratio difference between HM and LM (uncorrelated): " << sigmaDiffHM_LM_uncorr_PhiK0s << " sigma" << endl;
    cout << "Phi/K0s ratio between HM and LM: " << yHM_PhiK0s / yLM_PhiK0s << endl;

    // gPhiK0sRatio_UrQMDOffSmooth->SetLineColor(kGreen + 2);
    // gPhiK0sRatio_UrQMDOffSmooth->SetLineStyle(2);
    // gPhiK0sRatio_UrQMDOffSmooth->SetLineWidth(2);
    // gPhiK0sRatio_UrQMDOffSmooth->Draw("LX same");
    // gPhiK0sRatio_UrQMDOnSmooth->SetLineColor(kBlue - 6);
    // gPhiK0sRatio_UrQMDOnSmooth->SetLineStyle(9);
    // gPhiK0sRatio_UrQMDOnSmooth->SetLineWidth(2);
    // gPhiK0sRatio_UrQMDOnSmooth->Draw("LX same");

    SetGraphStyle(gPhiK0sRatio_EPOSNew_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gPhiK0sRatio_EPOSNew_UrQMDOff->Draw("LX same");
    SetGraphStyle(gPhiK0sRatio_EPOSNew_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gPhiK0sRatio_EPOSNew_UrQMDOn->Draw("LX same");

    gRatioPhiK0s_Thermal->SetLineColor(kMagenta + 2);
    gRatioPhiK0s_Thermal->SetLineStyle(5);
    gRatioPhiK0s_Thermal->SetLineWidth(2);
    gRatioPhiK0s_Thermal->Draw("LX same");

    for (auto model : modelsToPlot)
    {
        TGraphErrors *gRatioPhiKshortModel = MakeRatio(gMYield[model][kPhi], gMYield[model][kKshort], true);
        SetGraphStyle(gRatioPhiKshortModel, modelStyle[model].color, 1.0, 1, 1.0, 3, modelStyle[model].style);
        gRatioPhiKshortModel->Draw("l same");
    }

    TLegend *legendRatio5 = new TLegend(0.75, 0.83, 0.9, 0.9);
    SetLegend(legendRatio5);
    legendRatio5->SetTextSize(0.04);
    legendRatio5->AddEntry(gRatioPhiK0s[0], "#Phi / K^{0}_{S}", "p");
    legendRatio5->Draw();
    legendRatio2->Draw();

    if (isSavePlots)
        cRatioPhiK0s->SaveAs((plotOutputDirectory + "/Ratio_PhiK0s_Run3.pdf").c_str());

    //======================================================================
    // ========== Double Yield Ratio (Mult/LM): Kstar/K, Kstar/Pi ==========
    //======================================================================
    TCanvas *cYieldLMRatio = new TCanvas("cYieldLMRatio", "cYieldLMRatio", 720, 720);
    SetPad(cYieldLMRatio);

    TGraphErrors **gYieldKstarKaLMRatio = DivideByMult(gRatioKstarKaon, -1, 0.5, 3);
    TGraphErrors **gYieldKstarPiLMRatio = DivideByMult(gRatioKstarPion, -1, 0.5, 3);

    // Frame
    TH1D *hFrameYieldLMRatio = new TH1D("hFrameYieldLMRatio", "", 100, 1.8, 23.5);
    SetHistStyle(hFrameYieldLMRatio);
    hFrameYieldLMRatio->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    hFrameYieldLMRatio->GetYaxis()->SetTitle("Y / Y_{LM}");
    hFrameYieldLMRatio->GetYaxis()->SetRangeUser(0.75, 1.33);
    hFrameYieldLMRatio->GetXaxis()->SetTickLength(0.02);
    hFrameYieldLMRatio->GetYaxis()->SetTickLength(0.02);
    hFrameYieldLMRatio->Draw();

    // K*0 / K
    SetGraphStyle(gYieldKstarKaLMRatio[0], colors[5], 1.0, 20, 1.0, 1);
    gYieldKstarKaLMRatio[0]->Draw("PE SAME");

    SetGraphStyle(gYieldKstarKaLMRatio[1], colors[5], 1.0, 20, 1.0, 1);
    gYieldKstarKaLMRatio[1]->SetFillStyle(0);
    gYieldKstarKaLMRatio[1]->Draw("5 SAME");

    // K*0 / pi
    SetGraphStyle(gYieldKstarPiLMRatio[0], colors[4], 1.0, 21, 1.0, 1);
    gYieldKstarPiLMRatio[0]->Draw("PE SAME");

    SetGraphStyle(gYieldKstarPiLMRatio[1], colors[4], 1.0, 21, 1.0, 1);
    gYieldKstarPiLMRatio[1]->SetFillStyle(0);
    gYieldKstarPiLMRatio[1]->Draw("5 SAME");

    // EPOS model ratios
    TGraphErrors *gEPOSNew_KstarKaLMRatio_UrQMDOff = DivideByMultModel(gKstarKaRatio_UrQMDOffSmooth, 2.5);
    TGraphErrors *gEPOSNew_KstarKaLMRatio_UrQMDOn = DivideByMultModel(gKstarKaRatio_UrQMDOnSmooth, 2.5);
    TGraphErrors *gEPOSNew_KstarPiLMRatio_UrQMDOff = DivideByMultModel(gKstarPiRatio_UrQMDOffSmooth, 2.5);
    TGraphErrors *gEPOSNew_KstarPiLMRatio_UrQMDOn = DivideByMultModel(gKstarPiRatio_UrQMDOnSmooth, 2.5);

    SetGraphStyle(gEPOSNew_KstarKaLMRatio_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gEPOSNew_KstarKaLMRatio_UrQMDOff->Draw("LX SAME");

    SetGraphStyle(gEPOSNew_KstarKaLMRatio_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gEPOSNew_KstarKaLMRatio_UrQMDOn->Draw("LX SAME");

    SetGraphStyle(gEPOSNew_KstarPiLMRatio_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gEPOSNew_KstarPiLMRatio_UrQMDOff->Draw("LX SAME");

    SetGraphStyle(gEPOSNew_KstarPiLMRatio_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gEPOSNew_KstarPiLMRatio_UrQMDOn->Draw("LX SAME");

    // Data legend
    TLegend *legendYieldLMRatio = new TLegend(0.2, 0.89, 0.95, 0.98);
    SetLegend(legendYieldLMRatio);
    legendYieldLMRatio->SetTextSize(0.035);
    legendYieldLMRatio->SetNColumns(2);
    legendYieldLMRatio->AddEntry(gYieldKstarPiLMRatio[0], "K*^{0} / #pi", "p");
    legendYieldLMRatio->AddEntry(gYieldKstarKaLMRatio[0], "K*^{0} / K", "p");
    legendYieldLMRatio->Draw();

    // Model legend
    TLegend *legendYieldLMRatio2 = new TLegend(0.2, 0.78, 0.95, 0.92);
    SetLegend(legendYieldLMRatio2);
    legendYieldLMRatio2->SetTextSize(0.028);
    legendYieldLMRatio2->SetNColumns(2);
    legendYieldLMRatio2->AddEntry(gEPOSNew_KstarPiLMRatio_UrQMDOff, "EPOS4", "l");
    legendYieldLMRatio2->AddEntry(gEPOSNew_KstarKaLMRatio_UrQMDOff, "EPOS4", "l");
    legendYieldLMRatio2->AddEntry(gEPOSNew_KstarPiLMRatio_UrQMDOn, "EPOS4 + UrQMD", "l");
    legendYieldLMRatio2->AddEntry(gEPOSNew_KstarKaLMRatio_UrQMDOn, "EPOS4 + UrQMD", "l");
    legendYieldLMRatio2->Draw();

    if (isSavePlots)
    {
        cYieldLMRatio->SaveAs((plotOutputDirectory + "/Yield_LMRatio_KstarKaPi_Run3.pdf").c_str());
    }

    //======================================================================
    // ========== Double Yield Ratio (Mult/LM): Phi/K, Phi/Pi ==============
    //======================================================================
    TCanvas *cYieldLMRatioPhi = new TCanvas("cYieldLMRatioPhi", "cYieldLMRatioPhi", 720, 720);
    SetPad(cYieldLMRatioPhi);

    TGraphErrors **gYieldPhiKaLMRatio = DivideByMult(gRatioPhiK, -1, 0.5, 3);
    TGraphErrors **gYieldPhiPiLMRatio = DivideByMult(gRatioPhiPion, -1, 0.5, 3);

    // Frame
    TH1D *hFrameYieldLMRatioPhi = new TH1D("hFrameYieldLMRatioPhi", "", 100, 1.8, 23.5);
    SetHistStyle(hFrameYieldLMRatioPhi);
    hFrameYieldLMRatioPhi->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
    hFrameYieldLMRatioPhi->GetYaxis()->SetTitle("Y / Y_{LM}");
    hFrameYieldLMRatioPhi->GetYaxis()->SetRangeUser(0.95, 1.60);
    hFrameYieldLMRatioPhi->GetXaxis()->SetTickLength(0.02);
    hFrameYieldLMRatioPhi->GetYaxis()->SetTickLength(0.02);
    hFrameYieldLMRatioPhi->Draw();

    // Reference line at unity
    DrawUnityLine(1.8, 23.5, 2, 2, kBlack);

    // Phi / K
    SetGraphStyle(gYieldPhiKaLMRatio[0], colors[5], 1.0, 20, 1.0, 1);
    gYieldPhiKaLMRatio[0]->Draw("PE SAME");

    SetGraphStyle(gYieldPhiKaLMRatio[1], colors[5], 1.0, 20, 1.0, 1);
    gYieldPhiKaLMRatio[1]->SetFillStyle(0);
    gYieldPhiKaLMRatio[1]->Draw("5 SAME");

    // Phi / pi
    SetGraphStyle(gYieldPhiPiLMRatio[0], colors[4], 1.0, 21, 1.0, 1);
    gYieldPhiPiLMRatio[0]->Draw("PE SAME");

    SetGraphStyle(gYieldPhiPiLMRatio[1], colors[4], 1.0, 21, 1.0, 1);
    gYieldPhiPiLMRatio[1]->SetFillStyle(0);
    gYieldPhiPiLMRatio[1]->Draw("5 SAME");

    // EPOS model ratios
    TGraphErrors *gEPOS_PhiKaLMRatio_UrQMDOff = DivideByMultModel(gPhiKaRatio_UrQMDOffSmooth, 2.5);
    TGraphErrors *gEPOS_PhiKaLMRatio_UrQMDOn = DivideByMultModel(gPhiKaRatio_UrQMDOnSmooth, 2.5);
    TGraphErrors *gEPOS_PhiPiLMRatio_UrQMDOff = DivideByMultModel(gPhiPiRatio_UrQMDOffSmooth, 2.5);
    TGraphErrors *gEPOS_PhiPiLMRatio_UrQMDOn = DivideByMultModel(gPhiPiRatio_UrQMDOnSmooth, 2.5);

    SetGraphStyle(gEPOS_PhiKaLMRatio_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gEPOS_PhiKaLMRatio_UrQMDOff->Draw("LX SAME");

    SetGraphStyle(gEPOS_PhiKaLMRatio_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gEPOS_PhiKaLMRatio_UrQMDOn->Draw("LX SAME");

    SetGraphStyle(gEPOS_PhiPiLMRatio_UrQMDOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
    gEPOS_PhiPiLMRatio_UrQMDOff->Draw("LX SAME");

    SetGraphStyle(gEPOS_PhiPiLMRatio_UrQMDOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
    gEPOS_PhiPiLMRatio_UrQMDOn->Draw("LX SAME");

    // Legend
    TLegend *legendYieldLMRatioPhi = new TLegend(0.2, 0.78, 0.95, 0.95);
    SetLegend(legendYieldLMRatioPhi);
    legendYieldLMRatioPhi->SetTextSize(0.028);
    legendYieldLMRatioPhi->SetNColumns(2);

    legendYieldLMRatioPhi->AddEntry(gYieldPhiPiLMRatio[0], "#phi / #pi", "p");
    legendYieldLMRatioPhi->AddEntry(gYieldPhiKaLMRatio[0], "#phi / K", "p");
    legendYieldLMRatioPhi->AddEntry(gEPOS_PhiKaLMRatio_UrQMDOff, "#phi / K (UrQMD Off)", "l");
    legendYieldLMRatioPhi->AddEntry(gEPOS_PhiPiLMRatio_UrQMDOff, "#phi / #pi (UrQMD Off)", "l");
    legendYieldLMRatioPhi->AddEntry(gEPOS_PhiKaLMRatio_UrQMDOn, "#phi / K (UrQMD On)", "l");
    legendYieldLMRatioPhi->AddEntry(gEPOS_PhiPiLMRatio_UrQMDOn, "#phi / #pi (UrQMD On)", "l");

    legendYieldLMRatioPhi->Draw("SAME");

    if (isSavePlots)
    {
        cYieldLMRatioPhi->SaveAs((plotOutputDirectory + "/Yield_LMRatio_PhiKaPi_Run3.pdf").c_str());
    }

    // ======================================================
    //       FOUR-PANEL PARTICLE-RATIO COMPARISON
    //
    //       (a) K*0 / K       (b) Phi / K
    //       (c) K*0 / pi      (d) Phi / pi
    // ======================================================

    TCanvas *cParticleRatios = new TCanvas("cParticleRatios", "Particle ratios", 900, 900);
    SetPad(cParticleRatios);
    cParticleRatios->Divide(2, 2, 0.0, 0.0);

    // Pad 1: Kstar / K
    cParticleRatios->cd(1);

    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.005);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.0);
    gPad->SetTicks(1, 1);

    // Pad 2: Phi / K

    cParticleRatios->cd(2);

    gPad->SetLeftMargin(0.005);
    gPad->SetRightMargin(0.16);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.0);
    gPad->SetTicks(1, 1);

    // Pad 3: Kstar / pi

    cParticleRatios->cd(3);

    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.005);
    gPad->SetTopMargin(0.00);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicks(1, 1);

    // Pad 4: Phi / pi

    cParticleRatios->cd(4);

    gPad->SetLeftMargin(0.005);
    gPad->SetRightMargin(0.16);
    gPad->SetTopMargin(0.0);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicks(1, 1);

    // ------------------------------------------------------
    // Helper to draw one ratio panel
    // ------------------------------------------------------

    auto drawRatioPanel = [&](int pad,
                              const char *frameName,
                              const char *panelLabel,
                              TGraphErrors **data,
                              TGraphErrors *eposOff,
                              TGraphErrors *eposOn,
                              int numerator,
                              int denominator,
                              double modelScale,
                              double yMin,
                              double yMax,
                              bool rightAxis,
                              bool showXaxis,
                              bool showYaxis)
    {
        cParticleRatios->cd(pad);

        // Frame
        TH1D *hFrame = new TH1D(
            frameName, "", 100, 1.8, 23.5);

        SetHistStyle(hFrame);

        if (!showXaxis)
        {
            hFrame->GetXaxis()->SetTitleSize(0.0);
        }
        if (!showYaxis)
        {
            hFrame->GetYaxis()->SetTitleSize(0.0);
        }

        hFrame->GetXaxis()->SetTitle("#LT d#it{N_{#it{ch}}/#it{d}#eta} #GT_{|#it{#eta}|<0.5}");
        hFrame->GetYaxis()->SetTitle("d#it{N}/d#it{y} ratio");

        hFrame->GetXaxis()->SetLimits(1.8, 23.5);
        hFrame->GetYaxis()->SetRangeUser(yMin, yMax);

        hFrame->GetYaxis()->SetMaxDigits(3);

        hFrame->GetYaxis()->SetTitleOffset(
            rightAxis ? 1.25 : 1.35);

        // Draw the right-hand y axis for the right column
        if (rightAxis)
            hFrame->Draw("Y+");
        else
            hFrame->Draw();

        // --------------------------------------------------
        // Run 3 data
        // --------------------------------------------------

        SetGraphStyle(data[0], colors[5], 1.0, 20, 1.0, 1, 1);
        data[0]->Draw("PE SAME");

        // Total systematic uncertainty
        SetGraphStyle(data[1], colors[5], 1.0, 20, 1.0, 1, 1);
        data[1]->SetFillStyle(0);
        data[1]->Draw("5 SAME");

        // Uncorrelated systematic uncertainty
        SetGraphStyle(data[2], colors[5], 1.0, 20, 1.0, 1, 1);
        data[2]->SetFillStyle(1001);
        data[2]->SetFillColorAlpha(colors[5], 0.35);
        data[2]->Draw("E2 SAME");

        // --------------------------------------------------
        // EPOS4 models
        // --------------------------------------------------

        SetGraphStyle(eposOff, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 2, lineStyleEPOSWithoutUrQMD);
        eposOff->Draw("LX SAME");

        SetGraphStyle(eposOn, colorEPOSWithUrQMD, 1.0, 1, 1.0, 2, lineStyleEPOSWithUrQMD);
        eposOn->Draw("LX SAME");

        // --------------------------------------------------
        // PYTHIA models
        // --------------------------------------------------

        std::vector<TGraphErrors *> modelGraphs;

        for (auto model : modelsToPlot)
        {
            TGraphErrors *gRatioModel = MakeRatio(
                gMYield[model][numerator],
                gMYield[model][denominator],
                true);

            // Apply the same particle/antiparticle averaging
            // convention as in the individual ratio plots
            ScaleGraph(gRatioModel, modelScale);

            SetGraphStyle(gRatioModel, modelStyle[model].color, 1.0, 1, 1.0, 2, modelStyle[model].style);
            gRatioModel->Draw("L SAME");

            modelGraphs.push_back(gRatioModel);
        }

        return modelGraphs;
    };

    // ------------------------------------------------------
    // Draw all four panels
    // ------------------------------------------------------

    // (a) K*0 / K -- top left
    std::vector<TGraphErrors *> modelGraphsKstarKaon =
        drawRatioPanel(
            1,
            "hFrameFourRatioKstarKaon",
            "K*^{0} / K",
            gRatioKstarKaon,
            gKstarKaRatio_UrQMDOffSmooth,
            gKstarKaRatio_UrQMDOnSmooth,
            kKstar,
            kKaon,
            2.0,
            0.235, 0.44,
            false,
            false,
            true);

    // (b) Phi / K -- top right
    std::vector<TGraphErrors *> modelGraphsPhiKaon =
        drawRatioPanel(
            2,
            "hFrameFourRatioPhiKaon",
            "#phi / K",
            gRatioPhiK,
            gPhiKaRatio_UrQMDOffSmooth,
            gPhiKaRatio_UrQMDOnSmooth,
            kPhi,
            kKaon,
            2.0,
            0.07, 0.17,
            true,
            false,
            false);

    // (c) K*0 / pi -- bottom left
    std::vector<TGraphErrors *> modelGraphsKstarPion =
        drawRatioPanel(
            3,
            "hFrameFourRatioKstarPion",
            "K*^{0} / #pi",
            gRatioKstarPion,
            gKstarPiRatio_UrQMDOffSmooth,
            gKstarPiRatio_UrQMDOnSmooth,
            kKstar,
            kPion,
            2.0,
            0.031, 0.062,
            false,
            false,
            false);

    // (d) Phi / pi -- bottom right
    std::vector<TGraphErrors *> modelGraphsPhiPion =
        drawRatioPanel(
            4,
            "hFrameFourRatioPhiPion",
            "#phi / #pi",
            gRatioPhiPion,
            gPhiPiRatio_UrQMDOffSmooth,
            gPhiPiRatio_UrQMDOnSmooth,
            kPhi,
            kPion,
            2.0,
            0.007, 0.023,
            true,
            true,
            false);

    // --------------------------------------------------
    // Panel labels
    // --------------------------------------------------

    cParticleRatios->cd(1);

    TLatex panelLatex;
    panelLatex.SetNDC();
    panelLatex.SetTextFont(42);
    panelLatex.SetTextSize(0.045);
    panelLatex.DrawLatex(0.92, 0.85, "(a)");

    panelLatex.SetTextSize(0.05);
    panelLatex.DrawLatex(0.78, 0.845, "K*^{0} / K");
    // ALICE information

    panelLatex.SetTextSize(0.065);
    panelLatex.DrawLatex(0.22, 0.82, "#bf{ALICE}");

    TLegend *legendRun3 = new TLegend(0.2, 0.74, 0.54, 0.8);
    SetLegend(legendRun3);
    legendRun3->AddEntry(gRatioKstarKaon[0], "pp, #sqrt{s} = 13.6 TeV, |y| < 0.5", "p");
    legendRun3->Draw();

    // panelLatex.SetTextSize(0.055);
    // panelLatex.DrawLatex(0.22, 0.74, "pp, #sqrt{s} = 13.6 TeV");
    // panelLatex.DrawLatex(0.22, 0.67, "|y| < 0.5");

    // ------------------------------------------------------
    // Common model legend
    // ------------------------------------------------------

    cParticleRatios->cd(2);

    TLegend *legendFourRatios = new TLegend(0.05, 0.53, 0.54, 0.89);

    SetLegend(legendFourRatios);
    legendFourRatios->SetNColumns(1);

    legendFourRatios->AddEntry(gKstarKaRatio_UrQMDOffSmooth, "EPOS4", "l");
    legendFourRatios->AddEntry(gKstarKaRatio_UrQMDOnSmooth, "EPOS4 + UrQMD", "l");

    for (int i = 0; i < (int)modelGraphsKstarKaon.size(); i++)
    {
        auto model = modelsToPlot[i];
        legendFourRatios->AddEntry(modelGraphsKstarKaon[i], modelLabel[model], "l");
    }

    legendFourRatios->Draw();

    panelLatex.SetTextSize(0.045);
    panelLatex.DrawLatex(0.76, 0.85, "(b)");

    panelLatex.SetTextSize(0.05);
    panelLatex.DrawLatex(0.65, 0.85, "#phi / K");

    cParticleRatios->cd(3);
    panelLatex.SetTextSize(0.045);
    panelLatex.DrawLatex(0.92, 0.92, "(c)");

    panelLatex.SetTextSize(0.05);
    panelLatex.DrawLatex(0.78, 0.92, "K*^{0} / #pi");

    cParticleRatios->cd(4);
    panelLatex.SetTextSize(0.045);
    panelLatex.DrawLatex(0.76, 0.92, "(d)");

    panelLatex.SetTextSize(0.05);
    panelLatex.DrawLatex(0.65, 0.92, "#phi / #pi");

    // ------------------------------------------------------
    // Save
    // ------------------------------------------------------

    cParticleRatios->cd();

    if (isSavePlots)
    {
        cParticleRatios->SaveAs(
            (plotOutputDirectory + "/ParticleRatios_FourPanel_Run3.pdf").c_str());
    }

    //======================================================
    // COMBINED K*0 AND PHI OBSERVABLE PLOTS
    //======================================================
    // These panels reuse the measured, EPOS, and Pythia graphs already loaded
    // and drawn above. No additional physics scaling is applied here.

    gStyle->SetLineScalePS(1.7);

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

    auto setFourPanelMargins = [](int panel)
    {
        TVirtualPad *pad = gPad;
        const bool isLeft = (panel == 1 || panel == 3);
        const bool isTop = (panel == 1 || panel == 2);
        pad->SetLeftMargin(isLeft ? 0.16 : 0.005);
        pad->SetRightMargin(isLeft ? 0.005 : 0.16);
        pad->SetTopMargin(isTop ? 0.08 : 0.0);
        pad->SetBottomMargin(isTop ? 0.0 : 0.15);
        pad->SetTicks(1, 1);
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
                                   const char *panelLabel,
                                   const char *particleLabel,
                                   const PanelAnnotations &annotations)
    {
        pad->cd();

        TH1D *frame = new TH1D(frameName, "", 100, 2.0, 28.0);
        SetHistStyle(frame);
        frame->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
        frame->GetYaxis()->SetTitle(isMeanPt
                                        ? "#LT #it{p}_{T} #GT (GeV/#it{c})"
                                        : "dN/dy");
        frame->GetXaxis()->SetLimits(2.0, 28.0);
        frame->GetYaxis()->SetRangeUser(yMin, yMax);
        frame->GetYaxis()->SetMaxDigits(3);
        frame->GetYaxis()->SetTitleOffset(rightAxis ? 1.25 : 1.35);
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
        SetGraphStyle(run3Stat, kBlue + 1, 1.0, 20, 0.8, 1);
        run3Stat->Draw("PE SAME");
        SetGraphStyle(run3Sys, kBlue + 1, 1.0, 20, 0.8, 1);
        run3Sys->SetFillStyle(0);
        run3Sys->Draw("5 SAME");
        SetGraphStyle(run3Uncorr, kBlue + 1, 1.0, 20, 0.8, 1);
        run3Uncorr->SetFillStyle(1001);
        run3Uncorr->SetFillColorAlpha(kBlue + 1, 0.30);
        run3Uncorr->Draw("E2 SAME");

        // Run 2 data.
        SetGraphStyle(run2Stat, kRed + 1, 1.0, 24, 0.8, 1);
        run2Stat->Draw("PE SAME");
        SetGraphStyle(run2Sys, kRed + 1, 1.0, 24, 0.8, 1);
        run2Sys->SetFillStyle(0);
        run2Sys->Draw("5 SAME");
        SetGraphStyle(run2Uncorr, kRed + 1, 1.0, 24, 0.8, 1);
        run2Uncorr->SetFillStyle(1001);
        run2Uncorr->SetFillColorAlpha(kRed + 1, 0.30);
        run2Uncorr->Draw("E2 SAME");

        // EPOS with and without UrQMD.
        TGraphErrors *eposOffCopy = (TGraphErrors *)eposOff->Clone(Form("%s_EPOSOff", frameName));
        TGraphErrors *eposOnCopy = (TGraphErrors *)eposOn->Clone(Form("%s_EPOSOn", frameName));
        SetGraphStyle(eposOffCopy, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithoutUrQMD);
        eposOffCopy->Draw("L SAME");
        SetGraphStyle(eposOnCopy, colorEPOSWithUrQMD, 1.0, 1, 1.0, 3, lineStyleEPOSWithUrQMD);
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
            SetGraphStyle(modelGraph, modelStyle[model].color, 1.0, 1, 0.0, 3, modelStyle[model].style);
            modelGraph->Draw("LX SAME");
            drawnModelGraphs.push_back(modelGraph);
        }

        TLatex panelText;
        panelText.SetNDC();
        panelText.SetTextFont(42);
        panelText.SetTextSize(0.045);
        panelText.DrawLatex(annotations.panelLabelX, annotations.panelLabelY, panelLabel);
        panelText.SetTextSize(0.055);
        panelText.DrawLatex(annotations.particleLabelX, annotations.particleLabelY, particleLabel);

        if (annotations.legendKind == 1)
        {
            TLegend *dataLegend = new TLegend(annotations.legendX1, annotations.legendY1,
                                              annotations.legendX2, annotations.legendY2);
            SetLegend(dataLegend);
            dataLegend->SetTextAlign(12);
            dataLegend->AddEntry(run3Stat, "pp, #sqrt{s} = 13.6 TeV", "p");
            dataLegend->AddEntry(run2Stat, "pp, #sqrt{s} = 13 TeV", "p");
            dataLegend->Draw();
        }
        else if (annotations.legendKind == 2)
        {
            TLegend *modelLegend = new TLegend(annotations.legendX1, annotations.legendY1,
                                               annotations.legendX2, annotations.legendY2);
            SetLegend(modelLegend);
            modelLegend->AddEntry(eposOffCopy, "EPOS4", "l");
            modelLegend->AddEntry(eposOnCopy, "EPOS4 + UrQMD", "l");
            for (size_t iModel = 0; iModel < modelsToPlot.size(); ++iModel)
            {
                const auto model = modelsToPlot[iModel];
                modelLegend->AddEntry(drawnModelGraphs[iModel], modelLabel[model], "l");
            }
            modelLegend->Draw();
        }
    };

    // K*0 and phi side-by-side, matching the xT macro's two-pad layout.
    TCanvas *cMeanPtKstarPhi = new TCanvas("cMeanPtKstarPhi", "Mean pT: K*0 and phi", 850, 450);
    SetPad(cMeanPtKstarPhi);
    cMeanPtKstarPhi->cd();
    TPad *padMeanPtKstar = makeHorizontalPad("padMeanPtKstar", "K*0", 0.0, 0.5, true);
    TPad *padMeanPtPhi = makeHorizontalPad("padMeanPtPhi", "phi", 0.5, 1.0, false);
    const PanelAnnotations meanPtKstarAnnotations{
        0.27, // left x for particle label
        0.87, // left y for particle label
        0.2,  // left x for panel label
        0.88, // left y for panel label
        1,    // legend kind: 1 = data, 2 = models
        0.2,  // x1 for legend
        0.7,  // y1 for legend
        0.5,  // x2 for legend
        0.85  // y2 for legend
    };
    const PanelAnnotations meanPtPhiAnnotations{0.11, 0.88, 0.04, 0.88, 2, 0.4, 0.18, 0.92, 0.48};
    drawObservablePanel(padMeanPtKstar, "hFrameMeanPtKstarPhi_Kstar",
                        gMPtKstar, gMPtKstarRun2,
                        gEPOS_KstarMeanpT_UrQMDOff, gEPOS_KstarMeanpT_UrQMDOn,
                        kKstar, true, 0.58, 1.76, false, true, "(a)", "K*^{0}", meanPtKstarAnnotations);
    drawObservablePanel(padMeanPtPhi, "hFrameMeanPtKstarPhi_Phi",
                        gMPtPhiRun3, gMPtPhiRun2,
                        gEPOS_PhiMeanpT_UrQMDOff, gEPOS_PhiMeanpT_UrQMDOn,
                        kPhi, true, 0.58, 1.76, true, true, "(b)", "#phi", meanPtPhiAnnotations);

    // K*0 and phi dN/dy side-by-side.
    TCanvas *cYieldKstarPhi = new TCanvas("cYieldKstarPhi", "dN/dy: K*0 and phi", 850, 450);
    SetPad(cYieldKstarPhi);
    cYieldKstarPhi->cd();

    TPad *padYieldKstar = makeHorizontalPad("padYieldKstar", "K*0", 0.0, 0.5, true);
    TPad *padYieldPhi = makeHorizontalPad("padYieldPhi", "phi", 0.5, 1.0, false);
    const PanelAnnotations yieldKstarAnnotations{
        0.27, // left x for particle label
        0.87, // left y for particle label
        0.2,  // left x for panel label
        0.88, // left y for panel label
        1,    // legend kind: 1 = data, 2 = models
        0.2,  // x1 for legend
        0.7,  // y1 for legend
        0.5,  // x2 for legend
        0.85  // y2 for legend
    };
    const PanelAnnotations yieldPhiAnnotations{0.11, 0.88, 0.04, 0.88, 2, 0.04, 0.55, 0.5, 0.85};
    drawObservablePanel(padYieldKstar, "hFrameYieldKstarPhi_Kstar",
                        gMYieldKstar, gMYieldKstarRun2,
                        gEPOS_KstardNdy_UrQMDOff, gEPOS_KstardNdy_UrQMDOn,
                        kKstar, false, 0.0, 0.71, false, true, "(a)", "K*^{0}", yieldKstarAnnotations);
    drawObservablePanel(padYieldPhi, "hFrameYieldKstarPhi_Phi",
                        gMYieldPhiRun3, gMYieldPhiRun2,
                        gEPOS_PhidNdy_UrQMDOff, gEPOS_PhidNdy_UrQMDOn,
                        kPhi, false, 0.0, 0.31, true, true, "(b)", "#phi", yieldPhiAnnotations);

    // 2x2 summary: mean pT across the top row and dN/dy across the bottom.
    TCanvas *cMeanPtYieldKstarPhi = new TCanvas(
        "cMeanPtYieldKstarPhi", "Mean pT and dN/dy: K*0 and phi", 900, 900);
    SetPad(cMeanPtYieldKstarPhi);
    cMeanPtYieldKstarPhi->Divide(2, 2, 0.0, 0.0);
    for (int iPanel = 1; iPanel <= 4; ++iPanel)
    {
        cMeanPtYieldKstarPhi->cd(iPanel);
        setFourPanelMargins(iPanel);
    }

    const PanelAnnotations summaryKstarPtAnnotations{0.27, 0.85, 0.2, 0.86, 1, 0.2, 0.68, 0.5, 0.83};
    const PanelAnnotations summaryPhiPtAnnotations{0.11, 0.86, 0.04, 0.86, 2, 0.4, 0.04, 0.92, 0.35};
    const PanelAnnotations summaryKstarYieldAnnotations{0.27, 0.9, 0.2, 0.91, 0, 0.0, 0.0, 0.0, 0.0};
    const PanelAnnotations summaryPhiYieldAnnotations{0.11, 0.91, 0.04, 0.91, 0, 0.0, 0.0, 0.0, 0.0};
    cMeanPtYieldKstarPhi->cd(1);
    drawObservablePanel(gPad, "hFrameMeanPtYieldKstarPhi_KstarPt",
                        gMPtKstar, gMPtKstarRun2,
                        gEPOS_KstarMeanpT_UrQMDOff, gEPOS_KstarMeanpT_UrQMDOn,
                        kKstar, true, 0.58, 1.76, false, false, "(a)", "K*^{0}", summaryKstarPtAnnotations);

    TLatex panelLatex2;
    panelLatex2.SetNDC();
    panelLatex2.SetTextFont(42);
    panelLatex2.SetTextSize(0.065);
    panelLatex2.DrawLatex(0.78, 0.82, "#bf{ALICE}");

    cMeanPtYieldKstarPhi->cd(2);
    drawObservablePanel(gPad, "hFrameMeanPtYieldKstarPhi_PhiPt",
                        gMPtPhiRun3, gMPtPhiRun2,
                        gEPOS_PhiMeanpT_UrQMDOff, gEPOS_PhiMeanpT_UrQMDOn,
                        kPhi, true, 0.58, 1.76, true, false, "(b)", "#phi", summaryPhiPtAnnotations);
    cMeanPtYieldKstarPhi->cd(3);
    drawObservablePanel(gPad, "hFrameMeanPtYieldKstarPhi_KstarYield",
                        gMYieldKstar, gMYieldKstarRun2,
                        gEPOS_KstardNdy_UrQMDOff, gEPOS_KstardNdy_UrQMDOn,
                        kKstar, false, 0.0, 0.71, false, true, "(c)", "K*^{0}", summaryKstarYieldAnnotations);
    cMeanPtYieldKstarPhi->cd(4);
    drawObservablePanel(gPad, "hFrameMeanPtYieldKstarPhi_PhiYield",
                        gMYieldPhiRun3, gMYieldPhiRun2,
                        gEPOS_PhidNdy_UrQMDOff, gEPOS_PhidNdy_UrQMDOn,
                        kPhi, false, 0.0, 0.26, true, true, "(d)", "#phi", summaryPhiYieldAnnotations);

    // if (isSavePlots)
    {
        // cMeanPtKstarPhi->SaveAs((plotOutputDirectory + "/MeanPt_KstarPhi_Run3.pdf").c_str());
        // cYieldKstarPhi->SaveAs((plotOutputDirectory + "/dNdy_KstarPhi_Run3.pdf").c_str());
        cMeanPtYieldKstarPhi->SaveAs((plotOutputDirectory + "/MeanPtYield_KstarPhi_Run3.pdf").c_str());
    }

    //======================================================
    // 1. TWO-PANEL RATIO PLOT (2,1): K*0 / Phi and K*0 / K0s
    //======================================================
    TCanvas *cRatioKstarPhiK0s = new TCanvas("cRatioKstarPhiK0s", "Ratios: K*0/Phi and K*0/K0s", 850, 450);
    SetPad(cRatioKstarPhiK0s);
    cRatioKstarPhiK0s->cd();

    TPad *padRatioKstarPhi = makeHorizontalPad("padRatioKstarPhi", "K*0/Phi", 0.0, 0.51, true);
    TPad *padRatioKstarK0s = makeHorizontalPad("padRatioKstarK0s", "K*0/K0s", 0.51, 1.0, false);
    padRatioKstarK0s->SetRightMargin(0.15);

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
        if (!rightAxis)
            frame->GetYaxis()->SetTitle("dN / dy ratio");
        frame->GetXaxis()->SetLimits(1.8, 23.5);
        frame->GetYaxis()->SetRangeUser(yMin, yMax);
        frame->GetYaxis()->SetMaxDigits(3);
        frame->GetYaxis()->SetNdivisions(505);
        frame->GetYaxis()->SetTitleOffset(rightAxis ? 1.25 : 1.35);

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

        SetGraphStyle(eposOffCopy, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 1, lineStyleEPOSWithoutUrQMD);
        eposOffCopy->Draw("LX SAME");
        SetGraphStyle(eposOnCopy, colorEPOSWithUrQMD, 1.0, 1, 1.0, 1, lineStyleEPOSWithUrQMD);
        eposOnCopy->Draw("LX SAME");

        thermalCopy->SetLineColor(kMagenta + 2);
        thermalCopy->SetLineStyle(5);
        thermalCopy->SetLineWidth(2);
        thermalCopy->Draw("LX SAME");

        // Pythia models
        std::vector<TGraphErrors *> pythiaGraphs;
        for (auto model : modelsToPlot)
        {
            TGraphErrors *gMod = MakeRatio(gMYield[model][numParticle], gMYield[model][denParticle], true);
            ScaleGraph(gMod, modelScale);
            SetGraphStyle(gMod, modelStyle[model].color, 1.0, 1, 1.0, 1, modelStyle[model].style);
            gMod->Draw("L SAME");
            pythiaGraphs.push_back(gMod);
        }

        TLatex txt;
        txt.SetNDC();
        txt.SetTextFont(42);
        txt.SetTextSize(0.045);
        txt.DrawLatex(rightAxis ? 0.05 : 0.22, 0.87, panelLabel);
        txt.SetTextSize(0.052);
        txt.DrawLatex(rightAxis ? 0.13 : 0.3, 0.865, ratioLabel);

        if (showDataLegend)
        {
            txt.SetTextSize(0.052);
            txt.DrawLatex(0.65, 0.87, "#bf{ALICE}");
            TLegend *legData = new TLegend(0.57, 0.75, 0.90, 0.85);
            SetLegend(legData);
            legData->SetTextSize(0.038);
            legData->AddEntry((TObject *)0, "pp, #sqrt{s} = 13.6 TeV", "");
            legData->AddEntry((TObject *)0, "|y| < 0.5", "");
            legData->Draw();
        }

        if (showModelLegend)
        {
            // Create a legend for the model graphs
            TLegend *legMod = new TLegend(0.37, 0.59, 0.90, 0.92);
            SetLegend(legMod);
            legMod->SetTextSize(0.035);
            legMod->AddEntry(eposOffCopy, "EPOS4", "l");
            legMod->AddEntry(eposOnCopy, "EPOS4 + UrQMD", "l");
            legMod->AddEntry(thermalCopy, "Thermal-FIST", "l");
            for (size_t i = 0; i < modelsToPlot.size(); ++i)
            {
                legMod->AddEntry(pythiaGraphs[i], modelLabel[modelsToPlot[i]], "l");
            }
            legMod->Draw();
        }
    };

    // Panel (a): K*0 / Phi (already averaged in both numerator & denominator -> scale = 1.0)
    drawTwoPanelRatio(padRatioKstarPhi, "hFrame2Panel_KstarPhi",
                      gRatioKstarPhi,
                      gKstarPhiRatio_EPOSNew_UrQMDOff, gKstarPhiRatio_EPOSNew_UrQMDOn, gRatioKstarPhi_Thermal,
                      kKstar, kPhi, 1.0,
                      2.1, 5.15, false, "(a)", "K*^{0} / #phi", true, false);

    // Panel (b): K*0 / K0s (K*0 averaged, K0s not -> scale = 1.0)
    drawTwoPanelRatio(padRatioKstarK0s, "hFrame2Panel_KstarK0s",
                      gRatioKstarKshort,
                      gKstarK0sRatio_EPOSNew_UrQMDOff, gKstarK0sRatio_EPOSNew_UrQMDOn, gRatioKstarK0s_Thermal,
                      kKstar, kKshort, 1.0,
                      0.245, 0.54, true, "(b)", "K*^{0} / K_{S}^{0}", false, true);

    //======================================================
    // 2. THREE-PANEL RATIO PLOT (3,1) WITH SEPARATE Y-AXES:
    //    Phi / Pi, Phi / K0s, Phi / Xi
    //======================================================
    TCanvas *cRatioPhiThreePanel = new TCanvas("cRatioPhiThreePanel", "Ratios: Phi/Pi, Phi/K0s, Phi/Xi", 1350, 460);
    SetPad(cRatioPhiThreePanel);
    cRatioPhiThreePanel->cd();

    // Three equal-width pads with comfortable left margins so every y-axis label is visible
    TPad *padPhiPi = new TPad("padPhiPi", "Phi/Pi", 0.00, 0.0, 0.335, 1.0);
    TPad *padPhiK0s = new TPad("padPhiK0s", "Phi/K0s", 0.335, 0.0, 0.668, 1.0);
    TPad *padPhiXi = new TPad("padPhiXi", "Phi/Xi", 0.663, 0.0, 1.00, 1.0);

    for (TPad *p : {padPhiPi, padPhiK0s, padPhiXi})
    {
        p->SetLeftMargin(0.16);
        p->SetRightMargin(0.03);
        p->SetTopMargin(0.07);
        p->SetBottomMargin(0.15);
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
        gStyle->SetLineScalePS(1.5);

        TH1D *frame = new TH1D(frameName, "", 100, 1.8, 23.5);
        SetHistStyle(frame);
        frame->GetXaxis()->SetTitle("#LT dN_{ch}/d#eta #GT_{|#eta|<0.5}");
        frame->GetYaxis()->SetTitle(yTitle);
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

        SetGraphStyle(eposOffCopy, colorEPOSWithoutUrQMD, 1.0, 1, 1.0, 1, lineStyleEPOSWithoutUrQMD);
        eposOffCopy->Draw("LX SAME");
        SetGraphStyle(eposOnCopy, colorEPOSWithUrQMD, 1.0, 1, 1.0, 1, lineStyleEPOSWithUrQMD);
        eposOnCopy->Draw("LX SAME");

        thermalCopy->SetLineColor(kMagenta + 2);
        thermalCopy->SetLineStyle(5);
        thermalCopy->SetLineWidth(2);
        thermalCopy->Draw("LX SAME");

        std::vector<TGraphErrors *> pythiaGraphs;
        for (auto model : modelsToPlot)
        {
            TGraphErrors *gMod = MakeRatio(gMYield[model][numParticle], gMYield[model][denParticle], true);
            ScaleGraph(gMod, modelScale);
            SetGraphStyle(gMod, modelStyle[model].color, 1.0, 1, 1.0, 1, modelStyle[model].style);
            gMod->Draw("L SAME");
            pythiaGraphs.push_back(gMod);
        }

        TLatex txt;
        txt.SetNDC();
        txt.SetTextFont(42);
        txt.SetTextSize(0.048);
        txt.DrawLatex(0.2, 0.85, panelLabel);
        txt.SetTextSize(0.055);
        txt.DrawLatex(0.27, 0.85, ratioLabel);

        if (showInfo)
        {
            txt.SetTextSize(0.065);
            txt.DrawLatex(0.74, 0.84, "#bf{ALICE}");
        }

        if (showInfo2)
        {

            TLegend *legData = new TLegend(0.5, 0.77, 0.88, 0.90);
            SetLegend(legData);
            legData->SetTextSize(0.048);
            legData->AddEntry((TObject *)0, "pp, #sqrt{s} = 13.6 TeV", "");
            legData->AddEntry((TObject *)0, "|y| < 0.5", "");
            legData->Draw();
        }

        if (showModelLegend)
        {
            TLegend *legMod = new TLegend(0.55, 0.6, 0.90, 0.9);
            SetLegend(legMod);
            legMod->SetTextSize(0.038);
            legMod->AddEntry(eposOffCopy, "EPOS4", "l");
            legMod->AddEntry(eposOnCopy, "EPOS4 + UrQMD", "l");
            legMod->AddEntry(thermalCopy, "Thermal-FIST", "l");
            for (size_t i = 0; i < modelsToPlot.size(); ++i)
            {
                legMod->AddEntry(pythiaGraphs[i], modelLabel[modelsToPlot[i]], "l");
            }
            legMod->Draw();
        }
    };

    // Panel (a): Phi / Pi (pion is summed pi+ + pi- in Pythia -> scale = 2.0)
    drawSpacedRatioPanel(padPhiPi, "hFrame3Panel_PhiPi",
                         gRatioPhiPion,
                         gPhiPiRatio_EPOSNew_UrQMDOff, gPhiPiRatio_EPOSNew_UrQMDOn, gRatioPhiPi_Thermal,
                         kPhi, kPion, 2.0,
                         0.007, 0.028, "#phi / #pi", "(a)", "#phi / #pi", true, false, false, 1.25);

    // Panel (b): Phi / K0s
    drawSpacedRatioPanel(padPhiK0s, "hFrame3Panel_PhiK0s",
                         gRatioPhiK0s,
                         gPhiK0sRatio_EPOSNew_UrQMDOff, gPhiK0sRatio_EPOSNew_UrQMDOn, gRatioPhiK0s_Thermal,
                         kPhi, kKshort, 1.0,
                         0.07, 0.16, "#phi / K_{S}^{0}", "(b)", "#phi / K_{S}^{0}", false, true, false, 1.60);

    // Panel (c): Phi / Xi
    // Note: gPhiXiRatio_EPOSNew_UrQMDOff/On were already scaled by 0.5 earlier in the macro.
    // Pass nullptr for data while fXiRun3 / gMYieldXiRun3 are commented out;
    // once uncommented, pass MakeRatioUncorr(gMYieldPhiRun3, gMYieldXiRun3, 1.0, 3) instead.
    drawSpacedRatioPanel(padPhiXi, "hFrame3Panel_PhiXi",
                         gRatioPhiXi,
                         gPhiXiRatio_EPOSNew_UrQMDOff, gPhiXiRatio_EPOSNew_UrQMDOn, gRatioPhiXi_Thermal,
                         kPhi, kXi, 1.0,
                         1.1, 4.3, "#phi / #Xi", "(c)", "#phi / #Xi", false, false, true, 0.95);

    // if (isSavePlots)
    {
        cRatioKstarPhiK0s->SaveAs((plotOutputDirectory + "/Ratio_KstarPhi_KstarK0s_TwoPanel_Run3.pdf").c_str());
        cRatioPhiThreePanel->SaveAs((plotOutputDirectory + "/Ratio_PhiPi_PhiK0s_PhiXi_ThreePanel_Run3.pdf").c_str());
    }
}

// GRAPH HELPERS: ratios, input access, scaling, and multiplicity normalization
TGraphErrors *MakeRatio(const TGraphErrors *numerator, const TGraphErrors *denominator, bool isModel = false, float CorrFactorDen = 1)
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

        // if (std::round(xNumerator * 10.0) != std::round(xDenominator * 10.0))
        // {
        //     cout << "Warning!!!!! x values do not match for ratio at point "
        //          << i << ": " << xNumerator << " vs " << xDenominator
        //          << ", for graphs " << numerator->GetName()
        //          << " and " << denominator->GetName() << endl;
        // }

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

TGraphErrors **MakeRatioUncorr(TGraphErrors **numerator, TGraphErrors **denominator, float CorrFactorDen = 1, int nGraphs = 3)
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
    TGraphErrors *graph = (TGraphErrors *)f->Get(name.c_str());

    if (!graph || graph == nullptr)
    {
        cout << "Error: graph " << name
             << " not found in file " << f->GetName() << endl;
        return nullptr;
    }

    SetGraphStyle(graph);
    graph->SetTitle(0);
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

TGraphErrors **DivideByMult(TGraphErrors **gr, double WhichMultPoint, double tolerance = 0.5, int noOfGraphs = 3)
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

TGraphErrors *DivideByMultModel(TGraphErrors *gr, double WhichMultPoint, double tolerance = 0.5)
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
