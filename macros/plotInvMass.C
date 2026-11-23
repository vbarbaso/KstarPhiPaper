
#include <iostream>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "THnSparse.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"

#include "../include/StyleProposal.h"

void plotInvMass()
{
    TGaxis::SetMaxDigits(3);

    // =========================================================================
    // 1. Configuration
    // =========================================================================

    TString savepath = "";

    const std::string dataFilePath = "../inputs/Kstar0_136/749276.root";
    const std::string templFilePath = "../inputs/Kstar0_136/SignalMinusTrue_Norm1.root";

    double multLow = 20.0, multHigh = 30.0;
    double ptLow = 1.0, ptHigh = 1.2;

    double normLow = 1.20, normHigh = 1.30;
    double fitLow = 0.75, fitHigh = 1.06;

    int rebinVal = 1;

    double massPDG = 0.89555;
    double widthPDG = 0.0473;

    const std::string phiFilePath = "../inputs/Phi_136/inelgr0/PAR.root";
    double phiPtLow = 1.90, phiPtHigh = 2.00;

    // =========================================================================
    // 2. Input files and histograms
    // =========================================================================

    TFile *fInput = TFile::Open(dataFilePath.c_str(), "READ");
    TFile *fTempl = TFile::Open(templFilePath.c_str(), "READ");
    if (!fInput || fInput->IsZombie() || !fTempl || fTempl->IsZombie())
    {
        std::cerr << "Error opening data or template file!" << std::endl;
        return;
    }

    const std::string folderName = "kstarqa/hInvMass";

    THnSparseF *hSparseNum = (THnSparseF *)fInput->Get(Form("%s/h3KstarInvMassUnlikeSign", folderName.c_str()));
    THnSparseF *hSparseRot = (THnSparseF *)fInput->Get(Form("%s/h3KstarInvMassRotated", folderName.c_str()));

    TH1D *hReflRaw = (TH1D *)fTempl->Get(Form("%d-%d/hSigminusTrue_pt_%.1f_%.1f", (int)multLow, (int)multHigh, ptLow, ptHigh));
    if (!hSparseNum || !hSparseRot || !hReflRaw)
    {
        std::cerr << "Histograms or reflection template not found!" << std::endl;
        return;
    }

    TFile *fPhi = TFile::Open(phiFilePath.c_str(), "READ");
    if (!fPhi || fPhi->IsZombie())
    {
        std::cerr << "Error opening Phi results file!" << std::endl;
        return;
    }

    TString phiDir = Form("Data/%.2f_%.2f/%.2f_%.2f", multLow, multHigh, phiPtLow, phiPtHigh);

    TH1 *hPhiSignal = (TH1 *)fPhi->Get(phiDir + "/hSignal");
    TH1 *hPhiBgNorm = (TH1 *)fPhi->Get(phiDir + "/hBgNorm");
    TH1 *hPhiPeak = (TH1 *)fPhi->Get(phiDir + "/hPeak");

    if (!hPhiSignal || !hPhiBgNorm || !hPhiPeak)
    {
        std::cerr << "Phi histograms not found in " << phiDir << std::endl;
        return;
    }

    // =========================================================================
    // 3. Projections and rotational background subtraction
    // =========================================================================

    hSparseNum->GetAxis(0)->SetRangeUser(multLow, multHigh);
    hSparseNum->GetAxis(1)->SetRangeUser(ptLow, ptHigh);

    hSparseRot->GetAxis(0)->SetRangeUser(multLow, multHigh);
    hSparseRot->GetAxis(1)->SetRangeUser(ptLow, ptHigh);

    TH1D *hSigBkg = hSparseNum->Projection(2, "E");
    TH1D *hRotBkg = hSparseRot->Projection(2, "E");

    hSigBkg->SetName("hSigBkg_Kstar");
    hRotBkg->SetName("hRotBkg_Kstar");

    hSigBkg->Rebin(rebinVal);
    hRotBkg->Rebin(rebinVal);

    // Normalize rotational background

    int bNormLo = hSigBkg->FindBin(normLow + 1e-5);
    int bNormHi = hSigBkg->FindBin(normHigh - 1e-5);

    double normFactor = hSigBkg->Integral(bNormLo, bNormHi) / hRotBkg->Integral(bNormLo, bNormHi);

    TH1D *hBkgScaled = (TH1D *)hRotBkg->Clone("hBkgScaled_Kstar");
    hBkgScaled->Scale(normFactor);

    // Background-subtracted signal

    TH1D *hSubSignal = (TH1D *)hSigBkg->Clone("hSubSignal_Kstar");

    hSubSignal->Add(hBkgScaled, -1.0);

    // =========================================================================
    // 4. Reflection template preparation
    // =========================================================================

    TH1D *hReflection = (TH1D *)hReflRaw->Clone("hReflection_Kstar");
    hReflection->Rebin(rebinVal);

    for (int ib = 1; ib <= hReflection->GetNbinsX(); ib++)
    {
        if (hReflection->GetBinContent(ib) < 0)
        {
            hReflection->SetBinContent(ib, 0.0);
            hReflection->SetBinError(ib, 0.0);
        }
    }

    const double binWidthFit = hSigBkg->GetBinWidth(1);
    const double reflBinWidth = hReflection->GetBinWidth(1);

    int bReflLo = hReflection->FindBin(fitLow + 1e-6);
    int bReflHi = hReflection->FindBin(fitHigh - 1e-6);

    const double reflNorm = std::max(hReflection->Integral(bReflLo, bReflHi) * reflBinWidth, 1e-12);
    const double totalInFit = std::max(1.0, hSubSignal->Integral(hSubSignal->FindBin(fitLow), hSubSignal->FindBin(fitHigh)));

    // =========================================================================
    // 5. Residual background parametrization
    // =========================================================================

    auto evalBkgDensity = [fitLow, fitHigh](double x, const double *par) -> double
    {
        double t = 2.0 * (x - fitLow) / (fitHigh - fitLow) - 1.0;

        double T1 = t;
        double T2 = 2.0 * t * t - 1.0;
        double T3 = 4.0 * t * t * t - 3.0 * t;

        double chebRaw = 1.0 + par[0] * T1 + par[1] * T2 + par[2] * T3;
        double chebInt = ((fitHigh - fitLow) / 2.0) * (2.0 - (2.0 / 3.0) * par[1]);

        return chebRaw / std::max(chebInt, 1e-12);
    };

    // =========================================================================
    // 6. Total fit: Voigtian + reflection template + Chebyshev T3
    // =========================================================================

    TF1 *fTotal = new TF1("fTotal_Kstar", [hReflection, reflNorm, fitLow, fitHigh, binWidthFit, evalBkgDensity](double *xx, double *pp) -> double
                          {
            double x = xx[0];

            double Nsig = pp[0];
            double mass = pp[1];
            double width = pp[2];
            double sigma = pp[3];

            double Ncorr = pp[4];
            double Nbkg = pp[5];

            // Voigtian signal

            double voigtRaw =TMath::Voigt(x - mass, sigma, width);
            int nSteps = 40;
            double hStep = (fitHigh - fitLow) / nSteps;
            double vInt =TMath::Voigt(fitLow - mass, sigma, width) +TMath::Voigt(fitHigh - mass, sigma, width);

            for (int j = 1; j < nSteps; j++)
            {
                double xj = fitLow + j * hStep;

                vInt += TMath::Voigt(xj - mass, sigma, width) *((j % 2 == 0) ? 2.0 : 4.0);
            }

            vInt *= hStep / 3.0;

            double sigDensity =voigtRaw / std::max(vInt, 1e-12);

            // MC reflection template

            double corrDensity =(hReflection && hReflection->GetEntries() > 0)? hReflection->Interpolate(x) / reflNorm: 0.0;

            if (corrDensity < 0)
                corrDensity = 0.0;

            // Residual background

            double bkgDensity =evalBkgDensity(x, &pp[6]);

            return binWidthFit *(Nsig * sigDensity +Ncorr * corrDensity +Nbkg * bkgDensity); }, fitLow, fitHigh, 9);

    fTotal->SetParNames("Nsig", "mass", "width", "sigma", "Ncorr", "Nbkg", "c1", "c2", "c3");

    fTotal->SetParameters(
        totalInFit * 0.4,
        massPDG,
        widthPDG,
        0.005,
        totalInFit * 0.2,
        totalInFit * 0.3,
        -0.5, 0.1, -0.05);

    fTotal->FixParameter(2, widthPDG);
    fTotal->SetParLimits(1, massPDG - 0.01, massPDG + 0.01);
    fTotal->SetParLimits(3, 0.0005, 0.02);

    // =========================================================================
    // 7. Canvas and pad configuration
    // =========================================================================

    hSigBkg->GetYaxis()->SetNoExponent(kFALSE);
    hBkgScaled->GetYaxis()->SetNoExponent(kFALSE);
    hSubSignal->GetYaxis()->SetNoExponent(kFALSE);

    hPhiSignal->GetYaxis()->SetNoExponent(kFALSE);
    hPhiBgNorm->GetYaxis()->SetNoExponent(kFALSE);
    hPhiPeak->GetYaxis()->SetNoExponent(kFALSE);

    TCanvas *c1 = new TCanvas("c1_KstarPhi", "Kstar and Phi invariant mass", 900, 900);
    SetPad(c1);

    c1->Divide(2, 2, 0.0, 0.0);

    // Pad 1: Kstar unlike-sign pairs + background

    c1->cd(1);

    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.005);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.0);
    gPad->SetTicks(1, 1);

    // Pad 2: Phi signal + background normalization region

    c1->cd(2);

    gPad->SetLeftMargin(0.005);
    gPad->SetRightMargin(0.16);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.0);
    gPad->SetTicks(1, 1);

    // Pad 3: Kstar background-subtracted signal + fit

    c1->cd(3);

    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.005);
    gPad->SetTopMargin(0.00);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicks(1, 1);

    // Pad 4: Phi background-subtracted signal + fit

    c1->cd(4);

    gPad->SetLeftMargin(0.005);
    gPad->SetRightMargin(0.16);
    gPad->SetTopMargin(0.0);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicks(1, 1);

    // =========================================================================
    // 8. PAD 1: K*0 unlike-sign pairs + rotational background
    // =========================================================================

    c1->cd(1);

    SetHistStyle(hSigBkg, kBlack, 1.0, 20, 1.0, 1, 1);
    SetHistStyle(hBkgScaled, kRed + 1, 1.0, 24, 1.0, 1, 1);

    hSigBkg->SetTitle(Form("; ;Counts / (%.3f GeV/#it{c}^{2})", binWidthFit));

    hSigBkg->GetXaxis()->SetRangeUser(0.76, 1.048);
    hSigBkg->GetYaxis()->SetRangeUser(hSigBkg->GetMinimum() * 0.85, hSigBkg->GetMaximum() * 1.25);

    hSigBkg->GetYaxis()->SetTitleOffset(1.3);

    hSigBkg->Draw("E");
    hBkgScaled->Draw("E SAME");

    TLatex ltx;
    ltx.SetNDC();
    ltx.SetTextFont(42);
    ltx.SetTextSize(0.055);

    // ALICE labels
    ltx.DrawLatex(0.22, 0.83, "ALICE Preliminary");
    ltx.DrawLatex(0.22, 0.75, "pp #sqrt{#it{s}} = 13.6 TeV");
    ltx.DrawLatex(0.22, 0.67, "|#it{y}| < 0.5");

    // K*0 decay channel and multiplicity class
    ltx.SetTextSize(0.045);
    ltx.DrawLatex(0.22, 0.17, "K*^{0} (#bar{K}*^{0}) #rightarrow K^{#pm}#pi^{#mp}");
    ltx.DrawLatex(0.22, 0.11, Form("FT0M class: (%d - %d)%%", (int)multLow, (int)multHigh));
    ltx.DrawLatex(0.22, 0.05, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptLow, ptHigh));

    ltx.DrawLatex(0.92, 0.85, "(a)");

    // =========================================================================
    // 9. PAD 3: Background-subtracted signal + fit components
    // =========================================================================

    c1->cd(3);

    SetHistStyle(hSubSignal, kBlack, 1.0, 20, 1.0, 1, 1);

    hSubSignal->Fit(fTotal, "R Q S");
    hSubSignal->SetTitle(Form(";M_{K#pi} (GeV/#it{c}^{2});Counts / (%.3f GeV/#it{c}^{2})", binWidthFit));

    hSubSignal->GetXaxis()->SetRangeUser(0.76, 1.048);
    hSubSignal->GetYaxis()->SetRangeUser(0.0, hSubSignal->GetMaximum() * 1.5);
    hSubSignal->GetYaxis()->SetTitleOffset(1.3);

    hSubSignal->Draw("E");

    // -------------------------------------------------------------------------
    // Voigtian signal component
    // -------------------------------------------------------------------------

    TF1 *fSig = new TF1(
        "fSig_Kstar",
        [fitLow, fitHigh, binWidthFit](double *xx, double *pp) -> double
        {
            double x = xx[0];

            double voigtRaw =
                TMath::Voigt(x - pp[1], pp[3], pp[2]);

            int nSteps = 40;
            double hStep = (fitHigh - fitLow) / nSteps;

            double vInt =
                TMath::Voigt(fitLow - pp[1], pp[3], pp[2]) +
                TMath::Voigt(fitHigh - pp[1], pp[3], pp[2]);

            for (int j = 1; j < nSteps; j++)
            {
                vInt += TMath::Voigt(
                            fitLow + j * hStep - pp[1],
                            pp[3], pp[2]) *
                        ((j % 2 == 0) ? 2.0 : 4.0);
            }

            return binWidthFit * pp[0] * voigtRaw /
                   std::max(vInt * hStep / 3.0, 1e-12);
        },
        fitLow, fitHigh, 4);

    fSig->SetParameters(
        fTotal->GetParameter(0),
        fTotal->GetParameter(1),
        fTotal->GetParameter(2),
        fTotal->GetParameter(3));

    fSig->SetLineColor(kMagenta - 3);
    fSig->SetLineWidth(2);
    fSig->SetNpx(500);
    // fSig->Draw("SAME");

    // -------------------------------------------------------------------------
    // MC reflection template component
    // -------------------------------------------------------------------------

    TF1 *fCorr = new TF1(
        "fCorr_Kstar",
        [hReflection, reflNorm, binWidthFit](double *xx, double *pp) -> double
        {
            return binWidthFit * pp[0] *
                   (hReflection->Interpolate(xx[0]) / reflNorm);
        },
        fitLow, fitHigh, 1);

    fCorr->SetParameter(0, fTotal->GetParameter(4));

    SetFunctionStyle(fCorr, kGreen + 2, 1.0, 2, 1);
    fCorr->Draw("SAME");

    // -------------------------------------------------------------------------
    // Residual background: Chebyshev T3
    // -------------------------------------------------------------------------

    TF1 *fBkg = new TF1(
        "fBkg_Kstar",
        [fitLow, fitHigh, binWidthFit,
         evalBkgDensity](double *xx, double *pp) -> double
        {
            return binWidthFit * pp[0] *
                   evalBkgDensity(xx[0], &pp[1]);
        },
        fitLow, fitHigh, 4);

    fBkg->SetParameters(
        fTotal->GetParameter(5),
        fTotal->GetParameter(6),
        fTotal->GetParameter(7),
        fTotal->GetParameter(8));

    SetFunctionStyle(fBkg, kBlue, 1.0, 2, 1);
    fBkg->Draw("SAME");

    // -------------------------------------------------------------------------
    // Total background: MC reflection + Chebyshev T3
    // -------------------------------------------------------------------------

    TF1 *fTotalBkg = new TF1(
        "fTotalBkg_Kstar",
        [fCorr, fBkg](double *xx, double *pp) -> double
        {
            return fCorr->Eval(xx[0]) + fBkg->Eval(xx[0]);
        },
        fitLow, fitHigh, 0);

    SetFunctionStyle(fTotalBkg, kCyan + 2, 1.0, 2, 2);
    fTotalBkg->Draw("SAME");

    // -------------------------------------------------------------------------
    // Total signal + background fit
    // -------------------------------------------------------------------------

    SetFunctionStyle(fTotal, kRed, 1.0, 2, 1);
    fTotal->Draw("SAME");

    // -------------------------------------------------------------------------
    // Legend
    // -------------------------------------------------------------------------

    TLegend *leg2 = new TLegend(0.21, 0.8, 0.55, 0.96);
    // TLegend *leg2 = new TLegend(0.21, 0.75, 0.55, 0.96);
    SetLegend(leg2);
    // leg2->AddEntry(fSig, "Voigtian (Signal)", "l");
    leg2->AddEntry(fCorr, "MC Reflection Template", "l");
    leg2->AddEntry(fBkg, "Chebyshev T_{3}", "l");
    leg2->AddEntry(fTotalBkg, "Template + Chebyshev (Background)", "l");

    leg2->Draw();

    ltx.DrawLatex(0.92, 0.93, "(c)");

    // =========================================================================
    // 10. PAD 2: Phi unlike-sign pairs + background
    // =========================================================================

    c1->cd(2);

    hPhiSignal->GetYaxis()->SetMaxDigits(3);
    double phiBinWidth = hPhiSignal->GetBinWidth(1);

    SetHistStyle(hPhiSignal, kBlack, 1.0, 20, 1.0, 1, 1);
    SetHistStyle(hPhiBgNorm, kRed + 1, 1.0, 24, 1.0, 1, 1);

    hPhiSignal->SetTitle(Form("; ;Counts / (%.3f GeV/#it{c}^{2})", phiBinWidth));
    hPhiSignal->GetXaxis()->SetRangeUser(0.992, 1.058);
    hPhiSignal->GetYaxis()->SetRangeUser(0.001, hPhiSignal->GetMaximum() * 1.25);
    hPhiSignal->GetYaxis()->SetTitleOffset(1.3);

    hPhiSignal->Draw("E Y+");
    hPhiBgNorm->Draw("E SAME");

    // Legend

    TLegend *legPhi1 = new TLegend(0.38, 0.75, 0.73, 0.89);
    SetLegend(legPhi1);
    legPhi1->AddEntry(hPhiSignal, "Unlike-sign pairs", "p");
    legPhi1->AddEntry(hPhiBgNorm, "Rotational background", "p");

    legPhi1->Draw();

    ltx.DrawLatex(0.4, 0.17, "#phi #rightarrow K^{+}K^{-}");
    ltx.DrawLatex(0.4, 0.11, Form("FT0M class: (%d - %d)%%", (int)multLow, (int)multHigh));
    ltx.DrawLatex(0.4, 0.05, Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", phiPtLow, phiPtHigh));

    ltx.DrawLatex(0.04, 0.85, "(b)");

    // =========================================================================
    // PAD 4: Phi background-subtracted signal + fit components
    // =========================================================================

    c1->cd(4);

    hPhiPeak->GetYaxis()->SetMaxDigits(1);

    // Retrieve attached fit functions

    TF1 *fPhiFit = nullptr;
    TF1 *fPhiBkg = nullptr;

    TList *funcList = hPhiPeak->GetListOfFunctions();

    if (funcList)
    {
        for (int i = 0; i < funcList->GetSize(); i++)
        {
            TF1 *f = dynamic_cast<TF1 *>(funcList->At(i));
            if (!f)
                continue;

            TString funcName = f->GetName();
            funcName.ToLower();

            if (funcName.Contains("sigbgfit"))
                fPhiFit = f;

            else if (funcName.Contains("bgfit"))
                fPhiBkg = f;
        }
    }

    // -------------------------------------------------------------------------
    // Histogram styling
    // -------------------------------------------------------------------------

    SetHistStyle(hPhiPeak, kBlack, 1.0, 20, 1.0, 1, 1);

    hPhiPeak->SetTitle(Form(";M_{KK} (GeV/#it{c}^{2});Counts / (%.3f GeV/#it{c}^{2})", phiBinWidth));
    hPhiPeak->GetXaxis()->SetRangeUser(0.992, 1.058);
    hPhiPeak->GetYaxis()->SetRangeUser(0.0, hPhiPeak->GetMaximum() * 1.40);
    hPhiPeak->GetYaxis()->SetTitleOffset(1.3);

    SetFunctionStyle(fPhiFit, kRed, 1.0, 2, 1);
    SetFunctionStyle(fPhiBkg, kAzure + 2, 1.0, 2, 2);

    // Draw histogram and its attached functions
    hPhiPeak->Draw("E Y+");

    TF1 *fPhiSignal = new TF1(
        "fPhiSignal",
        [fPhiFit, fPhiBkg](double *xx, double *pp) -> double
        {
            return fPhiFit->Eval(xx[0]) - fPhiBkg->Eval(xx[0]);
        },
        0.992, 1.058, 0);

    fPhiSignal->SetRange(0.995, 1.05);

    SetFunctionStyle(fPhiSignal, kMagenta - 3, 1.0, 2, 1);
    fPhiSignal->SetNpx(500);
    // fPhiSignal->Draw("SAME");

    // -------------------------------------------------------------------------
    // Legend
    // -------------------------------------------------------------------------

    TLegend *legPhi2 = new TLegend(0.35, 0.8, 0.71, 0.96);
    SetLegend(legPhi2);
    legPhi2->AddEntry(hPhiPeak, "Data (stat. uncert.)", "p");
    legPhi2->AddEntry(fPhiFit, "Signal+Background fit", "l");
    legPhi2->AddEntry(fPhiBkg, "2nd-order Polynomial", "l");

    legPhi2->Draw();

    ltx.DrawLatex(0.04, 0.93, "(d)");

    // =========================================================================
    // 11. Save canvas
    // =========================================================================

    c1->cd();
    c1->SaveAs("../figures/KstarPhiInvariantMass.pdf");
}
