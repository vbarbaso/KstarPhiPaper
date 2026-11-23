// Draw phi and K*(892) multiplicity-differential spectra in the style of
// DrawFinal.C.  Run from the Phi-analysis directory so Functions.h and
// Style.h are available, e.g. .x outputs/DrawSpectra.C+
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <THnSparse.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TMath.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <string>
#include "../include/Functions.h"
#include "../include/StyleProposal.h"
#include "../include/YieldMean.C"

bool systematics = true;

namespace
{
    const char *kPhiFileDefault = "../inputs/Phi_136/inelgr0/PAR.root";
    const char *kKstarFileDefault = "../inputs/Kstar0_136/KstarSpectra.root";
    const char *kMult[10] = {"0 - 1 %", "1 - 5 %", "5 - 10 %", "10 - 15 %", "15 - 20 %", "20 - 30 %", "30 - 40 %", "40 - 50 %", "50 - 70 %", "70 - 100 %"};
    const double kPhiMass = 1.019461, kKstarMass = 0.895;

    TH1D *PhiProjection(THnSparseD *s, const char *variable, int cent = -1)
    {
        if (!s)
            return nullptr;
        int b = s->GetAxis(0)->FindBin(variable);
        if (b <= 0)
            return nullptr;
        s->GetAxis(0)->SetRange(b, b);
        if (cent > 0)
            s->GetAxis(2)->SetRange(cent, cent);
        TH1D *h = dynamic_cast<TH1D *>(s->Projection(1));
        s->GetAxis(0)->SetRange();
        if (cent > 0)
            s->GetAxis(2)->SetRange();
        if (h)
            h->SetDirectory(nullptr);
        return h;
    }

    TH1 *CloneHist(TFile *f, const TString &key)
    {
        TH1 *h = dynamic_cast<TH1 *>(f->Get(key));
        if (!h)
        {
            std::cerr << "ERROR: missing histogram '" << key << "'\n";
            return nullptr;
        }
        h = dynamic_cast<TH1 *>(h->Clone(Form("%s_clone", key.Data())));
        h->SetDirectory(nullptr);
        return h;
    }

    TGraphAsymmErrors *StatGraph(TH1 *h)
    {
        auto *g = new TGraphAsymmErrors(h->GetNbinsX());
        for (int b = 1; b <= h->GetNbinsX(); ++b)
        {
            int p = b - 1;
            double x = h->GetBinCenter(b), ex = h->GetBinWidth(b) / 2.;
            g->SetPoint(p, x, h->GetBinContent(b));
            g->SetPointError(p, ex, ex, h->GetBinError(b), h->GetBinError(b));
        }
        return g;
    }

    TGraphAsymmErrors *SysGraph(TH1 *central, TH1 *sys, bool relative = false)
    {
        auto *g = new TGraphAsymmErrors(central->GetNbinsX());
        for (int b = 1; b <= central->GetNbinsX(); ++b)
        {
            int p = b - 1;
            double y = central->GetBinContent(b), e = sys->GetBinError(b);
            if (relative)
                e *= y;
            double ex = central->GetBinWidth(b) / 2.;
            g->SetPoint(p, central->GetBinCenter(b), y);
            g->SetPointError(p, ex, ex, e, e);
        }
        return g;
    }

    TF1 *FitSpectrum(TH1 *h, TH1 *sys, const char *name, double mass,
                     int color, bool phi = false)
    {
        if (!h)
            return nullptr;

        double fitMin = 0.0;
        double fitMax = 20.0;
        double intMin = 0.0;
        double intMax = 20.0;
        Double_t loprecision = 0.01;
        Double_t hiprecision = 0.5;
        Option_t *opt = "RI0+";

        if (phi)
        {
            fitMin = 0.5;
            fitMax = 15.0;
            intMin = 0.;
            intMax = 15.0;
        }

        TF1 *f = nullptr;
        if (!phi)
        {
            f = new TF1(name, FuncLevy_Sawan, fitMin, fitMax, 4);
            f->SetParameter(0, 7.0);
            f->SetParameter(1, 0.4);
            f->FixParameter(2, mass);
            f->SetParameter(3, 0.25);
        }
        else
        {
            f = new TF1("Phi", FuncLevy_Sawan, fitMin, fitMax, 4);
            f->SetParameter(0, 8.0);
            f->SetParameter(1, 0.5);
            f->FixParameter(2, mass);
            f->SetParameter(3, 0.25);
        }

        TH1 *hout = YieldMean(h, sys, f, intMin, intMax, loprecision, hiprecision, opt, "log.root", fitMin, fitMax);

        f->SetLineColor(color);
        f->SetLineWidth(1);
        f->SetLineStyle(2);

        return f;
    }

    void ScaleGraph(TGraphAsymmErrors *g, double scale, bool omitLast)
    {
        for (int p = 0; p < g->GetN(); ++p)
        {
            double x, y;
            g->GetPoint(p, x, y);
            if (omitLast && p == g->GetN() - 1)
            {
                g->SetPoint(p, 0, 0);
                g->SetPointError(p, 0, 0, 0, 0);
                continue;
            }
            g->SetPoint(p, x, y * scale);
            g->SetPointError(p, g->GetErrorXlow(p), g->GetErrorXhigh(p), g->GetErrorYlow(p) * scale, g->GetErrorYhigh(p) * scale);
        }
    }

    struct Species
    {
        TString label;
        double mass;
        std::vector<TH1 *> stat, sys, uncorr;
    };

    void DrawOne(const Species &s, const char *out, bool phi, bool uncorrRelative, TH1 *ratioDen = nullptr)
    {
        gStyle->SetOptStat(0);

        TCanvas *c = new TCanvas(Form("c_%s", s.label.Data()), s.label, 700, 825);

        std::cout << "\nCanvas: " << c->GetName()
                  << "\nWidth:  " << c->GetWw()
                  << "\nHeight: " << c->GetWh()
                  << std::endl;

        TPad *top = new TPad(Form("top_%s", s.label.Data()), "", 0, .30, 1, 1);
        SetTopPad(top, true);
        TPad *bot = new TPad(Form("bot_%s", s.label.Data()), "", 0, 0, 1, .30);
        SetBottomPad(bot, true);

        top->cd();

        TH1 *frame = new TH1D(Form("frame_%s", s.label.Data()), "", 100, 0.0, s.stat[0]->GetXaxis()->GetXmax());
        frame->Reset();
        SetTopPadHistStyle(frame, kWhite, 0., 20, .5, 1, 1);
        frame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        frame->GetYaxis()->SetTitle("1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} [(GeV/#it{c})^{-1}]");
        frame->GetXaxis()->SetRangeUser(0.01, s.stat[0]->GetXaxis()->GetXmax());
        frame->GetYaxis()->SetRangeUser(.5e-10, 0.2e3);
        frame->Draw();

        TLegend *leg = new TLegend(.19, .03, .92, .2);
        SetLegend(leg);
        leg->SetTextSize(.025);
        leg->SetNColumns(4);
        leg->SetHeader("#bf{FT0M multiplicity event classes}");

        std::vector<TGraphAsymmErrors *> gs, gy;
        std::vector<TF1 *> fits;
        for (int i = 0; i < 10; ++i)
        {
            double scale = TMath::Power(2., 9 - i);

            auto *a = StatGraph(s.stat[i]);
            auto *b = SysGraph(s.stat[i], s.sys[i]);

            bool omit = phi && i == 9;

            ScaleGraph(a, scale, omit);
            ScaleGraph(b, scale, omit);
            SetGraphStyle(a, kColorsMult[i], 1., 20, .7, 1, 1);
            SetGraphStyle(b, kColorsMult[i], 1., 20, .7, 1, 1);
            b->SetFillStyle(0);
            b->SetLineWidth(1);
            a->SetFillStyle(0);
            a->SetLineWidth(1);

            if (systematics)
                b->Draw("2 SAME");
            a->Draw("P SAME");

            TF1 *fit = nullptr;

            fit = FitSpectrum(s.stat[i], s.sys[i], Form("fit_%s_%02d", s.label.Data(), i + 1), s.mass, kColorsMult[i], phi);
            fits.push_back(fit);

            top->cd();

            if (fit)
            {
                auto *line = new TGraph(400);
                for (int j = 0; j < 400; ++j)
                {
                    double x = 20. * j / 399.;
                    line->SetPoint(j, x, scale * fit->Eval(x));
                }
                line->SetLineColor(fit->GetLineColor());
                line->SetLineStyle(2);
                line->SetLineWidth(1);
                if (systematics)
                    line->Draw("L SAME");
            }
            gs.push_back(a);
            gy.push_back(b);
        }

        auto *gStat = StatGraph(s.stat[10]);
        auto *gSys = SysGraph(s.stat[10], s.sys[10]);
        ScaleGraph(gStat, 1e-2, false);
        ScaleGraph(gSys, 1e-2, false);
        SetGraphStyle(gStat, kBlack, 1., 20, .5, 1, 1);
        SetGraphStyle(gSys, kBlack, 1., 20, .5, 1, 1);
        gSys->SetFillStyle(0);
        gSys->SetLineWidth(1);
        if (systematics)
            gSys->Draw("2 SAME");
        gStat->Draw("P SAME");

        TF1 *fit = FitSpectrum(s.stat[10], s.sys[10], Form("fit_%s_INEL0", s.label.Data()), s.mass, kBlack, phi);
        fits.push_back(fit);

        top->cd();
        if (fit)
        {
            auto *line = new TGraph(400);
            for (int j = 0; j < 400; ++j)
            {
                double x = 20. * j / 399.;
                line->SetPoint(j, x, 1e-2 * fit->Eval(x));
            }
            line->SetLineColor(fit->GetLineColor());
            line->SetLineStyle(2);
            line->SetLineWidth(1);
            if (systematics)
                line->Draw("L SAME");
        }
        std::cout << "\nCanvas: " << c->GetName()
                  << "\nWidth:  " << c->GetWw()
                  << "\nHeight: " << c->GetWh()
                  << std::endl;

        for (int i = 0; i < 10; ++i)
        {
            int p = 9 - i;
            if (systematics)
                leg->AddEntry(gy[i], Form("%s%s", kMult[i], i < 9 ? Form(" (#times 2^{%d})", p) : ""), "Pfe");
            else
                leg->AddEntry(gy[i], Form("%s%s", kMult[i], i < 9 ? Form(" (#times 2^{%d})", p) : ""), "P");
        }
        if (systematics)
            leg->AddEntry(gSys, "INEL>0 (#times 10^{-2})", "Pfe");
        else
            leg->AddEntry(gSys, "INEL>0 (#times 10^{-2})", "P");

        if (fits[0])
            fits[0]->SetLineColor(kBlack); // black legend sample, matching the reference figure
        if (systematics)
            leg->AddEntry(fits[0], "Levy-Tsallis fit", "l");
        leg->Draw();

        TPaveText *txt = new TPaveText(.65, .7, .9, .9, "NDC");
        txt->SetFillStyle(0);
        txt->SetBorderSize(0);
        txt->SetTextAlign(31);
        txt->SetTextSize(.035);
        txt->SetTextFont(42);
        txt->AddText("ALICE Preliminary");
        txt->AddText("pp, 13.6 TeV");
        txt->AddText("|y| < 0.5");
        txt->AddText(s.label);

        txt->Draw();

        TPaveText *pt = new TPaveText(0.17, 0.2, 0.6, 0.24, "NDC");
        pt->SetFillStyle(0);
        pt->SetBorderSize(0);
        pt->SetTextAlign(11);
        pt->SetTextSize(0.025);
        pt->SetTextFont(42);

        pt->AddText("Uncertainties: stat. (bars), sys. (boxes)");
        pt->Draw();

        bot->cd();
        if (!ratioDen)
            ratioDen = s.stat[0];

        TH1 *rf = dynamic_cast<TH1 *>(frame->Clone(Form("ratio_frame_%s", s.label.Data())));
        rf->Reset();
        SetRatioHistStyle(rf, kWhite, 1., 20, .5, 1, 1);
        rf->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        rf->GetYaxis()->SetTitle("Ratio to INEL>0");
        rf->GetXaxis()->SetRangeUser(0, s.stat[0]->GetXaxis()->GetXmax());
        rf->GetYaxis()->SetRangeUser(.08, 14);
        rf->GetYaxis()->SetNdivisions(505);
        rf->Draw();

        std::cout << "\nCanvas: " << c->GetName()
                  << "\nWidth:  " << c->GetWw()
                  << "\nHeight: " << c->GetWh()
                  << std::endl;

        std::vector<TGraphAsymmErrors *> ratioSys;
        for (int i = 0; i < 10; ++i)
        {
            auto *r = new TGraphAsymmErrors(s.stat[i]->GetNbinsX());
            auto *q = new TGraphAsymmErrors(s.stat[i]->GetNbinsX());
            for (int b = 1; b <= s.stat[i]->GetNbinsX(); ++b)
            {
                int p = b - 1;
                double n = s.stat[i]->GetBinContent(b), d = ratioDen->GetBinContent(b), x = s.stat[i]->GetBinCenter(b), ex = s.stat[i]->GetBinWidth(b) / 2.;
                if (n <= 0 || d <= 0)
                    continue;
                double v = n / d, es = v * TMath::Sqrt(TMath::Power(s.stat[i]->GetBinError(b) / n, 2) + TMath::Power(ratioDen->GetBinError(b) / d, 2));
                r->SetPoint(p, x, v);
                r->SetPointError(p, ex, ex, es, es);

                double re = uncorrRelative ? s.uncorr[i]->GetBinContent(b) : s.uncorr[i]->GetBinError(b) / n;
                q->SetPoint(p, x, v);
                q->SetPointError(p, ex, ex, v * re, v * re);
            }
            if (phi && i == 9)
            {
                int p = r->GetN() - 1;
                if (p >= 0)
                {
                    r->SetPoint(p, 0, 0);
                    r->SetPointError(p, 0, 0, 0, 0);
                    q->SetPoint(p, 0, 0);
                    q->SetPointError(p, 0, 0, 0, 0);
                }
            }
            SetRatioGraphStyle(r, kColorsMult[i], 1., 20, .5, 1, 1);
            SetRatioGraphStyle(q, kColorsMult[i], 1., 20, .5, 1, 1);
            q->SetFillStyle(0);
            q->SetFillColorAlpha(kColorsMult[i], .3);
            q->SetLineWidth(1);
            if (systematics)
                q->Draw("2 SAME");
            r->Draw("P SAME");
            ratioSys.push_back(q);
        }
        TLine *one = new TLine(0, 1, s.stat[0]->GetXaxis()->GetXmax(), 1);
        one->SetLineStyle(2);
        one->Draw("SAME");

        c->SaveAs(Form("../figures/%s.pdf", out));
    }
}

void DrawCombined(const Species &ks, const Species &phi,
                  TH1 *kstarMB, TH1 *phiMB,
                  const char *out = "Combined_SpectraMultiplicity.pdf")
{
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0);

    // ------------------------------------------------------------
    // Canvas and pads
    // ------------------------------------------------------------

    TCanvas *c = new TCanvas("cCombined", "Combined spectra", 1300, 850);

    // K*0: left column
    // Phi: right column

    TPad *topK = new TPad("topK_combined", "", 0.0, 0.30, 0.5, 1.0);
    TPad *botK = new TPad("botK_combined", "", 0.0, 0.00, 0.5, 0.30);

    TPad *topPhi = new TPad("topPhi_combined", "", 0.5, 0.30, 1.0, 1.0);
    TPad *botPhi = new TPad("botPhi_combined", "", 0.5, 0.00, 1.0, 0.30);

    // Top pads
    topK->SetBottomMargin(0.00);
    topK->SetTopMargin(0.03);
    topK->SetLeftMargin(0.16);
    topK->SetRightMargin(0.005);
    topK->SetTickx();
    topK->SetTicky();
    topK->SetLogy();

    topPhi->SetBottomMargin(0.00);
    topPhi->SetTopMargin(0.03);
    topPhi->SetLeftMargin(0.009);
    topPhi->SetRightMargin(0.03);
    topPhi->SetTickx();
    topPhi->SetTicky();
    topPhi->SetLogy();

    // Bottom pads
    botK->SetTopMargin(0.00);
    botK->SetBottomMargin(0.35);
    botK->SetLeftMargin(0.16);
    botK->SetRightMargin(0.005);
    botK->SetTickx();
    botK->SetTicky();
    botK->SetLogy();

    botPhi->SetTopMargin(0.00);
    botPhi->SetBottomMargin(0.35);
    botPhi->SetLeftMargin(0.009);
    botPhi->SetRightMargin(0.03);
    botPhi->SetTickx();
    botPhi->SetTicky();
    botPhi->SetLogy();

    topK->Draw();
    botK->Draw();
    topPhi->Draw();
    botPhi->Draw();

    // ------------------------------------------------------------
    // Common drawing function for each particle
    // ------------------------------------------------------------

    auto DrawPanel = [&](const Species &s,
                         TPad *top, TPad *bot,
                         TH1 *ratioDen,
                         bool isPhi,
                         bool leftPanel,
                         bool uncorrRelative)
    {
        // const double xmax = s.stat[0]->GetXaxis()->GetXmax();
        const double xmin = isPhi ? 0.0 : 0.0;
        const double xmax = isPhi ? s.stat[0]->GetXaxis()->GetXmax() : 21.3;

        // --------------------------------------------------------
        // TOP: spectra
        // --------------------------------------------------------

        top->cd();

        // TH1D *frame = new TH1D(Form("frame_combined_%s", isPhi ? "phi" : "kstar"), "", 100, 0.0, xmax);
        TH1D *frame = new TH1D(Form("frame_combined_%s", isPhi ? "phi" : "kstar"), "", 100, xmin, xmax);

        frame->Reset();

        SetTopPadHistStyle(frame, kWhite, 0., 20, .7, 1, 1);

        frame->GetXaxis()->SetTitle("");
        frame->GetXaxis()->SetLabelSize(0.);
        frame->GetXaxis()->SetTickLength(0.02);

        frame->GetYaxis()->SetLabelOffset(0.005);
        frame->GetYaxis()->SetTitleOffset(1.15);
        frame->GetYaxis()->SetRangeUser(.4e-10, 0.2e3);

        if (leftPanel)
        {
            frame->GetYaxis()->SetTitle(
                "1/#it{N}_{ev} d^{2}#it{N}/d#it{p}_{T}d#it{y} [(GeV/#it{c})^{-1}]");
        }
        else
        {
            frame->GetYaxis()->SetTitle("");
            frame->GetYaxis()->SetLabelSize(0.);
            frame->GetYaxis()->SetTickLength(0.02);
        }

        frame->Draw();

        TLegend *leg = nullptr;

        if (leftPanel)
        {
            leg = new TLegend(.19, .03, .96, .2);
            SetLegend(leg);
            leg->SetTextSize(.025);
            leg->SetNColumns(4);
            leg->SetHeader("#bf{FT0M multiplicity event classes}");
        }

        std::vector<TGraphAsymmErrors *> statGraphs;
        std::vector<TGraphAsymmErrors *> sysGraphs;
        std::vector<TF1 *> fits;

        // Multiplicity classes
        for (int i = 0; i < 10; ++i)
        {
            double scale = TMath::Power(2., 9 - i);

            auto *a = StatGraph(s.stat[i]);
            auto *b = SysGraph(s.stat[i], s.sys[i]);

            bool omit = isPhi && i == 9;

            if (omit)
            {
                a->RemovePoint(a->GetN() - 1);
                b->RemovePoint(b->GetN() - 1);
            }

            ScaleGraph(a, scale, false);
            ScaleGraph(b, scale, false);

            SetGraphStyle(a, kColorsMult[i], 1., 20, .7, 1, 1);
            SetGraphStyle(b, kColorsMult[i], 1., 20, .7, 1, 1);

            a->SetFillStyle(0);
            b->SetFillStyle(0);
            b->SetLineWidth(1);

            if (systematics)
                b->Draw("2 SAME");

            a->Draw("P SAME");

            // Levy-Tsallis fit
            TF1 *fit = FitSpectrum(
                s.stat[i], s.sys[i],
                Form("fit_combined_%s_%02d",
                     isPhi ? "phi" : "kstar", i + 1),
                s.mass, kColorsMult[i], isPhi);

            top->cd();

            fits.push_back(fit);

            if (fit && systematics)
            {
                auto *line = new TGraph(400);

                for (int j = 0; j < 400; ++j)
                {
                    double x = 20. * j / 399.;
                    line->SetPoint(j, x, scale * fit->Eval(x));
                }

                line->SetLineColor(fit->GetLineColor());
                line->SetLineStyle(2);
                line->SetLineWidth(1);
                line->Draw("L SAME");
            }

            statGraphs.push_back(a);
            sysGraphs.push_back(b);
        }

        // INEL>0
        auto *gStat = StatGraph(s.stat[10]);
        auto *gSys = SysGraph(s.stat[10], s.sys[10]);

        ScaleGraph(gStat, 1e-2, false);
        ScaleGraph(gSys, 1e-2, false);

        SetGraphStyle(gStat, kBlack, 1., 20, .7, 1, 1);
        SetGraphStyle(gSys, kBlack, 1., 20, .7, 1, 1);

        gSys->SetFillStyle(0);
        gSys->SetLineWidth(1);

        if (systematics)
            gSys->Draw("2 SAME");

        gStat->Draw("P SAME");

        TF1 *fitMB = FitSpectrum(
            s.stat[10], s.sys[10],
            Form("fit_combined_%s_INEL0", isPhi ? "phi" : "kstar"),
            s.mass, kBlack, isPhi);

        top->cd();

        fits.push_back(fitMB);

        if (fitMB && systematics)
        {
            auto *line = new TGraph(400);

            for (int j = 0; j < 400; ++j)
            {
                double x = 20. * j / 399.;
                line->SetPoint(j, x, 1e-2 * fitMB->Eval(x));
            }

            line->SetLineColor(kBlack);
            line->SetLineStyle(2);
            line->SetLineWidth(1);
            line->Draw("L SAME");
        }

        // --------------------------------------------------------
        // Legend: only in the left panel (K*0)
        // --------------------------------------------------------

        if (leftPanel)
        {
            for (int i = 0; i < 10; ++i)
            {
                int p = 9 - i;

                TString label = Form("%s%s",
                                     kMult[i],
                                     i < 9 ? Form(" (#times 2^{%d})", p) : "");

                leg->AddEntry(
                    sysGraphs[i], label,
                    systematics ? "Pfe" : "P");
            }

            leg->AddEntry(
                gSys, "INEL>0 (#times 10^{-2})",
                systematics ? "Pfe" : "P");

            if (fits[0])
            {
                fits[0]->SetLineColor(kBlack);

                if (systematics)
                    leg->AddEntry(fits[0], "Levy-Tsallis fit", "l");
            }

            leg->Draw();

            // TPaveText *pt = new TPaveText(.17, .19, .65, .23, "NDC");
            // pt->SetFillStyle(0);
            // pt->SetBorderSize(0);
            // pt->SetTextAlign(11);
            // pt->SetTextSize(.025);
            // pt->SetTextFont(42);
            // pt->AddText("Uncertainties: stat. (bars), sys. (boxes)");
            // pt->Draw();
        }

        // --------------------------------------------------------
        // ALICE information: only in the right panel (Phi)
        // --------------------------------------------------------

        if (!leftPanel)
        {
            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(42);
            latex.SetTextSize(0.05);
            latex.DrawLatex(0.58, 0.89, "ALICE Preliminary");
            latex.SetTextSize(0.045);
            latex.DrawLatex(0.58, 0.84, "pp, #sqrt{#it{s}} = 13.6 TeV");
            latex.DrawLatex(0.68, 0.78, "|y| < 0.5");
            latex.DrawLatex(0.69, 0.72, s.label);

            latex.SetTextSize(0.03);
            latex.DrawLatex(0.05, 0.03, "Uncertainties: stat. (bars), sys. (boxes)");
        }
        else
        {
            TLatex KstarLabel;
            KstarLabel.SetNDC();
            KstarLabel.SetTextFont(42);
            KstarLabel.SetTextSize(0.045);
            KstarLabel.DrawLatex(0.82, 0.87, s.label);
        }

        // --------------------------------------------------------
        // BOTTOM: ratio to INEL>0
        // --------------------------------------------------------

        bot->cd();

        // TH1D *ratioFrame = new TH1D(Form("ratio_frame_combined_%s", isPhi ? "phi" : "kstar"), "", 100, 0.0, xmax);
        TH1D *ratioFrame = new TH1D(Form("ratio_frame_combined_%s", isPhi ? "phi" : "kstar"), "", 100, xmin, xmax);

        ratioFrame->Reset();

        SetRatioHistStyle(ratioFrame, kWhite, 1., 20, .7, 1, 1);

        ratioFrame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        ratioFrame->GetXaxis()->SetRangeUser(0., xmax);
        ratioFrame->GetYaxis()->SetRangeUser(.08, 14.);
        ratioFrame->GetYaxis()->SetNdivisions(505);

        if (!leftPanel)
        {
            ratioFrame->GetYaxis()->SetTitle("");
            ratioFrame->GetYaxis()->SetLabelSize(0.);
            ratioFrame->GetYaxis()->SetTickLength(.02);
        }
        else
        {
            ratioFrame->GetYaxis()->SetTitle("Ratio to INEL>0");
            ratioFrame->GetYaxis()->SetLabelOffset(0.005);
            ratioFrame->GetYaxis()->SetTitleOffset(0.75);
            ratioFrame->GetYaxis()->SetTitleSize(0.10);
            ratioFrame->GetYaxis()->CenterTitle(true);
        }

        ratioFrame->Draw();

        for (int i = 0; i < 10; ++i)
        {
            auto *r = new TGraphAsymmErrors(s.stat[i]->GetNbinsX());
            auto *q = new TGraphAsymmErrors(s.stat[i]->GetNbinsX());

            for (int b = 1; b <= s.stat[i]->GetNbinsX(); ++b)
            {
                // int p = b - 1;
                int p = r->GetN();

                double n = s.stat[i]->GetBinContent(b);
                double d = ratioDen->GetBinContent(b);
                double x = s.stat[i]->GetBinCenter(b);
                double ex = s.stat[i]->GetBinWidth(b) / 2.;

                if (n <= 0 || d <= 0)
                    continue;

                double v = n / d;

                double es = v * TMath::Sqrt(
                                    TMath::Power(s.stat[i]->GetBinError(b) / n, 2) +
                                    TMath::Power(ratioDen->GetBinError(b) / d, 2));

                r->SetPoint(p, x, v);
                r->SetPointError(p, ex, ex, es, es);

                double re = uncorrRelative
                                ? s.uncorr[i]->GetBinContent(b)
                                : s.uncorr[i]->GetBinError(b) / n;

                q->SetPoint(p, x, v);
                q->SetPointError(p, ex, ex, v * re, v * re);
            }

            // Preserve the original phi treatment of the last class
            if (isPhi && i == 9)
            {
                int p = r->GetN() - 1;

                if (p >= 0)
                {
                    r->RemovePoint(p);
                    q->RemovePoint(p);
                }
            }

            SetRatioGraphStyle(r, kColorsMult[i], 1., 20, .7, 1, 1);
            SetRatioGraphStyle(q, kColorsMult[i], 1., 20, .7, 1, 1);

            q->SetFillStyle(0);
            q->SetLineWidth(1);

            if (systematics)
                q->Draw("2 SAME");

            r->Draw("P SAME");
        }

        // TLine *one = new TLine(0, 1, xmax, 1);
        TLine *one = new TLine(xmin, 1, xmax, 1);
        one->SetLineStyle(2);
        one->Draw("SAME");
    };

    // ------------------------------------------------------------
    // Draw K*0 on the LEFT, phi on the RIGHT
    // ------------------------------------------------------------

    DrawPanel(ks, topK, botK, kstarMB,
              false, true, false);

    DrawPanel(phi, topPhi, botPhi, phiMB,
              true, false, true);

    c->cd();
    c->Modified();
    c->Update();

    gStyle->SetLineScalePS(1.0);

    c->SaveAs(Form("../figures/%s.pdf", out));
}

// Kstar uncertainty histogram convention: by default the bin CONTENT is the
// fractional uncorrelated systematic (the key is SpectraUncorrSys_*). Set
// kstarUncorrIsRelative=false if the bin ERROR stores an absolute uncertainty.

void DrawSpectra(const char *phiFileName = kPhiFileDefault,
                 const char *kstarFileName = kKstarFileDefault,
                 bool kstarUncorrIsRelative = false)
{
    TFile *fp = TFile::Open(phiFileName, "READ");
    if (!fp || fp->IsZombie())
    {
        std::cerr << "ERROR: cannot open phi file " << phiFileName << '\n';
        return;
    }

    auto *sm = dynamic_cast<THnSparseD *>(fp->Get("sResults"));
    auto *smb = dynamic_cast<THnSparseD *>(fp->Get("sResultsMB"));
    if (!sm || !smb)
    {
        std::cerr << "ERROR: phi sResults or sResultsMB missing\n";
        return;
    }

    Species phi;
    phi.label = "#phi(1020)";
    phi.mass = kPhiMass;

    for (int i = 0; i < 10; ++i)
    {
        phi.stat.push_back(PhiProjection(sm, "Yield_Corrected_Stat", i + 1));
        phi.sys.push_back(PhiProjection(sm, "Yield_Corrected_Sys_Mean", i + 1));
        phi.uncorr.push_back(PhiProjection(sm, "Sys_Uncorrelated_Mean", i + 1));
        if (!phi.stat.back() || !phi.sys.back() || !phi.uncorr.back())
        {
            std::cerr << "ERROR: missing phi stat/sys/uncorrelated spectrum for class " << i + 1 << '\n';
            return;
        }
    }
    phi.stat.push_back(PhiProjection(smb, "Yield_Corrected_Stat"));
    phi.sys.push_back(PhiProjection(smb, "Yield_Corrected_Sys_Mean"));

    if (!phi.stat.back() || !phi.sys.back() || !phi.uncorr.back())
    {
        std::cerr << "ERROR: missing phi stat/sys/uncorrelated spectrum for MB\n";
        return;
    }

    TFile *fk = TFile::Open(kstarFileName, "READ");
    if (!fk || fk->IsZombie())
    {
        std::cerr << "ERROR: cannot open Kstar file " << kstarFileName << '\n';
        return;
    }

    Species ks;
    ks.label = "K^{*0}(892)";
    ks.mass = kKstarMass;
    const char *classes[11] = {"0_1", "1_5", "5_10", "10_15", "15_20", "20_30", "30_40", "40_50", "50_70", "70_100", "0_100"};

    for (int i = 0; i < 11; ++i)
    {
        TString base = classes[i];
        ks.stat.push_back(CloneHist(fk, "SpectraStat_" + base));
        ks.sys.push_back(CloneHist(fk, "SpectraSys_" + base));
        if (i < 10)
        {
            ks.uncorr.push_back(CloneHist(fk, "SpectraUncorrSys_" + base));
        }
        if (!ks.stat.back() || !ks.sys.back() || !ks.uncorr.back())
        {
            std::cerr << "ERROR: Kstar file needs SpectraStat_*, SpectraSys_* and SpectraUncorrSys_* for all 10 classes. Ratios require the uncorrelated histogram; total systematics are not used as a substitute.\n";
            return;
        }
    }

    TH1D *phiMB = PhiProjection(smb, "Yield_Corrected_Stat");
    if (!phiMB)
    {
        std::cerr << "ERROR: phi MB spectrum missing\n";
        return;
    }

    DrawOne(phi, "Phi_SpectraMultiplicity", true, true, phiMB);

    TH1 *kstarMB = CloneHist(fk, "SpectraStat_0_100");
    if (!kstarMB)
    {
        std::cerr << "ERROR: Kstar INEL>0 reference SpectraStat_0_100 missing\n";
        return;
    }

    DrawOne(ks, "Kstar0_SpectraMultiplicity", false, false, kstarMB);

    DrawCombined(ks, phi, kstarMB, phiMB, "KStarPhi_SpectraMultiplicity");
}
