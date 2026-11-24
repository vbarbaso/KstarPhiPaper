// Kstar and Phi xT scaling in a common publication style.
// Run from the folder containing this macro and StyleProposal.h with:
//   root -l -q 'KstarPhiXTScaling.C()'

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGaxis.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <THnSparse.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSpline.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include "../include/StyleProposal.h"
#include "../include/Functions.h"
#include "../include/YieldMean.C"

namespace XTStyle
{
    // Same energy-to-colour and marker mapping for K*0 and phi
    constexpr int energyColors[] = {
        kMagenta - 3, kBlue, kCyan + 1, kYellow + 1, kAzure + 1, kBlack};

    constexpr int energyMarkers[] = {20, 20, 20, 20, 20, 20};

    constexpr const char *energyLabels[] = {"2.76 TeV", "5.02 TeV", "7 TeV", "8 TeV", "13 TeV", "13.6 TeV"};

    // Same exponent-ratio styles for both mesons
    constexpr int exponentColors[] = {
        kCyan + 1, kBlue, kMagenta - 3, kGreen + 2};

    constexpr int exponentMarkers[] = {20, 20, 20, 20};

    void SetEnergyGraph(TGraph *graph, int index)
    {
        SetGraphStyle(graph, energyColors[index], 1.0, energyMarkers[index], 0.8, 1);
    }

    void SetExponentGraph(TGraph *graph, int index)
    {
        SetGraphStyle(graph, exponentColors[index], 1.0, exponentMarkers[index], 0.8, 1);
    }
}

namespace KstarXT
{
    TFile *OpenFile(const std::string &path)
    {
        TFile *file = TFile::Open(path.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            std::cerr << "ERROR: Cannot open " << path << std::endl;
            return nullptr;
        }
        return file;
    }

    TH1 *GetHisto(TFile *file, const std::string &name)
    {
        if (!file)
            return nullptr;
        auto *hist = dynamic_cast<TH1 *>(file->Get(name.c_str()));
        if (!hist)
            std::cerr << "ERROR: Histogram " << name << " not found in " << file->GetName() << std::endl;
        return hist;
    }

    // In INEL_spectra.root the *_sys histogram stores absolute systematic
    // uncertainties in its bin errors; the histogram contents are the spectrum.
    // Convert those absolute errors to relative values for the existing helpers.
    TH1 *MakeRelativeSystematic(TH1 *stat, TH1 *sys, const char *name)
    {
        if (!stat || !sys || stat->GetNbinsX() != sys->GetNbinsX())
            return nullptr;

        for (int bin = 1; bin <= stat->GetNbinsX(); ++bin)
        {
            const double y = stat->GetBinContent(bin);
            const double sysY = sys->GetBinContent(bin);
            const double tolerance = 1.e-6 * std::max(1.0, std::abs(y));
            if (std::abs(y - sysY) > tolerance ||
                std::abs(stat->GetXaxis()->GetBinLowEdge(bin) -
                         sys->GetXaxis()->GetBinLowEdge(bin)) > 1.e-9 ||
                std::abs(stat->GetXaxis()->GetBinUpEdge(bin) -
                         sys->GetXaxis()->GetBinUpEdge(bin)) > 1.e-9)
            {
                std::cerr << "ERROR: " << name
                          << " does not share the spectrum binning/contents."
                          << std::endl;
                return nullptr;
            }
        }

        auto *relative = dynamic_cast<TH1 *>(sys->Clone(name));
        relative->SetDirectory(nullptr);
        for (int bin = 1; bin <= stat->GetNbinsX(); ++bin)
        {
            const double y = stat->GetBinContent(bin);
            relative->SetBinContent(
                bin, y > 0.0 ? sys->GetBinError(bin) / y : 0.0);
            relative->SetBinError(bin, 0.0);
        }
        return relative;
    }

    TGraphAsymmErrors *GetGraph(TFile *file, const std::string &name)
    {
        if (!file)
            return nullptr;
        auto *graph = dynamic_cast<TGraphAsymmErrors *>(file->Get(name.c_str()));
        if (!graph)
            std::cerr << "ERROR: Graph " << name << " not found in " << file->GetName() << std::endl;
        return graph;
    }

    TGraphErrors *HistToXTGraph(TH1 *stat, TH1 *relativeSys, double sqrtS, double sigmaInel, double energyFactor, bool alreadyNormalizedYield, double averageKstar)
    {
        auto *result = new TGraphErrors();
        for (int bin = 1; bin <= stat->GetNbinsX(); ++bin)
        {
            const double pt = stat->GetBinCenter(bin);
            const double dpt = stat->GetBinWidth(bin) / 2.0;
            const double y = stat->GetBinContent(bin);
            const double estat = stat->GetBinError(bin);
            const double esys = relativeSys ? relativeSys->GetBinContent(bin) * y : 0.0;
            const double total = std::sqrt(estat * estat + esys * esys);
            const double scale = sigmaInel * energyFactor * averageKstar;
            const double invY = alreadyNormalizedYield ? y * scale : y * scale / (2.0 * TMath::Pi() * pt);
            const double invErr = alreadyNormalizedYield ? total * scale : total * scale / (2.0 * TMath::Pi() * pt);
            const int point = result->GetN();
            result->SetPoint(point, 2.0 * pt / sqrtS, invY);
            result->SetPointError(point, 2.0 * dpt / sqrtS, invErr);
        }
        return result;
    }

    TGraphErrors *GraphToXTGraph(TGraphAsymmErrors *input, TH1 *stat, TH1 *sys,
                                 double sqrtS, double sigmaInel, double energyFactor,
                                 bool alreadyNormalizedYield, double averageKstar)
    {
        auto *result = new TGraphErrors();
        for (int i = 0; i < input->GetN(); ++i)
        {
            double pt, y;
            input->GetPoint(i, pt, y);
            const double statErr = stat->GetBinContent(i + 1);
            const double sysErr = sys->GetBinContent(i + 1);
            const double total = std::sqrt(statErr * statErr + sysErr * sysErr);
            const double scale = sigmaInel * energyFactor * averageKstar;
            const double invY = alreadyNormalizedYield ? y * scale : y * scale / (2.0 * TMath::Pi() * pt);
            const double invErr = alreadyNormalizedYield ? total * scale : total * scale / (2.0 * TMath::Pi() * pt);
            result->SetPoint(i, 2.0 * pt / sqrtS, invY);
            result->SetPointError(i, (input->GetErrorXlow(i) + input->GetErrorXhigh(i)) / sqrtS, invErr);
        }
        return result;
    }

    void GetGraphRange(TGraph *graph, double &xmin, double &xmax)
    {
        xmin = 1.e30;
        xmax = -1.e30;
        for (int i = 0; i < graph->GetN(); ++i)
        {
            double x, y;
            graph->GetPoint(i, x, y);
            xmin = std::min(xmin, x);
            xmax = std::max(xmax, x);
        }
    }

    struct InterpolatedInvariant
    {
        double value = 0.0;
        double error = 0.0;
        bool valid = false;
    };

    // Match PhiXT::EvalInvariant: interpolate the invariant yield in log(xT)
    // and log(y), and propagate endpoint relative errors with the same weights.
    InterpolatedInvariant EvalInvariant(TGraphErrors *graph, double x)
    {
        InterpolatedInvariant result;
        if (!graph || x <= 0.0 || graph->GetN() < 2)
            return result;

        for (int i = 0; i < graph->GetN() - 1; ++i)
        {
            double x0, y0, x1, y1;
            graph->GetPoint(i, x0, y0);
            graph->GetPoint(i + 1, x1, y1);
            if (x0 <= 0.0 || x1 <= 0.0 || y0 <= 0.0 || y1 <= 0.0 || x0 == x1)
                continue;
            if (!((x0 <= x && x <= x1) || (x1 <= x && x <= x0)))
                continue;

            const double t = std::log(x / x0) / std::log(x1 / x0);
            if (t < 0.0 || t > 1.0)
                continue;
            result.value = std::exp((1.0 - t) * std::log(y0) + t * std::log(y1));
            const double rel0 = graph->GetErrorY(i) / y0;
            const double rel1 = graph->GetErrorY(i + 1) / y1;
            result.error = result.value * std::sqrt(
                                              std::pow((1.0 - t) * rel0, 2) + std::pow(t * rel1, 2));
            result.valid = true;
            return result;
        }
        return result;
    }

    TGraphErrors *CalculateN(TGraphErrors *first, TGraphErrors *second, double sqrtS1, double sqrtS2)
    {
        auto *result = new TGraphErrors();
        const double denominator = std::log(sqrtS1 / sqrtS2);
        if (!first || !second || std::abs(denominator) < 1.e-12)
            return result;

        for (int i = 0; i < first->GetN(); ++i)
        {
            double x, y1;
            first->GetPoint(i, x, y1);
            if (x <= 0.0 || y1 <= 0.0)
                continue;

            const auto interp = EvalInvariant(second, x);
            if (!interp.valid || interp.value <= 0.0)
                continue;

            const double n = -std::log(y1 / interp.value) / denominator;
            const double en = std::sqrt(std::pow(first->GetErrorY(i) / y1, 2) + std::pow(interp.error / interp.value, 2)) / std::abs(denominator);
            const int point = result->GetN();
            result->SetPoint(point, x, n);
            result->SetPointError(point, 0.0, en);
        }
        return result;
    }

    TGraphErrors *CalculateN2(TGraphErrors *first, TGraphErrors *second, double sqrtS1, double sqrtS2)
    {
        // Use the same log-log interpolation and symmetric error propagation
        // as the other K*0 energy comparisons and the Phi calculation.
        return CalculateN(first, second, sqrtS1, sqrtS2);
    }

    void Draw(TPad *padN, TPad *padXT)
    {
        auto *fInel = OpenFile("../inputs/Kstar0_136/INEL_spectra.root");
        auto *f13 = OpenFile("../inputs/HEPData/HEPData-ins1797443-v1-root.root");
        auto *f7 = OpenFile("../inputs/HEPData/HEPData-ins1762364-v1-root.root");
        auto *f8 = OpenFile("../inputs/HEPData/HEPData-ins1762364-v1-root.root");
        auto *f276 = OpenFile("../inputs/HEPData/HEPData-ins1511864-v1-root.root");
        auto *f502 = OpenFile("../inputs/HEPData/HEPData-ins1762368-v1-root.root");

        if (!fInel || !f13 || !f7 || !f276 || !f502)
            return;

        auto *hStat136 = GetHisto(fInel, "hINEL_DataCorr_stat");
        auto *hSys136 = GetHisto(fInel, "hINEL_DataCorr_sys");
        auto *hRelSys136 = MakeRelativeSystematic(hStat136, hSys136, "hKstar136RelativeSys");

        auto *g13 = GetGraph(f13, "Table 4/Graph1D_y1");
        auto *hStat13 = GetHisto(f13, "Table 4/Hist1D_y1_e1");
        auto *hSys13 = GetHisto(f13, "Table 4/Hist1D_y1_e2");

        auto *g7 = GetGraph(f7, "Table 3/Graph1D_y1");
        auto *hStat7 = GetHisto(f7, "Table 3/Hist1D_y1_e1");
        auto *hSys7 = GetHisto(f7, "Table 3/Hist1D_y1_e2");

        auto *g8 = GetGraph(f8, "Table 1/Graph1D_y1");
        auto *hStat8 = GetHisto(f8, "Table 1/Hist1D_y1_e1");
        auto *hSys8 = GetHisto(f8, "Table 1/Hist1D_y1_e2");

        auto *g276 = GetGraph(f276, "Table 1/Graph1D_y1");
        auto *hStat276 = GetHisto(f276, "Table 1/Hist1D_y1_e1");
        auto *hSys276 = GetHisto(f276, "Table 1/Hist1D_y1_e2");

        auto *g502 = GetGraph(f502, "Table 2/Graph1D_y1");
        auto *hStat502 = GetHisto(f502, "Table 2/Hist1D_y1_e1");
        auto *hSys502 = GetHisto(f502, "Table 2/Hist1D_y1_e2");

        if (!hStat136 || !hRelSys136 || !g13 || !hStat13 || !hSys13 ||
            !g8 || !hStat8 || !hSys8 || !g7 || !hStat7 || !hSys7 || !g276 || !hStat276 || !hSys276 ||
            !g502 || !hStat502 || !hSys502)
            return;

        constexpr double s136 = 13600., s13 = 13000., s7 = 7000.;
        constexpr double s276 = 2760., s502 = 5020., s8 = 8000.;

        constexpr double sigma136 = 77.904, sigma13 = 77.6, sigma7 = 70.9;
        constexpr double sigma276 = 61.8, sigma502 = 77.6, sigma8 = 71.5;

        auto *n136 = HistToXTGraph(hStat136, hRelSys136, s136, sigma136, 1.0, false, 1.0);
        auto *n13 = GraphToXTGraph(g13, hStat13, hSys13, s13, sigma13, 1.0, false, 1.0);
        auto *n7 = GraphToXTGraph(g7, hStat7, hSys7, s7, sigma7, 1.0, false, 2.0);
        auto *n8 = GraphToXTGraph(g8, hStat8, hSys8, s8, sigma8, 1.0, false, 2.0);
        auto *n276 = GraphToXTGraph(g276, hStat276, hSys276, s276, sigma276, 1.0, true, 2.0);
        auto *n502 = GraphToXTGraph(g502, hStat502, hSys502, s502, sigma502, 1.0, false, 2.0);

        // ============================================================
        // Determine n(xT)
        // ============================================================

        auto *gN7 = CalculateN(n136, n7, s136, s7);
        auto *gN276 = CalculateN(n136, n276, s136, s276);
        auto *gN7_276 = CalculateN(n7, n276, s7, s276);
        auto *gN13_276 = CalculateN(n13, n276, s13, s276);
        auto *gN13_7 = CalculateN(n13, n7, s13, s7);
        auto *gN8_276 = CalculateN(n8, n276, s8, s276);
        auto *gN502 = CalculateN2(n136, n502, s136, s502);

        std::vector<TGraphErrors *> nGraphs = {gN7, gN502, gN276, gN7_276, gN8_276};

        const char *nLabels[] = {
            "n(#sigma_{13.6}/#sigma_{7})",
            "n(#sigma_{13.6}/#sigma_{5.02})",
            "n(#sigma_{13.6}/#sigma_{2.76})",
            "n(#sigma_{7}/#sigma_{2.76})"};

        auto *allN = new TGraphErrors();

        for (auto *graph : nGraphs)
        {
            for (int i = 0; i < graph->GetN(); ++i)
            {
                double x, y;
                graph->GetPoint(i, x, y);

                const int p = allN->GetN();
                allN->SetPoint(p, x, y);
                allN->SetPointError(p, graph->GetErrorX(i), graph->GetErrorY(i));
            }
        }

        auto *fitN = new TF1("kstar_n_default", "[0]", 1.6e-3, 4.1e-3);
        fitN->SetParameter(0, 4.7);
        allN->Fit(fitN, "REBMS0");

        auto *fitN1 = new TF1("kstar_n_sys1", "[0]", 1.3e-3, 3.4e-3);
        fitN1->SetParameter(0, 4.7);
        allN->Fit(fitN1, "REBMS0");

        auto *fitN2 = new TF1("kstar_n_sys2", "[0]", 1.9e-3, 4.4e-3);
        fitN2->SetParameter(0, 4.7);
        allN->Fit(fitN2, "REBMS0");

        const double nAverage = fitN->GetParameter(0);
        const double nSys1 = fitN1->GetParameter(0);
        const double nSys2 = fitN2->GetParameter(0);

        const double nSys = std::sqrt(
            (std::pow(nSys1 - nAverage, 2) + std::pow(nSys2 - nAverage, 2)) / 2.0);
        const double nTotalUncertainty =
            std::sqrt(std::pow(fitN->GetParError(0), 2) + std::pow(nSys, 2));

        std::cout << "Kstar n = " << nAverage << " +/- "
                  << fitN->GetParError(0) << " (stat) +/- "
                  << nSys << " (fit range), total +/- "
                  << nTotalUncertainty << std::endl;

        if (padN)
        {
            padN->cd();

            auto *frameN = padN->DrawFrame(0.0, -0.95, 4.7e-3, 7.4);
            frameN->SetTitle("");
            frameN->GetXaxis()->SetTitle("#it{x_{#it{T}}} (#times 10^{-3})");
            frameN->GetYaxis()->SetTitle("n(x_{T})");
            SetHistStyle(frameN);

            auto *legendN = new TLegend(0.5, 0.4, 0.9, 0.6);
            SetLegend(legendN);
            legendN->SetTextSize(0.027);
            legendN->SetNColumns(2);

            for (size_t i = 0; i < nGraphs.size() - 1; ++i)
            {
                XTStyle::SetExponentGraph(nGraphs[i], i);
                nGraphs[i]->Draw("P E SAME");
                legendN->AddEntry(nGraphs[i], nLabels[i], "p");
            }

            SetFunctionStyle(fitN, kBlack, 1.0, 1, 1);
            SetFunctionStyle(fitN1, kBlue + 2, 1.0, 1, 2);
            SetFunctionStyle(fitN2, kRed + 2, 1.0, 1, 2);

            fitN->Draw("SAME");
            fitN1->Draw("SAME");
            fitN2->Draw("SAME");
            // legendN->Draw();

            TLegend *legendFit = new TLegend(0.57, 0.2, 0.9, 0.3);
            SetLegend(legendFit);
            legendFit->SetTextSize(0.027);
            legendFit->AddEntry(fitN, Form("Default n = %.2f #pm %.2f", nAverage, nTotalUncertainty), "l");
            legendFit->AddEntry(fitN1, Form("Range 1 n = %.2f", nSys1), "l");
            legendFit->AddEntry(fitN2, Form("Range 2 n = %.2f", nSys2), "l");
            legendFit->Draw();

            TLatex label;
            label.SetNDC();
            label.SetTextFont(42);
            label.SetTextSize(0.065);
            label.DrawLatex(0.87, 0.87, "K*^{0}");

            label.SetTextSize(0.055);
            label.DrawLatex(0.23, 0.87, "#bf{ALICE}");
            label.SetTextSize(0.045);
            label.DrawLatex(0.23, 0.81, "pp, |y| < 0.5");
        }

        // ============================================================
        // xT distributions
        // ============================================================

        if (padXT)
        {
            auto *x136 = HistToXTGraph(hStat136, hRelSys136, s136, sigma136,
                                       std::pow(s136, nAverage), false, 1.0);
            auto *x13 = GraphToXTGraph(g13, hStat13, hSys13, s13, sigma13,
                                       std::pow(s13, nAverage), false, 1.0);
            auto *x7 = GraphToXTGraph(g7, hStat7, hSys7, s7, sigma7,
                                      std::pow(s7, nAverage), false, 2.0);
            auto *x8 = GraphToXTGraph(g8, hStat8, hSys8, s8, sigma8,
                                      std::pow(s8, nAverage), false, 2.0);
            auto *x276 = GraphToXTGraph(g276, hStat276, hSys276, s276, sigma276,
                                        std::pow(s276, nAverage), true, 2.0);
            auto *x502 = GraphToXTGraph(g502, hStat502, hSys502, s502, sigma502,
                                        std::pow(s502, nAverage), false, 2.0);

            TGraphErrors *xGraphs[] = {x276, x502, x7, x8, x13, x136};

            padXT->cd();
            padXT->SetLogx();
            padXT->SetLogy();

            auto *frameXT = padXT->DrawFrame(9.e-7, 8.e9, 3.e-2, 9.9e21);
            frameXT->SetTitle("");
            frameXT->GetXaxis()->SetTitle("#it{x_{#it{T}}}");
            frameXT->GetYaxis()->SetTitle(Form("#sqrt{s}^{%.2f} d^{3}#sigma/d^{3}p (mb GeV^{ -2}#it{c}^{ 3})", nAverage));

            SetHistStyle(frameXT);
            frameXT->GetYaxis()->SetNdivisions(505);

            auto *allXT = new TGraphErrors();

            for (int j = 0; j < 6; ++j)
            {
                XTStyle::SetEnergyGraph(xGraphs[j], j);
                xGraphs[j]->Draw("P E SAME");

                if (j == 0 || j == 1 || j == 2 || j == 4)
                {
                    for (int i = 0; i < xGraphs[j]->GetN(); ++i)
                    {
                        double x, y;
                        xGraphs[j]->GetPoint(i, x, y);

                        const int p = allXT->GetN();
                        allXT->SetPoint(p, x, y);
                        allXT->SetPointError(p, xGraphs[j]->GetErrorX(i),
                                             xGraphs[j]->GetErrorY(i));
                    }
                }
            }

            auto *fitXT = new TF1("kstar_xt_power",
                                  "[0]*pow(x,[1])*pow(1+x,[2])", 1.6e-3, 9.e-3);
            fitXT->SetParameters(1.8, -6.0, 5.0);
            fitXT->SetParLimits(0, -1000, 1000);
            fitXT->SetParLimits(1, -10, 10);
            fitXT->SetParLimits(2, -1000, 1000);

            SetFunctionStyle(fitXT, kRed, 1.0, 2, 1);
            allXT->Fit(fitXT, "REBMS0");
            fitXT->Draw("SAME");

            auto *legendXT = new TLegend(0.2, 0.18, 0.5, 0.5);
            SetLegend(legendXT);
            legendXT->SetTextSize(0.032);

            for (int j = 0; j < 6; ++j)
                legendXT->AddEntry(xGraphs[j], XTStyle::energyLabels[j], "p");

            legendXT->AddEntry(fitXT, "Power-law fit", "l");
            legendXT->Draw();

            TLatex xLabel;
            xLabel.SetNDC();
            xLabel.SetTextFont(42);

            // xLabel.SetTextSize(0.055);
            // xLabel.DrawLatex(0.8, 0.87, "#bf{ALICE}");
            // xLabel.SetTextSize(0.045);
            // xLabel.DrawLatex(0.74, 0.81, "pp, |y| < 0.5");

            xLabel.SetTextSize(0.065);
            xLabel.DrawLatex(0.87, 0.87, "K*^{0}");
        }
    }
} // namespace KstarXT

namespace PhiXT
{
    struct InvariantPoint
    {
        double value = 0.0;
        double errLow = 0.0;
        double errHigh = 0.0;
        bool valid = false;
    };

    InvariantPoint EvalInvariant(TGraphAsymmErrors *graph, double pt)
    {
        InvariantPoint result;
        if (!graph || pt <= 0.0 || graph->GetN() < 2)
            return result;
        const int n = graph->GetN();
        std::vector<double> pts(n), inv(n), relLow(n), relHigh(n);
        for (int i = 0; i < n; ++i)
        {
            double x, y;
            graph->GetPoint(i, x, y);
            if (x <= 0.0 || y <= 0.0)
                return result;
            pts[i] = x;
            inv[i] = y / (2.0 * TMath::Pi() * x);
            relLow[i] = graph->GetErrorYlow(i) / y;
            relHigh[i] = graph->GetErrorYhigh(i) / y;
        }
        int low = -1;
        for (int i = 0; i < n - 1; ++i)
            if ((pts[i] <= pt && pt <= pts[i + 1]) || (pts[i + 1] <= pt && pt <= pts[i]))
            {
                low = i;
                break;
            }
        if (low < 0)
            return result;
        const int high = low + 1;
        const double t = std::log(pt / pts[low]) / std::log(pts[high] / pts[low]);
        if (t < 0.0 || t > 1.0)
            return result;
        result.value = std::exp((1.0 - t) * std::log(inv[low]) + t * std::log(inv[high]));
        const double relErrLow = std::sqrt(std::pow((1.0 - t) * relLow[low], 2) + std::pow(t * relLow[high], 2));
        const double relErrHigh = std::sqrt(std::pow((1.0 - t) * relHigh[low], 2) + std::pow(t * relHigh[high], 2));
        result.errLow = result.value * relErrLow;
        result.errHigh = result.value * relErrHigh;
        result.valid = true;
        return result;
    }

    TGraphAsymmErrors *MakeExponentGraph(TGraphAsymmErrors *reference, double sqrtSRef,
                                         TGraphAsymmErrors *comparison, double sqrtSComp,
                                         const char *name)
    {
        auto *result = new TGraphAsymmErrors();
        result->SetName(name);
        if (!reference || !comparison || sqrtSRef <= sqrtSComp || sqrtSComp <= 0.0)
            return result;
        const double denominator = std::log(sqrtSRef / sqrtSComp);
        for (int i = 0; i < reference->GetN(); ++i)
        {
            double ptRef, sigmaRef;
            reference->GetPoint(i, ptRef, sigmaRef);
            if (ptRef <= 0.0 || sigmaRef <= 0.0)
                continue;
            const double xt = 2.0 * ptRef / sqrtSRef;
            const double ptComp = xt * sqrtSComp / 2.0;
            const auto invRef = EvalInvariant(reference, ptRef);
            const auto invComp = EvalInvariant(comparison, ptComp);
            if (!invRef.valid || !invComp.valid || invRef.value <= 0.0 || invComp.value <= 0.0)
                continue;
            const double n = std::log(invComp.value / invRef.value) / denominator;
            const double errLow = std::sqrt(std::pow(invComp.errLow / invComp.value, 2) + std::pow(invRef.errHigh / invRef.value, 2)) / denominator;
            const double errHigh = std::sqrt(std::pow(invComp.errHigh / invComp.value, 2) + std::pow(invRef.errLow / invRef.value, 2)) / denominator;
            const int p = result->GetN();
            result->SetPoint(p, 1.e3 * xt, n);
            // result->SetPointError(p, 1.e3 * 2.0 * reference->GetErrorXlow(i) / sqrtSRef,
            //                       1.e3 * 2.0 * reference->GetErrorXhigh(i) / sqrtSRef,
            //                       errLow, errHigh);
            result->SetPointError(p, 0.0, 0.0,
                                  errLow, errHigh);
        }
        return result;
    }

    TGraphAsymmErrors *MakeXTGraph(TGraphAsymmErrors *graph, double sqrtS, double nScaling, const char *name)
    {
        auto *result = new TGraphAsymmErrors();
        result->SetName(name);
        const double energyFactor = std::pow(sqrtS, nScaling);
        for (int i = 0; i < graph->GetN(); ++i)
        {
            double pt, crossSection;
            graph->GetPoint(i, pt, crossSection);
            if (pt <= 0.0 || crossSection <= 0.0)
                continue;
            const double xt = 2.0 * pt / sqrtS;
            const double inv = crossSection / (2.0 * TMath::Pi() * pt);
            const double scale = energyFactor / (2.0 * TMath::Pi() * pt);
            const int p = result->GetN();
            result->SetPoint(p, xt, energyFactor * inv);
            result->SetPointError(p, 2.0 * graph->GetErrorXlow(i) / sqrtS,
                                  2.0 * graph->GetErrorXhigh(i) / sqrtS,
                                  scale * graph->GetErrorYlow(i), scale * graph->GetErrorYhigh(i));
        }
        return result;
    }

    TGraphAsymmErrors *CombineExponentGraphs(const std::vector<TGraphAsymmErrors *> &graphs)
    {
        auto *combined = new TGraphAsymmErrors();
        for (auto *graph : graphs)
            for (int i = 0; i < graph->GetN(); ++i)
            {
                double x, y;
                graph->GetPoint(i, x, y);
                const int p = combined->GetN();
                combined->SetPoint(p, x, y);
                combined->SetPointError(p, graph->GetErrorXlow(i), graph->GetErrorXhigh(i),
                                        graph->GetErrorYlow(i), graph->GetErrorYhigh(i));
            }
        return combined;
    }

    TF1 *FitExponent(TGraphAsymmErrors *data, const char *name, double xmin, double xmax, int color, int style)
    {
        auto *fit = new TF1(name, "[0]", xmin, xmax);
        fit->SetParameter(0, 5.0);
        SetFunctionStyle(fit, color, 1.0, 2, style);
        TFitResultPtr result = data->Fit(fit, "QRSN");
        if (static_cast<int>(result) != 0 || result->Ndf() <= 0)
        {
            std::cerr << "Phi: constant fit failed: " << name << std::endl;
            delete fit;
            return nullptr;
        }
        return fit;
    }

    void ScaleGraphToCrossSection(TGraphAsymmErrors *graph, double sigmaInel)
    {
        for (int i = 0; i < graph->GetN(); ++i)
        {
            double x, y;
            graph->GetPoint(i, x, y);
            graph->SetPoint(i, x, y * sigmaInel);
            graph->SetPointError(i, graph->GetErrorXlow(i), graph->GetErrorXhigh(i),
                                 graph->GetErrorYlow(i) * sigmaInel,
                                 graph->GetErrorYhigh(i) * sigmaInel);
        }
    }

    void Draw(TPad *padN, TPad *padXT)
    {
        TFile *files[] = {
            TFile::Open("../inputs/Phi_136/inel/PAR.root", "READ"),
            TFile::Open("../inputs/HEPData/HEPData-ins1762364-v1-root.root", "READ"),
            TFile::Open("../inputs/HEPData/HEPData-ins1797443-v1-root.root", "READ"),
            TFile::Open("../inputs/HEPData/HEPData-ins1762368-v1-root.root", "READ"),
            TFile::Open("../inputs/HEPData/HEPData-ins1511864-v1-root.root", "READ")};

        for (auto *file : files)
        {
            if (!file || file->IsZombie())
            {
                std::cerr << "Phi: could not open one of the source ROOT files."
                          << std::endl;
                return;
            }
        }

        auto *g5 = dynamic_cast<TGraphAsymmErrors *>(files[3]->Get("Table 4/Graph1D_y1"));
        auto *g276 = dynamic_cast<TGraphAsymmErrors *>(files[4]->Get("Table 2/Graph1D_y1"));
        auto *g7 = dynamic_cast<TGraphAsymmErrors *>(files[1]->Get("Table 4/Graph1D_y1"));
        auto *g8 = dynamic_cast<TGraphAsymmErrors *>(files[1]->Get("Table 2/Graph1D_y1"));
        auto *g13 = dynamic_cast<TGraphAsymmErrors *>(files[2]->Get("Table 5/Graph1D_y1"));

        auto *sparse = dynamic_cast<THnSparseD *>(files[0]->Get("sResultsMB"));

        if (!g5 || !g276 || !g7 || !g8 || !g13 || !sparse)
        {
            std::cerr << "Phi: a source graph or sResultsMB is missing."
                      << std::endl;
            return;
        }

        const int statBin = sparse->GetAxis(0)->FindBin("Yield_Corrected_Stat");
        const int sysBin = sparse->GetAxis(0)->FindBin("Yield_Corrected_Sys_Mean");

        sparse->GetAxis(0)->SetRange(statBin, statBin);
        auto *h136Stat = dynamic_cast<TH1D *>(sparse->Projection(1));

        sparse->GetAxis(0)->SetRange(sysBin, sysBin);
        auto *h136Sys = dynamic_cast<TH1D *>(sparse->Projection(1));

        if (!h136Stat || !h136Sys)
            return;

        auto *g136Stat = new TGraphAsymmErrors(h136Stat);
        auto *g136Sys = new TGraphAsymmErrors(h136Sys);
        auto *g136 = new TGraphAsymmErrors(h136Stat);

        for (int i = 0; i < g136->GetN(); ++i)
        {
            double x, y;
            g136->GetPoint(i, x, y);

            const double stat = g136Stat->GetErrorY(i);
            const double sys = g136Sys->GetErrorY(i);
            const double total = std::sqrt(stat * stat + sys * sys);

            g136->SetPointError(i, g136Stat->GetErrorXlow(i),
                                g136Stat->GetErrorXhigh(i), total, total);
        }

        for (int i = 0; i < g276->GetN(); ++i)
        {
            double x, y;
            g276->GetPoint(i, x, y);
            g276->SetPoint(i, x, y * 2.0 * TMath::Pi() * x);
        }

        auto *xsec276 = new TGraphAsymmErrors(*g276);
        auto *xsec5 = new TGraphAsymmErrors(*g5);
        auto *xsec7 = new TGraphAsymmErrors(*g7);
        auto *xsec8 = new TGraphAsymmErrors(*g8);
        auto *xsec13 = new TGraphAsymmErrors(*g13);
        auto *xsec136 = new TGraphAsymmErrors(*g136);

        constexpr double sig276 = 61.8, sig5 = 77.6, sig7 = 70.9;
        constexpr double sig8 = 71.5, sig13 = 77.6, sig136 = 77.904;

        ScaleGraphToCrossSection(xsec276, sig276);
        ScaleGraphToCrossSection(xsec5, sig5);
        ScaleGraphToCrossSection(xsec7, sig7);
        ScaleGraphToCrossSection(xsec8, sig8);
        ScaleGraphToCrossSection(xsec13, sig13);
        ScaleGraphToCrossSection(xsec136, sig136);

        // ============================================================
        // Determine n(xT)
        // ============================================================

        auto *n136_5 = MakeExponentGraph(xsec136, 13600., xsec5, 5020., "phi_n_136_5");
        auto *n136_7 = MakeExponentGraph(xsec136, 13600., xsec7, 7000., "phi_n_136_7");
        auto *n136_8 = MakeExponentGraph(xsec136, 13600., xsec8, 8000., "phi_n_136_8");
        auto *n7_5 = MakeExponentGraph(xsec7, 7000., xsec5, 5020., "phi_n_7_5");
        auto *n8_5 = MakeExponentGraph(xsec8, 8000., xsec5, 5020., "phi_n_8_5");
        auto *n8_7 = MakeExponentGraph(xsec8, 8000., xsec7, 7000., "phi_n_8_7");
        auto *n136_276 = MakeExponentGraph(xsec136, 13600., xsec276, 2760., "phi_n_136_276");
        auto *n7_276 = MakeExponentGraph(xsec7, 7000., xsec276, 2760., "phi_n_7_276");

        auto *allN = CombineExponentGraphs({n136_5, n136_7, n136_8, n8_5, n136_276, n7_276});

        auto *fit1 = FitExponent(allN, "phi_n_range1", 1.9, 4.0, kRed + 1, 2);
        auto *fit2 = FitExponent(allN, "phi_n_range2", 2.5, 4.5, kBlue + 1, 2);

        if (!fit1 || !fit2)
            return;

        const double nCentral = 0.5 * (fit1->GetParameter(0) + fit2->GetParameter(0));
        const double nRangeSys = 0.5 * std::abs(fit2->GetParameter(0) - fit1->GetParameter(0));

        std::cout << "Phi n = " << nCentral << " +/- "
                  << nRangeSys << " (fit range)" << std::endl;

        if (padN)
        {
            padN->cd();

            auto *frameN = padN->DrawFrame(0.001, -0.95, 4.7, 7.4);
            frameN->SetTitle("");
            frameN->GetXaxis()->SetTitle("#it{x_{#it{T}}} (#times 10^{-3})");
            frameN->GetYaxis()->SetTitle("n(x_{T})");
            SetHistStyle(frameN);

            TGraphAsymmErrors *drawN[] = {
                n136_7, n136_5, n136_276, n7_276};

            const char *labelsN[] = {
                "n(#sigma_{13.6}/#sigma_{7})",
                "n(#sigma_{13.6}/#sigma_{5.02})",
                "n(#sigma_{13.6}/#sigma_{2.76})",
                "n(#sigma_{7}/#sigma_{2.76})"};

            auto *legendN = new TLegend(0.25, 0.35, 0.75, 0.55);
            SetLegend(legendN);
            legendN->SetTextSize(0.035);
            legendN->SetNColumns(2);

            for (int i = 0; i < 4; ++i)
            {
                XTStyle::SetExponentGraph(drawN[i], i);
                drawN[i]->Draw("P E SAME");
                legendN->AddEntry(drawN[i], labelsN[i], "p");
            }

            auto *fitCentral = new TF1("phi_n_central", "[0]", 1.9, 4.5);
            fitCentral->SetParameter(0, nCentral);

            SetFunctionStyle(fitCentral, kBlack, 1.0, 1, 1);
            SetFunctionStyle(fit1, kRed + 1, 1.0, 1, 2);
            SetFunctionStyle(fit2, kBlue + 1, 1.0, 1, 2);

            fit1->Draw("SAME");
            fit2->Draw("SAME");
            fitCentral->Draw("SAME");
            legendN->Draw();

            TLegend *legendFit = new TLegend(0.4, 0.2, 0.75, 0.3);
            SetLegend(legendFit);
            legendFit->SetTextSize(0.027);

            legendFit->AddEntry(fitCentral, Form("Default n = %.2f #pm %.2f", nCentral, nRangeSys), "l");
            legendFit->AddEntry(fit1, Form("Range 1: n = %.2f", fit1->GetParameter(0)), "l");
            legendFit->AddEntry(fit2, Form("Range 2: n = %.2f", fit2->GetParameter(0)), "l");
            legendFit->Draw();

            TLatex label;
            label.SetNDC();
            label.SetTextFont(42);
            label.SetTextSize(0.065);
            label.DrawLatex(0.74, 0.88, "#phi");
        }

        // ============================================================
        // xT distributions
        // ============================================================

        if (padXT)
        {
            auto *xt276 = MakeXTGraph(xsec276, 2760., nCentral, "phi_xt_276");
            auto *xt5 = MakeXTGraph(xsec5, 5020., nCentral, "phi_xt_5");
            auto *xt7 = MakeXTGraph(xsec7, 7000., nCentral, "phi_xt_7");
            auto *xt8 = MakeXTGraph(xsec8, 8000., nCentral, "phi_xt_8");
            auto *xt13 = MakeXTGraph(xsec13, 13000., nCentral, "phi_xt_13");
            auto *xt136 = MakeXTGraph(xsec136, 13600., nCentral, "phi_xt_136");

            XTStyle::SetEnergyGraph(xt136, 5);
            XTStyle::SetEnergyGraph(xt13, 4);
            XTStyle::SetEnergyGraph(xt8, 3);
            XTStyle::SetEnergyGraph(xt7, 2);
            XTStyle::SetEnergyGraph(xt5, 1);
            XTStyle::SetEnergyGraph(xt276, 0);

            padXT->cd();
            padXT->SetLogx();
            padXT->SetLogy();

            auto *frameXT = padXT->DrawFrame(3.e-5, 8.e9, 3.e-2, 9.9e21);
            frameXT->SetTitle("");
            frameXT->GetXaxis()->SetTitle("#it{x_{#it{T}}}");
            frameXT->GetYaxis()->SetTitle(Form("#sqrt{s}^{%.2f} d^{3}#sigma/d^{3}p (mb GeV^{ -2}#it{c}^{ 3})", nCentral));
            SetHistStyle(frameXT);
            frameXT->GetYaxis()->SetNdivisions(505);
            frameXT->Draw("E Y+");

            auto *mg = new TMultiGraph();
            mg->Add(xt276, "P");
            mg->Add(xt5, "P");
            mg->Add(xt7, "P");
            mg->Add(xt8, "P");
            mg->Add(xt13, "P");
            mg->Add(xt136, "P");

            mg->Draw("P");

            auto *fitPower = new TF1("phi_xt_power",
                                     "[0]*TMath::Power(x/1.e-3,-[1])",
                                     0.0015, 0.015);

            fitPower->SetParNames("F(xT = 0.001)", "n");
            fitPower->SetParameters(1.e17, 5.0);
            fitPower->SetParLimits(0, 0.0, 1.e25);
            fitPower->SetParLimits(1, 0.0, 20.0);

            SetFunctionStyle(fitPower, kRed, 1.0, 2, 1);

            TFitResultPtr fitResult = mg->Fit(fitPower, "RQS");

            if (static_cast<int>(fitResult) == 0)
                fitPower->Draw("SAME");
            else
                std::cerr << "Phi xT power-law fit status: "
                          << static_cast<int>(fitResult) << std::endl;

            auto *legendXT = new TLegend(0.19, 0.17, 0.48, 0.45);
            SetLegend(legendXT);
            legendXT->SetTextSize(0.032);

            legendXT->AddEntry(xt136, XTStyle::energyLabels[0], "p");
            legendXT->AddEntry(xt13, XTStyle::energyLabels[1], "p");
            legendXT->AddEntry(xt8, XTStyle::energyLabels[2], "p");
            legendXT->AddEntry(xt7, XTStyle::energyLabels[2], "p");
            legendXT->AddEntry(xt5, XTStyle::energyLabels[3], "p");
            legendXT->AddEntry(xt276, XTStyle::energyLabels[4], "p");
            legendXT->AddEntry(fitPower, "Power-law fit", "l");
            // legendXT->Draw();

            TLatex labelXT;
            labelXT.SetNDC();
            labelXT.SetTextFont(42);

            // labelXT.SetTextSize(0.055);
            // labelXT.DrawLatex(0.62, 0.80, "#bf{ALICE}");
            // labelXT.SetTextSize(0.045);
            // labelXT.DrawLatex(0.56, 0.74, "pp, |y| < 0.5");

            labelXT.SetTextSize(0.065);
            labelXT.DrawLatex(0.74, 0.88, "#phi");
        }
    }
} // namespace PhiXT

namespace CrossSectionPlot
{
    constexpr double sigmaINEL = 77.904; // mb, same normalization as xT macro
    constexpr double sqrtS = 13600.0;

    TF1 *MakeLevyTsallis(TH1 *hStat, TH1 *hSys,
                         double normalization,
                         const char *name,
                         double mass, int color,
                         bool isPhi = false)
    {
        if (!hStat || !hSys)
            return nullptr;

        const double fitMin = isPhi ? 0.4 : 0.0;
        const double fitMax = isPhi ? 15.0 : 20.0;

        // Work with a separate cross-section histogram.
        auto *hFit = dynamic_cast<TH1 *>(hStat->Clone(
            Form("hFit_%s", name)));

        hFit->SetDirectory(nullptr);
        hFit->Scale(normalization);

        int validBins = 0;

        for (int bin = 1; bin <= hFit->GetNbinsX(); ++bin)
        {
            const double pt = hFit->GetBinCenter(bin);

            const double statError = hStat->GetBinError(bin) * normalization * (1.0 / (2.0 * TMath::Pi() * pt));

            // Use the systematic uncertainty convention already
            // used to construct the systematic graph.
            const double sysError = hSys->GetBinContent(bin) * hStat->GetBinContent(bin) * normalization * (1.0 / (2.0 * TMath::Pi() * pt));

            const double totalError = std::sqrt(statError * statError + sysError * sysError);

            hFit->SetBinError(bin, totalError);

            if (pt >= fitMin && pt <= fitMax &&
                hFit->GetBinContent(bin) > 0.0 &&
                totalError > 0.0)
            {
                ++validBins;
            }
        }

        if (validBins < 4)
        {
            std::cerr << name << ": Not enough valid fit bins: "
                      << validBins << std::endl;
            delete hFit;
            return nullptr;
        }

        TF1 *fit = new TF1(name, FuncLevy_Sawan, fitMin, fitMax, 4);

        if (!isPhi)
        {
            fit->SetParameter(0, 7.0 * normalization);
            fit->SetParameter(1, 0.4);
        }
        else
        {
            fit->SetParameter(0, 8.0 * normalization);
            fit->SetParameter(1, 0.5);
        }

        fit->FixParameter(2, mass);
        fit->SetParameter(3, 0.25);

        fit->SetLineColor(color);
        fit->SetLineWidth(2);
        fit->SetLineStyle(7);

        std::cout << "\nFitting " << name
                  << "\nFit range: " << fitMin << " - " << fitMax
                  << "\nValid bins: " << validBins
                  << "\nNormalization: " << normalization
                  << std::endl;

        TFitResultPtr fitResult =
            hFit->Fit(fit, "RI0+S", "", fitMin, fitMax);

        if (static_cast<int>(fitResult) != 0)
        {
            std::cerr << "Fit failed for " << name
                      << ", status = "
                      << static_cast<int>(fitResult)
                      << std::endl;

            delete fit;
            delete hFit;
            return nullptr;
        }

        std::cout << "Fit successful: " << name
                  << "\nChi2/NDF = "
                  << fit->GetChisquare() << " / "
                  << fit->GetNDF()
                  << "\nT = " << fit->GetParameter(1)
                  << "\nn = " << fit->GetParameter(3)
                  << std::endl;

        // Keep hFit alive because it owns the fitted histogram.
        // ROOT retains the fitted function independently as well.
        return fit;
    }

    struct Spectrum
    {
        TGraphErrors *stat = nullptr;
        TGraphAsymmErrors *sys = nullptr;
        TGraphErrors *ratioStat = nullptr;
        TGraphAsymmErrors *ratioSys = nullptr;
        TF1 *fit = nullptr;
    };

    // Convert d2N/(dpT dy) per INEL event to d2sigma/(dpT dy).
    Spectrum MakeSpectrum(TH1 *hStat, TH1 *hSys, double mass, double normalization, const char *name)
    {
        Spectrum result;

        if (!hStat || !hSys)
            return result;

        result.stat = new TGraphErrors();
        result.sys = new TGraphAsymmErrors();

        result.stat->SetName(Form("g%sCrossSectionStat", name));
        result.sys->SetName(Form("g%sCrossSectionSys", name));

        for (int bin = 1; bin <= hStat->GetNbinsX(); ++bin)
        {
            const double pt = hStat->GetBinCenter(bin);
            const double width = hStat->GetBinWidth(bin);
            const double yield = hStat->GetBinContent(bin);

            if (pt <= 0.0 || yield <= 0.0)
                continue;

            const double crossSection = yield * normalization * (1.0 / (2.0 * TMath::Pi() * pt));
            const double statError = hStat->GetBinError(bin) * normalization * (1.0 / (2.0 * TMath::Pi() * pt));

            // hSys contains relative systematic uncertainties.
            const double sysError = hSys->GetBinContent(bin) * crossSection;
            const int p = result.stat->GetN();

            result.stat->SetPoint(p, pt, crossSection);
            result.stat->SetPointError(p, 0.0, statError);

            result.sys->SetPoint(p, pt, crossSection);
            result.sys->SetPointError(p, width / 2.0, width / 2.0, sysError, sysError);
        }

        if (result.stat->GetN() < 3)
            return result;

        double xmin, ymax, xmax, ymin;
        xmin = result.stat->GetX()[0];
        xmax = result.stat->GetX()[result.stat->GetN() - 1];

        const bool isPhi = TString(name) == "Phi";

        result.fit = MakeLevyTsallis(hStat, hSys, normalization, Form("levy%s", name), mass, kBlack, isPhi);

        result.ratioStat = new TGraphErrors();
        result.ratioSys = new TGraphAsymmErrors();

        for (int i = 0; i < result.stat->GetN(); ++i)
        {
            double pt, value;
            result.stat->GetPoint(i, pt, value);

            const double model = result.fit->Eval(pt);
            if (model <= 0.0)
                continue;

            const double statError = result.stat->GetErrorY(i);
            const double sysError =
                result.sys->GetErrorYhigh(i);

            const int p = result.ratioStat->GetN();

            result.ratioStat->SetPoint(p, pt, value / model);
            result.ratioStat->SetPointError(
                p, 0.0, statError / model);

            result.ratioSys->SetPoint(p, pt, value / model);
            result.ratioSys->SetPointError(
                p,
                result.sys->GetErrorXlow(i),
                result.sys->GetErrorXhigh(i),
                sysError / model,
                sysError / model);
        }

        return result;
    }

    void SetStyle(Spectrum &s, int color, int marker)
    {
        s.stat->SetMarkerStyle(marker);
        s.stat->SetMarkerColor(color);
        s.stat->SetLineColor(color);
        s.stat->SetMarkerSize(1.0);
        s.stat->SetFillStyle(0);

        s.sys->SetMarkerColor(color);
        s.sys->SetLineColor(color);
        s.sys->SetFillStyle(0);
        s.sys->SetFillColorAlpha(color, 0.20);

        s.ratioStat->SetMarkerStyle(marker);
        s.ratioStat->SetMarkerColor(color);
        s.ratioStat->SetLineColor(color);
        s.ratioStat->SetMarkerSize(0.8);

        s.ratioSys->SetLineColor(color);
        s.ratioSys->SetFillColorAlpha(color, 0.20);
        s.ratioSys->SetFillStyle(0);

        s.fit->SetLineColor(color);
        s.fit->SetLineWidth(2);
    }

    void DrawInvariant(const Spectrum &kstar, const Spectrum &phi)
    {
        // Convert differential cross-section graphs to invariant cross sections.
        auto ConvertStat = [](TGraphErrors *input, const char *name)
        {
            auto *output = new TGraphErrors();
            output->SetName(name);

            for (int i = 0; i < input->GetN(); ++i)
            {
                double pt, value;
                input->GetPoint(i, pt, value);

                if (pt <= 0.0)
                    continue;

                const double factor = 1.0 / (2.0 * TMath::Pi() * pt);

                const int p = output->GetN();
                output->SetPoint(p, pt, value * factor);
                output->SetPointError(
                    p,
                    input->GetErrorX(i),
                    input->GetErrorY(i) * factor);
            }

            return output;
        };

        auto ConvertSys = [](TGraphAsymmErrors *input, const char *name)
        {
            auto *output = new TGraphAsymmErrors();
            output->SetName(name);

            for (int i = 0; i < input->GetN(); ++i)
            {
                double pt, value;
                input->GetPoint(i, pt, value);

                if (pt <= 0.0)
                    continue;

                const double factor = 1.0 / (2.0 * TMath::Pi() * pt);

                const int p = output->GetN();
                output->SetPoint(p, pt, value * factor);

                output->SetPointError(
                    p,
                    input->GetErrorXlow(i),
                    input->GetErrorXhigh(i),
                    input->GetErrorYlow(i) * factor,
                    input->GetErrorYhigh(i) * factor);
            }

            return output;
        };

        // Transform the fitted differential curves as well.
        // Use sampled TGraphs to preserve the original fit functions.
        auto ConvertFit = [](TF1 *input, const char *name,
                             double ptMin, double ptMax, int color)
        {
            auto *output = new TGraph();
            output->SetName(name);

            constexpr int nPoints = 1000;

            for (int i = 0; i < nPoints; ++i)
            {
                const double pt = ptMin +
                                  (ptMax - ptMin) * i / (nPoints - 1.0);

                if (pt <= 0.0)
                    continue;

                const double value =
                    input->Eval(pt) / (2.0 * TMath::Pi() * pt);

                if (value <= 0.0)
                    continue;

                output->SetPoint(output->GetN(), pt, value);
            }

            output->SetLineColor(color);
            output->SetLineWidth(2);
            output->SetLineStyle(7);

            return output;
        };

        // ------------------------------------------------------------
        // Convert data
        // ------------------------------------------------------------

        auto *kstarStat = ConvertStat(
            kstar.stat, "gKstarInvariantStat");

        auto *kstarSys = ConvertSys(
            kstar.sys, "gKstarInvariantSys");

        auto *phiStat = ConvertStat(
            phi.stat, "gPhiInvariantStat");

        auto *phiSys = ConvertSys(
            phi.sys, "gPhiInvariantSys");

        // ------------------------------------------------------------
        // Convert Levy-Tsallis fits
        // ------------------------------------------------------------

        auto *kstarFit = ConvertFit(
            kstar.fit, "gKstarInvariantFit",
            0.05, 30.0, kOrange + 7);

        auto *phiFit = ConvertFit(
            phi.fit, "gPhiInvariantFit",
            0.05, 20.0, kCyan + 2);

        // ------------------------------------------------------------
        // Styles
        // ------------------------------------------------------------

        SetGraphStyle(kstarStat, kOrange + 7, 1.0, 20, 0.8, 1);
        SetGraphStyle(phiStat, kCyan + 2, 1.0, 20, 0.8, 1);

        kstarSys->SetLineColor(kOrange + 7);
        kstarSys->SetFillColorAlpha(kOrange + 7, 0.20);
        kstarSys->SetFillStyle(0);

        phiSys->SetLineColor(kCyan + 2);
        phiSys->SetFillColorAlpha(kCyan + 2, 0.20);
        phiSys->SetFillStyle(0);

        // ------------------------------------------------------------
        // Canvas
        // ------------------------------------------------------------

        auto *c = new TCanvas(
            "cKstarPhiInvariantCrossSection",
            "Kstar and Phi invariant cross sections",
            700, 700);

        SetPad(c);

        c->SetBottomMargin(0.15);
        c->SetLeftMargin(0.21);
        c->SetRightMargin(0.03);
        c->SetLogy();

        auto *frame = c->DrawFrame(
            0.0, 1.e-8, 30.0, 15.0);

        SetHistStyle(frame);

        frame->GetXaxis()->SetTitle(
            "#it{p}_{T} (GeV/#it{c})");

        frame->GetYaxis()->SetTitle(
            "#it{E} #frac{d^{3}#sigma}{d^{3}#it{p}} "
            "(mb GeV^{ -2}#it{c}^{ 3})");

        frame->GetYaxis()->SetTitleOffset(1.7);

        // Systematic uncertainties
        kstarSys->Draw("2 SAME");
        phiSys->Draw("2 SAME");

        // Levy-Tsallis curves
        kstarFit->Draw("L SAME");
        phiFit->Draw("L SAME");

        // Statistical uncertainties and points
        kstarStat->Draw("P E SAME");
        phiStat->Draw("P E SAME");

        TF1 *fitLegend = new TF1("fitLegend", "1.0", 0.0, 1.0);
        fitLegend->SetLineColor(kBlack);
        fitLegend->SetLineWidth(2);
        fitLegend->SetLineStyle(7);

        // ------------------------------------------------------------
        // Legend
        // ------------------------------------------------------------

        auto *legend = new TLegend(0.55, 0.55, 0.85, 0.75);
        SetLegend(legend);
        legend->SetTextSize(0.045);
        legend->AddEntry(kstarStat, "K*^{0} + #bar{K*}^{0}", "p");
        legend->AddEntry(phiStat, "#phi", "p");
        legend->AddEntry(fitLegend, "Levy-Tsallis fit", "l");
        legend->Draw();

        // ------------------------------------------------------------
        // ALICE labels
        // ------------------------------------------------------------

        TLatex label;
        label.SetNDC();
        label.SetTextFont(42);

        label.SetTextSize(0.055);
        label.DrawLatex(0.45, 0.87, "#bf{ALICE}");

        label.SetTextSize(0.045);
        label.DrawLatex(0.45, 0.81,"pp, #sqrt{s} = 13.6 TeV, |y| < 0.5");

        label.SetTextSize(0.038);
        label.DrawLatex(0.25, 0.19,"Global norm. unc. 2.6% (not shown)");

        c->SaveAs("../figures/KstarPhi_InvariantCrossSection_13p6TeV.pdf");
    }

    void Draw()
    {
        const char *root = "/Users/veronikabarbasova/git/Phi-analysis/";

        // ------------------------------------------------------------
        // K*0: Run 3
        // ------------------------------------------------------------

        TFile *fKstar = TFile::Open("../inputs/Kstar0_136/INEL_spectra.root", "READ");

        if (!fKstar || fKstar->IsZombie())
        {
            std::cerr << "CrossSection: Cannot open Kstar input file."
                      << std::endl;
            return;
        }

        auto *hKstat = dynamic_cast<TH1 *>(
            fKstar->Get("hINEL_DataCorr_stat"));
        auto *hKsysAbsolute = dynamic_cast<TH1 *>(
            fKstar->Get("hINEL_DataCorr_sys"));
        auto *hKsys = KstarXT::MakeRelativeSystematic(
            hKstat, hKsysAbsolute, "hKstarCrossSectionRelativeSys");

        // ------------------------------------------------------------
        // Phi: Run 3
        // ------------------------------------------------------------

        TFile *fPhi = TFile::Open(Form("%soutputs/CBT_hadronPID/inel/PAR.root", root), "READ");

        if (!fPhi || fPhi->IsZombie())
        {
            std::cerr << "CrossSection: Cannot open Phi file."
                      << std::endl;
            return;
        }

        auto *sparse = dynamic_cast<THnSparseD *>(
            fPhi->Get("sResultsMB"));

        if (!sparse || !hKstat || !hKsys)
        {
            std::cerr << "CrossSection: Missing input histograms."
                      << std::endl;
            return;
        }

        const int statBin = sparse->GetAxis(0)->FindBin("Yield_Corrected_Stat");
        const int sysBin = sparse->GetAxis(0)->FindBin("Yield_Corrected_Sys_Mean");

        sparse->GetAxis(0)->SetRange(statBin, statBin);
        auto *hPhiStat = dynamic_cast<TH1D *>(
            sparse->Projection(1));

        sparse->GetAxis(0)->SetRange(sysBin, sysBin);
        auto *hPhiSys = dynamic_cast<TH1D *>(
            sparse->Projection(1));

        sparse->GetAxis(0)->SetRange(0, -1);

        if (!hPhiStat || !hPhiSys)
        {
            std::cerr << "CrossSection: Phi projections failed."
                      << std::endl;
            return;
        }

        // ------------------------------------------------------------
        // Convert to differential cross sections
        // ------------------------------------------------------------

        constexpr double kstarFactor = 1.0; // K*0 + K*0bar

        auto kstar = MakeSpectrum(hKstat, hKsys, 0.89555, sigmaINEL * kstarFactor, "Kstar");

        // For Phi, the systematic histogram stores absolute
        // uncertainties, so convert it to relative uncertainties first.
        auto *hPhiRelSys = dynamic_cast<TH1D *>(hPhiSys->Clone("hPhiRelativeSys"));

        for (int bin = 1; bin <= hPhiRelSys->GetNbinsX(); ++bin)
        {
            const double y = hPhiStat->GetBinContent(bin);
            const double error = hPhiSys->GetBinError(bin);

            hPhiRelSys->SetBinContent(bin, y > 0.0 ? error / y : 0.0);
        }

        auto phi = MakeSpectrum(hPhiStat, hPhiRelSys, 1.01946, sigmaINEL, "Phi");

        if (!kstar.stat || !phi.stat ||
            !kstar.fit || !phi.fit)
        {
            std::cerr << "CrossSection: Spectrum construction failed."
                      << std::endl;
            return;
        }

        // ------------------------------------------------------------
        // Styles
        // ------------------------------------------------------------

        SetStyle(phi, kCyan + 2, 20);
        SetStyle(kstar, kOrange + 7, 20);

        // ------------------------------------------------------------
        // Canvas and pads
        // ------------------------------------------------------------

        auto *c = new TCanvas("cKstarPhiCrossSection", "Kstar and Phi cross sections", 700, 700);
        SetPad(c);
        c->SetBottomMargin(0.15);
        c->SetLeftMargin(0.21);
        c->SetRightMargin(0.03);
        c->SetLogy();

        auto *frame = c->DrawFrame(0.0, 0.15 * 1.e-5, 30.0, 15.0);

        SetHistStyle(frame);
        frame->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        frame->GetYaxis()->SetTitle("d^{2}#sigma/(d#it{p}_{T} d#it{y}) (mb /GeV/#it{c})");
        // frame->GetYaxis()->SetTitle("#it{E} #frac{d^{3}#sigma}{d^{3}#it{p}} (mb GeV^{ -2}#it{c}^{ 3})");
        frame->GetYaxis()->SetTitleOffset(1.7);

        frame->GetXaxis()->SetRangeUser(0.0, 30.0);

        kstar.fit->SetLineColor(kOrange + 7);
        phi.fit->SetLineColor(kCyan + 2);

        TF1 *fitClone = dynamic_cast<TF1 *>(kstar.fit->Clone("fitClone"));
        fitClone->SetLineColor(kBlack);

        kstar.sys->Draw("2 SAME");
        phi.sys->Draw("2 SAME");

        kstar.fit->SetRange(0.0, 30.0);
        phi.fit->SetRange(0.0, 20.0);
        kstar.fit->Draw("SAME");
        phi.fit->Draw("SAME");

        kstar.stat->Draw("P E SAME");
        phi.stat->Draw("P E SAME");

        auto *legend = new TLegend(0.55, 0.55, 0.85, 0.75);
        SetLegend(legend);
        legend->SetTextSize(0.045);

        legend->AddEntry(kstar.stat, "K*^{0} + #bar{K*}^{0}", "p");
        legend->AddEntry(phi.stat, "#phi", "p");
        legend->AddEntry(fitClone, "Levy-Tsallis fit", "l");

        legend->Draw();

        TLatex label;
        label.SetNDC();
        label.SetTextFont(42);
        label.SetTextSize(0.055);
        label.DrawLatex(0.45, 0.87, "#bf{ALICE}");
        label.SetTextSize(0.045);
        label.DrawLatex(0.45, 0.81, "pp, #sqrt{s} = 13.6 TeV, |y| < 0.5");

        label.SetTextSize(0.038);
        label.DrawLatex(0.25, 0.19, "Global norm. unc. 2.6% (not shown)");

        c->SaveAs("../figures/KstarPhi_CrossSection_13p6TeV.pdf");

        DrawInvariant(kstar, phi);
    }

}

void xTScaling()
{
    TGaxis::SetMaxDigits(3);
    TH1::AddDirectory(kFALSE);

    auto *cN = new TCanvas("cKstarPhiN", "xT exponent", 850, 450);
    auto *padKN = new TPad("padKstarN", "Kstar", 0.0, 0.0, 0.5, 1.0);
    padKN->SetTickx(1);
    padKN->SetTicky(1);
    padKN->SetLeftMargin(0.18);
    padKN->SetRightMargin(0.0);
    padKN->SetBottomMargin(0.15);
    padKN->SetTopMargin(0.05);
    padKN->Draw();

    auto *padPN = new TPad("padPhiN", "Phi", 0.5, 0.0, 1.0, 1.0);
    padPN->SetTickx(1);
    padPN->SetTicky(1);
    padPN->SetLeftMargin(0.00);
    padPN->SetRightMargin(0.18);
    padPN->SetBottomMargin(0.15);
    padPN->SetTopMargin(0.05);
    padPN->Draw();

    auto *cXT = new TCanvas("cKstarPhiXT", "xT distributions", 850, 450);
    auto *padKXT = new TPad("padKstarXT", "Kstar", 0.0, 0.0, 0.5, 1.0);
    padKXT->SetTickx(1);
    padKXT->SetTicky(1);
    padKXT->SetLeftMargin(0.18);
    padKXT->SetRightMargin(0.0);
    padKXT->SetBottomMargin(0.15);
    padKXT->SetTopMargin(0.05);
    padKXT->Draw();

    auto *padPXT = new TPad("padPhiXT", "Phi", 0.5, 0.0, 1.0, 1.0);
    padPXT->SetTickx(1);
    padPXT->SetTicky(1);
    padPXT->SetLeftMargin(0.00);
    padPXT->SetRightMargin(0.18);
    padPXT->SetBottomMargin(0.15);
    padPXT->SetTopMargin(0.05);
    padPXT->Draw();

    // Each species is calculated once; its n(xT) and scaled distribution are
    // drawn into the corresponding pads on the two canvases.
    KstarXT::Draw(padKN, padKXT);
    PhiXT::Draw(padPN, padPXT);
    CrossSectionPlot::Draw();

    cN->Modified();
    cN->Update();
    cXT->Modified();
    cXT->Update();
    cN->SaveAs("../figures/KstarPhi_n_xT.pdf");
    cXT->SaveAs("../figures/KstarPhi_xT_scaling.pdf");
}
