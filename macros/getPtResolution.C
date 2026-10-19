// Author: Pavel Larionov
// Email:  pavel.larionov@cern.ch
// Event:  ACTS Tutorial, June '25

#include "TreeReader.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TMath.h>
#include "TError.h"
#include "TSystem.h"
#include <iostream>
#include <string>

void FillQopResidualVsEta(TTree *inTree, TH2F *h_ptres, TFile *out);
TH1F *GetResolutionHisto(TH2F *inHist, const float &fitRangeSigma, const std::string &paramStr, const std::string &axisStr, const std::string &plotsDir, bool useCustomBinning = false);
void plotHist(TH1F *input, const std::string &canvTitle, const std::string &yAxisTitle, const std::string &xAxisLabel, const std::string &plotsDir, bool xlog, bool ylog, float yAxisMin, float yAxisMax, float rngmin, float rngmax);

std::string plotFileName(const std::string &base, const std::string &suffix = "")
{
    return "Plots/" + base + "_" + suffix + ".pdf";
}

int getPtResolution(std::string strFolder = "reco_output_pythia")
{

    gErrorIgnoreLevel = kWarning;

    /*
    Open the tree with track summary information
    */

    std::string targetDir = "../" + strFolder;

    TFile *inFile = TFile::Open(Form("%s/tracksummary_ambi.root", targetDir.c_str()));

    if (inFile == nullptr)
    {
        std::cerr << "Error: Could not open input file." << std::endl;
        return 0;
    }

    TTree *tree = (TTree *)inFile->Get("tracksummary");

    if (tree == nullptr)
    {
        std::cerr << "Error: Could not retrieve TTree from input file." << std::endl;
        inFile->Close();
        return 0;
    }

    const std::string treeOutputDir = targetDir + "/treeoutput";
    const std::string plotsOutputDir = targetDir + "/Plots";

    // Check and create the directory if it doesn't exist
    if (gSystem->AccessPathName(treeOutputDir.c_str()))
    {
        gSystem->mkdir(treeOutputDir.c_str(), kTRUE); // kTRUE for recursive creation
    }

    if (gSystem->AccessPathName(plotsOutputDir.c_str()))
    {
        gSystem->mkdir(plotsOutputDir.c_str(), kTRUE); // kTRUE for recursive creation
    }

    TFile *outFile = TFile::Open(Form("%s/OutputHistos.root", treeOutputDir.c_str()), "RECREATE");

    /*
    TH2F histograms with the residual.
    Here we set the desired binning for eta and pT.
    The eta and pT ranges are accessed later by checking the first/last bin edges.
    */

    TH2F *h_qop_residual_vs_eta = new TH2F("h_qop_residual_vs_eta", "normalised qop residual vs eta", 80, -4., 4., 60, -0.3, 0.3);

    /*
    Fill the residuals for qop.
    */

    FillQopResidualVsEta(tree, h_qop_residual_vs_eta, outFile);

    /*
    Get the qop resolution histograms.
    The histograms are filled with the resolution of the qop parameter.
    The resolution is defined as the difference between the generated and reconstructed qop divided by the generated qop.
    */

    auto h_pT_resolution_vs_eta = GetResolutionHisto(h_qop_residual_vs_eta, 2.5, "qop", "eta", plotsOutputDir);

    plotHist(h_pT_resolution_vs_eta, "pT_resolution_vs_eta", "#Delta#it{p}_{T}/#it{p}_{T}", "#it{#eta}", plotsOutputDir, false, false, 0., 0.2, 0., 4.);
    outFile->WriteTObject(h_pT_resolution_vs_eta, h_pT_resolution_vs_eta->GetName());

    /*
    Close files, delete pointers
    */

    inFile->Close();
    outFile->Close();

    return 1;
}

void FillQopResidualVsEta(TTree *inTree, TH2F *h_residual, TFile *out)
{
    //
    // Here we process the tracksummary tree
    //

    TrackSummaryReader tsReader(inTree, false);
    h_residual->Sumw2();

    int nbins = h_residual->GetNbinsX();

    for (unsigned iEntry = 0; iEntry < inTree->GetEntries(); iEntry++)
    {
        tsReader.getEntry(iEntry);

        // Select only fitted tracks
        if (not tsReader.hasFittedParams)
            continue;

        // Loop over tracks in the event
        for (unsigned ip = 0; ip < tsReader.eQOP_fit->size(); ip++)
        {
            if (tsReader.t_eta->at(ip) < h_residual->GetXaxis()->GetBinLowEdge(1) || tsReader.t_eta->at(ip) > h_residual->GetXaxis()->GetBinUpEdge(nbins))
                continue;

            // Fill 2D histogram with normalized QOP residual vs eta
            auto momGen = tsReader.t_p->at(ip) * tsReader.t_charge->at(ip);
            auto momReco = 1 / tsReader.eQOP_fit->at(ip);
            h_residual->Fill(tsReader.t_eta->at(ip), (momGen - momReco) / momGen);
        }
    }

    out->cd();
    out->WriteTObject(h_residual, h_residual->GetName());
}

TH1F *GetResolutionHisto(TH2F *inHist, const float &fitRangeSigma, const std::string &paramStr, const std::string &axisStr, const std::string &plotsDir, bool useCustomBinning)
{
    //
    // This method gets the qop resolution
    //

    if (inHist == nullptr)
    {
        std::cerr << "[GetResolutionHisto] Error: empty histogram!" << std::endl;
        return nullptr;
    }

    /// Output histogram
    auto nbins = inHist->GetNbinsX();
    auto lowedge = inHist->GetXaxis()->GetBinLowEdge(1);
    auto upedge = inHist->GetXaxis()->GetBinUpEdge(inHist->GetNbinsX());
    TH1F *outHist = nullptr;
    if (useCustomBinning)
    {
        // pT binning
        const int nBinsPt = 35;
        const double xBinsPt[nBinsPt + 1] = {
            0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.125, 0.15, 0.175, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
            0.8, 0.9, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};

        outHist = new TH1F((paramStr + "_resolution_vs_" + axisStr).c_str(), (paramStr + " resolution vs " + axisStr).c_str(), nBinsPt, xBinsPt);
    }
    else
        outHist = new TH1F((paramStr + "_resolution_vs_" + axisStr).c_str(), (paramStr + " resolution vs " + axisStr).c_str(), nbins, lowedge, upedge);
    outHist->Sumw2();

    /// Gaussian to fit the signal
    TF1 *gfnc = new TF1("gfnc", "gaus", -1., 1.);
    gfnc->SetLineColor(kRed);
    gfnc->SetParameter(0, 100);
    gfnc->SetParameter(1, 0.);
    gfnc->SetParameter(2, 0.05);

    /// Fitting options
    std::string opt = "LI Q";

    unsigned cnt = 0;

    TCanvas *canvfit = new TCanvas("", "", 800, 600);

    gStyle->SetOptFit(1111);

    for (int ib = 1; ib < inHist->GetNbinsX() + 1; ib++)
    {
        canvfit->SetName(Form("canvfit_%s_%.1f", axisStr.c_str(), inHist->GetXaxis()->GetBinLowEdge(ib)));
        canvfit->SetTitle(Form("canvfit %s = %.1f", axisStr.c_str(), inHist->GetXaxis()->GetBinLowEdge(ib)));

        auto hProj = inHist->ProjectionY(Form("_py_bin%d", ib), ib, ib);
        hProj->SetTitle(Form("%s = %.1f", axisStr.c_str(), inHist->GetXaxis()->GetBinLowEdge(ib)));

        auto fitRange = hProj->GetRMS() * fitRangeSigma;

        if (hProj->GetEntries() < 100)
            continue;

        /// Legend
        TLegend *leg = new TLegend(0.12, 0.6, 0.45, 0.8);

        leg->SetBorderSize(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.03);
        leg->SetFillStyle(0);

        canvfit->cd();

        hProj->SetMarkerColor(kBlue - 4);
        hProj->SetLineColor(kBlue - 4);
        hProj->GetXaxis()->SetTitle(Form("(%s_(MC) - %s_(fit))", paramStr.c_str(), paramStr.c_str()));
        hProj->GetXaxis()->SetNdivisions(505);
        hProj->GetXaxis()->SetTitleOffset(0.85);
        hProj->GetXaxis()->SetTitleSize(0.05);
        hProj->GetXaxis()->SetLabelSize(0.05);
        hProj->GetYaxis()->SetLabelSize(0.05);
        hProj->DrawCopy();

        hProj->Fit(gfnc, Form("%s", opt.c_str()), "", -1 * fitRange, fitRange);

        auto mean = gfnc->GetParameter(1);
        auto sigma = gfnc->GetParameter(2);
        auto sigmaerr = gfnc->GetParError(2);

        outHist->SetBinContent(ib, sigma);
        outHist->SetBinError(ib, sigmaerr);

        canvfit->SaveAs(Form("%s/fit_%s_res_vs_%s.pdf[", plotsDir.c_str(), paramStr.c_str(), axisStr.c_str()));

        hProj->Clear();
        canvfit->Clear();
        cnt++;
        delete leg;
    }

    canvfit->SaveAs(Form("%s/fit_%s_res_vs_%s.pdf]", plotsDir.c_str(), paramStr.c_str(), axisStr.c_str()));

    TH1F *result = outHist;

    delete gfnc;
    delete canvfit;

    return result;
}

void plotHist(TH1F *input, const std::string &canvTitle, const std::string &yAxisTitle, const std::string &xAxisLabel, const std::string &plotsDir, bool xlog, bool ylog, float yAxisMin, float yAxisMax, float rngmin, float rngmax)
{
    // Plot the histogram

    if (input == nullptr)
    {
        std::cerr << "[plotHist] Error: empty histogram!" << std::endl;
        return;
    }

    TCanvas *canv = new TCanvas(Form("canv_%s", canvTitle.c_str()), Form("%s", canvTitle.c_str()), 800, 600);
    TLegend *leg = new TLegend(0.2, 0.7, 0.8, 0.8);
    leg->AddEntry((TObject *)0, "", "");

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    canv->SetLogx(xlog);
    canv->SetLogy(ylog);
    canv->SetTicks(1, 1);

    TH1F *dummy = new TH1F("dummy", input->GetTitle(), 1, rngmin, rngmax);
    dummy->SetMinimum(yAxisMin);
    dummy->SetMaximum(yAxisMax);
    dummy->GetYaxis()->SetTitle(yAxisTitle.c_str());
    dummy->GetXaxis()->SetTitle(xAxisLabel.c_str());
    dummy->DrawCopy();

    input->SetMarkerStyle(kFullCircle);
    input->SetMarkerSize(1.2);
    input->SetMarkerColor(9);
    input->SetLineColor(9);
    input->SetLineWidth(2);

    input->DrawCopy("Same");

    leg->AddEntry(input, input->GetTitle(), "lp");

    leg->DrawClone("Same");

    {
        canv->SaveAs(Form("%s/%s%s.pdf", plotsDir.c_str(), canv->GetName(), (ylog) ? "_log" : ""));
    }

    delete dummy;
    delete leg;
    delete canv;
}