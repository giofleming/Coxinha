#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TString.h>
#include <TTree.h>

#include <cmath>
#include <iostream>

void plotLutRoot(const char* inputFile,
                 const char* outputPrefix = "lut_root_qa",
                 float etaMin = -1, float etaMax = 1,
                 int nchBin = 0,
                 int radBin = 0)
{
  TFile* fin = new TFile(inputFile, "READ");
  if (fin->IsZombie()) {
    std::cerr << "Cannot open input ROOT file: " << inputFile << std::endl;
    return;
  }

  auto* header = dynamic_cast<TTree*>(fin->Get("lutHeader"));
  auto* entries = dynamic_cast<TTree*>(fin->Get("lutEntries"));
  if (header == nullptr || entries == nullptr) {
    std::cerr << "Missing lutHeader or lutEntries tree in file: " << inputFile << std::endl;
    return;
  }

  int pdg = 0;
  float field = 0.f;
  int pt_nbins = 0;
  float pt_min = 0.f;
  float pt_max = 0.f;
  bool pt_log = false;
  bool covarianceStoresSigma = false;

  header->SetBranchAddress("pdg", &pdg);
  header->SetBranchAddress("field", &field);
  header->SetBranchAddress("pt_nbins", &pt_nbins);
  header->SetBranchAddress("pt_min", &pt_min);
  header->SetBranchAddress("pt_max", &pt_max);
  header->SetBranchAddress("pt_log", &pt_log);
  header->SetBranchAddress("covarianceStoresSigma", &covarianceStoresSigma);
  header->GetEntry(0);

  std::cout << "LUT header: pdg=" << pdg
            << " field=" << field << " T"
            << " pt bins=" << pt_nbins
            << " mode=" << (covarianceStoresSigma ? "sigma" : "mean")
            << std::endl;

  auto* hEffVsPt = new TH2F("hEffVsPt", "Efficiency vs p_{T};p_{T} [GeV/c];eff", 200, 0., 100., 100, 0., 1.2);
  hEffVsPt->SetDirectory(0);
  auto* hCov44VsPt = new TH2F("hCov44VsPt", "covm[14] vs p_{T};p_{T} [GeV/c];covm[14]", 200, 0., 100., 200, 0., 0.1);
  hCov44VsPt->SetDirectory(0);
  auto* hCov00 = new TH1F("hCov00", "Distribution of covm[0];covm[0];entries", 200, 0., 0.02);
  hCov00->SetDirectory(0);
  auto* hEtaVsPt = new TH2F("hEtaVsPt", "#eta vs p_{T};p_{T} [GeV/c];#eta", 200, 0., 100., 200, -4., 4.);
  hEtaVsPt->SetDirectory(0);
  auto* hPtResVsPt = new TH2F("hPtResVsPt", "p_{T} resolution vs p_{T};p_{T} [GeV/c];#sigma(p_{T})/p_{T}", 200, 0., 100., 200, 0., 0.5);
  hPtResVsPt->SetDirectory(0);
  auto* pAvgPtResVsPt = new TProfile("pAvgPtResVsPt", "Average p_{T} resolution vs p_{T};p_{T} [GeV/c];<#sigma(p_{T})/p_{T}>", 200, 0., 100.);
  pAvgPtResVsPt->SetDirectory(0);
  auto* pAvgEffVsPt = new TProfile("pAvgEffVsPt", "Average efficiency vs p_{T};p_{T} [GeV/c];<eff>", 200, 0., 100.);
  pAvgEffVsPt->SetDirectory(0);

  int entryNchBin = -1;
  int entryRadBin = -1;
  int entryEtaBin = -1;
  float entryPt = 0.f;
  float entryEta = 0.f;
  float entryEff = 0.f;
  bool entryValid = false;
  float entryCovm[15] = {0.f};

  entries->SetBranchAddress("nchBin", &entryNchBin);
  entries->SetBranchAddress("radBin", &entryRadBin);
  entries->SetBranchAddress("etaBin", &entryEtaBin);
  entries->SetBranchAddress("pt", &entryPt);
  entries->SetBranchAddress("eta", &entryEta);
  entries->SetBranchAddress("eff", &entryEff);
  entries->SetBranchAddress("valid", &entryValid);
  entries->SetBranchAddress("covm", entryCovm);

  const Long64_t nEntries = entries->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    entries->GetEntry(i);

    if (!entryValid) {
      continue;
    }
    if (entryNchBin != nchBin || entryRadBin != radBin) {
      continue;
    }
    if (etaMin < etaMax && (entryEta < etaMin || entryEta > etaMax)) {
      continue;
    }

    hEffVsPt->Fill(entryPt, entryEff);
    pAvgEffVsPt->Fill(entryPt, entryEff);
    hCov44VsPt->Fill(entryPt, entryCovm[14]);
    hCov00->Fill(entryCovm[0]);
    hEtaVsPt->Fill(entryPt, entryEta);

    if (entryPt > 0.f) {
      float sigmaQOverPt = 0.f;
      if (covarianceStoresSigma) {
        sigmaQOverPt = std::abs(entryCovm[14]);
      } else if (entryCovm[14] > 0.f) {
        sigmaQOverPt = std::sqrt(entryCovm[14]);
      }

      if (sigmaQOverPt > 0.f) {
        const float relPtResolution = entryPt * sigmaQOverPt;
        hPtResVsPt->Fill(entryPt, relPtResolution);
        pAvgPtResVsPt->Fill(entryPt, relPtResolution);
      }
    }
  }

  auto* c1 = new TCanvas("c1", "LUT QA", 1500, 900);
  c1->Divide(3, 2);
  c1->cd(1);
  hEffVsPt->Draw("COLZ");
  c1->cd(2);
  pAvgEffVsPt->Draw();
  c1->cd(3);
  hCov44VsPt->Draw("COLZ");
  c1->cd(4);
  hCov00->Draw();

  c1->cd(5);
  hEtaVsPt->Draw("COLZ");
  c1->cd(6);
  hPtResVsPt->Draw("COLZ");

  auto* c2 = new TCanvas("c2", "LUT 1D averages", 900, 700);
  c2->Divide(1, 2);
  c2->cd(1);
  pAvgEffVsPt->Draw();
  c2->cd(2);
  pAvgPtResVsPt->Draw();

  if (0) {
    TString suffix = TString::Format("pdg%d", pdg);
    if (etaMin < etaMax) {
      suffix += TString::Format("_eta%.1f_%+.1f", etaMin, etaMax);
    }

    const TString outPng = TString(outputPrefix) + "_" + suffix + ".png";
    c1->SaveAs(outPng);

    auto* fout = TFile::Open((TString(outputPrefix) + "_" + suffix + ".root").Data(), "RECREATE");
    if (fout != nullptr && !fout->IsZombie()) {
      hEffVsPt->Write();
      hCov44VsPt->Write();
      hCov00->Write();
      hEtaVsPt->Write();
      hPtResVsPt->Write();
      pAvgPtResVsPt->Write();
      pAvgEffVsPt->Write();
      c1->Write("canvas");
      c2->Write("averages_canvas");
      fout->Close();
    }

    std::cout << "Wrote QA outputs with prefix: " << outputPrefix << "_" << suffix << std::endl;
  }
}
