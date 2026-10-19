// drawCumulativeEfficienciesVsEta.C
//
// Companion to drawIterativeEfficiencies.C: draws the CUMULATIVE (merged)
// seeding and tracking efficiencies as a function of pseudorapidity for the
// ACTS iterative tracking chain. Only the cumulative curves are drawn - the
// per-iteration breakdown is what drawIterativeEfficiencies.C shows vs pT.
//
// Everything is computed directly from the output TTrees, exactly as in
// drawIterativeEfficiencies.C:
//   - denominator : particles_digitized_selected.root   (tree "particles")
//   - seeding eff : union over all iterations of the particles that are the
//                   (pure) truth match of at least one seed
//                   (branch "truthMatched" of tree "estimatedparams")
//                   iter 0 in <dir>/iter_0/, iter N in <dir>/iter_N/
//   - tracking eff: tracksummary_ambi-tracks-merged.root (tree "tracksummary")
//                   numerator = particles with at least one track whose
//                   trackClassification == Matched (== 1)
//
// A particle is identified by (event number, barcode components), and its
// pT/eta are always taken from the particles tree, so numerator and
// denominator are guaranteed to be consistent.
//
// Outputs (written next to the input dir, name encodes the pT range):
//   efficiencies_vs_eta<range>_pt<...>.png / .pdf  - the drawn canvas
//     (plus _seedingOnly / _trackingOnly variants)
//   efficiencies_vs_eta<range>_pt<...>.root        - the two efficiency graphs
//                                              (gSeedEff_merged, gTrkEff_merged)
//                                              + the canvas
//
// Usage:
//   root -l -q 'drawCumulativeEfficienciesVsEta.C("/path/to/output_dir")'
//   root -l -q 'drawCumulativeEfficienciesVsEta.C("/path/to/out", -4, 4, 0.1, 100., 4, 40)'
//
// Claude, September 2026

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TEfficiency.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"

#include <array>
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <cmath>

// ---------------------------------------------------------------- settings
static const int kNEtaBinsDefault = 40;

// eta bin edges, filled at startup
static std::vector<double> gEtaEdges;

// annotation text on the plot (edit to match your setup); font 72 is
// helvetica-bold-oblique, i.e. the usual bold-italic experiment tag
static const char* kEtaInfoLine1 = "#font[72]{ALICE3-ACTS,} Pythia8+G4";
static const char* kEtaInfoLine2 = "pp #LT#mu#GT = 30, N^{truth}_{hits} #geq 7";

static const int kCumulColor = kBlack;

// marker / axis sizes
static const double kMarkerSize = 1.3;
static const double kAxisTitleSize = 0.05;
static const double kAxisLabelSize = 0.043;

// the three flavours of the plot: both, seeding only, tracking only
static const char* kModeSuffix[3] = {"", "_seedingOnly", "_trackingOnly"};

// ------------------------------------------------------------------ helpers

// unique particle key: {event, vertex_primary, vertex_secondary, particle,
//                       generation, sub_particle}
typedef std::array<UInt_t, 6> EtaPKey;

struct EtaParticleInfo {
  float pt;
  float eta;
};

TH1D* makeEtaHist(const char* name)
{
  return new TH1D(name, name, gEtaEdges.size() - 1, gEtaEdges.data());
}

// common canvas / axis styling
void styleCanvasEta(TCanvas* c)
{
  c->SetTicks(1, 1);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.04);
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.13);
}

void styleAxesEta(TH1D* h, const char* xTitle, const char* yTitle)
{
  h->SetTitle("");
  h->GetXaxis()->SetTitle(xTitle);
  h->GetXaxis()->SetTitleSize(kAxisTitleSize);
  h->GetXaxis()->SetLabelSize(kAxisLabelSize);
  h->GetXaxis()->SetTitleOffset(1.15);
  h->GetYaxis()->SetTitle(yTitle);
  h->GetYaxis()->SetTitleSize(kAxisTitleSize);
  h->GetYaxis()->SetLabelSize(kAxisLabelSize);
  h->GetYaxis()->SetTitleOffset(1.20);
}

// a coloured line/marker used only as a legend key
TGraph* legendKeyEta(int color, int markerStyle)
{
  TGraph* g = new TGraph();
  g->SetLineColor(color);
  g->SetLineWidth(2);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(markerStyle);
  g->SetMarkerSize(kMarkerSize);
  return g;
}

TFile* openFileEta(const TString& path)
{
  if (gSystem->AccessPathName(path)) { // true if NOT accessible
    std::cout << "WARNING: file not found: " << path << std::endl;
    return nullptr;
  }
  TFile* f = TFile::Open(path);
  if (!f || f->IsZombie()) {
    std::cout << "WARNING: cannot open: " << path << std::endl;
    return nullptr;
  }
  return f;
}

// ---- read the selected-particles tree -> map key->(pt,eta), fill hAll -----
// The eta window is SIGNED here (it is the plotted axis); the pT window is the
// selection applied to the denominator.
bool readParticlesEta(const TString& path, double etaMin, double etaMax,
                      double ptMin, double ptMax,
                      std::map<EtaPKey, EtaParticleInfo>& partMap, TH1D* hAll)
{
  TFile* f = openFileEta(path);
  if (!f)
    return false;
  TTree* t = (TTree*)f->Get("particles");
  if (!t) {
    std::cout << "ERROR: no 'particles' tree in " << path << std::endl;
    return false;
  }

  UInt_t eventId = 0;
  std::vector<float>* pt = nullptr;
  std::vector<float>* eta = nullptr;
  std::vector<UInt_t>* vp = nullptr;
  std::vector<UInt_t>* vs = nullptr;
  std::vector<UInt_t>* par = nullptr;
  std::vector<UInt_t>* gen = nullptr;
  std::vector<UInt_t>* sub = nullptr;

  t->SetBranchStatus("*", 0);
  for (const char* b : {"event_id", "pt", "eta", "vertex_primary",
                        "vertex_secondary", "particle", "generation",
                        "sub_particle"})
    t->SetBranchStatus(b, 1);

  t->SetBranchAddress("event_id", &eventId);
  t->SetBranchAddress("pt", &pt);
  t->SetBranchAddress("eta", &eta);
  t->SetBranchAddress("vertex_primary", &vp);
  t->SetBranchAddress("vertex_secondary", &vs);
  t->SetBranchAddress("particle", &par);
  t->SetBranchAddress("generation", &gen);
  t->SetBranchAddress("sub_particle", &sub);

  Long64_t nEnt = t->GetEntries();
  for (Long64_t i = 0; i < nEnt; ++i) {
    t->GetEntry(i);
    for (size_t j = 0; j < pt->size(); ++j) {
      if ((*eta)[j] < etaMin || (*eta)[j] > etaMax)
        continue;
      if ((*pt)[j] < ptMin || (*pt)[j] > ptMax)
        continue;
      EtaPKey key = {eventId, (*vp)[j], (*vs)[j], (*par)[j], (*gen)[j], (*sub)[j]};
      if (partMap.count(key))
        continue; // safety against duplicates
      partMap[key] = {(*pt)[j], (*eta)[j]};
      hAll->Fill((*eta)[j]);
    }
  }
  f->Close();
  std::cout << "  particles passing eta/pT selection: " << partMap.size()
            << std::endl;
  return true;
}

// ---- collect seed-matched particles from an estimatedparams tree ----------
void collectSeedMatchesEta(const TString& path,
                           const std::map<EtaPKey, EtaParticleInfo>& partMap,
                           std::set<EtaPKey>& matched)
{
  TFile* f = openFileEta(path);
  if (!f)
    return;
  TTree* t = (TTree*)f->Get("estimatedparams");
  if (!t) {
    std::cout << "ERROR: no 'estimatedparams' tree in " << path << std::endl;
    return;
  }

  Int_t eventNr = 0;
  Bool_t isMatched = false;
  UInt_t vp = 0, vs = 0, par = 0, gen = 0, sub = 0;

  t->SetBranchStatus("*", 0);
  for (const char* b : {"event_nr", "truthMatched",
                        "particleId_vertex_primary",
                        "particleId_vertex_secondary",
                        "particleId_particle", "particleId_generation",
                        "particleId_sub_particle"})
    t->SetBranchStatus(b, 1);

  t->SetBranchAddress("event_nr", &eventNr);
  t->SetBranchAddress("truthMatched", &isMatched);
  t->SetBranchAddress("particleId_vertex_primary", &vp);
  t->SetBranchAddress("particleId_vertex_secondary", &vs);
  t->SetBranchAddress("particleId_particle", &par);
  t->SetBranchAddress("particleId_generation", &gen);
  t->SetBranchAddress("particleId_sub_particle", &sub);

  Long64_t nEnt = t->GetEntries();
  for (Long64_t i = 0; i < nEnt; ++i) {
    t->GetEntry(i);
    if (!isMatched)
      continue;
    EtaPKey key = {(UInt_t)eventNr, vp, vs, par, gen, sub};
    if (!partMap.count(key))
      continue; // not in selected denominator
    matched.insert(key);
  }
  f->Close();
}

// ---- collect track-matched particles from a tracksummary tree -------------
void collectTrackMatchesEta(const TString& path,
                            const std::map<EtaPKey, EtaParticleInfo>& partMap,
                            std::set<EtaPKey>& matched)
{
  TFile* f = openFileEta(path);
  if (!f)
    return;
  TTree* t = (TTree*)f->Get("tracksummary");
  if (!t) {
    std::cout << "ERROR: no 'tracksummary' tree in " << path << std::endl;
    return;
  }

  UInt_t eventNr = 0;
  std::vector<int>* cls = nullptr;
  std::vector<UInt_t>* vp = nullptr;
  std::vector<UInt_t>* vs = nullptr;
  std::vector<UInt_t>* par = nullptr;
  std::vector<UInt_t>* gen = nullptr;
  std::vector<UInt_t>* sub = nullptr;

  t->SetBranchStatus("*", 0);
  for (const char* b : {"event_nr", "trackClassification",
                        "majorityParticleId_vertex_primary",
                        "majorityParticleId_vertex_secondary",
                        "majorityParticleId_particle",
                        "majorityParticleId_generation",
                        "majorityParticleId_sub_particle"})
    t->SetBranchStatus(b, 1);

  t->SetBranchAddress("event_nr", &eventNr);
  t->SetBranchAddress("trackClassification", &cls);
  t->SetBranchAddress("majorityParticleId_vertex_primary", &vp);
  t->SetBranchAddress("majorityParticleId_vertex_secondary", &vs);
  t->SetBranchAddress("majorityParticleId_particle", &par);
  t->SetBranchAddress("majorityParticleId_generation", &gen);
  t->SetBranchAddress("majorityParticleId_sub_particle", &sub);

  // TrackMatchClassification: 0=Unknown, 1=Matched, 2=Duplicate, 3=Fake
  const int kMatched = 1;

  Long64_t nEnt = t->GetEntries();
  for (Long64_t i = 0; i < nEnt; ++i) {
    t->GetEntry(i);
    for (size_t j = 0; j < cls->size(); ++j) {
      if ((*cls)[j] != kMatched)
        continue;
      EtaPKey key = {eventNr, (*vp)[j], (*vs)[j], (*par)[j], (*gen)[j], (*sub)[j]};
      if (!partMap.count(key))
        continue;
      matched.insert(key);
    }
  }
  f->Close();
}

// ---- matched set -> numerator histogram ------------------------------------
TH1D* makeNumHistEta(const char* name, const std::set<EtaPKey>& matched,
                     const std::map<EtaPKey, EtaParticleInfo>& partMap)
{
  // NB: no Sumw2 - unweighted counts keep TEfficiency in Clopper-Pearson mode
  TH1D* h = makeEtaHist(name);
  for (const auto& key : matched) {
    auto it = partMap.find(key);
    if (it != partMap.end())
      h->Fill(it->second.eta);
  }
  return h;
}

// ---- TEfficiency -> styled graph, drawn on current pad ---------------------
TGraphAsymmErrors* drawEffEta(TH1D* hNum, TH1D* hDen, int color,
                              int markerStyle, const char* name)
{
  if (!TEfficiency::CheckConsistency(*hNum, *hDen)) {
    std::cout << "ERROR: inconsistent histograms for " << name << std::endl;
    return nullptr;
  }
  TEfficiency eff(*hNum, *hDen);
  eff.SetStatisticOption(TEfficiency::kFCP); // Clopper-Pearson
  TGraphAsymmErrors* g = eff.CreateGraph();
  g->SetName(name);
  g->SetMarkerStyle(markerStyle);
  g->SetMarkerSize(kMarkerSize);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetLineWidth(2);
  g->Draw("PZ same");
  return g;
}

// ============================================================== entry point
void drawCumulativeEfficienciesVsEta(
  const char* dir = "path_to_reco_task_output",
  double etaMin = -4.0, double etaMax = 4.0,
  double ptMin = 0.0, double ptMax = 1.e9,
  int nIterations = 4, int nEtaBins = kNEtaBinsDefault)
{
  TString d(dir);
  if (d.EndsWith("/"))
    d.Remove(d.Length() - 1);

  // ----------------------------------------------------------- eta binning
  gEtaEdges.clear();
  for (int i = 0; i <= nEtaBins; ++i)
    gEtaEdges.push_back(etaMin + (etaMax - etaMin) * i / nEtaBins);

  // ---------------------------------------------------------- denominator
  std::map<EtaPKey, EtaParticleInfo> partMap;
  TH1D* hAll = makeEtaHist("hAllEta");
  std::cout << "Reading particles..." << std::endl;
  if (!readParticlesEta(d + "/particles_digitized_selected.root", etaMin, etaMax,
                        ptMin, ptMax, partMap, hAll))
    return;

  // ------------------------------- cumulative seeding: union over the passes
  std::set<EtaPKey> seedMergedSet;
  for (int it = 0; it < nIterations; ++it) {
    TString fSeed;
    if (nIterations <= 1)
      fSeed = d + "/estimatedparameters.root";
    else if (it == 0)
      fSeed = d + Form("/iter_%d/estimatedparameters.root", it);
    else
      fSeed = d + Form("/iter_%d/iter_%d_estimatedparameters.root", it, it);
    std::set<EtaPKey> seedMatched;
    collectSeedMatchesEta(fSeed, partMap, seedMatched);
    std::cout << "Iteration " << it << ": seed-matched particles : "
              << seedMatched.size() << std::endl;
    seedMergedSet.insert(seedMatched.begin(), seedMatched.end());
  }

  // ---------------------------------- cumulative tracking: merged collection
  std::set<EtaPKey> trkMergedSet;
  std::cout << "Merged tracks:" << std::endl;
  collectTrackMatchesEta(d + "/tracksummary_ambi-tracks-merged.root", partMap,
                         trkMergedSet);
  std::cout << "  seed-matched (union)   : " << seedMergedSet.size()
            << "\n  track-matched particles: " << trkMergedSet.size()
            << std::endl;

  // ------------------------------------------------- numerator histograms
  TH1D* hSeedMergedNum = makeNumHistEta("hSeedNumEta_merged", seedMergedSet, partMap);
  TH1D* hTrkMergedNum = makeNumHistEta("hTrkNumEta_merged", trkMergedSet, partMap);

  // -------------------------------------------------------------- canvases
  // One canvas per mode: 0 = seeding + tracking, 1 = seeding only,
  // 2 = tracking only. Only the graphs of mode 0 are written to the .root
  // file, so that their names stay gSeedEff_merged / gTrkEff_merged.
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  std::vector<TGraphAsymmErrors*> allGraphs;
  TCanvas* cEff[3] = {nullptr, nullptr, nullptr};

  for (int mode = 0; mode < 3; ++mode) {
    const char* sfx = kModeSuffix[mode];
    TCanvas* cm = new TCanvas(Form("cEffEta%s", sfx), "efficiencies vs eta",
                              1100, 750);
    styleCanvasEta(cm);
    cm->cd();

    TH1D* hFrame = makeEtaHist(Form("hFrameEta%s", sfx));
    hFrame->Reset();
    styleAxesEta(hFrame, "#eta", "Technical Efficiency");
    hFrame->GetYaxis()->SetRangeUser(0.0, 1.1);
    hFrame->Draw("axis");

    if (mode != 2) { // open markers: seeding
      TGraphAsymmErrors* g = drawEffEta(hSeedMergedNum, hAll, kCumulColor, 24,
                                        Form("gSeedEff_merged%s", sfx));
      if (mode == 0)
        allGraphs.push_back(g);
    }
    if (mode != 1) { // full markers: tracking
      TGraphAsymmErrors* g = drawEffEta(hTrkMergedNum, hAll, kCumulColor, 20,
                                        Form("gTrkEff_merged%s", sfx));
      if (mode == 0)
        allGraphs.push_back(g);
    }

    // ------------------------------------------------------- annotations
    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.038);
    tex.DrawLatex(0.38, 0.36, kEtaInfoLine1);
    tex.DrawLatex(0.38, 0.30, kEtaInfoLine2);

    // marker convention: shown as markers, not described in words
    TLegend* legMark = new TLegend(0.40, 0.15, 0.64, 0.27);
    legMark->SetBorderSize(0);
    legMark->SetFillStyle(0);
    legMark->SetTextFont(42);
    legMark->SetTextSize(0.038);
    if (mode != 2)
      legMark->AddEntry(legendKeyEta(kCumulColor, 24), "seeding", "P");
    if (mode != 1)
      legMark->AddEntry(legendKeyEta(kCumulColor, 20), "tracking", "P");
    legMark->Draw();

    cm->RedrawAxis();
    cEff[mode] = cm;
  }
  TCanvas* c = cEff[0];

  // ------------------------------------------------------------- output
  // the eta range is part of the name: restricted-range versions of the same
  // sample must not overwrite each other
  TString out = (ptMax > 1.e8)
                  ? Form("%s/efficiencies_vs_eta%g-%g_ptGt%g", dir, etaMin,
                         etaMax, ptMin)
                  : Form("%s/efficiencies_vs_eta%g-%g_pt%g-%g", dir, etaMin,
                         etaMax, ptMin, ptMax);
  out.ReplaceAll(".", "p");
  for (int mode = 0; mode < 3; ++mode) {
    cEff[mode]->SaveAs(out + kModeSuffix[mode] + ".png");
    cEff[mode]->SaveAs(out + kModeSuffix[mode] + ".pdf");
  }

  TFile* fOut = TFile::Open(out + ".root", "RECREATE");
  if (fOut && !fOut->IsZombie()) {
    for (TGraphAsymmErrors* g : allGraphs)
      if (g)
        g->Write();
    c->Write("cEffEta");
    fOut->Close();
    std::cout << "Saved " << out << ".png / .pdf / .root" << std::endl;
  } else {
    std::cout << "WARNING: could not open " << out << ".root for writing"
              << std::endl;
    std::cout << "Saved " << out << ".png / .pdf" << std::endl;
  }
}
