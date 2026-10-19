// draw_iterative_efficiencies.C
//
// Draws seeding and tracking efficiencies vs pT for all tracking iterations
// (passes) of the ACTS iterative tracking chain, plus the efficiencies of the
// merged track collection - like the reference "Tracking efficiency with
// |eta| in [x,y]" figure.
//
// Everything is computed directly from the output TTrees:
//   - denominator : particles_digitized_selected.root   (tree "particles")
//                   == the "particles_selected" collection used by the ACTS
//                   performance writers
//   - seeding eff : estimatedparameters.root            (tree "estimatedparams")
//                   iter 0 in <dir>/, iter N in <dir>/iter_N/iter_N_estimatedparameters.root
//                   numerator = particles that are the (pure) truth match of
//                   at least one seed (branch "truthMatched")
//   - tracking eff: tracksummary_ambi.root              (tree "tracksummary")
//                   iter 0 in <dir>/, iter N in <dir>/iter_N/tracksummary_iter_N_ambi.root
//                   numerator = particles that have at least one track with
//                   trackClassification == Matched (== 1)
//   - merged      : tracking -> tracksummary_ambi-tracks-merged.root
//                   seeding  -> union of the per-iteration seed-matched particles
//
// A particle is identified by (event number, barcode components), and its
// pT/eta are always taken from the particles tree, so numerator and
// denominator are guaranteed to be consistent.
//
// Note on definitions: the seed numerator uses "truthMatched" from
// estimatedparams, which requires ALL seed space points to come from one
// particle (purity = 1). The ACTS TrackFinderPerformanceWriter instead uses a
// double-majority criterion, so the seeding efficiency here can be very
// slightly lower than in performance_seeding.root.
//
// Fake and duplicate rates are per-candidate fractions (NOT per-particle): the
// denominator is the number of reconstructed candidates in the pT bin, so
// within one collection matched + duplicate + fake = 1 - except for tracks left
// "Unknown" (classification 0): those with no contributing particle at all, or
// whose majority particle is absent from the input particle collection, get no
// entry in the matching map and are therefore counted as neither.
//
//   - FAKE RATE. Seeds: fraction of seeds whose space points do not all come
//     from one and the same particle ("truthMatched == false" in
//     estimatedparams is literally "more than one contributing particle", i.e.
//     purity < 1), a combinatorial seed. Tracks: fraction of tracks with
//     trackClassification == Fake, which in this chain means
//     nMajorityHits / nMeasurements < matchingRatio = 0.5, i.e. no single
//     particle owns the majority of the track's measurements. Note that
//     writers.py sets doubleMatching = false, so only this reco-side leg is
//     applied; the particle-side leg (the track holding >= 50% of that
//     particle's hits, and the looper protection patching it) is computed but
//     unused.
//   - DUPLICATE RATE. Seeds: fraction of seeds that are truth-matched to a
//     particle already seeded by another seed of the same collection - the
//     first seed of a particle counts as matched, every further one as a
//     duplicate. Tracks: fraction of tracks with trackClassification ==
//     Duplicate, i.e. several tracks pass the match above for the same
//     particle; TrackTruthMatcher keeps the one with more measurements or the
//     smaller chi2 as Matched and flags the others. For the merged/union curves
//     the seed bookkeeping runs across all iterations, so a particle seeded in
//     two different passes contributes one match and one duplicate.
//
// Outputs (next to the input dir, name encodes eta range / logPt):
//   efficiencies_vs_pT_eta<...>.png / .pdf   - the efficiency canvas
//   fakerate_seeds_vs_pT_eta<...>.png / .pdf
//   fakerate_tracks_vs_pT_eta<...>.png / .pdf
//   duprate_seeds_vs_pT_eta<...>.png / .pdf
//   duprate_tracks_vs_pT_eta<...>.png / .pdf
//   duprate_seeds_and_tracks_vs_pT_eta<...>.png / .pdf  - both on one canvas,
//                                              open markers seeds, closed tracks
//   efficiencies_vs_pT_eta<...>.root         - all efficiency / fake /
//                                              duplicate graphs + the canvases
//
// Usage:
//   drawIterativeEfficiencies(dir, etaMin, etaMax, ptMin, ptMax, nIterations,
//                             logPt, nPtBins, yMax, nameSuffix)
//
//   root -l 'drawIterativeEfficiencies.C("/path/to/output_dir", 0.0, 0.8, 0.05, 2.0, 3)'
//   logarithmic pT axis (log-spaced bins):
//   root -l 'drawIterativeEfficiencies.C("/path/to/output_dir", 0.0, 0.8, 0.05, 2.0, 3, true)'
//   |eta| < 1, 0.04-1 GeV in 0.02 GeV bins, y up to 1.5, tagged file names:
//   root -l 'drawIterativeEfficiencies.C("/path/to/output_dir", 0.0, 1.0, 0.04, 1.0, 4, false, 48, 1.5, "nhits_comparison")'
//
// yMax sets the upper edge of the y axis on every canvas; nameSuffix is
// appended to every output file name (a leading "_" is added if missing), which
// is how several runs of the same range are kept side by side.
//
// Igor / Claude, June 2026

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TEfficiency.h"
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
static const int kNPtBins = 60;
// static const double kPtMin = 0.0;     // GeV/c (linear axis)
// static const double kPtMinLog = 0.05; // GeV/c (log axis: lowest bin edge, must be > 0)
// static const double kPtMax = 2.0;     // GeV/c

// pT bin edges, filled at startup (linear or log-spaced)
static std::vector<double> gPtEdges;

// annotation text on the plot (edit to match your setup)
// experiment tag, drawn in Helvetica bold italic (ROOT font 72) on every plot
static const char* kExperimentLabel = "ALICE3-ACTS";
static const char* kInfoLine1 = "Pythia8 pp, pileup 30";
static const char* kInfoLine2 = "Geant4, %d-pass reco";

// colors per iteration (extend if you run more passes)
static const int kIterColors[] = {kRed + 1, kBlue, kGreen + 2,
                                  kMagenta + 1, kOrange + 7, kCyan + 2};
static const int kMergedColor = kGray + 1;

// ------------------------------------------------------------------ helpers

// unique particle key: {event, vertex_primary, vertex_secondary, particle,
//                       generation, sub_particle}
typedef std::array<UInt_t, 6> PKey;

struct ParticleInfo {
  float pt;
  float eta;
};

// create a TH1D with the global pT binning
TH1D* makePtHist(const char* name)
{
  return new TH1D(name, name, gPtEdges.size() - 1, gPtEdges.data());
}

// open a file, return nullptr (with warning) if missing
TFile* openFile(const TString& path)
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
bool readParticles(const TString& path, double etaMin, double etaMax,
                   std::map<PKey, ParticleInfo>& partMap, TH1D* hAll)
{
  TFile* f = openFile(path);
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
      double aeta = TMath::Abs((*eta)[j]);
      if (aeta < etaMin || aeta > etaMax)
        continue;
      PKey key = {eventId, (*vp)[j], (*vs)[j], (*par)[j], (*gen)[j], (*sub)[j]};
      if (partMap.count(key))
        continue; // safety against duplicates
      partMap[key] = {(*pt)[j], (*eta)[j]};
      hAll->Fill((*pt)[j]);
    }
  }
  f->Close();
  std::cout << "  particles passing |eta| cut: " << partMap.size() << std::endl;
  return true;
}

// ---- collect seed-matched particles from an estimatedparams tree ----------
void collectSeedMatches(const TString& path,
                        const std::map<PKey, ParticleInfo>& partMap,
                        std::set<PKey>& matched)
{
  TFile* f = openFile(path);
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
    PKey key = {(UInt_t)eventNr, vp, vs, par, gen, sub};
    if (!partMap.count(key))
      continue; // not in selected denominator
    matched.insert(key);
  }
  f->Close();
}

// ---- collect track-matched particles from a tracksummary tree -------------
void collectTrackMatches(const TString& path,
                         const std::map<PKey, ParticleInfo>& partMap,
                         std::set<PKey>& matched)
{
  TFile* f = openFile(path);
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
      PKey key = {eventNr, (*vp)[j], (*vs)[j], (*par)[j], (*gen)[j], (*sub)[j]};
      if (!partMap.count(key))
        continue;
      matched.insert(key);
    }
  }
  f->Close();
}

// ---- fill (total, fake, duplicate) SEED histograms from estimatedparams ----
// One entry per seed; binned by the seed's ESTIMATED pT, eta cut on the
// estimated eta.
//   fake      = seed that is not a pure truth match (truthMatched == false)
//   duplicate = truth-matched seed whose particle was already seeded by an
//               earlier seed of this file (the first seed of a particle is the
//               match, every further one a duplicate)
// If seenMerged / hDupMerged are given, the same duplicate bookkeeping is
// repeated against a set shared between files, which is what makes the merged
// curve count a particle re-seeded in a later iteration as a duplicate.
void fillSeedCounts(const TString& path, double etaMin, double etaMax,
                    TH1D* hTot, TH1D* hFake, TH1D* hDup,
                    std::set<PKey>* seenMerged = nullptr,
                    TH1D* hDupMerged = nullptr)
{
  TFile* f = openFile(path);
  if (!f)
    return;
  TTree* t = (TTree*)f->Get("estimatedparams");
  if (!t) {
    std::cout << "ERROR: no 'estimatedparams' tree in " << path << std::endl;
    return;
  }
  Float_t pt = 0.f, eta = 0.f;
  Bool_t isMatched = false;
  Int_t eventNr = 0;
  UInt_t vp = 0, vs = 0, par = 0, gen = 0, sub = 0;
  t->SetBranchStatus("*", 0);
  for (const char* b : {"pt", "eta", "truthMatched", "event_nr",
                        "particleId_vertex_primary",
                        "particleId_vertex_secondary",
                        "particleId_particle", "particleId_generation",
                        "particleId_sub_particle"})
    t->SetBranchStatus(b, 1);
  t->SetBranchAddress("pt", &pt);
  t->SetBranchAddress("eta", &eta);
  t->SetBranchAddress("truthMatched", &isMatched);
  t->SetBranchAddress("event_nr", &eventNr);
  t->SetBranchAddress("particleId_vertex_primary", &vp);
  t->SetBranchAddress("particleId_vertex_secondary", &vs);
  t->SetBranchAddress("particleId_particle", &par);
  t->SetBranchAddress("particleId_generation", &gen);
  t->SetBranchAddress("particleId_sub_particle", &sub);

  // particles already seeded within this file
  std::set<PKey> seen;

  Long64_t nEnt = t->GetEntries();
  for (Long64_t i = 0; i < nEnt; ++i) {
    t->GetEntry(i);
    double aeta = TMath::Abs(eta);
    if (aeta < etaMin || aeta > etaMax)
      continue;
    hTot->Fill(pt);
    if (!isMatched) {
      hFake->Fill(pt);
      continue;
    }
    PKey key = {(UInt_t)eventNr, vp, vs, par, gen, sub};
    if (!seen.insert(key).second) // already seeded in this file
      hDup->Fill(pt);
    if (seenMerged && hDupMerged && !seenMerged->insert(key).second)
      hDupMerged->Fill(pt);
  }
  f->Close();
}

// ---- fill (total, fake, duplicate) TRACK histograms from tracksummary ------
// fake = trackClassification == 3 (Fake), duplicate == 2 (Duplicate). Binned by
// the fitted pT computed from (theta, q/p); eta cut on the fitted eta. Tracks
// with an unusable fit (q/p == 0 or non-finite) are skipped.
void fillTrackCounts(const TString& path, double etaMin, double etaMax,
                     TH1D* hTot, TH1D* hFake, TH1D* hDup)
{
  TFile* f = openFile(path);
  if (!f)
    return;
  TTree* t = (TTree*)f->Get("tracksummary");
  if (!t) {
    std::cout << "ERROR: no 'tracksummary' tree in " << path << std::endl;
    return;
  }
  std::vector<int>* cls = nullptr;
  std::vector<float>* theta = nullptr;
  std::vector<float>* qop = nullptr;
  t->SetBranchStatus("*", 0);
  for (const char* b : {"trackClassification", "eTHETA_fit", "eQOP_fit"})
    t->SetBranchStatus(b, 1);
  t->SetBranchAddress("trackClassification", &cls);
  t->SetBranchAddress("eTHETA_fit", &theta);
  t->SetBranchAddress("eQOP_fit", &qop);

  // 0=Unknown, 1=Matched, 2=Duplicate, 3=Fake
  const int kDuplicate = 2;
  const int kFake = 3;
  Long64_t nEnt = t->GetEntries();
  for (Long64_t i = 0; i < nEnt; ++i) {
    t->GetEntry(i);
    for (size_t j = 0; j < cls->size(); ++j) {
      double th = (*theta)[j];
      double q = (*qop)[j];
      if (q == 0. || !std::isfinite(q) || !std::isfinite(th))
        continue;
      double pt = std::sin(th) / std::abs(q);
      double eta = -std::log(std::tan(th / 2.));
      double aeta = TMath::Abs(eta);
      if (aeta < etaMin || aeta > etaMax)
        continue;
      hTot->Fill(pt);
      if ((*cls)[j] == kFake)
        hFake->Fill(pt);
      else if ((*cls)[j] == kDuplicate)
        hDup->Fill(pt);
    }
  }
  f->Close();
}

// ---- matched set -> numerator histogram ------------------------------------
TH1D* makeNumHist(const char* name, const std::set<PKey>& matched,
                  const std::map<PKey, ParticleInfo>& partMap)
{
  // NB: no Sumw2 - unweighted counts keep TEfficiency in Clopper-Pearson mode
  TH1D* h = makePtHist(name);
  for (const auto& key : matched) {
    auto it = partMap.find(key);
    if (it != partMap.end())
      h->Fill(it->second.pt);
  }
  return h;
}

// ---- TEfficiency -> styled graph, drawn on current pad ---------------------
TGraphAsymmErrors* drawEff(TH1D* hNum, TH1D* hDen, int color, int markerStyle,
                           const char* name)
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
  g->SetMarkerSize(0.9);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetLineWidth(2);
  g->Draw("PZ same");
  return g;
}

// ---- per-iteration input files --------------------------------------------
void iterationFiles(const TString& d, int it, int nIterations, TString& fSeed,
                    TString& fTrk)
{
  if (nIterations <= 1) {
    fSeed = d + "/estimatedparameters.root";
    fTrk = d + "/tracksummary_ambi_tracks_writer.root";
  } else if (it == 0) {
    fSeed = d + Form("/iter_%d/estimatedparameters.root", it);
    fTrk = d + Form("/iter_%d/tracksummary_ambi_tracks_writer.root", it);
  } else {
    fSeed = d + Form("/iter_%d/iter_%d_estimatedparameters.root", it, it);
    fTrk = d + Form("/iter_%d/tracksummary_iter_%d_ambi_tracks_writer.root", it,
                    it);
  }
}

// ---- experiment tag -------------------------------------------------------
// "ALICE3-ACTS" in Helvetica bold italic (font 72), identical on every canvas.
void drawExperimentLabel(double x, double y, double size = 0.045)
{
  TLatex t;
  t.SetNDC();
  t.SetTextFont(72); // Helvetica bold italic
  t.SetTextSize(size);
  t.SetTextColor(kBlack);
  t.DrawLatex(x, y, kExperimentLabel);
}

// ---- canvas + empty frame for a rate plot (fake / duplicate) --------------
TCanvas* makeRateCanvas(const char* cname, const char* title,
                        const char* yTitle, bool logPt, double yMax)
{
  TCanvas* c = new TCanvas(cname, title, 1100, 750);
  c->SetTicks(1, 1);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.09);
  if (logPt)
    c->SetLogx();
  c->cd();

  TH1D* h = makePtHist(Form("hFrame_%s", cname));
  h->Reset();
  h->SetTitle(title);
  h->GetXaxis()->SetTitle("pT [GeV/c]");
  if (logPt)
    h->GetXaxis()->SetMoreLogLabels();
  h->GetYaxis()->SetTitle(yTitle);
  h->GetYaxis()->SetRangeUser(0.0, yMax);
  h->Draw("axis");
  return c;
}

// ---- info lines + per-iteration colour legend on a rate canvas -------------
void annotateRateCanvas(int nIterations, const char* openNote = nullptr,
                        const char* closedNote = nullptr)
{
  const int nColors = sizeof(kIterColors) / sizeof(kIterColors[0]);
  // with a marker key on the right the info block moves left, out of its way
  const bool hasNotes = (openNote != nullptr || closedNote != nullptr);
  const double xInfo = hasNotes ? 0.16 : 0.42;
  drawExperimentLabel(xInfo, 0.855);
  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(42);
  tex.SetTextSize(0.038);
  tex.DrawLatex(xInfo, 0.805, kInfoLine1);
  tex.DrawLatex(xInfo, 0.755, Form(kInfoLine2, nIterations));
  if (hasNotes) {
    tex.SetTextAlign(31); // right-aligned, so the longer note cannot run off
    if (openNote != nullptr)
      tex.DrawLatex(0.955, 0.855, openNote);
    if (closedNote != nullptr)
      tex.DrawLatex(0.955, 0.805, closedNote);
    tex.SetTextAlign(11);
  }
  double y = 0.685;
  for (int it = 0; it < nIterations; ++it) {
    tex.SetTextColor(kIterColors[it % nColors]);
    tex.DrawLatex(xInfo, y, Form("iteration %d", it));
    y -= 0.05;
  }
  tex.SetTextColor(kMergedColor);
  tex.DrawLatex(xInfo, y, "merged");
  tex.SetTextColor(kBlack);
  gPad->RedrawAxis();
}

// ---- save a canvas as .png / .pdf, return the (extension-less) path --------
TString saveCanvas(TCanvas* c, const char* dir, const char* stem, double etaMin,
                   double etaMax, bool logPt, const TString& nameSuffix)
{
  TString name = Form("%s_vs_pT_eta%g-%g%s%s", stem, etaMin, etaMax,
                      logPt ? "_logPt" : "", nameSuffix.Data());
  name.ReplaceAll(".", "p"); // 4.2 -> 4p2; the directory is prepended after
  TString out = TString(dir) + "/" + name;
  c->SaveAs(out + ".png");
  c->SaveAs(out + ".pdf");
  return out;
}

// ============================================================== entry point
void drawIterativeEfficiencies(
  const char* dir = "path_to_reco_task_output",
  double etaMin = -1.0, double etaMax = 1.0, double kPtMin = 0.05, double kPtMax = 2.0, int nIterations = 3, bool logPt = false,
  int nPtBins = kNPtBins,  // nPtBins=1 -> single integrated efficiency over [kPtMin,kPtMax]
  double yMax = 1.1,       // upper edge of the y axis of every canvas
  const char* nameSuffix = "")  // appended to every output file name
{
  TString d(dir);
  if (d.EndsWith("/"))
    d.Remove(d.Length() - 1);

  // "nhits_comparison" and "_nhits_comparison" both give ..._nhits_comparison
  TString sfx(nameSuffix);
  if (!sfx.IsNull() && !sfx.BeginsWith("_"))
    sfx.Prepend("_");

  // ------------------------------------------------------------ pT binning
  gPtEdges.clear();
  if (logPt) { // log-spaced bins from kPtMin to kPtMax
    const double lmin = TMath::Log10(kPtMin);
    const double lmax = TMath::Log10(kPtMax);
    for (int i = 0; i <= nPtBins; ++i)
      gPtEdges.push_back(TMath::Power(10., lmin + (lmax - lmin) * i / nPtBins));
  } else { // equidistant bins
    for (int i = 0; i <= nPtBins; ++i)
      gPtEdges.push_back(kPtMin + (kPtMax - kPtMin) * i / nPtBins);
  }

  // ---------------------------------------------------------- denominator
  std::map<PKey, ParticleInfo> partMap;
  TH1D* hAll = makePtHist("hAll");
  std::cout << "Reading particles..." << std::endl;
  if (!readParticles(d + "/particles_digitized_selected.root", etaMin, etaMax,
                     partMap, hAll))
    return;

  // ------------------------------------------------- per-iteration matches
  std::vector<std::set<PKey>> seedMatched(nIterations), trkMatched(nIterations);
  std::set<PKey> seedMergedSet;

  for (int it = 0; it < nIterations; ++it) {
    TString fSeed, fTrk;
    iterationFiles(d, it, nIterations, fSeed, fTrk);
    std::cout << "Iteration " << it << ": " << std::endl;
    collectSeedMatches(fSeed, partMap, seedMatched[it]);
    collectTrackMatches(fTrk, partMap, trkMatched[it]);
    std::cout << "  seed-matched particles : " << seedMatched[it].size()
              << "\n  track-matched particles: " << trkMatched[it].size()
              << std::endl;
    seedMergedSet.insert(seedMatched[it].begin(), seedMatched[it].end());
  }

  // ------------------------------------------------------------- merged
  std::set<PKey> trkMergedSet;
  std::cout << "Merged tracks:" << std::endl;
  collectTrackMatches(d + "/tracksummary_ambi-tracks-merged.root", partMap,
                      trkMergedSet);
  std::cout << "  seed-matched (union)   : " << seedMergedSet.size()
            << "\n  track-matched particles: " << trkMergedSet.size()
            << std::endl;

  // -------------------------------------------------------------- canvas
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas* c = new TCanvas("cEff", "efficiencies", 1100, 750);
  c->SetTicks(1, 1);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.09);

  if (logPt)
    c->SetLogx();

  TH1D* hFrame = makePtHist("hFrame");
  hFrame->Reset();
  hFrame->SetTitle(Form("Tracking efficiency with |#eta| in [%g, %g]",
                        etaMin, etaMax));
  hFrame->GetXaxis()->SetTitle("pT [GeV/c]");
  if (logPt)
    hFrame->GetXaxis()->SetMoreLogLabels();
  hFrame->GetYaxis()->SetTitle("Efficiency");
  hFrame->GetYaxis()->SetRangeUser(0.0, yMax);
  hFrame->Draw("axis");

  // ------------------------------------------------------- draw all curves
  // Keep every efficiency graph so it can be written to the .root file below.
  std::vector<TGraphAsymmErrors*> allGraphs;
  const int nColors = sizeof(kIterColors) / sizeof(kIterColors[0]);
  for (int it = 0; it < nIterations; ++it) {
    int color = kIterColors[it % nColors];
    TH1D* hSeedNum = makeNumHist(Form("hSeedNum_%d", it), seedMatched[it], partMap);
    TH1D* hTrkNum = makeNumHist(Form("hTrkNum_%d", it), trkMatched[it], partMap);
    allGraphs.push_back(drawEff(hSeedNum, hAll, color, 24, Form("gSeedEff_%d", it))); // open
    allGraphs.push_back(drawEff(hTrkNum, hAll, color, 20, Form("gTrkEff_%d", it)));   // closed
  }
  TH1D* hSeedMergedNum = makeNumHist("hSeedNum_merged", seedMergedSet, partMap);
  TH1D* hTrkMergedNum = makeNumHist("hTrkNum_merged", trkMergedSet, partMap);
  allGraphs.push_back(drawEff(hSeedMergedNum, hAll, kMergedColor, 24, "gSeedEff_merged"));
  allGraphs.push_back(drawEff(hTrkMergedNum, hAll, kMergedColor, 20, "gTrkEff_merged"));

  // --------------------------------------------------------- annotations
  // An efficiency never exceeds 1, so when the y axis is stretched past ~1.15
  // the band above 1 is empty and the annotations go there, in two columns.
  // With the default yMax they keep their historical mid-pad position.
  const bool headroom = (yMax > 1.15);

  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(42);
  tex.SetTextSize(headroom ? 0.034 : 0.038);

  if (headroom) {
    drawExperimentLabel(0.16, 0.875, 0.040);
    tex.DrawLatex(0.16, 0.830, kInfoLine1);
    tex.DrawLatex(0.16, 0.790, Form(kInfoLine2, nIterations));
    tex.DrawLatex(0.16, 0.735, "open markers #minus seeding eff");
    tex.DrawLatex(0.16, 0.695, "closed markers #minus tracking eff");

    double y = 0.875;
    for (int it = 0; it < nIterations; ++it) {
      tex.SetTextColor(kIterColors[it % nColors]);
      tex.DrawLatex(0.62, y, Form("iteration %d", it));
      y -= 0.045;
    }
    tex.SetTextColor(kMergedColor);
    tex.DrawLatex(0.62, y, "merged");
    tex.SetTextColor(kBlack);
  } else {
    drawExperimentLabel(0.42, 0.655);
    tex.DrawLatex(0.42, 0.60, kInfoLine1);
    tex.DrawLatex(0.42, 0.55, Form(kInfoLine2, nIterations));

    double y = 0.47;
    for (int it = 0; it < nIterations; ++it) {
      tex.SetTextColor(kIterColors[it % nColors]);
      tex.DrawLatex(0.42, y, Form("iteration %d", it));
      y -= 0.05;
    }
    tex.SetTextColor(kMergedColor);
    tex.DrawLatex(0.42, y, "merged");

    tex.SetTextColor(kBlack);
    tex.DrawLatex(0.60, 0.38, "open markers #minus seeding eff");
    tex.DrawLatex(0.60, 0.33, "closed markers #minus tracking eff");
  }

  c->RedrawAxis();

  // ================================== further canvases: fake / duplicate rate
  // Per-candidate fractions vs pT (nFake / nTotal and nDuplicate / nTotal),
  // binned by the reconstructed pT, with Clopper-Pearson errors (reuse
  // drawEff). Seeds and tracks get their own canvas: the two collections have
  // different denominators and very different magnitudes, so overlaying them
  // only made the plot ambiguous.
  std::vector<TGraphAsymmErrors*> rateGraphs;

  TCanvas* cFakeSeed = makeRateCanvas(
    "cFakeSeed", Form("Seed fake rate with |#eta| in [%g, %g]", etaMin, etaMax),
    "Fake rate", logPt, yMax);
  TCanvas* cFakeTrk = makeRateCanvas(
    "cFakeTrk", Form("Track fake rate with |#eta| in [%g, %g]", etaMin, etaMax),
    "Fake rate", logPt, yMax);
  TCanvas* cDupSeed = makeRateCanvas(
    "cDupSeed",
    Form("Seed duplicate rate with |#eta| in [%g, %g]", etaMin, etaMax),
    "Duplicate rate", logPt, yMax);
  TCanvas* cDupTrk = makeRateCanvas(
    "cDupTrk",
    Form("Track duplicate rate with |#eta| in [%g, %g]", etaMin, etaMax),
    "Duplicate rate", logPt, yMax);
  // seeds and tracks on one canvas as well: open markers seeds, closed tracks
  TCanvas* cDupBoth = makeRateCanvas(
    "cDupBoth",
    Form("Duplicate rate with |#eta| in [%g, %g]", etaMin, etaMax),
    "Duplicate rate", logPt, yMax);

  // merged-seed counts = accumulation over the per-iteration seed collections;
  // the merged duplicates use a set shared across iterations (see
  // fillSeedCounts), so a particle re-seeded in a later pass counts as one.
  TH1D* hSeedTotMerged = makePtHist("hSeedTot_merged");
  TH1D* hSeedFakeMerged = makePtHist("hSeedFake_merged");
  TH1D* hSeedDupMerged = makePtHist("hSeedDup_merged");
  std::set<PKey> seedSeenMerged;

  for (int it = 0; it < nIterations; ++it) {
    int color = kIterColors[it % nColors];
    TString fSeed, fTrk;
    iterationFiles(d, it, nIterations, fSeed, fTrk);

    TH1D* hSeedTot = makePtHist(Form("hSeedTot_%d", it));
    TH1D* hSeedFake = makePtHist(Form("hSeedFake_%d", it));
    TH1D* hSeedDup = makePtHist(Form("hSeedDup_%d", it));
    fillSeedCounts(fSeed, etaMin, etaMax, hSeedTot, hSeedFake, hSeedDup,
                   &seedSeenMerged, hSeedDupMerged);
    hSeedTotMerged->Add(hSeedTot);
    hSeedFakeMerged->Add(hSeedFake);
    cFakeSeed->cd();
    rateGraphs.push_back(
      drawEff(hSeedFake, hSeedTot, color, 24, Form("gSeedFake_%d", it)));
    cDupSeed->cd();
    rateGraphs.push_back(
      drawEff(hSeedDup, hSeedTot, color, 24, Form("gSeedDup_%d", it)));
    cDupBoth->cd(); // same numbers, kept with the canvas only
    drawEff(hSeedDup, hSeedTot, color, 24, Form("gSeedDup_%d_both", it));

    TH1D* hTrkTot = makePtHist(Form("hTrkTot_%d", it));
    TH1D* hTrkFake = makePtHist(Form("hTrkFake_%d", it));
    TH1D* hTrkDup = makePtHist(Form("hTrkDup_%d", it));
    fillTrackCounts(fTrk, etaMin, etaMax, hTrkTot, hTrkFake, hTrkDup);
    cFakeTrk->cd();
    rateGraphs.push_back(
      drawEff(hTrkFake, hTrkTot, color, 20, Form("gTrkFake_%d", it)));
    cDupTrk->cd();
    rateGraphs.push_back(
      drawEff(hTrkDup, hTrkTot, color, 20, Form("gTrkDup_%d", it)));
    cDupBoth->cd();
    drawEff(hTrkDup, hTrkTot, color, 20, Form("gTrkDup_%d_both", it));
  }

  TH1D* hTrkTotMerged = makePtHist("hTrkTot_merged");
  TH1D* hTrkFakeMerged = makePtHist("hTrkFake_merged");
  TH1D* hTrkDupMerged = makePtHist("hTrkDup_merged");
  fillTrackCounts(d + "/tracksummary_ambi-tracks-merged.root", etaMin, etaMax,
                  hTrkTotMerged, hTrkFakeMerged, hTrkDupMerged);

  cFakeSeed->cd();
  rateGraphs.push_back(drawEff(hSeedFakeMerged, hSeedTotMerged, kMergedColor,
                               24, "gSeedFake_merged"));
  annotateRateCanvas(nIterations);

  cDupSeed->cd();
  rateGraphs.push_back(drawEff(hSeedDupMerged, hSeedTotMerged, kMergedColor, 24,
                               "gSeedDup_merged"));
  annotateRateCanvas(nIterations);

  cFakeTrk->cd();
  rateGraphs.push_back(drawEff(hTrkFakeMerged, hTrkTotMerged, kMergedColor, 20,
                               "gTrkFake_merged"));
  annotateRateCanvas(nIterations);

  cDupTrk->cd();
  rateGraphs.push_back(drawEff(hTrkDupMerged, hTrkTotMerged, kMergedColor, 20,
                               "gTrkDup_merged"));
  annotateRateCanvas(nIterations);

  cDupBoth->cd();
  drawEff(hSeedDupMerged, hSeedTotMerged, kMergedColor, 24, "gSeedDup_merged_both");
  drawEff(hTrkDupMerged, hTrkTotMerged, kMergedColor, 20, "gTrkDup_merged_both");
  annotateRateCanvas(nIterations, "open markers #minus seed duplicates",
                     "closed markers #minus track duplicates");

  // ------------------------------------------------------------- output
  TString out = saveCanvas(c, dir, "efficiencies", etaMin, etaMax, logPt, sfx);
  TString outFS = saveCanvas(cFakeSeed, dir, "fakerate_seeds", etaMin, etaMax, logPt, sfx);
  TString outFT = saveCanvas(cFakeTrk, dir, "fakerate_tracks", etaMin, etaMax, logPt, sfx);
  TString outDS = saveCanvas(cDupSeed, dir, "duprate_seeds", etaMin, etaMax, logPt, sfx);
  TString outDT = saveCanvas(cDupTrk, dir, "duprate_tracks", etaMin, etaMax, logPt, sfx);
  TString outDB = saveCanvas(cDupBoth, dir, "duprate_seeds_and_tracks", etaMin,
                             etaMax, logPt, sfx);

  // Save all efficiency graphs (per-iteration seeding + tracking, and merged),
  // the fake / duplicate rate graphs and the canvases to a ROOT file. Each
  // graph keeps its name set in drawEff (gSeedEff_<it> / gTrkEff_<it> /
  // gSeedFake_<it> / gTrkFake_<it> / gSeedDup_<it> / gTrkDup_<it>, and the
  // corresponding _merged).
  TFile* fOut = TFile::Open(out + ".root", "RECREATE");
  if (fOut && !fOut->IsZombie()) {
    for (TGraphAsymmErrors* g : allGraphs)
      if (g)
        g->Write();
    for (TGraphAsymmErrors* g : rateGraphs)
      if (g)
        g->Write();
    c->Write("cEff");
    cFakeSeed->Write("cFakeSeed");
    cFakeTrk->Write("cFakeTrk");
    cDupSeed->Write("cDupSeed");
    cDupTrk->Write("cDupTrk");
    cDupBoth->Write("cDupBoth");
    fOut->Close();
    std::cout << "Saved " << out << ".png / .pdf / .root" << std::endl;
  } else {
    std::cout << "WARNING: could not open " << out << ".root for writing"
              << std::endl;
    std::cout << "Saved " << out << ".png / .pdf" << std::endl;
  }
  for (const TString& o : {outFS, outFT, outDS, outDT, outDB})
    std::cout << "Saved " << o << ".png / .pdf" << std::endl;
}
