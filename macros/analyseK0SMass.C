// analyseK0SMass.C
//
// K0S -> pi+ pi- analysis for the ALICE 3 sim output, at TRUTH level and using
// RECONSTRUCTED tracks.
//
//   particles_simulation.root    (tree "particles") : per-particle truth info.
//   vertices_gen_and_geant.root  (tree "vertices")   : mother (incoming) /
//        daughters (outgoing) barcodes; K0S decay vertices identified here.
//        (Needs the ACTS parentParticleId fix + a re-run of the simulation.)
//   <reco>/tracksummary_*.root    (tree "tracksummary") : reconstructed tracks
//        with fitted params and the truth majority particle.
//
// Truth part:  m(pi+ pi-) from the truth daughter momenta.
// Reco part:   each daughter pion may be reconstructed as SEVERAL track
//        "segments" (loopers, esp. low pT). First attempt: build the K0S mass
//        from ALL segment combinations (every reco pi+ segment x every reco pi-
//        segment). Also:
//          - hNSeg : number of reco segments per K0S daughter pion;
//          - hPrel : for pions with >1 segment, p_segment / max_p_segment.
//
// Usage:  root -l 'analyseK0SMass.C("/path/to/sim_dir")'
//         root -l 'analyseK0SMass.C("/path/to/sim_dir","/path/to/reco_dir")'

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TString.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <vector>
#include <iostream>

namespace {
using Key = std::array<std::uint32_t, 5>;  // barcode (vp, vs, particle, gen, sub)
Key keyOf(std::uint32_t vp, std::uint32_t vs, std::uint32_t pa,
          std::uint32_t ge, std::uint32_t sp) { return Key{vp, vs, pa, ge, sp}; }
struct PInfo { int pdg; float px, py, pz, m; };
struct Seg { double px, py, pz, p; };
}  // namespace

void analyseK0SMass(
    const char* dir ="folder_with_sim_output",
    const char* recoDir = "folder_with_rec_output",
    const char* recoFile = "tracksummary_ambi-tracks-merged.root",
    int nMinMeas = 5,                // min reco measurements per pion segment
    double effRMax = 300.0,          // max K0S decay radius [mm] for the eff plot
    int effRBins = 30,
    int kNPhiBins = 32,              // bins for the delta-phi distributions / eff
    double ptMax = 5.0,              // K0S pT axis range [GeV] for eff-vs-pT
    double etaMax = 2.0,             // K0S |eta| axis range for eff-vs-eta
    int nIter = 4)                   // # seeding iterations to union (iter_0..iter_{n-1})
{
  const int kPdgK0S = 310, kPdgPiP = 211, kPdgPiM = -211;
  const double kMassK0SPdg = 0.497611;  // GeV/c^2
  const double kMassPi = 0.13957039;    // GeV/c^2
  if (recoDir == nullptr) recoDir = dir;

  // ---- open files / trees ----
  TFile* fPart = TFile::Open(Form("%s/particles_simulation.root", dir));
  TFile* fVtx  = TFile::Open(Form("%s/vertices_gen_and_geant.root", dir));
  TFile* fReco = TFile::Open(Form("%s/%s", recoDir, recoFile));
  if (!fPart || fPart->IsZombie() || !fVtx || fVtx->IsZombie()) {
    std::cerr << "ERROR: cannot open particles/vertices files in " << dir << "\n";
    return;
  }
  const bool haveReco = (fReco && !fReco->IsZombie());
  if (!haveReco)
    std::cerr << "WARNING: no reco tracksummary at " << recoDir << "/" << recoFile
              << " -- doing truth part only\n";

  TTree* tP = (TTree*)fPart->Get("particles");
  TTree* tV = (TTree*)fVtx->Get("vertices");
  TTree* tR = haveReco ? (TTree*)fReco->Get("tracksummary") : nullptr;
  if (!tP || !tV || (haveReco && !tR)) { std::cerr << "ERROR: missing tree\n"; return; }

  // ---- particle branches ----
  std::uint32_t pEvt = 0;
  std::vector<int>*           pPdg = nullptr;
  std::vector<float>*         pPx = nullptr; std::vector<float>* pPy = nullptr;
  std::vector<float>*         pPz = nullptr; std::vector<float>* pM  = nullptr;
  std::vector<std::uint32_t>* pVp = nullptr; std::vector<std::uint32_t>* pVs = nullptr;
  std::vector<std::uint32_t>* pPa = nullptr; std::vector<std::uint32_t>* pGe = nullptr;
  std::vector<std::uint32_t>* pSp = nullptr;
  tP->SetBranchAddress("event_id", &pEvt);
  tP->SetBranchAddress("particle_type", &pPdg);
  tP->SetBranchAddress("px", &pPx); tP->SetBranchAddress("py", &pPy);
  tP->SetBranchAddress("pz", &pPz); tP->SetBranchAddress("m", &pM);
  tP->SetBranchAddress("vertex_primary", &pVp);
  tP->SetBranchAddress("vertex_secondary", &pVs);
  tP->SetBranchAddress("particle", &pPa);
  tP->SetBranchAddress("generation", &pGe);
  tP->SetBranchAddress("sub_particle", &pSp);

  // ---- vertex branches ----
  std::vector<std::vector<std::uint32_t>> *inVp=nullptr,*inVs=nullptr,*inPa=nullptr,*inGe=nullptr,*inSp=nullptr;
  std::vector<std::vector<std::uint32_t>> *ouVp=nullptr,*ouVs=nullptr,*ouPa=nullptr,*ouGe=nullptr,*ouSp=nullptr;
  tV->SetBranchAddress("incoming_particles_vertex_primary", &inVp);
  tV->SetBranchAddress("incoming_particles_vertex_secondary", &inVs);
  tV->SetBranchAddress("incoming_particles_particle", &inPa);
  tV->SetBranchAddress("incoming_particles_generation", &inGe);
  tV->SetBranchAddress("incoming_particles_sub_particle", &inSp);
  tV->SetBranchAddress("outgoing_particles_vertex_primary", &ouVp);
  tV->SetBranchAddress("outgoing_particles_vertex_secondary", &ouVs);
  tV->SetBranchAddress("outgoing_particles_particle", &ouPa);
  tV->SetBranchAddress("outgoing_particles_generation", &ouGe);
  tV->SetBranchAddress("outgoing_particles_sub_particle", &ouSp);
  std::vector<float> *vVx = nullptr, *vVy = nullptr, *vVz = nullptr;  // vertex position
  tV->SetBranchAddress("vx", &vVx);
  tV->SetBranchAddress("vy", &vVy);
  tV->SetBranchAddress("vz", &vVz);

  // ---- reco (tracksummary) branches ----
  std::uint32_t rEvt = 0;
  std::vector<float>*         rPhi = nullptr; std::vector<float>* rTheta = nullptr;
  std::vector<float>*         rQop = nullptr;
  std::vector<unsigned int>*  rNMeas = nullptr;
  std::vector<std::uint32_t>* rVp = nullptr; std::vector<std::uint32_t>* rVs = nullptr;
  std::vector<std::uint32_t>* rPa = nullptr; std::vector<std::uint32_t>* rGe = nullptr;
  std::vector<std::uint32_t>* rSp = nullptr;
  std::map<std::uint32_t, Long64_t> tsEntryByEvent;
  if (haveReco) {
    tR->SetBranchAddress("event_nr", &rEvt);
    tR->SetBranchAddress("ePHI_fit", &rPhi);
    tR->SetBranchAddress("eTHETA_fit", &rTheta);
    tR->SetBranchAddress("eQOP_fit", &rQop);
    tR->SetBranchAddress("nMeasurements", &rNMeas);
    tR->SetBranchAddress("majorityParticleId_vertex_primary", &rVp);
    tR->SetBranchAddress("majorityParticleId_vertex_secondary", &rVs);
    tR->SetBranchAddress("majorityParticleId_particle", &rPa);
    tR->SetBranchAddress("majorityParticleId_generation", &rGe);
    tR->SetBranchAddress("majorityParticleId_sub_particle", &rSp);
    for (Long64_t i = 0; i < tR->GetEntries(); ++i) {
      tR->GetEntry(i);
      tsEntryByEvent[rEvt] = i;
    }
  }

  // ---- seeding: union of truth-matched seeds over iterations ----
  // seededByEvent[event] = set of particle barcodes with >=1 truth-matched seed
  // (from iter_0/estimatedparameters.root and iter_N/iter_N_estimatedparameters.root).
  std::map<std::uint32_t, std::set<Key>> seededByEvent;
  for (int it = 0; it < nIter; ++it) {
    TString sp = (it == 0)
        ? Form("%s/iter_0/estimatedparameters.root", recoDir)
        : Form("%s/iter_%d/iter_%d_estimatedparameters.root", recoDir, it, it);
    TFile* fs = TFile::Open(sp);
    if (!fs || fs->IsZombie()) { std::cerr << "  (no seed file " << sp << ")\n"; continue; }
    TTree* ts = (TTree*)fs->Get("estimatedparams");
    if (!ts) { fs->Close(); continue; }
    Int_t sEvt = 0; Bool_t sMatched = false;
    UInt_t svp = 0, svs = 0, spa = 0, sge = 0, ssp = 0;
    ts->SetBranchAddress("event_nr", &sEvt);
    ts->SetBranchAddress("truthMatched", &sMatched);
    ts->SetBranchAddress("particleId_vertex_primary", &svp);
    ts->SetBranchAddress("particleId_vertex_secondary", &svs);
    ts->SetBranchAddress("particleId_particle", &spa);
    ts->SetBranchAddress("particleId_generation", &sge);
    ts->SetBranchAddress("particleId_sub_particle", &ssp);
    for (Long64_t i = 0; i < ts->GetEntries(); ++i) {
      ts->GetEntry(i);
      if (!sMatched) continue;
      seededByEvent[(std::uint32_t)sEvt].insert(keyOf(svp, svs, spa, sge, ssp));
    }
    fs->Close();
  }
  const bool haveSeeds = !seededByEvent.empty();

  // ---- histograms ----
  gStyle->SetOptStat(0);
  TH1F* hMassTruth = new TH1F("hMassTruth", ";m(#pi^{+}#pi^{-}) truth [GeV/c^{2}];K^{0}_{S}",
                              200, 0.40, 0.60);
  TH1F* hMassReco  = new TH1F("hMassReco", ";m(#pi^{+}#pi^{-}) reco [GeV/c^{2}];segment pairs",
                              200, 0.40, 0.60);
  TH1F* hMassRecoBest = new TH1F("hMassRecoBest",
                                 ";m(#pi^{+}#pi^{-}) reco [GeV/c^{2}];K^{0}_{S}",
                                 200, 0.40, 0.60);  // highest-pT segment per pion
  TH1F* hNSeg      = new TH1F("hNSeg", ";# reco segments per K^{0}_{S} daughter #pi;daughter pions",
                              11, -0.5, 10.5);
  TH1F* hPrel      = new TH1F("hPrel", ";p_{segment} / p_{max} (pions with >1 segment);segments",
                              55, 0.0, 1.1);
  // K0S reco efficiency vs decay radius: denominator = all K0S->pi+pi-,
  // numerator = both daughter pions have >=1 reconstructed segment (>= nMinMeas).
  TH1F* hDenR = new TH1F("hDenR", ";K^{0}_{S} decay radius [mm];K^{0}_{S}", effRBins, 0.0, effRMax);
  TH1F* hNumR = new TH1F("hNumR", "", effRBins, 0.0, effRMax);
  // single daughter-pion reco efficiency vs decay radius
  TH1F* hDenPiR = new TH1F("hDenPiR", "", effRBins, 0.0, effRMax);
  TH1F* hNumPiR = new TH1F("hNumPiR", "", effRBins, 0.0, effRMax);
  for (auto* h : {hMassTruth, hMassReco, hMassRecoBest, hNSeg, hPrel, hDenR, hNumR, hDenPiR, hNumPiR})
    h->SetDirectory(nullptr);

  // pion reco efficiency vs delta-phi(pion, K0S) in decay-radius bins
  const int NRB = 6;
  const double rbLo[NRB] = {0, 10, 20, 40, 60, 100};
  const double rbHi[NRB] = {10, 20, 40, 60, 100, 200};
  const int    rbCol[NRB] = {kBlue + 1, kCyan + 2, kGreen + 2, kOrange + 7, kMagenta + 1, kRed + 1};
  auto radBin = [&](double r) -> int {
    for (int b = 0; b < NRB; ++b) if (r >= rbLo[b] && r < rbHi[b]) return b;
    return -1;
  };
  // hDenPhi = truth, hNumPhi = reconstructed (track), hSeedPhi = seeded
  std::vector<TH1F*> hDenPhi(NRB), hNumPhi(NRB), hSeedPhi(NRB);
  for (int b = 0; b < NRB; ++b) {
    hDenPhi[b]  = new TH1F(Form("hDenPhi_%d", b),  "", kNPhiBins, -TMath::Pi(), TMath::Pi());
    hNumPhi[b]  = new TH1F(Form("hNumPhi_%d", b),  "", kNPhiBins, -TMath::Pi(), TMath::Pi());
    hSeedPhi[b] = new TH1F(Form("hSeedPhi_%d", b), "", kNPhiBins, -TMath::Pi(), TMath::Pi());
    hDenPhi[b]->SetDirectory(nullptr); hNumPhi[b]->SetDirectory(nullptr); hSeedPhi[b]->SetDirectory(nullptr);
  }
  // pion seeding-efficiency numerator vs decay radius (has a truth-matched seed)
  TH1F* hSeedPiR = new TH1F("hSeedPiR", "", effRBins, 0.0, effRMax);
  hSeedPiR->SetDirectory(nullptr);

  // delta-phi(pion, K0S) distributions themselves (integrated over radius):
  // truth = all daughter pions, reco = reconstructed daughter pions.
  TH1F* hDphiTruth = new TH1F("hDphiTruth", ";#Delta#phi(#pi, K^{0}_{S}) [rad];daughter pions",
                              kNPhiBins, -TMath::Pi(), TMath::Pi());
  TH1F* hDphiReco  = new TH1F("hDphiReco", "", kNPhiBins, -TMath::Pi(), TMath::Pi());
  // K0S truth-level distribution vs decay radius (all K0S decays)
  TH1F* hK0SRadius = new TH1F("hK0SRadius", ";K^{0}_{S} decay radius [mm];K^{0}_{S}",
                              effRBins, 0.0, effRMax);
  // efficiency inputs vs K0S pT and eta (K0S = both pions, Pi = single daughter)
  TH1F* hDenPt   = new TH1F("hDenPt", "", effRBins, 0.0, ptMax);
  TH1F* hNumPt   = new TH1F("hNumPt", "", effRBins, 0.0, ptMax);
  TH1F* hDenPiPt = new TH1F("hDenPiPt", "", effRBins, 0.0, ptMax);
  TH1F* hNumPiPt = new TH1F("hNumPiPt", "", effRBins, 0.0, ptMax);
  TH1F* hDenEta   = new TH1F("hDenEta", "", effRBins, -etaMax, etaMax);
  TH1F* hNumEta   = new TH1F("hNumEta", "", effRBins, -etaMax, etaMax);
  TH1F* hDenPiEta = new TH1F("hDenPiEta", "", effRBins, -etaMax, etaMax);
  TH1F* hNumPiEta = new TH1F("hNumPiEta", "", effRBins, -etaMax, etaMax);
  for (auto* h : {hDphiTruth, hDphiReco, hK0SRadius, hDenPt, hNumPt, hDenPiPt, hNumPiPt,
                  hDenEta, hNumEta, hDenPiEta, hNumPiEta}) h->SetDirectory(nullptr);

  long nK0S = 0, nWith2Pi = 0, nRecoPairs = 0;
  const Long64_t nEnt = std::min(tP->GetEntries(), tV->GetEntries());
  for (Long64_t i = 0; i < nEnt; ++i) {
    tP->GetEntry(i);
    tV->GetEntry(i);

    // barcode -> truth particle info (this event)
    std::map<Key, PInfo> pmap;
    for (size_t k = 0; k < pPdg->size(); ++k)
      pmap[keyOf((*pVp)[k], (*pVs)[k], (*pPa)[k], (*pGe)[k], (*pSp)[k])] =
          {(*pPdg)[k], (*pPx)[k], (*pPy)[k], (*pPz)[k], (*pM)[k]};

    // truth-matched seeds for this event
    const std::set<Key>* seededSet = nullptr;
    { auto sit = seededByEvent.find(pEvt); if (sit != seededByEvent.end()) seededSet = &sit->second; }
    auto isSeeded = [&](const Key& k) { return seededSet && seededSet->count(k) > 0; };

    // reco segments for this event: majority-particle barcode -> list of segments
    std::map<Key, std::vector<Seg>> recoByPart;
    if (haveReco) {
      auto eit = tsEntryByEvent.find(pEvt);
      if (eit != tsEntryByEvent.end()) {
        tR->GetEntry(eit->second);
        for (size_t k = 0; k < rQop->size(); ++k) {
          if ((*rQop)[k] == 0) continue;
          if (rNMeas && (*rNMeas)[k] < (unsigned)nMinMeas) continue;  // min-measurements cut
          const double p = 1.0 / std::abs((*rQop)[k]);
          const double th = (*rTheta)[k], ph = (*rPhi)[k];
          Seg s{p * std::sin(th) * std::cos(ph), p * std::sin(th) * std::sin(ph),
                p * std::cos(th), p};
          recoByPart[keyOf((*rVp)[k], (*rVs)[k], (*rPa)[k], (*rGe)[k], (*rSp)[k])]
              .push_back(s);
        }
      }
    }

    // loop vertices; keep those whose mother is a K0S
    for (size_t j = 0; j < inPa->size(); ++j) {
      bool motherK0S = false;
      const PInfo* momK0S = nullptr;
      for (size_t m = 0; m < (*inPa)[j].size(); ++m) {
        auto it = pmap.find(keyOf((*inVp)[j][m], (*inVs)[j][m], (*inPa)[j][m],
                                  (*inGe)[j][m], (*inSp)[j][m]));
        if (it != pmap.end() && it->second.pdg == kPdgK0S) {
          motherK0S = true; momK0S = &it->second; break;
        }
      }
      if (!motherK0S) continue;
      ++nK0S;

      // K0S decay radius (this vertex' position) and mother kinematics
      const double rDecay = std::sqrt((*vVx)[j]*(*vVx)[j] + (*vVy)[j]*(*vVy)[j]);
      const double pK  = std::sqrt(momK0S->px*momK0S->px + momK0S->py*momK0S->py + momK0S->pz*momK0S->pz);
      const double ptK = std::sqrt(momK0S->px*momK0S->px + momK0S->py*momK0S->py);
      const double etaK = (pK > std::abs(momK0S->pz))
                          ? 0.5 * std::log((pK + momK0S->pz) / (pK - momK0S->pz)) : 0.0;
      hK0SRadius->Fill(rDecay);  // truth-level K0S vs decay radius (all K0S)

      // find the pi+ / pi- daughters (barcode + truth info)
      const PInfo *tPiP = nullptr, *tPiM = nullptr;
      Key kPiP{}, kPiM{};
      for (size_t d = 0; d < (*ouPa)[j].size(); ++d) {
        Key kd = keyOf((*ouVp)[j][d], (*ouVs)[j][d], (*ouPa)[j][d],
                       (*ouGe)[j][d], (*ouSp)[j][d]);
        auto it = pmap.find(kd);
        if (it == pmap.end()) continue;
        if (it->second.pdg == kPdgPiP && !tPiP) { tPiP = &it->second; kPiP = kd; }
        if (it->second.pdg == kPdgPiM && !tPiM) { tPiM = &it->second; kPiM = kd; }
      }
      if (!tPiP || !tPiM) continue;  // not the pi+ pi- channel
      ++nWith2Pi;

      // ---- TRUTH mass ----
      auto Etruth = [&](const PInfo* p) {
        return std::sqrt(p->px*p->px + p->py*p->py + p->pz*p->pz + p->m*p->m);
      };
      {
        const double E = Etruth(tPiP) + Etruth(tPiM);
        const double Px = tPiP->px + tPiM->px, Py = tPiP->py + tPiM->py, Pz = tPiP->pz + tPiM->pz;
        const double m2 = E*E - (Px*Px + Py*Py + Pz*Pz);
        if (m2 > 0) hMassTruth->Fill(std::sqrt(m2));
      }

      // efficiency denominators for the K0S (both pions) case
      hDenR->Fill(rDecay);
      hDenPt->Fill(ptK);
      hDenEta->Fill(etaK);

      // ---- RECO: segments of each daughter pion ----
      if (!haveReco) continue;
      const std::vector<Seg> emptySeg;
      auto sit1 = recoByPart.find(kPiP); const auto& segP = (sit1 != recoByPart.end()) ? sit1->second : emptySeg;
      auto sit2 = recoByPart.find(kPiM); const auto& segM = (sit2 != recoByPart.end()) ? sit2->second : emptySeg;

      // K0S reconstructed if BOTH daughter pions have >=1 reco segment
      if (!segP.empty() && !segM.empty()) {
        hNumR->Fill(rDecay); hNumPt->Fill(ptK); hNumEta->Fill(etaK);
      }

      // per-daughter efficiency vs decay radius / K0S pT / K0S eta, and the
      // delta-phi(pion, K0S) distributions.
      const double phiMom = std::atan2(momK0S->py, momK0S->px);
      const int rbin = radBin(rDecay);
      struct Dau { const PInfo* t; bool reco; bool seeded; };
      const Dau daus[2] = {{tPiP, !segP.empty(), isSeeded(kPiP)},
                           {tPiM, !segM.empty(), isSeeded(kPiM)}};
      for (const auto& d : daus) {
        hDenPiR->Fill(rDecay); hDenPiPt->Fill(ptK); hDenPiEta->Fill(etaK);
        if (d.reco)   { hNumPiR->Fill(rDecay); hNumPiPt->Fill(ptK); hNumPiEta->Fill(etaK); }
        if (d.seeded) hSeedPiR->Fill(rDecay);
        double dphi = std::atan2(d.t->py, d.t->px) - phiMom;
        while (dphi >  TMath::Pi()) dphi -= 2 * TMath::Pi();
        while (dphi <= -TMath::Pi()) dphi += 2 * TMath::Pi();
        hDphiTruth->Fill(dphi);
        if (d.reco) hDphiReco->Fill(dphi);
        if (rbin >= 0) {
          hDenPhi[rbin]->Fill(dphi);
          if (d.reco)   hNumPhi[rbin]->Fill(dphi);
          if (d.seeded) hSeedPhi[rbin]->Fill(dphi);
        }
      }

      // number of segments per daughter pion
      hNSeg->Fill((double)segP.size());
      hNSeg->Fill((double)segM.size());

      // p_segment / max_p for pions with >1 segment
      auto fillPrel = [&](const std::vector<Seg>& segs) {
        if (segs.size() < 2) return;
        double pmax = 0.0;
        for (const auto& s : segs) pmax = std::max(pmax, s.p);
        if (pmax <= 0) return;
        for (const auto& s : segs) hPrel->Fill(s.p / pmax);
      };
      fillPrel(segP);
      fillPrel(segM);

      // K0S mass from ALL segment combinations (every pi+ seg x every pi- seg)
      for (const auto& a : segP) {
        const double Ea = std::sqrt(a.p*a.p + kMassPi*kMassPi);
        for (const auto& b : segM) {
          const double Eb = std::sqrt(b.p*b.p + kMassPi*kMassPi);
          const double Px = a.px + b.px, Py = a.py + b.py, Pz = a.pz + b.pz;
          const double m2 = (Ea+Eb)*(Ea+Eb) - (Px*Px + Py*Py + Pz*Pz);
          if (m2 > 0) { hMassReco->Fill(std::sqrt(m2)); ++nRecoPairs; }
        }
      }

      // K0S mass using only the highest-pT segment of each daughter pion
      if (!segP.empty() && !segM.empty()) {
        auto bestByPt = [](const std::vector<Seg>& segs) -> const Seg* {
          const Seg* best = &segs.front();
          double ptmax2 = best->px*best->px + best->py*best->py;
          for (const auto& s : segs) {
            const double pt2 = s.px*s.px + s.py*s.py;
            if (pt2 > ptmax2) { ptmax2 = pt2; best = &s; }
          }
          return best;
        };
        const Seg* a = bestByPt(segP);
        const Seg* b = bestByPt(segM);
        const double Ea = std::sqrt(a->p*a->p + kMassPi*kMassPi);
        const double Eb = std::sqrt(b->p*b->p + kMassPi*kMassPi);
        const double Px = a->px + b->px, Py = a->py + b->py, Pz = a->pz + b->pz;
        const double m2 = (Ea+Eb)*(Ea+Eb) - (Px*Px + Py*Py + Pz*Pz);
        if (m2 > 0) hMassRecoBest->Fill(std::sqrt(m2));
      }
    }
  }

  std::cout << "K0S decay vertices found          : " << nK0S << "\n"
            << "  ... with pi+ pi- daughters      : " << nWith2Pi << "\n"
            << "reco segment pairs (all combos)   : " << nRecoPairs << "\n";

  // ---- fit the truth peak (for reference) ----
  double fitMean = kMassK0SPdg, fitSigma = 0;
  if (hMassTruth->GetEntries() > 20) {
    TF1 g("g", "gaus", 0.40, 0.60);
    hMassTruth->Fit(&g, "QN0", "", 0.46, 0.54);
    double m = g.GetParameter(1), s = g.GetParameter(2);
    hMassTruth->Fit(&g, "QN0", "", m - 3*s, m + 3*s);
    fitMean = g.GetParameter(1); fitSigma = g.GetParameter(2);
  }
  std::cout.setf(std::ios::fixed); std::cout.precision(3);
  std::cout << "truth m(pipi) peak (fit)          : " << 1e3*fitMean
            << " +/- " << 1e3*fitSigma << " MeV\n"
            << "K0S PDG mass                      : " << 1e3*kMassK0SPdg << " MeV\n";

  // ---- draw ----
  auto drawPdg = [&](TH1F* h, int col) {
    h->SetLineColor(col); h->SetLineWidth(2); h->SetMinimum(0.0);
    h->Draw("HIST E");
  };
  TCanvas* c = new TCanvas("cK0S", "K0S analysis", 1200, 850);
  c->Divide(2, 2);

  c->cd(1); gPad->SetGridx(); gPad->SetGridy();
  drawPdg(hMassTruth, kBlue + 1);
  { TLine* l = new TLine(kMassK0SPdg, 0, kMassK0SPdg, 1.05*hMassTruth->GetMaximum());
    l->SetLineColor(kGreen + 2); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw(); }
  { TLatex t; t.SetNDC(); t.SetTextSize(0.04);
    t.DrawLatex(0.15, 0.85, "truth m(#pi^{+}#pi^{-})"); }

  c->cd(2); gPad->SetGridx(); gPad->SetGridy();
  { const double ymax2 = std::max(hMassReco->GetMaximum(), hMassRecoBest->GetMaximum());
    hMassReco->SetLineColor(kRed + 1); hMassReco->SetLineWidth(2);
    hMassReco->SetMinimum(0.0); hMassReco->SetMaximum(1.15 * std::max(1.0, ymax2));
    hMassReco->Draw("HIST E");
    hMassRecoBest->SetLineColor(kBlue + 1); hMassRecoBest->SetLineWidth(2);
    hMassRecoBest->Draw("HIST E SAME");
    TLine* l = new TLine(kMassK0SPdg, 0, kMassK0SPdg, 1.05*std::max(1.0, ymax2));
    l->SetLineColor(kGreen + 2); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw();
    TLegend* leg = new TLegend(0.14, 0.73, 0.58, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(hMassReco, "all segment pairs", "l");
    leg->AddEntry(hMassRecoBest, "highest-p_{T} segment / pion", "l");
    leg->Draw(); }

  c->cd(3); gPad->SetGridx(); gPad->SetGridy(); gPad->SetLogy();
  hNSeg->SetLineColor(kBlack); hNSeg->SetLineWidth(2); hNSeg->Draw("HIST");

  c->cd(4); gPad->SetGridx(); gPad->SetGridy();
  hPrel->SetLineColor(kMagenta + 1); hPrel->SetLineWidth(2); hPrel->SetMinimum(0.0);
  hPrel->Draw("HIST");

  c->SaveAs(Form("%s/K0S_analysis.pdf", recoDir));
  c->SaveAs(Form("%s/K0S_analysis.png", recoDir));

  // ---- canvas: efficiency vs decay radius (K0S track, pion track, pion seed) ----
  TH1F* hEffR = (TH1F*)hNumR->Clone("hEffR"); hEffR->SetDirectory(nullptr);
  hEffR->Divide(hNumR, hDenR, 1.0, 1.0, "B");
  hEffR->SetTitle(";decay radius [mm];efficiency");
  hEffR->SetMinimum(0.0); hEffR->SetMaximum(1.05);
  hEffR->SetLineColor(kBlue + 1); hEffR->SetMarkerColor(kBlue + 1); hEffR->SetMarkerStyle(20); hEffR->SetLineWidth(2);
  TH1F* hEffPiR = (TH1F*)hNumPiR->Clone("hEffPiR"); hEffPiR->SetDirectory(nullptr);
  hEffPiR->Divide(hNumPiR, hDenPiR, 1.0, 1.0, "B");
  hEffPiR->SetLineColor(kRed + 1); hEffPiR->SetMarkerColor(kRed + 1); hEffPiR->SetMarkerStyle(24); hEffPiR->SetLineWidth(2);
  TH1F* hEffSeedPiR = (TH1F*)hSeedPiR->Clone("hEffSeedPiR"); hEffSeedPiR->SetDirectory(nullptr);
  hEffSeedPiR->Divide(hSeedPiR, hDenPiR, 1.0, 1.0, "B");
  hEffSeedPiR->SetLineColor(kGreen + 2); hEffSeedPiR->SetMarkerColor(kGreen + 2); hEffSeedPiR->SetMarkerStyle(25); hEffSeedPiR->SetLineWidth(2);

  TCanvas* ce = new TCanvas("cK0SEff", "efficiency vs radius", 900, 650);
  ce->SetGridx(); ce->SetGridy();
  hEffR->Draw("E1"); hEffSeedPiR->Draw("E1 SAME"); hEffPiR->Draw("E1 SAME");
  { TLegend* leg = new TLegend(0.44, 0.18, 0.88, 0.38);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(hEffR, Form("K^{0}_{S} (both #pi track, #geq%d meas)", nMinMeas), "lp");
    leg->AddEntry(hEffPiR, Form("single #pi track (#geq%d meas)", nMinMeas), "lp");
    leg->AddEntry(hEffSeedPiR, "single #pi seed", "lp");
    leg->Draw(); }
  ce->SaveAs(Form("%s/K0S_eff_vs_r.pdf", recoDir));
  ce->SaveAs(Form("%s/K0S_eff_vs_r.png", recoDir));

  // ---- canvas: pion efficiency vs delta-phi(pion, K0S), one pad per radius bin ----
  //      per pad: track reco efficiency (red full) + seeding efficiency (green open).
  const int ncol = 3, nrow = (NRB + ncol - 1) / ncol;
  std::vector<TH1F*> hEffPhi(NRB, nullptr), hEffSeedPhi(NRB, nullptr);
  TCanvas* cphi = new TCanvas("cPiEffPhi", "pion eff vs dphi in r bins", 1200, 700);
  cphi->Divide(ncol, nrow);
  for (int b = 0; b < NRB; ++b) {
    cphi->cd(b + 1); gPad->SetGridx(); gPad->SetGridy();
    hEffPhi[b] = (TH1F*)hNumPhi[b]->Clone(Form("hEffPhi_%d", b)); hEffPhi[b]->SetDirectory(nullptr);
    hEffPhi[b]->Divide(hNumPhi[b], hDenPhi[b], 1.0, 1.0, "B");
    hEffPhi[b]->SetTitle(Form("r #in [%g, %g] mm;#Delta#phi(#pi, K^{0}_{S}) [rad];#pi efficiency", rbLo[b], rbHi[b]));
    hEffPhi[b]->SetMinimum(0.0); hEffPhi[b]->SetMaximum(1.05);
    hEffPhi[b]->SetLineColor(kRed + 1); hEffPhi[b]->SetMarkerColor(kRed + 1); hEffPhi[b]->SetMarkerStyle(20); hEffPhi[b]->SetLineWidth(2);
    hEffPhi[b]->Draw("E1");
    hEffSeedPhi[b] = (TH1F*)hSeedPhi[b]->Clone(Form("hEffSeedPhi_%d", b)); hEffSeedPhi[b]->SetDirectory(nullptr);
    hEffSeedPhi[b]->Divide(hSeedPhi[b], hDenPhi[b], 1.0, 1.0, "B");
    hEffSeedPhi[b]->SetLineColor(kGreen + 2); hEffSeedPhi[b]->SetMarkerColor(kGreen + 2); hEffSeedPhi[b]->SetMarkerStyle(25); hEffSeedPhi[b]->SetLineWidth(2);
    hEffSeedPhi[b]->Draw("E1 SAME");
    if (b == 0) { TLegend* lg = new TLegend(0.15, 0.15, 0.55, 0.32);
      lg->SetBorderSize(0); lg->SetFillStyle(0);
      lg->AddEntry(hEffPhi[b], "track", "lp"); lg->AddEntry(hEffSeedPhi[b], "seed", "lp"); lg->Draw(); }
  }
  cphi->SaveAs(Form("%s/K0S_pi_eff_vs_dphi.pdf", recoDir));
  cphi->SaveAs(Form("%s/K0S_pi_eff_vs_dphi.png", recoDir));

  // ---- canvas: delta-phi(pion, K0S) distributions, one pad per radius bin ----
  //      truth (black open) and reconstructed (red full).
  TCanvas* cdphi = new TCanvas("cDphi", "dphi distributions", 1200, 700);
  cdphi->Divide(ncol, nrow);
  for (int b = 0; b < NRB; ++b) {
    cdphi->cd(b + 1); gPad->SetGridx(); gPad->SetGridy();
    hDenPhi[b]->SetTitle(Form("r #in [%g, %g] mm;#Delta#phi(#pi, K^{0}_{S}) [rad];daughter pions", rbLo[b], rbHi[b]));
    hDenPhi[b]->SetMinimum(0.0); hDenPhi[b]->SetMaximum(1.2 * std::max(1.0, hDenPhi[b]->GetMaximum()));
    hDenPhi[b]->SetLineColor(kBlack); hDenPhi[b]->SetMarkerColor(kBlack); hDenPhi[b]->SetMarkerStyle(24); hDenPhi[b]->SetLineWidth(2);
    hDenPhi[b]->Draw("E1");
    hNumPhi[b]->SetLineColor(kRed + 1); hNumPhi[b]->SetMarkerColor(kRed + 1); hNumPhi[b]->SetMarkerStyle(20); hNumPhi[b]->SetLineWidth(2);
    hNumPhi[b]->Draw("E1 SAME");
    if (b == 0) { TLegend* lg = new TLegend(0.15, 0.72, 0.55, 0.88);
      lg->SetBorderSize(0); lg->SetFillStyle(0);
      lg->AddEntry(hDenPhi[b], "truth", "lp"); lg->AddEntry(hNumPhi[b], "reco", "lp"); lg->Draw(); }
  }
  cdphi->SaveAs(Form("%s/K0S_dphi_dist.pdf", recoDir));
  cdphi->SaveAs(Form("%s/K0S_dphi_dist.png", recoDir));

  // ---- overlay canvas: pion TRACK efficiency vs dphi, one graph per radius bin ----
  std::vector<TH1F*> hEffPhiOv(NRB), hEffSeedPhiOv(NRB);
  TCanvas* cphiOv = new TCanvas("cPiEffPhiOverlay", "pion track eff vs dphi (r bins)", 900, 650);
  cphiOv->SetGridx(); cphiOv->SetGridy();
  TLegend* legTr = new TLegend(0.44, 0.16, 0.88, 0.40);
  legTr->SetBorderSize(0); legTr->SetFillStyle(0); legTr->SetHeader("decay radius (track)");
  for (int b = 0; b < NRB; ++b) {
    hEffPhiOv[b] = (TH1F*)hEffPhi[b]->Clone(Form("hEffPhiOv_%d", b)); hEffPhiOv[b]->SetDirectory(nullptr);
    hEffPhiOv[b]->SetTitle(";#Delta#phi(#pi, K^{0}_{S}) [rad];#pi track efficiency");
    hEffPhiOv[b]->SetMinimum(0.0); hEffPhiOv[b]->SetMaximum(1.05);
    hEffPhiOv[b]->SetLineColor(rbCol[b]); hEffPhiOv[b]->SetMarkerColor(rbCol[b]); hEffPhiOv[b]->SetMarkerStyle(20);
    hEffPhiOv[b]->Draw(b == 0 ? "E1" : "E1 SAME");
    legTr->AddEntry(hEffPhiOv[b], Form("r #in [%g, %g] mm", rbLo[b], rbHi[b]), "lp");
  }
  legTr->Draw();
  cphiOv->SaveAs(Form("%s/K0S_pi_eff_vs_dphi_overlay.pdf", recoDir));
  cphiOv->SaveAs(Form("%s/K0S_pi_eff_vs_dphi_overlay.png", recoDir));

  // ---- overlay canvas: pion SEEDING efficiency vs dphi, one graph per radius bin ----
  TCanvas* cseedOv = new TCanvas("cPiSeedEffPhiOverlay", "pion seed eff vs dphi (r bins)", 900, 650);
  cseedOv->SetGridx(); cseedOv->SetGridy();
  TLegend* legSd = new TLegend(0.44, 0.16, 0.88, 0.40);
  legSd->SetBorderSize(0); legSd->SetFillStyle(0); legSd->SetHeader("decay radius (seed)");
  for (int b = 0; b < NRB; ++b) {
    hEffSeedPhiOv[b] = (TH1F*)hEffSeedPhi[b]->Clone(Form("hEffSeedPhiOv_%d", b)); hEffSeedPhiOv[b]->SetDirectory(nullptr);
    hEffSeedPhiOv[b]->SetTitle(";#Delta#phi(#pi, K^{0}_{S}) [rad];#pi seeding efficiency");
    hEffSeedPhiOv[b]->SetMinimum(0.0); hEffSeedPhiOv[b]->SetMaximum(1.05);
    hEffSeedPhiOv[b]->SetLineColor(rbCol[b]); hEffSeedPhiOv[b]->SetMarkerColor(rbCol[b]); hEffSeedPhiOv[b]->SetMarkerStyle(20);
    hEffSeedPhiOv[b]->Draw(b == 0 ? "E1" : "E1 SAME");
    legSd->AddEntry(hEffSeedPhiOv[b], Form("r #in [%g, %g] mm", rbLo[b], rbHi[b]), "lp");
  }
  legSd->Draw();
  cseedOv->SaveAs(Form("%s/K0S_pi_seedeff_vs_dphi_overlay.pdf", recoDir));
  cseedOv->SaveAs(Form("%s/K0S_pi_seedeff_vs_dphi_overlay.png", recoDir));

  // ---- canvas: truth-level K0S vs decay radius ----
  TCanvas* crad = new TCanvas("cK0SRadius", "K0S vs decay radius", 900, 650);
  crad->SetGridx(); crad->SetGridy();
  hK0SRadius->SetLineColor(kBlue + 1); hK0SRadius->SetLineWidth(2); hK0SRadius->SetMinimum(0.0);
  hK0SRadius->Draw("HIST E");
  crad->SaveAs(Form("%s/K0S_vs_radius.pdf", recoDir));
  crad->SaveAs(Form("%s/K0S_vs_radius.png", recoDir));

  // ---- canvases: efficiency vs K0S pT and vs eta (K0S + single pion) ----
  std::vector<TObject*> gExtra;  // extra eff hists + canvases to write
  auto makeEffCanvas = [&](TH1F* numK, TH1F* denK, TH1F* numPi, TH1F* denPi,
                           const char* cname, const char* xt, const char* fbase) {
    TH1F* eK = (TH1F*)numK->Clone(Form("%s_effK", cname)); eK->SetDirectory(nullptr);
    eK->Divide(numK, denK, 1.0, 1.0, "B");
    eK->SetTitle(Form(";%s;reconstruction efficiency", xt));
    eK->SetMinimum(0.0); eK->SetMaximum(1.05);
    eK->SetLineColor(kBlue + 1); eK->SetMarkerColor(kBlue + 1); eK->SetMarkerStyle(20); eK->SetLineWidth(2);
    TH1F* ePi = (TH1F*)numPi->Clone(Form("%s_effPi", cname)); ePi->SetDirectory(nullptr);
    ePi->Divide(numPi, denPi, 1.0, 1.0, "B");
    ePi->SetLineColor(kRed + 1); ePi->SetMarkerColor(kRed + 1); ePi->SetMarkerStyle(24); ePi->SetLineWidth(2);
    TCanvas* cc = new TCanvas(cname, cname, 900, 650);
    cc->SetGridx(); cc->SetGridy();
    eK->Draw("E1"); ePi->Draw("E1 SAME");
    TLegend* lg = new TLegend(0.50, 0.20, 0.88, 0.34); lg->SetBorderSize(0); lg->SetFillStyle(0);
    lg->AddEntry(eK, Form("K^{0}_{S} (both #pi, #geq%d meas)", nMinMeas), "lp");
    lg->AddEntry(ePi, Form("single #pi (#geq%d meas)", nMinMeas), "lp");
    lg->Draw();
    cc->SaveAs(Form("%s/%s.pdf", recoDir, fbase));
    cc->SaveAs(Form("%s/%s.png", recoDir, fbase));
    gExtra.push_back(eK); gExtra.push_back(ePi); gExtra.push_back(cc);
  };
  makeEffCanvas(hNumPt, hDenPt, hNumPiPt, hDenPiPt, "cK0SEffPt", "K^{0}_{S} p_{T} [GeV]", "K0S_eff_vs_pt");
  makeEffCanvas(hNumEta, hDenEta, hNumPiEta, hDenPiEta, "cK0SEffEta", "K^{0}_{S} #eta", "K0S_eff_vs_eta");

  TFile* fout = TFile::Open(Form("%s/K0S_analysis.root", recoDir), "RECREATE");
  if (fout && !fout->IsZombie()) {
    hMassTruth->Write(); hMassReco->Write(); hMassRecoBest->Write();
    hNSeg->Write(); hPrel->Write();
    hDenR->Write(); hNumR->Write(); hEffR->Write();
    hDenPiR->Write(); hNumPiR->Write(); hEffPiR->Write();
    hSeedPiR->Write(); hEffSeedPiR->Write();
    for (int b = 0; b < NRB; ++b) {
      hDenPhi[b]->Write(); hNumPhi[b]->Write(); hEffPhi[b]->Write();
      hSeedPhi[b]->Write(); hEffSeedPhi[b]->Write();
    }
    hDphiTruth->Write(); hDphiReco->Write(); hK0SRadius->Write();
    hDenPt->Write(); hNumPt->Write(); hDenPiPt->Write(); hNumPiPt->Write();
    hDenEta->Write(); hNumEta->Write(); hDenPiEta->Write(); hNumPiEta->Write();
    for (auto* o : gExtra) o->Write();
    c->Write("cK0S"); ce->Write("cK0SEff"); cphi->Write("cPiEffPhi");
    cphiOv->Write("cPiEffPhiOverlay"); cseedOv->Write("cPiSeedEffPhiOverlay");
    cdphi->Write("cDphi"); crad->Write("cK0SRadius");
    fout->Close();
  }
  std::cout << "Saved K0S_analysis.*, K0S_eff_vs_r/pt/eta.*, K0S_pi_eff_vs_dphi.*, "
               "K0S_dphi_dist.*, K0S_vs_radius.* in " << recoDir << std::endl;
}
