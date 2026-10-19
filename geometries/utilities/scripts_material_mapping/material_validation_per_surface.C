// ---------------------------------------------------------------------------
//  material_validation_per_surface.C
//
//  Per-surface material-validation plots for the ALICE 3 tracking geometry
//  (Gen1 or Gen3): binning fidelity of the material map. Compares, per surface
//  and globally,
//    A = binned map   (propagation-material*.root: the assigner reading the map
//                       applied to the tracking geometry), against
//    B = assigned map (material_*_mapped.root: the real geant material assigned
//                       to surfaces by the mapping step),
//  i.e. "does the binned map reproduce the material that was assigned to it?"
//  (ratio A/B). The track samples differ (validation gun vs recording
//  geantinos), so the comparison is valid in the mean / per surface.
//
//  Self-contained (no geometry-map.json): surface name/position come from the
//  GeometryIdentifier and a volume->region map selected by `geomGeneration`
//  ("gen1"/"gen3"); cyl/disc is data-driven from the hit-position spread.
//
//  Args: valid_file, mapped_file, outdir="material_mapping_validation",
//        geomGeneration="gen3",
//        drawTotalX=false (also draw total-X/X0 vs eta/phi),
//        storeHists=false (also write outdir/surfaces/hists.root),
//        zBins=800, rBins=200, phiBins=72, etaBins=84 (bins per axis in the
//        1D/2D histograms), nbprocess, treeName.
//
//  All output goes to <outdir>/surfaces/: map_rz_A/B/ratio.pdf,
//  map_etaphi_A/B/ratio.pdf (2D), per_surface_X0.pdf, [compare_vs_eta/phi.pdf if
//  drawTotalX], and per surface <name>_<id>_1D.pdf (vs z|r, +ratio), _phi.pdf
//  (vs phi, +ratio), _eta.pdf (vs eta, +ratio), _2D_A/B/ratio.pdf (2D),
//  _1Dratio.pdf / _phiratio.pdf
//  (per-cell 2D ratio averaged over cells passing a relative-error cut).
//
//  Usage:
//    root -l -b 'material_validation_per_surface.C(
//      "propagation-material.root","material_mapped.root",
//      "material_mapping_validation","gen3")'
// ---------------------------------------------------------------------------
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>

namespace {

const double kPi = 3.14159265358979323846;

// 2D-ratio cells whose relative statistical uncertainty exceeds this are dropped
// from the 1D ratio projections; the surviving cells are plain-averaged. Tune to
// taste (smaller = stricter).
const double kRatioMaxRelErr = 0.2;

// ACTS GeometryIdentifier bit layout.
std::uint64_t gidVolume(std::uint64_t id) { return (id >> 56) & 0xffULL; }
std::uint64_t gidLayer(std::uint64_t id) { return (id >> 36) & 0xfffULL; }

// Set from the geomGeneration argument; selects the volume->region name map.
bool gIsGen1 = false;

// Number of bins along each axis, set from the zBins/rBins/phiBins arguments.
//   z  : cylinder coord axis + the r-z map z-axis
//   r  : disc coord axis + the r-z map r-axis
//   phi: per-surface 1D/2D phi + the eta-phi map phi-axis
int gZBins = 800;
int gRBins = 200;
int gPhiBins = 72;
int gEtaBins = 84;  // per-surface X/X0 vs eta

// Volume-id -> region name. Gen3: kVol* in GeometryGen3_*.h. Gen1: volume IDs of
// the buildALICE3Geometry TrackingGeometry (from the Gen1 digi NAME fields).
std::string regionName(std::uint64_t v) {
  if (gIsGen1) {
    switch (v) {
      case 5:  return "OuterMostEndCapNeg";
      case 14: return "InnerPixelsBarrel";
      case 17: return "OuterPixelsNeg";
      case 18: return "OuterPixelsBarrel";
      case 19: return "OuterPixelsPos";
      case 21: return "OuterMostBarrel";
      case 22: return "OuterMostEndCapPos";
      case 24: return "oTOF";
      default: return "vol" + std::to_string(v);
    }
  }
  switch (v) {  // gen3
    case 20: return "OuterDisksNeg";
    case 30: return "MiddleDisksNeg";
    case 40: return "InnerBarrel";
    case 50: return "OuterTrackerBarrel";
    case 60: return "MiddleDisksPos";
    case 70: return "OuterDisksPos";
    case 90: return "OTOF";
    default: return "vol" + std::to_string(v);
  }
}
// Fallback classification only (data-driven z/r spread overrides it): barrels.
bool isCylVolume(std::uint64_t v) {
  if (gIsGen1) return v == 14 || v == 18 || v == 21 || v == 24;
  return v == 40 || v == 50 || v == 90;
}
std::string surfaceLabel(std::uint64_t id) {
  return regionName(gidVolume(id)) + "_vol" + std::to_string(gidVolume(id)) +
         "_L" + std::to_string(gidLayer(id));
}

struct SurfHist {
  bool isCyl = true;
  TProfile* prof = nullptr;      // <X/X0> vs z (cyl) or r (disc)
  TProfile* profPhi = nullptr;   // <X/X0> vs phi
  TProfile* profEta = nullptr;   // <X/X0> vs eta (of the crossing point)
  TProfile2D* prof2 = nullptr;   // <X/X0> vs (z|r, phi)  -- 2D map of the surface
  double sumX0 = 0.;             // per-track-normalised later
  double sumPos = 0.;            // mean fixed coordinate (r for cyl, z for disc)
  long nhit = 0;
};

struct Dataset {
  std::map<std::uint64_t, SurfHist> surf;
  TProfile* pEta = nullptr;
  TProfile* pPhi = nullptr;
  TProfile2D* pEtaPhi = nullptr;  // <X/X0> per track vs (eta, phi) -- 2D map
  TH2D* hZR = nullptr;
  long nTracks = 0;
  double meanTot = 0.;  // mean total X/X0 per track
  bool ok = false;
};

// Load one material-track file, filling per-track eta/phi profiles, an r-z map,
// and (if the file is surface-resolved) per-surface <X/X0> profiles. Histograms
// are created in `dir` (unique names via `tag`); `dir` must be a writable file.
Dataset loadDataset(const std::string& fname, const std::string& tag,
                    const std::string& treeName, long nbproc, TDirectory* dir) {
  Dataset d;
  TFile* f = TFile::Open(fname.c_str());
  if (!f || f->IsZombie()) {
    std::cerr << "WARN: cannot open " << fname << std::endl;
    return d;
  }
  TTree* t = dynamic_cast<TTree*>(f->Get(treeName.c_str()));
  if (!t) {
    std::cerr << "WARN: no tree '" << treeName << "' in " << fname << std::endl;
    f->Close();
    return d;
  }
  float v_eta = 0.f, v_phi = 0.f;
  std::vector<float>*mx = nullptr, *my = nullptr, *mz = nullptr, *mlen = nullptr,
                    *mx0 = nullptr;
  std::vector<std::uint64_t>* sid = nullptr;
  std::vector<float>*sx = nullptr, *sy = nullptr, *sz = nullptr;
  t->SetBranchAddress("v_eta", &v_eta);
  t->SetBranchAddress("v_phi", &v_phi);
  t->SetBranchAddress("mat_x", &mx);
  t->SetBranchAddress("mat_y", &my);
  t->SetBranchAddress("mat_z", &mz);
  t->SetBranchAddress("mat_step_length", &mlen);
  t->SetBranchAddress("mat_X0", &mx0);
  t->SetBranchAddress("sur_id", &sid);
  t->SetBranchAddress("sur_x", &sx);
  t->SetBranchAddress("sur_y", &sy);
  t->SetBranchAddress("sur_z", &sz);

  dir->cd();  // histograms below belong to the output file, not the input
  d.pEta = new TProfile((tag + "_pEta").c_str(),
                        ";#eta;total X/X_{0} per track", 84, -4.2, 4.2);
  d.pPhi = new TProfile((tag + "_pPhi").c_str(),
                        ";#varphi;total X/X_{0} per track", 72, -kPi, kPi);
  d.hZR = new TH2D((tag + "_hZR").c_str(), ";z [mm];r [mm];<X/X_{0}> per track",
                   gZBins, -4000., 4000., gRBins, 0., 900.);
  d.pEtaPhi = new TProfile2D((tag + "_pEtaPhi").c_str(),
                             ";#eta;#varphi;<X/X_{0}> per track", 84, -4.2, 4.2,
                             gPhiBins, -kPi, kPi);

  Long64_t n = t->GetEntries();
  if (nbproc > 0 && nbproc < n) n = nbproc;
  std::cout << "  reading " << n << " tracks from " << fname << std::endl;

  // Pass 1 (surface branches only): classify each surface as cylinder (z
  // spreads) or disc (r spreads) from the actual hit positions, so surfaces that
  // don't follow their volume's majority type get the right axis. In particular
  // the forward SERVICE SHELLS (TRK_*_Services) are cylinders living in the
  // forward disc volumes (20/70); a volume-based guess would mis-plot them vs r.
  std::map<std::uint64_t, bool> isCylMap;
  if (sid) {
    t->SetBranchStatus("*", 0);
    for (const char* bn : {"sur_id", "sur_x", "sur_y", "sur_z"})
      if (t->GetBranch(bn)) t->SetBranchStatus(bn, 1);
    std::map<std::uint64_t, std::array<double, 5>> cls;  // id->{sz,sz2,sr,sr2,n}
    for (Long64_t i = 0; i < n; ++i) {
      t->GetEntry(i);
      if (!sid || !sx || !sy || !sz) continue;
      for (std::size_t j = 0; j < sid->size(); ++j) {
        const double zz = sz->at(j), rr = std::hypot(sx->at(j), sy->at(j));
        auto& c = cls[sid->at(j)];
        c[0] += zz;
        c[1] += zz * zz;
        c[2] += rr;
        c[3] += rr * rr;
        c[4] += 1.0;
      }
    }
    for (auto& kv : cls) {
      const double nn = kv.second[4];
      const double vz = kv.second[1] / nn - std::pow(kv.second[0] / nn, 2);
      const double vr = kv.second[3] / nn - std::pow(kv.second[2] / nn, 2);
      isCylMap[kv.first] = (vz >= vr);  // cylinder if z spreads more than r
    }
    t->SetBranchStatus("*", 1);  // re-enable all branches for pass 2
  }

  double sumTot = 0.;
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    const std::size_t nStep = mx0->size();
    const bool haveSur = sid && sid->size() >= nStep;
    double tot = 0.;

    // Aggregate steps per surface CROSSING. In the assigned (mapped) file a
    // crossing is split into many fine Geant4 sub-steps sharing one sur_id; the
    // assigner (validation) file lumps them into one. Summing X/X0 per
    // (track, surface) and filling the profile once per crossing makes the two
    // comparable (tracks originate at the IP, so each surface is crossed once).
    std::map<std::uint64_t, std::array<double, 5>> cross;  // id->{sumt,z,r,phi,n}

    for (std::size_t j = 0; j < nStep; ++j) {
      const double X0 = mx0->at(j);
      if (!(X0 > 0.)) continue;
      const double tX0 = mlen->at(j) / X0;  // path-length material (incl. 1/cos)
      tot += tX0;
      d.hZR->Fill(mz->at(j), std::hypot(mx->at(j), my->at(j)), tX0);
      if (haveSur) {
        auto& c = cross[sid->at(j)];
        c[0] += tX0;
        c[1] += sz->at(j);
        c[2] += std::hypot(sx->at(j), sy->at(j));
        c[3] += std::atan2(sy->at(j), sx->at(j));
        c[4] += 1.0;
      }
    }

    for (auto& kv : cross) {
      const std::uint64_t id = kv.first;
      const double sumt = kv.second[0], nsub = kv.second[4];
      const double zc = kv.second[1] / nsub, rc = kv.second[2] / nsub,
                   phic = kv.second[3] / nsub;
      // pseudorapidity of the crossing point
      const double etac =
          (rc > 1e-6) ? std::asinh(zc / rc) : (zc >= 0. ? 10. : -10.);
      SurfHist& s = d.surf[id];
      if (s.prof == nullptr) {
        auto itc = isCylMap.find(id);
        s.isCyl = (itc != isCylMap.end()) ? itc->second
                                          : isCylVolume(gidVolume(id));
        const std::string cax = s.isCyl ? "z [mm]" : "r [mm]";
        const double clo = s.isCyl ? -4000. : 0., chi = s.isCyl ? 4000. : 900.;
        // coord axis: z (cylinder) uses gZBins, r (disc) uses gRBins.
        const int nCoord = s.isCyl ? gZBins : gRBins;
        const std::string nm = tag + "_p_" + std::to_string(id);
        s.prof = new TProfile(
            nm.c_str(),
            (surfaceLabel(id) + ";" + cax + ";X/X_{0} per crossing").c_str(),
            nCoord, clo, chi);
        const std::string npf = tag + "_pphi_" + std::to_string(id);
        s.profPhi = new TProfile(
            npf.c_str(),
            (surfaceLabel(id) + ";#varphi;X/X_{0} per crossing").c_str(),
            gPhiBins, -kPi, kPi);
        const std::string npe = tag + "_peta_" + std::to_string(id);
        s.profEta = new TProfile(
            npe.c_str(),
            (surfaceLabel(id) + ";#eta;X/X_{0} per crossing").c_str(), gEtaBins,
            -4.2, 4.2);
        const std::string n2 = tag + "_p2_" + std::to_string(id);
        s.prof2 = new TProfile2D(
            n2.c_str(),
            (surfaceLabel(id) + ";" + cax + ";#varphi;X/X_{0} per crossing")
                .c_str(),
            nCoord, clo, chi, gPhiBins, -kPi, kPi);
      }
      s.prof->Fill(s.isCyl ? zc : rc, sumt);          // one entry per crossing
      s.profPhi->Fill(phic, sumt);                    // vs phi
      s.profEta->Fill(etac, sumt);                    // vs eta
      s.prof2->Fill(s.isCyl ? zc : rc, phic, sumt);   // 2D map of the surface
      s.sumX0 += sumt;
      s.sumPos += s.isCyl ? rc : zc;  // the ~constant coordinate
      ++s.nhit;                       // now counts crossings, not sub-steps
    }

    d.pEta->Fill(v_eta, tot);
    d.pPhi->Fill(v_phi, tot);
    d.pEtaPhi->Fill(v_eta, v_phi, tot);
    sumTot += tot;
  }
  d.nTracks = n;
  d.meanTot = (n > 0) ? sumTot / n : 0.;
  if (n > 0) d.hZR->Scale(1.0 / n);  // -> <X/X0> per track per cell
  d.ok = true;
  f->Close();
  return d;
}

// Overlay two profiles (A over B) with a ratio A/B pad. profB may be null.
void drawRatio(TProfile* A, TProfile* B, const std::string& xlabel,
               const std::string& ylabel, const std::string& legA,
               const std::string& legB, const std::string& outpath) {
  if (!A) return;
  TCanvas c("c", "cmp", 1000, 850);
  const bool hasB = (B != nullptr);
  TPad* p1 = new TPad("p1", "", 0, hasB ? 0.32 : 0.0, 1, 1);
  p1->SetGridx();
  p1->SetGridy();
  if (hasB) p1->SetBottomMargin(0.02);
  p1->Draw();
  p1->cd();
  A->SetStats(0);
  A->SetMarkerStyle(20);
  A->SetMarkerColor(kRed + 1);
  A->SetLineColor(kRed + 1);
  A->SetTitle((";" + (hasB ? std::string("") : xlabel) + ";" + ylabel).c_str());
  double ymax = A->GetMaximum();
  if (hasB) ymax = std::max(ymax, B->GetMaximum());
  A->SetMaximum(1.30 * ymax);
  A->SetMinimum(0.);
  A->Draw("P HIST");
  if (hasB) {
    B->SetStats(0);
    B->SetLineColor(kBlue + 1);
    B->SetLineWidth(2);
    B->SetMarkerColor(kBlue + 1);
    B->Draw("HIST SAME");
  }
  TLegend* leg = new TLegend(0.32, 0.78, 0.68, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(A, legA.c_str(), "lp");
  if (hasB) leg->AddEntry(B, legB.c_str(), "l");
  leg->Draw();
  if (hasB) {
    c.cd();
    TPad* p2 = new TPad("p2", "", 0, 0.0, 1, 0.32);
    p2->SetGridx();
    p2->SetGridy();
    p2->SetTopMargin(0.03);
    p2->SetBottomMargin(0.30);
    p2->Draw();
    p2->cd();
    TH1D* rm = A->ProjectionX((std::string(A->GetName()) + "_rm").c_str());
    TH1D* rg = B->ProjectionX((std::string(B->GetName()) + "_rg").c_str());
    rm->Divide(rg);
    rm->SetStats(0);
    rm->SetTitle((";" + xlabel + ";A / B").c_str());
    rm->SetMarkerStyle(20);
    rm->SetMarkerColor(kBlack);
    rm->SetLineColor(kBlack);
    rm->GetYaxis()->SetRangeUser(0.0, 2.0);
    rm->GetYaxis()->SetNdivisions(505);
    rm->GetYaxis()->SetLabelSize(0.09);
    rm->GetYaxis()->SetTitleSize(0.10);
    rm->GetYaxis()->SetTitleOffset(0.45);
    rm->GetXaxis()->SetLabelSize(0.09);
    rm->GetXaxis()->SetTitleSize(0.10);
    rm->Draw("P");
  }
  c.Print(outpath.c_str());
}

void drawH2(TH2D* h, const std::string& outpath, double zmax = -1) {
  if (!h) return;
  TCanvas c("c", "h2", 1200, 700);
  c.SetGridx();
  c.SetGridy();
  c.SetRightMargin(0.15);
  if (zmax > 0) h->SetMaximum(zmax);
  h->Draw("COLZ");
  c.Print(outpath.c_str());
}

void drawProf2(TProfile2D* h, const std::string& outpath,
               const std::string& label = "") {
  if (!h) return;
  TCanvas c("c", "p2", 1100, 700);
  c.SetGridx();
  c.SetGridy();
  c.SetRightMargin(0.15);
  h->Draw("COLZ");
  if (!label.empty()) {
    TLatex* t = new TLatex();
    t->SetNDC();
    t->SetTextFont(42);
    t->SetTextSize(0.035);
    t->DrawLatex(0.12, 0.945, label.c_str());
  }
  c.Print(outpath.c_str());
}

// 2D ratio A/B of two profiles: project each to a TH2D of per-bin means, then
// divide bin by bin (analogous to the 1D ProjectionX + Divide). Returns the
// ratio TH2D (0 where B is empty) for further 1D projection.
TH2D* drawProf2Ratio(TProfile2D* A, TProfile2D* B, const std::string& name,
                     const std::string& outpath) {
  if (!A || !B) return nullptr;
  TH2D* a = A->ProjectionXY((name + "_rA").c_str());
  TH2D* b = B->ProjectionXY((name + "_rB").c_str());
  a->Divide(b);
  a->GetZaxis()->SetTitle("A / B");
  a->SetMinimum(0.0);
  a->SetMaximum(2.0);
  TCanvas c("c", "p2r", 1100, 700);
  c.SetGridx();
  c.SetGridy();
  c.SetRightMargin(0.15);
  a->Draw("COLZ");
  c.Print(outpath.c_str());
  return a;
}

// Collapse a 2D ratio map onto one axis by plain-averaging its cells, after
// dropping cells with no entries or with relative statistical error (from
// ProjectionXY + Divide) above kRatioMaxRelErr. Averaging the per-cell A/B ratio
// (each ~1 for a faithful map) is not diluted by partial-phi/partial-z coverage.
// alongX=true -> vs coord; false -> vs phi.
TH1D* ratioProjection(TH2D* ratio, bool alongX, const std::string& name,
                      const std::string& title) {
  if (!ratio) return nullptr;
  const int nx = ratio->GetNbinsX(), ny = ratio->GetNbinsY();
  const int nOut = alongX ? nx : ny, nIn = alongX ? ny : nx;
  TAxis* ax = alongX ? ratio->GetXaxis() : ratio->GetYaxis();
  TH1D* h = new TH1D(name.c_str(), title.c_str(), nOut, ax->GetXmin(),
                     ax->GetXmax());
  for (int a = 1; a <= nOut; ++a) {
    double sum = 0., sum2 = 0.;
    int n = 0;
    for (int b = 1; b <= nIn; ++b) {
      const int ix = alongX ? a : b, iy = alongX ? b : a;
      const int gbin = ratio->GetBin(ix, iy);
      const double r = ratio->GetBinContent(gbin);
      const double e = ratio->GetBinError(gbin);
      if (r == 0. || e <= 0.) continue;       // empty / no error estimate
      if (e / r > kRatioMaxRelErr) continue;  // statistical-significance cut
      sum += r;
      sum2 += r * r;
      ++n;
    }
    if (n > 0) {
      const double mean = sum / n;
      h->SetBinContent(a, mean);
      const double var = std::max(0.0, sum2 / n - mean * mean);
      h->SetBinError(a, n > 1 ? std::sqrt(var / n) : 0.0);
    }
  }
  return h;
}

// Draw a single 1D ratio histogram with a dashed reference line at 1.
void drawRatio1D(TH1D* h, const std::string& xlabel,
                 const std::string& outpath) {
  if (!h) return;
  TCanvas c("c", "r1", 1000, 600);
  c.SetGridx();
  c.SetGridy();
  h->SetStats(0);
  h->SetMarkerStyle(20);
  h->SetMarkerColor(kBlack);
  h->SetLineColor(kBlack);
  h->SetTitle((";" + xlabel + ";A / B (mean over good cells)").c_str());
  h->GetYaxis()->SetRangeUser(0.0, 2.0);
  h->Draw("P");
  TLine* l = new TLine(h->GetXaxis()->GetXmin(), 1.0, h->GetXaxis()->GetXmax(),
                       1.0);
  l->SetLineColor(kRed + 1);
  l->SetLineStyle(2);
  l->Draw();
  c.Print(outpath.c_str());
}

// One validation step: load A (numerator) and B (reference), draw comparisons.
void compareStep(const std::string& outsub, TFile* fout,
                 const std::string& stepTag, const std::string& fileA,
                 const std::string& legA, const std::string& fileB,
                 const std::string& legB, const std::string& treeName,
                 long nbproc, bool doPerSurface, bool drawTotalX) {
  std::cout << "\n===== " << stepTag << "  (" << legA << " / " << legB
            << ") =====" << std::endl;
  gSystem->mkdir(outsub.c_str(), true);

  // Histograms live in the output file (if kept) or the in-memory directory.
  TDirectory* dir = fout ? fout->mkdir(stepTag.c_str()) : gDirectory;
  Dataset A = loadDataset(fileA, stepTag + "_A", treeName, nbproc, dir);
  Dataset B = loadDataset(fileB, stepTag + "_B", treeName, nbproc, dir);
  if (!A.ok) return;

  printf("  tracks: A=%ld  B=%ld   mean X/X0 per track: A=%.4f  B=%.4f  "
         "A/B=%.3f\n",
         A.nTracks, B.nTracks, A.meanTot, B.meanTot,
         B.meanTot > 0 ? A.meanTot / B.meanTot : 0.);

  // total-X/X0-per-track comparison vs eta / phi (optional)
  if (drawTotalX) {
    drawRatio(A.pEta, B.ok ? B.pEta : nullptr, "#eta",
              "total X/X_{0} per track", legA, legB,
              outsub + "/compare_vs_eta.pdf");
    drawRatio(A.pPhi, B.ok ? B.pPhi : nullptr, "#varphi",
              "total X/X_{0} per track", legA, legB,
              outsub + "/compare_vs_phi.pdf");
  }

  // r-z maps + ratio (localises where material differs)
  drawH2(A.hZR, outsub + "/map_rz_A.pdf");
  if (B.ok) {
    drawH2(B.hZR, outsub + "/map_rz_B.pdf");
    TH2D* rz = static_cast<TH2D*>(A.hZR->Clone((stepTag + "_rzRatio").c_str()));
    rz->SetTitle(";z [mm];r [mm];A / B");
    rz->Divide(B.hZR);
    drawH2(rz, outsub + "/map_rz_ratio.pdf", 2.0);
  }

  // 2D eta-phi material maps (per dataset) + ratio
  drawProf2(A.pEtaPhi, outsub + "/map_etaphi_A.pdf", legA);
  if (B.ok) {
    drawProf2(B.pEtaPhi, outsub + "/map_etaphi_B.pdf", legB);
    drawProf2Ratio(A.pEtaPhi, B.pEtaPhi, stepTag + "_etaphi",
                   outsub + "/map_etaphi_ratio.pdf");
  }

  if (!doPerSurface) return;

  // per-surface: A vs B profile + ratio, and per-surface totals bar chart
  TH1D* hA = new TH1D((stepTag + "_surfA").c_str(), ";surface;<X/X_{0}> per track",
                      std::max<int>(1, A.surf.size()), 0,
                      std::max<int>(1, A.surf.size()));
  TH1D* hB = new TH1D((stepTag + "_surfB").c_str(), ";surface;<X/X_{0}> per track",
                      std::max<int>(1, A.surf.size()), 0,
                      std::max<int>(1, A.surf.size()));
  int ib = 1;
  std::cout << "  name                              type    pos[mm]   A<X/X0>  "
               "B<X/X0>  A/B\n";
  for (auto& kv : A.surf) {
    const std::uint64_t id = kv.first;
    SurfHist& sa = kv.second;
    const std::string label = surfaceLabel(id);
    const double posA = sa.nhit ? sa.sumPos / sa.nhit : 0.;
    const double aAvg = sa.nhit ? sa.sumX0 / sa.nhit : 0.;

    TProfile* pb = nullptr;
    double bAvg = 0.;
    auto it = B.surf.find(id);
    if (B.ok && it != B.surf.end()) {
      pb = it->second.prof;
      bAvg = it->second.nhit ? it->second.sumX0 / it->second.nhit : 0.;
    }
    hA->SetBinContent(ib, A.nTracks > 0 ? sa.sumX0 / A.nTracks : 0.);
    hA->GetXaxis()->SetBinLabel(ib, label.c_str());
    if (pb) hB->SetBinContent(ib, B.nTracks > 0 ? it->second.sumX0 / B.nTracks
                                                : 0.);
    ++ib;
    printf("  %-32s  %-4s  %8.1f  %7.4f  %7.4f  %5.2f\n", label.c_str(),
           sa.isCyl ? "cyl" : "disc", posA, aAvg, bAvg,
           bAvg > 0 ? aAvg / bAvg : 0.);

    const std::string cax = sa.isCyl ? "z [mm]" : "r [mm]";
    const std::string sbase = outsub + "/" + label + "_" + std::to_string(id);

    // per-crossing 1D profiles A vs B (+ratio pad)
    drawRatio(sa.prof, pb, cax, "X/X_{0} per crossing", legA, legB,
              sbase + "_1D.pdf");
    TProfile* pbPhi = (B.ok && it != B.surf.end()) ? it->second.profPhi : nullptr;
    drawRatio(sa.profPhi, pbPhi, "#varphi", "X/X_{0} per crossing", legA, legB,
              sbase + "_phi.pdf");
    TProfile* pbEta = (B.ok && it != B.surf.end()) ? it->second.profEta : nullptr;
    drawRatio(sa.profEta, pbEta, "#eta", "X/X_{0} per crossing", legA, legB,
              sbase + "_eta.pdf");

    // per-surface 2D maps of A and B + ratio A/B
    drawProf2(sa.prof2, sbase + "_2D_A.pdf", legA);
    if (B.ok && it != B.surf.end()) {
      drawProf2(it->second.prof2, sbase + "_2D_B.pdf", legB);
      const std::string tg = stepTag + "_s" + std::to_string(id);
      TH2D* r2 = drawProf2Ratio(sa.prof2, it->second.prof2, tg,
                                sbase + "_2D_ratio.pdf");
      // 1D means of the per-cell 2D ratio, averaged over non-empty cells only
      // (not diluted by partial-phi/partial-z: empty cells are excluded).
      drawRatio1D(ratioProjection(r2, true, tg + "_pX",
                                  label + ";" + cax + ";A / B"),
                  cax, sbase + "_1Dratio.pdf");
      drawRatio1D(ratioProjection(r2, false, tg + "_pY",
                                  label + ";#varphi;A / B"),
                  "#varphi", sbase + "_phiratio.pdf");
    }
  }

  TCanvas c("c", "surf", 1500, 750);
  c.SetGridx();
  c.SetGridy();
  c.SetBottomMargin(0.24);
  hA->SetFillColor(kRed - 9);
  hA->SetLineColor(kRed + 1);
  hA->LabelsOption("v");
  hA->SetMaximum(1.25 * std::max(hA->GetMaximum(), hB->GetMaximum()));
  hA->Draw("HIST");
  hB->SetLineColor(kBlue + 1);
  hB->SetLineWidth(2);
  hB->SetFillStyle(0);
  hB->Draw("HIST SAME");
  TLegend* leg = new TLegend(0.5, 0.80, 0.88, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hA, legA.c_str(), "f");
  leg->AddEntry(hB, legB.c_str(), "l");
  leg->Draw();
  c.Print((outsub + "/per_surface_X0.pdf").c_str());
}

}  // namespace

void material_validation_per_surface(
    std::string valid_file = "propagation-material-gen3.root",  // A: binned map
    std::string mapped_file = "material_gen3_mapped.root",      // B: assigned map
    std::string outdir = "material_mapping_validation",
    std::string geomGeneration = "gen3",  // "gen1" or "gen3": volume->region map
    bool drawTotalX = false,   // also draw the total-X/X0 vs eta/phi comparison
    bool storeHists = false,   // also write outdir/surfaces/hists.root
    int zBins = 800,           // number of z bins in the 1D/2D histograms
    int rBins = 200,           // number of r bins
    int phiBins = 72,          // number of phi bins
    int etaBins = 84,          // number of eta bins (per-surface vs eta)
    long nbprocess = -1, std::string treeName = "material_tracks") {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kBird);
  gIsGen1 = (geomGeneration == "gen1");
  gZBins = std::max(1, zBins);
  gRBins = std::max(1, rBins);
  gPhiBins = std::max(1, phiBins);
  gEtaBins = std::max(1, etaBins);

  gSystem->mkdir((outdir + "/surfaces").c_str(), true);
  TFile* fout =
      storeHists
          ? new TFile((outdir + "/surfaces/hists.root").c_str(), "RECREATE")
          : nullptr;

  // Binning fidelity: A/B = binned map / assigned map. Does the binned map
  // reproduce the material assigned to it? Different track samples -> valid in
  // the mean / per surface.
  compareStep(outdir + "/surfaces", fout, "validation", valid_file,
              "binned map (assigner)", mapped_file, "assigned (mapped)",
              treeName, nbprocess, /*doPerSurface=*/true, drawTotalX);

  if (fout) {
    fout->cd();
    fout->Write();
    fout->Close();
  }
  std::cout << "\nDone. Output in '" << outdir << "/surfaces/'." << std::endl;
}
