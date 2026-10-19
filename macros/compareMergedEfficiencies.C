// ============================================================================
//  compareMergedEfficiencies.C
//
//  Overlay the MERGED (total) seeding and tracking efficiencies from several
//  draw_iterative_efficiencies.C ROOT outputs (one per directory) on a single
//  canvas.
//
//  Each ROOT file (efficiencies_vs_pT_eta<...>.root) is expected to contain the
//  TGraphAsymmErrors:
//     gSeedEff_merged   - total seeding efficiency  -> drawn with OPEN   markers
//     gTrkEff_merged    - total tracking efficiency -> drawn with FILLED markers
//  Every directory gets its own colour; open + filled of the same colour are the
//  seeding / tracking efficiency of that directory.
//
//  Usage:   root -l compareMergedEfficiencies.C
//           root -l 'compareMergedEfficiencies.C(true)'   // log pT axis
//  (edit the dir[][2] array below)
// ============================================================================

#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TMath.h>

void compareMergedEfficiencies(bool logx = true)
{
  gStyle->SetOptStat(0);

  // ---- directories (full path to the efficiencies_vs_pT_*.root) + label ----
  TString dir[][2] = {
    {"task_rec_output_A/efficiencies_vs_pT_eta-1-1_logPt.root", "file A"},
    {"task_rec_output_B/efficiencies_vs_pT_eta-1-1_logPt.root", "file B"},
  };
  const int ndir = sizeof(dir) / sizeof(dir[0]);
  const int col[] = {kBlue + 1, kMagenta + 1, kRed, kGreen + 2, kOrange + 7, kAzure + 2, kBlack};

  // Get a graph by name from a file; clone + detach so it survives closing.
  static int gClone = 0;
  auto getG = [&](const TString& path, const char* name) -> TGraphAsymmErrors* {
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "ERROR: cannot open " << path << std::endl;
      return nullptr;
    }
    TGraphAsymmErrors* g = dynamic_cast<TGraphAsymmErrors*>(f->Get(name));
    if (!g) {
      std::cerr << "WARNING: '" << name << "' not in " << path << std::endl;
      f->Close();
      return nullptr;
    }
    g = (TGraphAsymmErrors*)g->Clone(Form("%s_c%d", name, gClone++));
    f->Close();
    return g;
  };

  // ---- read the merged seeding (open) and tracking (filled) graphs ----------
  std::vector<TGraphAsymmErrors*> gSeed(ndir, nullptr), gTrk(ndir, nullptr);
  double xlo = 1e30, xhi = -1e30;
  auto scanX = [&](TGraphAsymmErrors* g) {
    if (!g)
      return;
    for (int i = 0; i < g->GetN(); ++i) {
      double x, y;
      g->GetPoint(i, x, y);
      xlo = std::min(xlo, x - g->GetErrorXlow(i));
      xhi = std::max(xhi, x + g->GetErrorXhigh(i));
    }
  };
  for (int i = 0; i < ndir; ++i) {
    gSeed[i] = getG(dir[i][0], "gSeedEff_merged");
    gTrk[i] = getG(dir[i][0], "gTrkEff_merged");
    scanX(gSeed[i]);
    scanX(gTrk[i]);
  }
  if (xhi <= xlo) {
    std::cerr << "ERROR: no graphs found - nothing to draw" << std::endl;
    return;
  }
  if (logx && xlo <= 0)
    xlo = 0.05; // safe lower edge for a log axis

  // -------------------------------------------------------------- canvas
  TCanvas* c = new TCanvas("cmp_merged_eff", "merged efficiencies", 900, 650);
  c->SetGridx();
  c->SetGridy();
  if (logx)
    c->SetLogx();
  c->SetTicks(1, 1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridColor(kGray + 1);

  TH1F* frame = new TH1F("frame_merged_eff", "", 100, xlo, xhi);
  frame->SetDirectory(nullptr);
  frame->SetTitle("Total (merged) efficiency vs p_{T};p_{T} [GeV/c];Efficiency");
  frame->GetYaxis()->SetRangeUser(0.0, 1.15);
  if (logx) {
    frame->GetXaxis()->SetMoreLogLabels();
    frame->GetXaxis()->SetNoExponent();
  }
  frame->Draw("AXIS");

  TLegend* leg = new TLegend(0.58, 0.16, 0.89, 0.40);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  for (int i = 0; i < ndir; ++i) {
    const int cc = col[i % 7];
    if (gTrk[i]) {
      gTrk[i]->SetLineColor(cc);
      gTrk[i]->SetMarkerColor(cc);
      gTrk[i]->SetMarkerStyle(20);
      gTrk[i]->SetMarkerSize(0.9);
      gTrk[i]->SetLineWidth(2);
      gTrk[i]->Draw("PZ SAME");
    }
    if (gSeed[i]) {
      gSeed[i]->SetLineColor(cc);
      gSeed[i]->SetMarkerColor(cc);
      gSeed[i]->SetMarkerStyle(24);
      gSeed[i]->SetMarkerSize(0.9);
      gSeed[i]->SetLineWidth(2);
      gSeed[i]->Draw("PZ SAME");
    }
    // one legend entry per directory (colour); use the tracking graph if
    // present, otherwise the seeding one.
    TGraphAsymmErrors* rep = gTrk[i] ? gTrk[i] : gSeed[i];
    if (rep)
      leg->AddEntry(rep, dir[i][1], "lp");
  }
  leg->Draw();

  TLatex tl;
  tl.SetNDC();
  tl.SetTextSize(0.03);
  tl.DrawLatex(0.16, 0.85, "filled: tracking eff    open: seeding eff");

  c->RedrawAxis();
  gPad->RedrawAxis("g"); // force the grid lines to be (re)painted on top
  c->Update();
  const char* suf = logx ? "_logPt" : "";
  c->SaveAs(Form("cmp_merged_efficiency%s.pdf", suf));
  c->SaveAs(Form("cmp_merged_efficiency%s.png", suf));
  std::cout << "Saved cmp_merged_efficiency" << suf << ".pdf / .png" << std::endl;
}
