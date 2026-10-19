// ============================================================================
//  compareEffResFakes.C
//
//  Overlay results from several EffResFakes_out.root files (one per directory)
//  on separate canvases. NOTE: the histograms are in EffResFakes_out.root
//  (the "_plots" file is the PDF), so point the paths below at the *_out.root.
//
//  For each comparison, every directory gets its own colour; where two series
//  are shown, the first ("all") is drawn with FILLED markers and the second
//  (physics / clean) with OPEN markers of the same colour.
//
//  Usage:   root -l compareEffResFakes.C
//  (edit the dir[][2] array below)
// ============================================================================

#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>

void compareEffResFakes()
{
  gStyle->SetOptStat(0);

  // ---- directories (full path to the EffResFakes_out.root) + legend label ----
  TString dir[][2] = {

    {"EffResFakes_out_fileA.root", "file A"},
    {"EffResFakes_out_fileB.root", "file B"},
    // ... add more files as needed, each with a label for the legend

  };
  const int ndir = sizeof(dir) / sizeof(dir[0]);
  const int col[] = {kBlue + 1, kMagenta + 1, kRed, kGreen + 1, kOrange + 7, kCyan + 2, kBlack};

  // Get a histogram by name from a file; clone + detach so it survives closing.
  static int gClone = 0;
  auto getH = [&](const TString& path, const char* name) -> TH1* {
    TFile* f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "ERROR: cannot open " << path << std::endl;
      return nullptr;
    }
    TH1* h = dynamic_cast<TH1*>(f->Get(name));
    if (!h) {
      std::cerr << "WARNING: '" << name << "' not in " << path << std::endl;
      f->Close();
      return nullptr;
    }
    h = (TH1*)h->Clone(Form("%s_c%d", name, gClone++));
    h->SetDirectory(nullptr);
    f->Close();
    return h;
  };

  // Collected outputs: every canvas and every drawn histogram, for the
  // multi-page PDF and the ROOT file written at the end.
  std::vector<TCanvas*> gCanvases;
  std::vector<TH1*> gHists;

  // Generic overlay: hMain (filled markers) + optional hOpen (open markers),
  // one colour per directory. perLayer => vertical bin labels + wide canvas.
  auto draw = [&](const char* cname, const char* title, const char* hMain, const char* hOpen,
                  const char* xt, const char* yt, bool logx, double forceMax, bool legBottom,
                  const char* mainTag, const char* openTag, bool perLayer) {
    TCanvas* c = new TCanvas(cname, title, perLayer ? 1100 : 850, 650);
    c->SetGrid();
    if (logx)
      c->SetLogx();
    if (perLayer) {
      c->SetBottomMargin(0.40);
      c->SetGridy();
    }

    std::vector<TH1*> M(ndir, nullptr), O(ndir, nullptr);
    double ymax = 0.0;
    for (int i = 0; i < ndir; ++i) {
      M[i] = getH(dir[i][0], hMain);
      if (hOpen && hOpen[0])
        O[i] = getH(dir[i][0], hOpen);
      if (M[i])
        ymax = std::max(ymax, M[i]->GetMaximum());
      if (O[i])
        ymax = std::max(ymax, O[i]->GetMaximum());
    }
    const double ymaxUse = (forceMax > 0) ? forceMax : 1.25 * ymax;

    TLegend* leg = legBottom ? new TLegend(0.58, 0.16, 0.89, 0.40)
                             : new TLegend(0.58, 0.66, 0.89, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    bool first = true;
    for (int i = 0; i < ndir; ++i) {
      if (!M[i])
        continue;
      const int cc = col[i % 7];
      M[i]->SetLineColor(cc);
      M[i]->SetMarkerColor(cc);
      M[i]->SetMarkerStyle(i == 0 ? 21 : 20);
      M[i]->SetMarkerSize(0.8);
      M[i]->SetLineWidth(2);
      M[i]->SetMinimum(0.0);
      M[i]->SetMaximum(ymaxUse);
      M[i]->SetTitle(Form("%s;%s;%s", title, xt, yt));
      if (logx) {
        M[i]->GetXaxis()->SetMoreLogLabels();
        M[i]->GetXaxis()->SetNoExponent();
      }
      if (perLayer)
        M[i]->LabelsOption("v");
      M[i]->Draw(first ? "E1" : "E1 SAME");
      first = false;
      leg->AddEntry(M[i], dir[i][1], "lp");
      if (O[i]) {
        O[i]->SetLineColor(cc);
        O[i]->SetMarkerColor(cc);
        O[i]->SetMarkerStyle(i == 0 ? 25 : 24);
        O[i]->SetMarkerSize(0.8);
        O[i]->SetLineWidth(2);
        O[i]->Draw("E1 SAME");
      }
    }
    leg->Draw();

    // Common info legend, drawn on every canvas: middle-right, above the
    // per-graph legend.
    {
      TLatex ti;
      ti.SetNDC();
      ti.SetTextSize(0.035);
      ti.DrawLatex(0.58, 0.56, "ACTS, Pythia8, Geant4");
      ti.DrawLatex(0.58, 0.51, "pp, pileup 30, |#eta| < 0.8");
    }

    if (hOpen && hOpen[0]) {
      TLatex tl;
      tl.SetNDC();
      tl.SetTextSize(0.035);
      // inside the plot frame (top-left), not in the top margin
      tl.DrawLatex(0.16, 0.85, Form("filled: %s    open: %s", mainTag, openTag));
    }
    c->Update();
    //    c->SaveAs(Form("%s.pdf", cname));

    // Collect the canvas and the individual histograms (renamed so they are
    // identifiable in the ROOT file: <canvas>_dir<i>[ _open ]).
    for (int i = 0; i < ndir; ++i) {
      if (M[i]) { M[i]->SetName(Form("%s_dir%d", cname, i)); gHists.push_back(M[i]); }
      if (O[i]) { O[i]->SetName(Form("%s_dir%d_open", cname, i)); gHists.push_back(O[i]); }
    }
    gCanvases.push_back(c);
  };

  // 1) efficiency vs pT: all reconstructed (filled) + physics reconstructed (open)
  draw("cmp_eff_vs_pT", "Efficiency vs p_{T}",
       "eff_all_vs_pT", "physeff_all_vs_pT",
       "p_{T} [GeV]", "efficiency", true, 1.15, true,
       "all reconstructed", "physics reconstructed", false);

  // 2) efficiency vs pT: all reconstructed (filled) + clean 0-fake tracks (open)
  draw("cmp_eff_clean_vs_pT", "Efficiency vs p_{T}",
       "eff_all_vs_pT", "eff_exact0_vs_pT",
       "p_{T} [GeV]", "efficiency", true, 1.15, true,
       "all reconstructed", "clean (0 fake)", false);

  // 2b) efficiency vs eta: all reconstructed (filled) + physics reconstructed (open)
  draw("cmp_eff_vs_eta", "Efficiency vs #eta",
       "eff_all_vs_eta", "physeff_all_vs_eta",
       "#eta", "efficiency", false, 1.15, true,
       "all reconstructed", "physics reconstructed", false);

  // 2c) efficiency vs eta: all reconstructed (filled) + clean 0-fake tracks (open)
  draw("cmp_eff_clean_vs_eta", "Efficiency vs #eta",
       "eff_all_vs_eta", "eff_exact0_vs_eta",
       "#eta", "efficiency", false, 1.15, true,
       "all reconstructed", "clean (0 fake)", false);

  // 3) fake-track rate (>=1 fake) vs pT
  draw("cmp_fakeRate_vs_pT", "Fake-track rate (#geq1 fake) vs p_{T}",
       "fakeRate_atLeast1_vs_pT", "",
       "p_{T} [GeV]", "fraction of tracks", true, -1, false, "", "", false);

  // 4) fake-track rate (>=1 fake) vs eta
  draw("cmp_fakeRate_vs_eta", "Fake-track rate (#geq1 fake) vs #eta",
       "fakeRate_atLeast1_vs_eta", "",
       "#eta", "fraction of tracks", false, -1, false, "", "", false);

  // 5) per-layer wrong-hit rate
  draw("cmp_fakeRate_perLayer", "Per-layer wrong-hit rate",
       "fakeRate_perLayer", "",
       "", "wrong / crossing tracks", false, -1, false, "", "", true);

  // 6) pT resolution vs pT: all tracks (filled) + clean (open)
  draw("cmp_resPt_vs_pT", "#sigma(p_{T})/p_{T} vs p_{T}",
       //     "res_ptrel_pt_c0_sigma", "res_ptrel_pt_c1_sigma",
       "res_ptrel_pt_c0_sigma", "",
       //     "res_ptrel_pt_c1_sigma", "",
       "p_{T} [GeV]", "#sigma(p_{T})/p_{T}", true, -1, false,
       "all tracks", "clean (0 fake)", false);

  // 7) pT resolution vs eta: all tracks (filled) + clean (open)
  draw("cmp_resPt_vs_eta", "#sigma(p_{T})/p_{T} vs #eta",
       "res_ptrel_eta_c0_sigma", "res_ptrel_eta_c1_sigma",
       "#eta", "#sigma(p_{T})/p_{T}", false, -1, false,
       "all tracks", "clean (0 fake)", false);

  // 8) DCAxy resolution vs pT: all tracks (filled) + clean (open)
  draw("cmp_resDcaXy_vs_pT", "#sigma(DCA_{xy}) vs p_{T}",
       "res_dcaxy_pt_c0_sigma", "res_dcaxy_pt_c1_sigma",
       "p_{T} [GeV]", "#sigma(DCA_{xy}) [#mum]", true, -1, false,
       "all tracks", "clean (0 fake)", false);

  // 9) DCAz resolution vs pT: all tracks (filled) + clean (open)
  draw("cmp_resDcaZ_vs_pT", "#sigma(DCA_{z}) vs p_{T}",
       "res_dcaz_pt_c0_sigma", "res_dcaz_pt_c1_sigma",
       "p_{T} [GeV]", "#sigma(DCA_{z}) [#mum]", true, -1, false,
       "all tracks", "clean (0 fake)", false);

  // ---- save all canvases to a multi-page PDF ----
  const TString pdfName = "compareEffResFakes.pdf";
  const size_t nc = gCanvases.size();
  for (size_t i = 0; i < nc; ++i) {
    TString suf = (nc == 1) ? "" : (i == 0 ? "(" : (i + 1 == nc ? ")" : ""));
    gCanvases[i]->Print(pdfName + suf);
  }
  if (nc) std::cout << "Wrote " << pdfName << std::endl;

  // ---- save canvases and individual histograms to a ROOT file ----
  TFile* fout = TFile::Open("compareEffResFakes.root", "RECREATE");
  if (fout && !fout->IsZombie()) {
    for (auto* h : gHists) if (h) h->Write();
    for (auto* c : gCanvases) if (c) c->Write();
    fout->Close();
    std::cout << "Wrote compareEffResFakes.root" << std::endl;
  } else {
    std::cerr << "ERROR: cannot open compareEffResFakes.root for writing" << std::endl;
  }
}
