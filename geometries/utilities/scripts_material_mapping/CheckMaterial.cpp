#include <TString.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TTree.h>
#include <TH2F.h>
#include <TProfile.h>
#include <vector>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TProfile2D.h>
#include <string>
#include <iostream>
#include <sstream>
#include <TMath.h>
#include <TGraph.h>
#include <TEfficiency.h>
#include <map>
#include <fstream>
#include <TLatex.h>
#include <TPave.h>
#include "TApplication.h"
#include "TROOT.h"

#include <algorithm>
#include <filesystem>
#include <regex>

namespace fs = std::filesystem;

bool MakeComparisonPlots(std::string geant_filename, std::string trk_filename, std::string outputdir,
                         std::string layout, std::string config, int entries);
void Draw_ratio(TCanvas* c, TProfile* h1, TProfile* h2, TLegend* leg, std::string name,
                std::string ALICE_label, std::string ALICE_text, std::string layout,
                std::string config, std::string outputdir);

// ./CheckMaterial.exe --G4Input MaterialTracks_mapping-ATLAS-P2-RUN4-01-00-00.root --TrkInput MaterialTracks-ATLAS-P2-RUN4-01-00-00.root --config ACTS
// ./CheckMaterial.exe --G4Input AtlasGeant4Geometry.root --TrkInput AtlasTrackingGeometry.root

static void show_usage(std::string name)
{
  std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "Usage: " << name << " <option(s)> SOURCES" << std::endl;
  printf("Produce comparison plots for validation of material maps:\n");
  printf("--help\t\t\tShow this help message\n");
  printf("Necessary arguments\n");
  printf("--G4Input\t\t\tGeant4 material (produced in material mapping process)\n");
  printf("--TrkInput\t\t\tMaterial collected on tracking layers (produced in validation job)\n");
  printf("--outdir\t\t\tfollowed by the name of the desirded output directory: [default : material_mapping_validation]\n");
  printf("--geoTag\t\t\tfollowed by the geoTag: default will be [default : ATLAS-P2-RUN4-01-00-00]\n");
  printf("--entries\t\t\tnumber of entries to process [default : -1]\n");
  printf("--config\t\t\tdefines if ATLAS or ACTS maps [default : ATLAS]\n");
  std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
  printf("For additional information on how to produce the input files refer to:\n");
  printf("https://twiki.cern.ch/twiki/bin/viewauth/Atlas/ProducingITkMaterialMapsMaster\n");
  std::cout << "---------------------------------------------------------------------------------------------------" << std::endl;
}

void myText(double x, double y, int color,
            std::string text, double text_size = 0.03)
{
  TLatex l;
  l.SetTextAlign(10);
  l.SetTextSize(text_size);
  l.SetNDC();
  l.SetTextFont(42); // Helvetica regular (TLatex default 62 is bold)
  l.SetTextColor(color);
  l.DrawLatex(x, y, text.c_str());
}

void ALICE_LABEL(double x, double y, std::string label = "ALICE3",
                 int color = 1, double text_size = 0.05)
{
  TLatex l; // l.SetTextAlign(12);
  l.SetTextSize(text_size);
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  l.DrawLatex(x, y, label.c_str());
}

void makeEtaLines(TCanvas* canvas, double zMin, double zMax, double rMax)
{

  canvas->cd();

  // Draw pseudorapidity lines
  std::vector<double> etas = {};
  for (unsigned int bin = 0; bin < 41; bin++)
    etas.push_back(-4.0 + bin * 0.2);

  for (auto& eta : etas) {
    double z_max = eta >= 0. ? zMax : zMin;
    double r_max = rMax;

    double theta = 2 * atan(exp(-eta));
    r_max = eta != 0. ? z_max * tan(theta) : rMax;

    if (r_max >= rMax) {
      r_max = rMax;
      z_max = eta != 0. ? r_max / tan(theta) : 0.;
    }

    // Draw the lines a bit shorter than the full frame so the eta labels at
    // their ends stay inside the canvas frame.
    const double frac = 0.85;
    const double z_line = z_max * frac;
    const double r_line = r_max * frac;

    TLine* track = new TLine(0., 0., z_line, r_line);

    if (abs(static_cast<int>(std::round(eta * 10))) % 10 == 0) {
      track->SetLineColor(kBlack);

      // Label the eta value just above the (shortened) end of each black line.
      TLatex* etaLabel = new TLatex();
      etaLabel->SetTextFont(42);
      etaLabel->SetTextSize(0.028);
      etaLabel->SetTextColor(kBlack);
      etaLabel->SetTextAlign(21);  // centered horizontally, anchored above point
      etaLabel->DrawLatex(z_line, r_line + rMax * 0.02, Form("#eta=%.1f", eta));
    } else {
      track->SetLineColor(17);
      track->SetLineStyle(7);
    }

    track->Draw("same");
  }
  return;
}

int main(int argc, char* argv[])
{

  gROOT->SetBatch(true);
  std::cout << "Checking your command..." << std::endl;
  for (int ar = 0; ar < argc; ar++)
    std::cout << std::string(argv[ar]) << " ";
  std::cout << std::endl;
  if (argc < 3) {
    show_usage(argv[0]);
    return 0;
  }
  std::cout << std::endl;

  // std::string layout = "ALICE3-23-12-26-peacock";
  std::string layout = "ALICE3-2026_06_10";
  std::string outputdir = "./material_mapping_validation";
  std::string config = "ACTS";
  std::string geant_filename = "";
  std::string trk_filename = "";
  int entries = -1;

  bool inputG4Set = false;
  bool inputTrkSet = false;
  bool outdirSet = false;

  for (int ar = 1; ar < argc; ++ar) {
    std::string arg = argv[ar];
    if (arg == "--help") {
      show_usage(argv[0]);
      return 0;
    } else if (arg == "--G4Input") { // getting input files
      if ((ar + 1 < argc) && (std::string(argv[ar + 1]).find("--") == std::string::npos)) {
        // Make sure we aren't at the end of argv!
        geant_filename = argv[ar + 1];
        if (geant_filename == "") {
          std::cerr << "Define an inputfile..." << std::endl;
          return 1;
        }
        if (!fs::exists(geant_filename)) {
          std::cerr << "File doesn't exist:  " << geant_filename << std::endl;
          return 1;
        }
        fs::path p(geant_filename);
        ar += 1; // Increment arg once to get next command
        inputG4Set = true;
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--G4Input argument necessary to specify the input file..." << std::endl;
        return 1;
      }
    } else if (arg == "--TrkInput") { // getting input files
      if ((ar + 1 < argc) && (std::string(argv[ar + 1]).find("--") == std::string::npos)) {
        // Make sure we aren't at the end of argv!
        trk_filename = argv[ar + 1];
        if (trk_filename == "") {
          std::cerr << "Define an inputfile..." << std::endl;
          return 1;
        }
        if (!fs::exists(trk_filename)) {
          std::cerr << "File doesn't exist:  " << trk_filename << std::endl;
          return 1;
        }
        fs::path p(trk_filename);
        ar += 1; // Increment arg once to get next command
        inputTrkSet = true;
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--TrkInput argument necessary to specify the input file..." << std::endl;
        return 1;
      }
    } else if (arg == "--outdir") {                                                         // getting output directory
      if ((ar + 1 < argc) && (std::string(argv[ar + 1]).find("--") == std::string::npos)) { // Make sure we aren't at the end of argv!
        outputdir = (argv[ar + 1] != "" && argv[ar + 1] != " ") ? argv[ar + 1] : outputdir;
        if (argv[ar + 1] == "" || argv[ar + 1] != " ") {
          std::cout << "output directory command used but without an arguement, setting to default value: " << outputdir << std::endl;
        }
        if (!fs::exists(outputdir)) {
          if (fs::create_directory(outputdir))
            std::cout << "Made directory for output: " << outputdir << std::endl;
          else {
            std::cerr << "Failed to make directory for output: " << outputdir << std::endl;
            return 1;
          }
        }
        ar += 1; // Increment ar once to get next command
        outdirSet = true;
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--outdir option requires one argument." << std::endl;
        return 1;
      }
    } else if (arg == "--geoTag") {                                                         // getting input directoris
      if ((ar + 1 < argc) && (std::string(argv[ar + 1]).find("--") == std::string::npos)) { // Make sure we aren't at the end of argv!
        layout = argv[ar + 1] != "" ? argv[ar + 1] : layout;
        if (argv[ar + 1] == "") {
          std::cout << "geoTag command used but without an argument, setting to default value:  " << layout << std::endl;
        }
        ar += 1; // Increment arg once to get next command
      } else {   // Uh-oh, there was no argument to the destination option.
        std::cerr << "--geoTag option requires an argument." << std::endl;
        return 1;
      }
    } else if (arg == "--config") {                                                         // getting input directoris
      if ((ar + 1 < argc) && (std::string(argv[ar + 1]).find("--") == std::string::npos)) { // Make sure we aren't at the end of argv!
        config = argv[ar + 1] != "" ? argv[ar + 1] : config;
        if (argv[ar + 1] == "") {
          std::cout << "config command used but without an argument, setting to default value:  " << config << std::endl;
        }
        ar += 1; // Increment arg once to get next command
        if (config != "ATLAS" and config != "ACTS") {
          std::cerr << "invalid argument for --config. Use ATLAS or ACTS." << std::endl;
          return 1;
        }
      } else { // Uh-oh, there was no argument to the destination option.
        std::cerr << "--config option requires an argument." << std::endl;
        return 1;
      }
    } else if (arg == "--entries") {                                                        // getting input directoris
      if ((ar + 1 < argc) && (std::string(argv[ar + 1]).find("--") == std::string::npos)) { // Make sure we aren't at the end of argv!
        entries = argv[ar + 1] != "" ? std::stoi(argv[ar + 1]) : entries;
        if (argv[ar + 1] == "") {
          std::cout << "entries command used but without an argument, setting to default value:  " << entries << std::endl;
        }
        ar += 1; // Increment arg once to get next command
      } else {   // Uh-oh, there was no argument to the destination option.
        std::cerr << "--entries option requires an argument." << std::endl;
        return 1;
      }
    }
  }

  if (not inputG4Set) {
    std::cout << "Define --G4Input path/to/dir1/file.root" << std::endl;
    return 1;
  }

  if (not inputTrkSet) {
    std::cout << "Define --TrkInput path/to/dir1/file.root" << std::endl;
    return 1;
  }

  if (not outdirSet) {
    if (!fs::exists(outputdir)) {
      if (fs::create_directory(outputdir))
        std::cout << "Made directory for output: " << outputdir << std::endl;
      else {
        std::cerr << "Failed to make directory for output: " << outputdir << std::endl;
        return 1;
      }
    }
  }

  TApplication theApp("Analysis", &argc, argv);
  if (!MakeComparisonPlots(geant_filename, trk_filename, outputdir, layout, config, entries))
    return 1;

  return 0;
}

bool MakeComparisonPlots(std::string geant_filename, std::string trk_filename, std::string outputdir,
                         std::string layout, std::string config, int entries)
{

  std::string ALICE_label = "ALICE3";
  std::string ALICE_text = "Simulation Internal";

  double zMin = -4200.;
  double zMax = 4200.;
  int zBin = 4200;

  double rMin = 0.;
  double rMax = 1500.;
  int rBin = 1500;

  double etaMin = -4.2;
  double etaMax = 4.2;
  int etaBin = 105;

  double phiMin = -4.;
  double phiMax = 4.;
  int phiBin = 150;

  float eta, phi, t_x0, t_L0;
  std::vector<float>* mat_x = 0;
  std::vector<float>* mat_y = 0;
  std::vector<float>* mat_z = 0;

  int steps = 0;
  float mat_x_steps[8000];
  float mat_y_steps[8000];
  float mat_z_steps[8000];

  std::cout << "Reading file " << geant_filename << " and comparing to " << trk_filename << std::endl;

  TFile* geant_file = TFile::Open(geant_filename.c_str());
  TFile* trk_file = TFile::Open(trk_filename.c_str());

  TProfile* geant_t_x0_vs_v_eta = new TProfile("geant_t_x0_vs_v_eta", ";#eta; <t/X_{0}>", etaBin, etaMin, etaMax);
  geant_t_x0_vs_v_eta->SetMarkerSize(1.2);
  geant_t_x0_vs_v_eta->SetMarkerStyle(21);
  geant_t_x0_vs_v_eta->SetMarkerColor(kRed + 1);
  geant_t_x0_vs_v_eta->SetLineColor(kRed + 1);
  geant_t_x0_vs_v_eta->SetStats(0);
  geant_t_x0_vs_v_eta->SetTitle(0);

  TProfile* geant_t_x0_vs_v_phi = new TProfile("geant_t_x0_vs_v_phi", ";#phi; <t/X_{0}>", phiBin, phiMin, phiMax);
  geant_t_x0_vs_v_phi->SetMarkerSize(1.2);
  geant_t_x0_vs_v_phi->SetMarkerStyle(21);
  geant_t_x0_vs_v_phi->SetMarkerColor(kRed + 1);
  geant_t_x0_vs_v_phi->SetLineColor(kRed + 1);
  geant_t_x0_vs_v_phi->SetStats(0);
  geant_t_x0_vs_v_phi->SetTitle(0);

  TProfile* geant_t_l0_vs_v_eta = new TProfile("geant_t_l0_vs_v_eta", ";#eta; <t/L_{0}>", etaBin, etaMin, etaMax);
  geant_t_l0_vs_v_eta->SetMarkerSize(1.2);
  geant_t_l0_vs_v_eta->SetMarkerStyle(21);
  geant_t_l0_vs_v_eta->SetMarkerColor(kRed + 1);
  geant_t_l0_vs_v_eta->SetLineColor(kRed + 1);
  geant_t_l0_vs_v_eta->SetStats(0);
  geant_t_l0_vs_v_eta->SetTitle(0);

  TProfile* geant_t_l0_vs_v_phi = new TProfile("geant_t_l0_vs_v_phi", ";#phi; <t/L_{0}>", phiBin, phiMin, phiMax);
  geant_t_l0_vs_v_phi->SetMarkerSize(1.2);
  geant_t_l0_vs_v_phi->SetMarkerStyle(21);
  geant_t_l0_vs_v_phi->SetMarkerColor(kRed + 1);
  geant_t_l0_vs_v_phi->SetLineColor(kRed + 1);
  geant_t_l0_vs_v_phi->SetStats(0);
  geant_t_l0_vs_v_phi->SetTitle(0);

  TH2D* geant_mat_z_vs_mat_r = new TH2D("geant_mat_z_vs_mat_r", ";z [mm]; r [mm]", zBin, zMin, zMax, rBin, rMin, rMax);
  geant_mat_z_vs_mat_r->SetMarkerSize(0.004);
  geant_mat_z_vs_mat_r->SetMarkerColor(kRed + 1);
  geant_mat_z_vs_mat_r->SetLineColor(kRed + 1);
  geant_mat_z_vs_mat_r->SetStats(0);
  geant_mat_z_vs_mat_r->SetTitle(0);

  TProfile* trk_t_x0_vs_v_eta = new TProfile("trk_t_x0_vs_v_eta", ";#eta; <t/X_{0}>", etaBin, etaMin, etaMax);
  trk_t_x0_vs_v_eta->SetMarkerSize(0.8);
  trk_t_x0_vs_v_eta->SetMarkerStyle(20);
  trk_t_x0_vs_v_eta->SetMarkerColor(kBlack);
  trk_t_x0_vs_v_eta->SetLineColor(kBlack);
  trk_t_x0_vs_v_eta->SetStats(0);
  trk_t_x0_vs_v_eta->SetTitle(0);

  TProfile* trk_t_x0_vs_v_phi = new TProfile("trk_t_x0_vs_v_phi", ";#phi; <t/X_{0}>", phiBin, phiMin, phiMax);
  trk_t_x0_vs_v_phi->SetMarkerSize(0.8);
  trk_t_x0_vs_v_phi->SetMarkerStyle(20);
  trk_t_x0_vs_v_phi->SetMarkerColor(kBlack);
  trk_t_x0_vs_v_phi->SetLineColor(kBlack);
  trk_t_x0_vs_v_phi->SetStats(0);
  trk_t_x0_vs_v_phi->SetTitle(0);

  TProfile* trk_t_l0_vs_v_eta = new TProfile("trk_t_l0_vs_v_eta", ";#eta; <t/L_{0}>", etaBin, etaMin, etaMax);
  trk_t_l0_vs_v_eta->SetMarkerSize(0.8);
  trk_t_l0_vs_v_eta->SetMarkerStyle(20);
  trk_t_l0_vs_v_eta->SetMarkerColor(kBlack);
  trk_t_l0_vs_v_eta->SetLineColor(kBlack);
  trk_t_l0_vs_v_eta->SetStats(0);
  trk_t_l0_vs_v_eta->SetTitle(0);

  TProfile* trk_t_l0_vs_v_phi = new TProfile("trk_t_l0_vs_v_phi", ";#phi; <t/L_{0}>", phiBin, phiMin, phiMax);
  trk_t_l0_vs_v_phi->SetMarkerSize(0.8);
  trk_t_l0_vs_v_phi->SetMarkerStyle(20);
  trk_t_l0_vs_v_phi->SetMarkerColor(kBlack);
  trk_t_l0_vs_v_phi->SetLineColor(kBlack);
  trk_t_l0_vs_v_phi->SetStats(0);
  trk_t_l0_vs_v_phi->SetTitle(0);

  TH2D* trk_mat_z_vs_mat_r = new TH2D("trk_mat_z_vs_mat_r", ";z [mm]; r [mm]", zBin, zMin, zMax, rBin, rMin, rMax);
  trk_mat_z_vs_mat_r->SetMarkerSize(0.25);
  trk_mat_z_vs_mat_r->SetMarkerStyle(2);
  trk_mat_z_vs_mat_r->SetMarkerColor(kBlack);
  trk_mat_z_vs_mat_r->SetLineColor(kBlack);
  trk_mat_z_vs_mat_r->SetStats(0);
  trk_mat_z_vs_mat_r->SetTitle(0);

  TTree* geant_tree = nullptr;

  if (config == "ACTS") {
    geant_tree = (TTree*)geant_file->Get("material_tracks");
    geant_tree->SetBranchAddress("v_eta", &eta);
    geant_tree->SetBranchAddress("v_phi", &phi);
    geant_tree->SetBranchAddress("t_X0", &t_x0);
    geant_tree->SetBranchAddress("t_L0", &t_L0);
    geant_tree->SetBranchAddress("mat_x", &mat_x);
    geant_tree->SetBranchAddress("mat_y", &mat_y);
    geant_tree->SetBranchAddress("mat_z", &mat_z);
  } else {
    mat_x = new std::vector<float>();
    mat_y = new std::vector<float>();
    mat_z = new std::vector<float>();

    geant_tree = (TTree*)geant_file->Get("MaterialMapper");
    geant_tree->SetBranchAddress("Eta", &eta);
    geant_tree->SetBranchAddress("Phi", &phi);
    geant_tree->SetBranchAddress("PathInX0", &t_x0);
    geant_tree->SetBranchAddress("PathInL0", &t_L0);
    geant_tree->SetBranchAddress("MaterialSteps", &steps);
    geant_tree->SetBranchAddress("MaterialStepPositionX", &mat_x_steps);
    geant_tree->SetBranchAddress("MaterialStepPositionY", &mat_y_steps);
    geant_tree->SetBranchAddress("MaterialStepPositionZ", &mat_z_steps);
  }

  // read all entries and fill the histograms
  Long64_t nentries = geant_tree->GetEntries();
  if (entries > 0 and nentries > entries)
    nentries = entries;

  int maxMatEntries = 100000;

  std::cout << "---- From file " << geant_filename << " getting tree with " << nentries << " entries" << std::endl;

  for (Long64_t i = 0; i < nentries; i++) {

    geant_tree->GetEntry(i);

    geant_t_x0_vs_v_eta->Fill(eta, t_x0);
    geant_t_l0_vs_v_eta->Fill(eta, t_L0);

    if (abs(eta) < 4) {
      geant_t_x0_vs_v_phi->Fill(phi, t_x0);
      geant_t_l0_vs_v_phi->Fill(phi, t_L0);
    }

    if (i > maxMatEntries)
      continue;

    // std::cout << eta << ", " << phi << ", " << t_x0 << ", " << t_L0 << std::endl;

    if (config == "ATLAS") {
      mat_x->clear();
      mat_y->clear();
      mat_z->clear();

      for (unsigned int step = 0; step < steps; step++) {
        mat_x->push_back(mat_x_steps[step]);
        mat_y->push_back(mat_y_steps[step]);
        mat_z->push_back(mat_z_steps[step]);
      }
    }

    // std::cout << "---------------------------------" << std::endl;
    for (unsigned int hit = 0; hit < mat_x->size(); hit++) {
      float mat_r = sqrt(mat_x->at(hit) * mat_x->at(hit) + mat_y->at(hit) * mat_y->at(hit));
      geant_mat_z_vs_mat_r->Fill(mat_z->at(hit), mat_r);
      // std::cout << mat_x->at(hit) << ", " << mat_y->at(hit) << ", " << mat_z->at(hit) << std::endl;
    }
    // std::cout << "---------------------------------" << std::endl;
  }

  TTree* trk_tree = nullptr;

  if (config == "ACTS") {
    trk_tree = (TTree*)trk_file->Get("material_tracks");
    trk_tree->SetBranchAddress("v_eta", &eta);
    trk_tree->SetBranchAddress("v_phi", &phi);
    trk_tree->SetBranchAddress("t_X0", &t_x0);
    trk_tree->SetBranchAddress("t_L0", &t_L0);
    trk_tree->SetBranchAddress("mat_x", &mat_x);
    trk_tree->SetBranchAddress("mat_y", &mat_y);
    trk_tree->SetBranchAddress("mat_z", &mat_z);
  } else {
    trk_tree = (TTree*)trk_file->Get("MaterialMapper");
    trk_tree->SetBranchAddress("Eta", &eta);
    trk_tree->SetBranchAddress("Phi", &phi);
    trk_tree->SetBranchAddress("PathInX0", &t_x0);
    trk_tree->SetBranchAddress("PathInL0", &t_L0);
    trk_tree->SetBranchAddress("MaterialSteps", &steps);
    trk_tree->SetBranchAddress("MaterialStepPositionX", &mat_x_steps);
    trk_tree->SetBranchAddress("MaterialStepPositionY", &mat_y_steps);
    trk_tree->SetBranchAddress("MaterialStepPositionZ", &mat_z_steps);
  }

  // read all entries and fill the histograms
  nentries = trk_tree->GetEntries();
  if (entries > 0 and nentries > entries)
    nentries = entries;

  std::cout << "---- From file " << trk_filename << " getting tree with " << nentries << " entries" << std::endl;

  for (Long64_t i = 0; i < nentries; i++) {

    trk_tree->GetEntry(i);

    trk_t_x0_vs_v_eta->Fill(eta, t_x0);
    trk_t_l0_vs_v_eta->Fill(eta, t_L0);

    if (abs(eta) < 4.0) {
      trk_t_x0_vs_v_phi->Fill(phi, t_x0);
      trk_t_l0_vs_v_phi->Fill(phi, t_L0);
    }

    if (i > maxMatEntries)
      continue;

    if (config == "ATLAS") {
      mat_x->clear();
      mat_y->clear();
      mat_z->clear();

      for (unsigned int step = 0; step < steps; step++) {
        mat_x->push_back(mat_x_steps[step]);
        mat_y->push_back(mat_y_steps[step]);
        mat_z->push_back(mat_z_steps[step]);
      }
    }

    for (unsigned int hit = 0; hit < mat_x->size(); hit++) {
      float mat_r = sqrt(mat_x->at(hit) * mat_x->at(hit) + mat_y->at(hit) * mat_y->at(hit));
      trk_mat_z_vs_mat_r->Fill(mat_z->at(hit), mat_r);
    }
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
  canvas->SetTicks(1, 1); // ticks also on top and right

  geant_mat_z_vs_mat_r->Draw("SCAT");
  trk_mat_z_vs_mat_r->Draw("SCAT same");

  makeEtaLines(canvas, zMin, zMax, rMax - 200);

  ALICE_LABEL(0.14, 0.82, ALICE_label);
  myText(0.28, 0.818, 1, "Simulation", 0.05);
  // the scatter markers are too small to show in a legend: use proxy lines
  TLine* geant_line = new TLine();
  geant_line->SetLineColor(kRed + 1);
  geant_line->SetLineWidth(3);
  TLine* trk_line = new TLine();
  trk_line->SetLineColor(kBlack);
  trk_line->SetLineWidth(3);

  TLegend* map_legend = new TLegend(0.66, 0.76, 0.87, 0.86);
  map_legend->SetMargin(0.25); // fraction of the width used for the line column
  map_legend->SetBorderSize(0);
  map_legend->SetFillStyle(0);
  map_legend->SetTextFont(42);
  map_legend->SetTextSize(0.03);
  map_legend->AddEntry(geant_line, "Geant4 Material", "l");
  map_legend->AddEntry(trk_line, "ACTS Material", "l");
  map_legend->Draw();

  canvas->SaveAs((outputdir + "/material_comparison.pdf").c_str());

  TLegend* legend = new TLegend(0.66, 0.72, 0.88, 0.84);
  legend->SetBorderSize(0);
  legend->SetTextAlign(12);
  legend->SetTextSize(0.04);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->AddEntry(geant_t_x0_vs_v_eta, "Geant4 Material", "epl");
  legend->AddEntry(trk_t_x0_vs_v_eta, "ACTS Material", "epl");

  Draw_ratio(canvas, geant_t_x0_vs_v_eta, trk_t_x0_vs_v_eta, legend, "t_x0_vs_v_eta", ALICE_label, ALICE_text, layout, config, outputdir);
  Draw_ratio(canvas, geant_t_l0_vs_v_eta, trk_t_l0_vs_v_eta, legend, "t_l0_vs_v_eta", ALICE_label, ALICE_text, layout, config, outputdir);
  Draw_ratio(canvas, geant_t_x0_vs_v_phi, trk_t_x0_vs_v_phi, legend, "t_x0_vs_v_phi", ALICE_label, ALICE_text, layout, config, outputdir);
  Draw_ratio(canvas, geant_t_l0_vs_v_phi, trk_t_l0_vs_v_phi, legend, "t_l0_vs_v_phi", ALICE_label, ALICE_text, layout, config, outputdir);

  return 0;
}

void Draw_ratio(TCanvas* c, TProfile* h1, TProfile* h2, TLegend* leg, std::string name,
                std::string ALICE_label, std::string ALICE_text,
                std::string layout, std::string config, std::string outputdir)
{

  c->Clear();

  // Upper plot will be in pad1
  TPad* pad1 = new TPad("pad1", "pad1", 0., 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.07); // Upper and lower plot are joined
  pad1->SetLeftMargin(0.1);    // Upper and lower plot are joined
  pad1->SetTicks(1, 1);        // ticks also on top and right
  pad1->Draw();                // Draw the upper pad: pad1

  TPad* pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.3);
  pad2->SetLeftMargin(0.1); // Upper and lower plot are joined
  pad2->SetTicks(1, 1);
  pad2->Draw();

  pad1->cd();
  pad1->Clear();

  TH1D* hnum = h1->ProjectionX();
  hnum->SetLineColor(h1->GetLineColor());
  hnum->SetMarkerStyle(h1->GetMarkerStyle());
  hnum->SetMarkerColor(h1->GetMarkerColor());
  hnum->SetMarkerSize(h1->GetMarkerSize());
  hnum->SetStats(0);
  hnum->SetTitle(0);
  TH1D* hden = h2->ProjectionX();
  hden->SetLineColor(h2->GetLineColor());
  hden->SetMarkerStyle(h2->GetMarkerStyle());
  hden->SetMarkerColor(h2->GetMarkerColor());
  hden->SetMarkerSize(h2->GetMarkerSize());
  hden->SetStats(0);
  hden->SetTitle(0);

  double max_hist[4];
  max_hist[0] = hnum->GetMaximum() * 1.45; // headroom for label and legend
  max_hist[1] = hden->GetMaximum() * 1.45;

  hnum->GetYaxis()->SetTitleOffset(0.8);
  hnum->GetYaxis()->SetRangeUser(0, *std::max_element(std::begin(max_hist), std::end(max_hist)));

  hnum->Draw("E");
  hden->Draw("E SAME");
  leg->Draw("SAME");

  ALICE_LABEL(0.14, 0.78, ALICE_label);
  myText(0.24, 0.778, 1, "Simulation", 0.05);

  // go to bottom
  pad2->cd();
  pad2->Clear();

  // Define the ratio plot

  TLine line = TLine(hnum->GetXaxis()->GetXmin(), 1, hnum->GetXaxis()->GetXmax(), 1);
  line.SetLineColor(kRed + 1);
  line.SetLineWidth(1);

  TProfile* ratio = (TProfile*)hnum->Clone((name + "_ratio").c_str());
  ratio->SetLineColor(kBlack);

  ratio->SetMarkerStyle(20);
  ratio->SetMarkerSize(0.8);
  ratio->SetStats(0); // No statistics on lower plot
  ratio->Divide(hden);

  // ratio->SetMarkerStyle(7);
  ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
  ratio->GetYaxis()->SetTitle("#frac{Geant4}{ACTS}");
  ratio->GetYaxis()->CenterTitle(1);
  ratio->GetYaxis()->SetNdivisions(502);
  ratio->GetYaxis()->SetTitleSize(0.12);
  ratio->GetYaxis()->SetTitleFont(42);
  ratio->GetYaxis()->SetTitleOffset(0.25);
  ratio->GetYaxis()->SetLabelSize(0.12);
  ratio->GetXaxis()->SetNdivisions(510);
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetTitleFont(42);
  ratio->GetXaxis()->SetTitleOffset(1.);
  ratio->GetXaxis()->SetLabelSize(0.12);
  ratio->GetXaxis()->SetTickSize(0.07);
  ratio->GetXaxis()->SetTickSize(0.07);

  ratio->Draw("ep");
  line.Draw("same");

  c->SaveAs((outputdir + "/ratio_" + name + ".pdf").c_str());

  delete ratio;
}
