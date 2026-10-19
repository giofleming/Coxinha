// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "LutMaker.hpp"
#include "LutRootData.hpp"

#include "ActsToO2Track.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TTree.h>
#include <TSystem.h>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <stdexcept>
#include <vector>

using namespace ActsExamples;

namespace AliceActsTrk
{

namespace
{

std::vector<double> makeBinEdges(const o2::delphes::map_t& map)
{
  std::vector<double> edges(static_cast<std::size_t>(map.nbins) + 1u);
  const double width = static_cast<double>(map.max - map.min) /
                       static_cast<double>(map.nbins);
  for (int i = 0; i <= map.nbins; ++i) {
    const double coordinate = static_cast<double>(map.min) +
                              static_cast<double>(i) * width;
    edges[static_cast<std::size_t>(i)] =
      map.log ? std::pow(10.0, coordinate) : coordinate;
  }
  return edges;
}

std::vector<double> makeLinearBinEdges(int nbins, double min, double max)
{
  std::vector<double> edges(static_cast<std::size_t>(nbins) + 1u);
  const double width = (max - min) / static_cast<double>(nbins);
  for (int i = 0; i <= nbins; ++i) {
    edges[static_cast<std::size_t>(i)] = min + static_cast<double>(i) * width;
  }
  return edges;
}

/// Compute the truth local coordinates on exactly the same reference surface
/// used by the fitted ACTS track. This intentionally mirrors the definition
/// used by the ACTS ROOT track-summary writer for t_d0/t_z0.
std::optional<Acts::Vector2> truthLocalAtReference(const Acts::BoundTrackParameters& track,
                                                   const ActsExamples::SimParticle& particle,
                                                   const Acts::GeometryContext& geometryContext)
{
  const Acts::Surface& referenceSurface = track.referenceSurface();

  const auto intersection =
    referenceSurface
      .intersect(geometryContext,
                 particle.position(),
                 particle.direction(),
                 Acts::BoundaryTolerance::Infinite())
      .closest();

  const auto localResult =
    referenceSurface.globalToLocal(geometryContext,
                                   intersection.position(),
                                   particle.direction());

  if (!localResult.ok()) {
    return std::nullopt;
  }

  return localResult.value();
}

} // namespace

std::map<int, LutContainer> LutMaker::mLutMap;
std::map<int, std::vector<std::uint64_t>> LutMaker::mPtTruthCounts;
std::map<int, std::vector<std::uint64_t>> LutMaker::mPtMatchedCounts;
std::map<int, std::vector<std::array<double, 15>>> LutMaker::mCovSums;
std::map<int, std::vector<std::array<double, 15>>> LutMaker::mCov2Sums;
std::map<int, std::vector<ResidualAccumulator>> LutMaker::mResiduals;
std::map<int, std::vector<std::uint64_t>> LutMaker::mEntryCounts;
std::map<int, TH3F*> LutMaker::res_acts_d0_vs_pT;
std::map<int, TH3F*> LutMaker::res_o2_y_vs_pT;
std::map<int, TH3F*> LutMaker::res_qopt_vs_pT;
std::map<int, TH3F*> LutMaker::res_pt_vs_pT;
TH1D* LutMaker::h_stat = nullptr;
TTree* LutMaker::mTree = nullptr;
LutMaker::TreeEntry LutMaker::mTreeContainer;

LutMaker::LutMaker(const Config& config,
                   std::unique_ptr<const Acts::Logger> _logger) : IAlgorithm("LutMaker", std::move(_logger)),
                                                                  m_cfg(config)
{
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles");
  }
  if (m_cfg.inputTrackParticleMatching.empty()) {
    throw std::invalid_argument("LutMaker requires an explicit inputTrackParticleMatching collection");
  }

  if (m_cfg.produceOutputTree) {
    mTree = new TTree("ParticleTrackTree", "ParticleTrackTree");
    mTree->Branch("pdg", &mTreeContainer.pdg, "pdg/I");
    mTree->Branch("px", &mTreeContainer.px, "px/F");
    mTree->Branch("py", &mTreeContainer.py, "py/F");
    mTree->Branch("pz", &mTreeContainer.pz, "pz/F");
    mTree->Branch("pxTrue", &mTreeContainer.pxTrue, "pxTrue/F");
    mTree->Branch("pyTrue", &mTreeContainer.pyTrue, "pyTrue/F");
    mTree->Branch("pzTrue", &mTreeContainer.pzTrue, "pzTrue/F");
    mTree->SetDirectory(0);
  }

  std::vector<int> pdgCodes = {11, 13, 211, 321, 2212}; // example PDG codes for electron, muon, pion, kaon, proton
  std::map<int, float> pdgMasses = {
    {11, 0.000511f}, // electron mass in GeV/c^2
    {13, 0.10566f},  // muon mass in GeV/c^2
    {211, 0.13957f}, // pion mass in GeV/c^2
    {321, 0.49367f}, // kaon mass in GeV/c^2
    {2212, 0.93827f} // proton mass in GeV/c^2
  };
  for (int pdg : pdgCodes) {
    ACTS_DEBUG("Initializing LUT for PDG " << pdg);
    o2::delphes::lutHeader_t header;
    header.version = LUTCOVM_VERSION;
    header.pdg = pdg;
    header.mass = pdgMasses.at(pdg);
    // Magnetic field written into the LUT header, in Tesla. Taken from the
    // configuration - it used to be hardcoded to 0.5, which silently disagreed
    // with the field the tracks were actually reconstructed in (ChainConfig
    // defaults to 2.0 T) and with the value put into the LUT file name below.
    header.field = m_cfg.magneticField;
    header.nchmap = {1, 0.f, 1000.f, false};
    header.radmap = {1, 0.f, 100.0f, false};
    header.etamap = {10, -4.0f, 4.0f, false};
    // Log10 pT binning from 1e-2 to 1e2 GeV/c.
    header.ptmap = {200, -2.0f, 2.0f, true};
    mLutMap[pdg].header = header;

    ACTS_DEBUG("Initializing lutmap");
    o2::delphes::lutEntry_t entry;
    // Now filling with empty LUT entries
    for (int nch_bin = 0; nch_bin < header.nchmap.nbins; ++nch_bin) {
      for (int rad_bin = 0; rad_bin < header.radmap.nbins; ++rad_bin) {
        for (int eta_bin = 0; eta_bin < header.etamap.nbins; ++eta_bin) {
          for (int pt_bin = 0; pt_bin < header.ptmap.nbins; ++pt_bin) {
            entry.nch = header.nchmap.eval(nch_bin);
            entry.eta = header.etamap.eval(eta_bin);
            entry.pt = header.ptmap.eval(pt_bin);
            mLutMap[pdg].entries.push_back(entry);
          }
        }
      }
    }

    mPtTruthCounts[pdg] = std::vector<std::uint64_t>(header.ptmap.nbins, 0u);
    mPtMatchedCounts[pdg] = std::vector<std::uint64_t>(header.ptmap.nbins, 0u);

    const std::size_t nEntries = static_cast<std::size_t>(header.nchmap.nbins) * header.radmap.nbins * header.etamap.nbins * header.ptmap.nbins;
    mCovSums[pdg] = std::vector<std::array<double, 15>>(nEntries, {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.});
    mCov2Sums[pdg] = std::vector<std::array<double, 15>>(nEntries, {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.});
    mResiduals[pdg] = std::vector<ResidualAccumulator>(nEntries);
    mEntryCounts[pdg] = std::vector<std::uint64_t>(nEntries, 0u);
    // QA histograms use the exact physical pT edges of the LUT. For a
    // logarithmic map the header stores log10(pT) limits, so using min/max
    // directly as ROOT axis limits would be wrong.
    const auto ptEdges = makeBinEdges(header.ptmap);
    const auto etaEdges = makeBinEdges(header.etamap);
    const auto diffEdges = makeLinearBinEdges(2000, -1, 1);

    // Both position residual histograms below are stored in cm, matching the
    // O2 track parametrization and therefore the units of covm[0].
    res_acts_d0_vs_pT[pdg] = new TH3F(
      Form("res_acts_d0_vs_pT_pdg%d", pdg),
      Form("ACTS d_{0} residual vs p_{T} for PDG %d;p_{T}^{true} [GeV/c];d_{0}^{fit}-d_{0}^{true} [cm];#eta^{true}", pdg),
      header.ptmap.nbins, ptEdges.data(),
      diffEdges.size() - 1, diffEdges.data(),
      etaEdges.size() - 1, etaEdges.data());
    res_acts_d0_vs_pT[pdg]->SetDirectory(0);

    res_o2_y_vs_pT[pdg] = new TH3F(
      Form("res_o2_y_vs_pT_pdg%d", pdg),
      Form("O2 local-Y residual vs p_{T} for PDG %d;p_{T}^{true} [GeV/c];Y^{reco}-Y^{truth(prod)} [cm];#eta^{true}", pdg),
      header.ptmap.nbins, ptEdges.data(),
      diffEdges.size() - 1, diffEdges.data(),
      etaEdges.size() - 1, etaEdges.data());
    res_o2_y_vs_pT[pdg]->SetDirectory(0);

    res_qopt_vs_pT[pdg] = new TH3F(
      Form("res_qopT_res_vs_pT_pdg%d", pdg),
      Form("q/pT residuals vs pT for PDG %d;p_{T}^{true} [GeV/c];q/p_{T} residual;#eta^{true}", pdg),
      header.ptmap.nbins, ptEdges.data(),
      diffEdges.size() - 1, diffEdges.data(),
      etaEdges.size() - 1, etaEdges.data());
    res_qopt_vs_pT[pdg]->SetDirectory(0);

    res_pt_vs_pT[pdg] = new TH3F(
      Form("res_pT_res_vs_pT_pdg%d", pdg),
      Form("pT residuals vs pT for PDG %d;p_{T}^{true} [GeV/c];p_{T}^{reco}-p_{T}^{true} [GeV/c];#eta^{true}", pdg),
      header.ptmap.nbins, ptEdges.data(),
      diffEdges.size() - 1, diffEdges.data(),
      etaEdges.size() - 1, etaEdges.data());
    res_pt_vs_pT[pdg]->SetDirectory(0);
  }
  h_stat = new TH1D("h_stat", "LUT Maker Statistics;statistic;count", 10, -0.5, 9.5);
  h_stat->SetDirectory(0);

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);
}

LutMaker::~LutMaker()
{
  ACTS_DEBUG("Called LutMaker destructor, with lutMap of size: " << mLutMap.size());

  const std::string outputDir = m_cfg.outputDirectory + "/";
  if (gSystem->AccessPathName(outputDir.c_str())) {
    gSystem->mkdir(outputDir.c_str(), true);
  }

  for (auto& [pdg, lut] : mLutMap) {
    // First fill the efficiency values into the LUT
    const auto& truthCounts = mPtTruthCounts[pdg];
    const auto& matchedCounts = mPtMatchedCounts[pdg];

    const int nnch = lut.header.nchmap.nbins;
    const int nrad = lut.header.radmap.nbins;
    const int neta = lut.header.etamap.nbins;
    const int npt = lut.header.ptmap.nbins;

    // Fill efficiency and average covariance + eigen decomposition per entry
    const auto& covSums = mCovSums[pdg];
    const auto& cov2Sums = mCov2Sums[pdg];
    const auto& entryCounts = mEntryCounts[pdg];

    for (int nch_bin = 0; nch_bin < nnch; ++nch_bin) {
      for (int rad_bin = 0; rad_bin < nrad; ++rad_bin) {
        for (int eta_bin = 0; eta_bin < neta; ++eta_bin) {
          TH2F* h_pull0 = new TH2F(Form("h_pull0_pdg%d_nch%d_rad%d_eta%d", pdg, nch_bin, rad_bin, eta_bin),
                                   Form("Pull distribution for PDG %d, nch %d, rad %d, eta %d;pull[0];entries", pdg, nch_bin, rad_bin, eta_bin),
                                   npt, -0.5, -0.5 + npt,
                                   100, -5., 5.);
          TH2F* h_pull1 = new TH2F(Form("h_pull1_pdg%d_nch%d_rad%d_eta%d", pdg, nch_bin, rad_bin, eta_bin),
                                   Form("Pull distribution for PDG %d, nch %d, rad %d, eta %d;pull[1];entries", pdg, nch_bin, rad_bin, eta_bin),
                                   npt, -0.5, -0.5 + npt,
                                   100, -5., 5.);
          TH2F* h_pull2 = new TH2F(Form("h_pull2_pdg%d_nch%d_rad%d_eta%d", pdg, nch_bin, rad_bin, eta_bin),
                                   Form("Pull distribution for PDG %d, nch %d, rad %d, eta %d;pull[2];entries", pdg, nch_bin, rad_bin, eta_bin),
                                   npt, -0.5, -0.5 + npt,
                                   100, -5., 5.);
          TH2F* h_pull3 = new TH2F(Form("h_pull3_pdg%d_nch%d_rad%d_eta%d", pdg, nch_bin, rad_bin, eta_bin),
                                   Form("Pull distribution for PDG %d, nch %d, rad %d, eta %d;pull[3];entries", pdg, nch_bin, rad_bin, eta_bin),
                                   npt, -0.5, -0.5 + npt,
                                   100, -5., 5.);
          TH2F* h_pull4 = new TH2F(Form("h_pull4_pdg%d_nch%d_rad%d_eta%d", pdg, nch_bin, rad_bin, eta_bin),
                                   Form("Pull distribution for PDG %d, nch %d, rad %d, eta %d;pull[4];entries", pdg, nch_bin, rad_bin, eta_bin),
                                   npt, -0.5, -0.5 + npt,
                                   100, -5., 5.);

          for (int pt_bin = 0; pt_bin < npt; ++pt_bin) {
            const std::size_t entryIdx = static_cast<std::size_t>(nch_bin) * (nrad * neta * npt) + rad_bin * (neta * npt) + eta_bin * npt + pt_bin;
            auto& entry = lut.entries[entryIdx];

            // Efficiency
            const auto denom = truthCounts[pt_bin];
            const auto numer = matchedCounts[pt_bin];
            entry.eff = denom > 0u ? static_cast<float>(numer) / static_cast<float>(denom) : 0.f;

            // Average covariance and eigen decomposition
            Eigen::Matrix<float, 5, 5> cov5 = Eigen::Matrix<float, 5, 5>::Zero();
            if (m_cfg.useResiduals) { // In this case we fill the LUT from the average residuals
              const auto& acc = mResiduals[pdg][entryIdx];
              if (acc.entries() >= 2) {
                const auto covResidual = acc.coreCovariance();
                for (int i = 0; i < 5; ++i) {
                  for (int j = 0; j < 5; ++j) {
                    cov5(i, j) = static_cast<float>(covResidual(i, j));
                  }
                }
                int covIdx = 0;
                for (int i = 0; i < 5; ++i) {
                  for (int j = 0; j <= i; ++j) {
                    const float v = cov5(i, j);
                    entry.covm[covIdx++] = v;
                  }
                }
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix<float, 5, 5>> eigSolver(cov5);
                if (eigSolver.info() == Eigen::Success) {
                  const auto& eigValues = eigSolver.eigenvalues();
                  const auto& eigVectors = eigSolver.eigenvectors();
                  const auto eigVectorsInv = eigVectors.inverse();
                  for (int i = 0; i < 5; ++i) {
                    entry.eigval[i] = eigValues(i);
                    for (int j = 0; j < 5; ++j) {
                      entry.eigvec[i][j] = eigVectors(i, j);
                      entry.eiginv[i][j] = eigVectorsInv(i, j);
                    }
                  }
                } else {
                  for (int i = 0; i < 5; ++i) {
                    entry.eigval[i] = cov5(i, i);
                    for (int j = 0; j < 5; ++j) {
                      entry.eigvec[i][j] = (i == j) ? 1.f : 0.f;
                      entry.eiginv[i][j] = (i == j) ? 1.f : 0.f;
                    }
                  }
                }
                entry.valid = true;

                acc.computePulls(h_pull0, pt_bin, 0);
                acc.computePulls(h_pull1, pt_bin, 1);
                acc.computePulls(h_pull2, pt_bin, 2);
                acc.computePulls(h_pull3, pt_bin, 3);
                acc.computePulls(h_pull4, pt_bin, 4);
              }
            } else { // In this case we fill the LUT from the average covariance matrix
              const auto count = entryCounts[entryIdx];
              if (count > 0u) {
                int covIdx = 0;
                for (int i = 0; i < 5; ++i) {
                  for (int j = 0; j <= i; ++j) {
                    const double mean = covSums[entryIdx][covIdx] / static_cast<double>(count);
                    double value = mean;
                    if (m_cfg.storeCovarianceSigma) {
                      if (count > 1u) {
                        const double variance = cov2Sums[entryIdx][covIdx] / static_cast<double>(count) - mean * mean;
                        value = std::sqrt(std::max(variance, 0.0));
                      } else {
                        value = 0.0;
                      }
                    }
                    const float storedValue = static_cast<float>(value);
                    entry.covm[covIdx] = storedValue;
                    cov5(i, j) = storedValue;
                    cov5(j, i) = storedValue;
                    ++covIdx;
                  }
                }
              }

              Eigen::SelfAdjointEigenSolver<Eigen::Matrix<float, 5, 5>> eigSolver(cov5);
              if (eigSolver.info() == Eigen::Success) {
                const auto& eigValues = eigSolver.eigenvalues();
                const auto& eigVectors = eigSolver.eigenvectors();
                const auto eigVectorsInv = eigVectors.inverse();
                for (int i = 0; i < 5; ++i) {
                  entry.eigval[i] = eigValues(i);
                  for (int j = 0; j < 5; ++j) {
                    entry.eigvec[i][j] = eigVectors(i, j);
                    entry.eiginv[i][j] = eigVectorsInv(i, j);
                  }
                }
              } else {
                for (int i = 0; i < 5; ++i) {
                  entry.eigval[i] = cov5(i, i);
                  for (int j = 0; j < 5; ++j) {
                    entry.eigvec[i][j] = (i == j) ? 1.f : 0.f;
                    entry.eiginv[i][j] = (i == j) ? 1.f : 0.f;
                  }
                }
              }
              entry.valid = true;
            }
          }
          if (h_pull0->GetEntries() > 0) {
            ACTS_DEBUG("Saving pull histogram for PDG " << pdg << ", nch " << nch_bin << ", rad " << rad_bin << ", eta " << eta_bin);
            h_pull0->SaveAs(Form("/tmp/%s.root", h_pull0->GetName()));
          }
          if (h_pull1->GetEntries() > 0) {
            ACTS_DEBUG("Saving pull histogram for PDG " << pdg << ", nch " << nch_bin << ", rad " << rad_bin << ", eta " << eta_bin);
            h_pull1->SaveAs(Form("/tmp/%s.root", h_pull1->GetName()));
          }
          if (h_pull2->GetEntries() > 0) {
            ACTS_DEBUG("Saving pull histogram for PDG " << pdg << ", nch " << nch_bin << ", rad " << rad_bin << ", eta " << eta_bin);
            h_pull2->SaveAs(Form("/tmp/%s.root", h_pull2->GetName()));
          }
          if (h_pull3->GetEntries() > 0) {
            ACTS_DEBUG("Saving pull histogram for PDG " << pdg << ", nch " << nch_bin << ", rad " << rad_bin << ", eta " << eta_bin);
            h_pull3->SaveAs(Form("/tmp/%s.root", h_pull3->GetName()));
          }
          if (h_pull4->GetEntries() > 0) {
            ACTS_DEBUG("Saving pull histogram for PDG " << pdg << ", nch " << nch_bin << ", rad " << rad_bin << ", eta " << eta_bin);
            h_pull4->SaveAs(Form("/tmp/%s.root", h_pull4->GetName()));
          }
        }
      }
    }

    LutRootData rootLut(lut.header, m_cfg.storeCovarianceSigma);
    rootLut.reserve(lut.entries.size());

    const std::string lutFileTag = Form("%s%s_%.2fT", m_cfg.lutTag.empty() ? "" : "_", m_cfg.lutTag.c_str(), m_cfg.magneticField);
    const std::string lutFileName = outputDir + std::string("lut_pdg") + std::to_string(pdg) + lutFileTag;
    ACTS_DEBUG("Saving LUT for PDG " << pdg << " to file " << lutFileName << ".dat");
    // output file
    std::ofstream lutFile(lutFileName + ".dat", std::ofstream::binary);
    if (!lutFile.is_open()) {
      ACTS_ERROR("Did not manage to open output file " << lutFileName << ".dat");
      return;
    }

    o2::delphes::lutHeader_t lutHeader = lut.header;
    lutFile.write(reinterpret_cast<char*>(&lutHeader), sizeof(o2::delphes::lutHeader_t));

    ACTS_DEBUG("LUT dimensions: nch_bins = " << nnch << ", rad_bins = " << nrad
                                             << ", eta_bins = " << neta << ", pt_bins = " << npt);

    for (int nch_bin = 0; nch_bin < nnch; ++nch_bin) {
      for (int rad_bin = 0; rad_bin < nrad; ++rad_bin) {
        for (int eta_bin = 0; eta_bin < neta; ++eta_bin) {
          for (int pt_bin = 0; pt_bin < npt; ++pt_bin) {
            const int entryIndex = nch_bin * (nrad * neta * npt) + rad_bin * (neta * npt) + eta_bin * npt + pt_bin;
            auto* entry = &lut.entries[entryIndex];
            if (entry->valid) {
              entry->print();
            }
            lutFile.write(reinterpret_cast<const char*>(entry), sizeof(o2::delphes::lutEntry_t));
            rootLut.addEntry(entryIndex, nch_bin, rad_bin, eta_bin, pt_bin, *entry);
          }
        }
      }
    }

    lutFile.close();

    ACTS_DEBUG("Saving ROOT LUT for PDG " << pdg << " to file " << lutFileName << ".root");
    if (!rootLut.writeToFile(lutFileName + ".root")) {
      ACTS_ERROR("Did not manage to write ROOT LUT file " << lutFileName << ".root");
      return;
    }
  }

  TFile output((outputDir + "lutMakerQA.root").c_str(), "RECREATE");
  h_stat->Write();
  auto writeResidualProfiles = [](const std::map<int, TH3F*>& histMap) {
    for (const auto& [pdg, hist] : histMap) {
      hist->Write();
      const char* name = hist->GetName();

      // Store the residual-vs-pT projection integrated over eta.
      TH2D* res2D = static_cast<TH2D*>(hist->Project3D("yx"));
      res2D->SetName(Form("%s_residualVsPt", name));
      res2D->Write();

      // Explicitly compute mean and full RMS in each physical pT bin,
      // integrated over all eta. This avoids relying on TProfile error
      // semantics and preserves the logarithmic LUT bin edges.
      TH1D* resWidth = hist->ProjectionX(Form("%s_width", name));
      TH1D* resMean = hist->ProjectionX(Form("%s_mean", name));
      resWidth->Reset("ICES");
      resMean->Reset("ICES");

      for (int iPt = 1; iPt <= hist->GetNbinsX(); ++iPt) {
        TH1D* residualProjection = hist->ProjectionY(
          Form("%s_py_%d", name, iPt),
          iPt, iPt,
          1, hist->GetNbinsZ());

        if (residualProjection->GetEntries() > 0.) {
          resWidth->SetBinContent(iPt, residualProjection->GetRMS());
          resWidth->SetBinError(iPt, residualProjection->GetRMSError());
          resMean->SetBinContent(iPt, residualProjection->GetMean());
          resMean->SetBinError(iPt, residualProjection->GetMeanError());
        }
        delete residualProjection;
      }

      resWidth->Write();
      resMean->Write();
      delete res2D;
      delete resWidth;
      delete resMean;
    }
  };

  writeResidualProfiles(res_qopt_vs_pT);
  writeResidualProfiles(res_pt_vs_pT);
  writeResidualProfiles(res_acts_d0_vs_pT);
  writeResidualProfiles(res_o2_y_vs_pT);

  output.Close();

  if (mTree) {
    mTree->SaveAs((outputDir + "tree.root").c_str());
  }
}

ActsExamples::ProcessCode LutMaker::execute(const ActsExamples::AlgorithmContext& ctx) const
{
  // Read the final reconstructed tracks and their explicit truth matching.
  const ConstTrackContainer& tracks = m_inputTracks(ctx);
  const auto& particles = m_inputParticles(ctx);
  const auto& trackParticleMatching = m_inputTrackParticleMatching(ctx);

  h_stat->Fill(0);
  ACTS_DEBUG("Number of input tracks: " << tracks.size());
  h_stat->Fill(1, static_cast<double>(tracks.size()));
  ACTS_DEBUG("Number of input particles: " << particles.size());
  h_stat->Fill(2, static_cast<double>(particles.size()));

  for (const auto& particle : particles) {

    const int pdg = std::abs(particle.pdg());
    auto lutIt = mLutMap.find(pdg);
    if (lutIt == mLutMap.end()) {
      continue;
    }

    const auto& p = particle.momentum();
    const double ptTrue = std::hypot(p.x(), p.y());
    const int ptBin = lutIt->second.header.ptmap.find(static_cast<float>(ptTrue));
    ++mPtTruthCounts[pdg][ptBin];
  }

  for (const auto& trackProxy : tracks) {
    const auto trackIndex = trackProxy.index();
    ACTS_VERBOSE("Processing track index " << static_cast<std::size_t>(trackIndex));

    // A fitted reference state is required to build the LUT entry.
    if (!trackProxy.hasReferenceSurface()) {
      continue;
    }

    const auto matchIt = trackParticleMatching.find(trackIndex);
    if (matchIt == trackParticleMatching.end()) {
      continue;
    }

    const auto& match = matchIt->second;
    // Reject fake and duplicate tracks. The matcher assigns exactly one
    // TrackMatchClassification::Matched track to each reconstructed particle.
    if (match.classification != TrackMatchClassification::Matched ||
        !match.particle.has_value()) {
      continue;
    }

    const auto particleIt = particles.find(match.particle.value());
    if (particleIt == particles.end()) {
      continue;
    }
    const SimParticle* particlePtr = &(*particleIt);
    h_stat->Fill(3, 1.0);

    // Preserve the original track index for truth lookup, but use the fitted
    // parameters at the reference surface for the existing LUT calculations.
    const Acts::BoundTrackParameters track = trackProxy.createParametersAtReference();

    const int pdg = std::abs(particlePtr->pdg());
    LutMaker::mTreeContainer.import(particlePtr, &track);
    if (mLutMap.find(pdg) == mLutMap.end()) {
      ACTS_DEBUG("Skipping particle with unsupported PDG " << particlePtr->pdg());
      continue;
    } else {
      ACTS_DEBUG(" + Processing matched particle with PDG " << particlePtr->pdg());
    }

    const double ptReco = track.transverseMomentum();
    const auto& p = particlePtr->momentum();
    const double ptTrue = std::hypot(p.x(), p.y());

    const double relRes = (ptReco - ptTrue) / ptTrue;
    ACTS_DEBUG("matched particle: ptTrue=" << ptTrue << ", ptReco=" << ptReco
                                           << ", relRes=" << relRes);

    ACTS_VERBOSE("  pxReco = " << track.momentum().x() << ", pxTrue = " << p.x());
    ACTS_VERBOSE("  pyReco = " << track.momentum().y() << ", pyTrue = " << p.y());
    ACTS_VERBOSE("  pzReco = " << track.momentum().z() << ", pzTrue = " << p.z());
    const float mod = std::sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z());
    const float eta = 0.5f * std::log((mod + p.z()) / (mod - p.z()));

    // Fill the LUT information
    auto& lut = mLutMap[pdg];
    const int nch_bin = 0;
    const int rad_bin = 0; // example with only one radial bin
    const int eta_bin = lut.header.etamap.find(eta);
    const int pt_bin = lut.header.ptmap.find(ptTrue);

    ++mPtMatchedCounts[pdg][pt_bin];

    const std::size_t entryIdx = static_cast<std::size_t>(nch_bin) * (lut.header.radmap.nbins * lut.header.etamap.nbins * lut.header.ptmap.nbins) + rad_bin * (lut.header.etamap.nbins * lut.header.ptmap.nbins) + eta_bin * lut.header.ptmap.nbins + pt_bin;

    // Convert the reconstructed track to the standalone O2 five-parameter
    // representation (Y, Z, snp, tgl, q/pT) instead of re-deriving the
    // q/p -> q/pT covariance rotation by hand.
    std::optional<standalone_o2::TrackParametrizationWithError> recoO2;
    try {
      recoO2 = standalone_o2::convertActsToO2TrackWithCovariance(track, ctx.simGeoContext);
    } catch (const std::exception& e) {
      ACTS_VERBOSE("Skipping track: could not convert to O2 parametrization: " << e.what());
      continue;
    }

    const int charge = recoO2->getCharge();

    // ----------------------------------------------------------------------
    // ACTS-native d0 residual, using exactly the fitted track reference
    // surface and the same truth construction as RootTrackSummaryWriter.
    // The ACTS bound local coordinate is converted to cm for direct comparison
    // with O2 Y and with sqrt(covm[0]).
    // ----------------------------------------------------------------------
    const auto truthLocal =
      truthLocalAtReference(track, *particlePtr, ctx.simGeoContext);

    if (truthLocal.has_value()) {
      const double d0Reco = static_cast<double>(track.parameters()[Acts::eBoundLoc0]);
      const double d0Truth = static_cast<double>((*truthLocal)[Acts::eBoundLoc0]);
      const double actsD0ResidualCm = (d0Reco - d0Truth) / Acts::UnitConstants::cm;
      res_acts_d0_vs_pT[pdg]->Fill(ptTrue, actsD0ResidualCm, eta);
    } else {
      ACTS_VERBOSE("Could not determine truth local coordinates on the track reference surface");
    }

    if (charge != 0 && ptTrue > 0.0) {
      const double qOverPtResidual = static_cast<double>(recoO2->getQ2Pt()) - static_cast<double>(charge) / ptTrue;
      res_qopt_vs_pT[pdg]->Fill(ptTrue, qOverPtResidual, eta);
      res_pt_vs_pT[pdg]->Fill(ptTrue, ptReco - ptTrue, eta);
    }

    // Accumulate covariance lower-triangle statistics for later mean/sigma
    // extraction, and keep the CURRENT O2 residual definition used by the LUT.
    // This is intentionally not changed yet: the new QA histograms let us
    // compare it directly with the ACTS-native d0 residual first.
    {
      auto& covSum = mCovSums[pdg][entryIdx];
      auto& cov2Sum = mCov2Sums[pdg][entryIdx];
      auto& residual = mResiduals[pdg][entryIdx];

      const standalone_o2::TrackParametrization::dim3_t truthPositionCm{static_cast<float>(particlePtr->position().x() / Acts::UnitConstants::cm),
                                                                        static_cast<float>(particlePtr->position().y() / Acts::UnitConstants::cm),
                                                                        static_cast<float>(particlePtr->position().z() / Acts::UnitConstants::cm)};
      const standalone_o2::TrackParametrization::dim3_t truthMomentumGeV{static_cast<float>(p.x() / Acts::UnitConstants::GeV),
                                                                         static_cast<float>(p.y() / Acts::UnitConstants::GeV),
                                                                         static_cast<float>(p.z() / Acts::UnitConstants::GeV)};

      // IMPORTANT: this is the old LUT residual definition. The reconstructed
      // and truth O2 objects choose their own local alpha, and the truth object
      // is built at the generator production point. We retain it here only to
      // diagnose how much this differs from the ACTS-native d0 residual.
      const standalone_o2::TrackParametrization truthO2{
        truthPositionCm, truthMomentumGeV, charge};

      ResidualAccumulator::Vector residualVector;
      for (int i = 0; i < ResidualAccumulator::kSize; ++i) {
        residualVector(i) = static_cast<double>(recoO2->getParams()[i]) -
                            static_cast<double>(truthO2.getParams()[i]);
      }

      // Parameter 0 is O2 local Y, not ACTS d0. Store it explicitly as QA.
      res_o2_y_vs_pT[pdg]->Fill(ptTrue, residualVector(0), eta);

      residual.add(residualVector);

      const auto& covm = recoO2->getCov();
      for (std::size_t covIdx = 0; covIdx < covm.size(); ++covIdx) {
        const double covElement = static_cast<double>(covm[covIdx]);
        covSum[covIdx] += covElement;
        cov2Sum[covIdx] += covElement * covElement;
      }
      ++mEntryCounts[pdg][entryIdx];
    }
  }

  return ProcessCode::SUCCESS;
}

} // namespace AliceActsTrk
