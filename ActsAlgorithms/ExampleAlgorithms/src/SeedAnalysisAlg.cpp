// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "SeedAnalysisAlg.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsPlugins/Root/HistogramConverter.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <TEfficiency.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TTree.h>


using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::phi;
using ActsPlugins::toRoot;

using namespace ActsExamples;

namespace AliceActsTrk
{
    SeedAnalysisAlg::SeedAnalysisAlg(SeedAnalysisAlg::Config config, Acts::Logging::Level level)
    : WriterT(config.inputSeeds, "SeedAnalysisAlg", level),
    m_cfg(std::move(config)),
    m_seedAnalysisPlotTool(m_cfg.seedAnalysisPlotToolConfig, level) {


        // the output file can not be given externally since TFile accesses to the
        // same file from multiple threads are unsafe.
        // must always be opened internally
        auto path = m_cfg.filePath;
        m_outputFile = TFile::Open(path.c_str(), "RECREATE");
        if (m_outputFile == nullptr) {
            throw std::invalid_argument("Could not open '" + path + "'");
        }


        m_inputTracks.initialize(m_cfg.inputTracks);

        // Optional truth inputs -- each drives one of the two fake flags.
        if (!m_cfg.inputTrackParticleMatching.empty()) {
            m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);
            m_haveMatching = true;
        }
        if (!m_cfg.inputMeasurementParticlesMap.empty()) {
            m_inputMeasParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
            m_haveMeasMap = true;
        }

        // Per-seed output tree.
        m_outputFile->cd();
        m_tree = new TTree("seeds", "per-seed kinematics, quality and truth flags");
        m_tree->Branch("event", &m_bEvent);
        // Local hit coords of the bottom SP on its sensor (what the seed-track
        // params are actually referenced to -- these are NOT impact parameters).
        m_tree->Branch("bottomLoc0", &m_bBotLoc0);
        m_tree->Branch("bottomLoc1", &m_bBotLoc1);
        // Seed-finder linear transverse impact estimate that impactMax cuts [mm].
        m_tree->Branch("im", &m_bIm);
        m_tree->Branch("pt", &m_bPt);
        m_tree->Branch("eta", &m_bEta);
        m_tree->Branch("phi", &m_bPhi);
        m_tree->Branch("seedQuality", &m_bSeedQuality);
        m_tree->Branch("vertexZ", &m_bVertexZ);
        m_tree->Branch("nSeedsPerSpM", &m_bNSeedsPerSpM);
        m_tree->Branch("botR", &m_bBotR);
        m_tree->Branch("botZ", &m_bBotZ);
        m_tree->Branch("midR", &m_bMidR);
        m_tree->Branch("midZ", &m_bMidZ);
        m_tree->Branch("topR", &m_bTopR);
        m_tree->Branch("topZ", &m_bTopZ);
        // 1 = fake, 0 = truth, -1 = flag unavailable (input not configured)
        m_tree->Branch("isFake_seedTrack", &m_bIsFakeSeedTrack);
        m_tree->Branch("isFake_sharedParticle", &m_bIsFakeShared);
    }//


    SeedAnalysisAlg::~SeedAnalysisAlg()
    {
    delete m_outputFile;  
    }

    ProcessCode SeedAnalysisAlg::finalize()
    {
        if (m_outputFile == nullptr) {
            return ProcessCode::SUCCESS;
        }   
        
        m_outputFile->cd();

        if (m_tree != nullptr) {
            m_tree->Write();
        }

        for (const auto& [name, hist] : m_seedAnalysisPlotTool.RvsZ()) {
            std::cout<<"Writing down "<< name<<std::endl;
            toRoot(hist)->Write();
        }

        if (m_outputFile != nullptr) {
            m_outputFile->Close();
        }

    
    return ProcessCode::SUCCESS;        
    }

    ProcessCode SeedAnalysisAlg::writeT(
        const AlgorithmContext& ctx, const SeedContainer& seeds) {

            // Read seed-track collection
            const auto& seed_tracks = m_inputTracks(ctx);    

            // Exclusive access to the tree while writing
            std::lock_guard<std::mutex> lock(m_writeMutex);

	    if (seeds.size() != seed_tracks.size()) {  
                ACTS_WARNING("Seed/seed-track count mismatch: "
                    << seeds.size() << " vs " << seed_tracks.size());
}
            const auto n = std::min(seeds.size(), seed_tracks.size());

            // Optional truth inputs for the two fake flags.
            const TrackParticleMatching* matching =
                m_haveMatching ? &m_inputTrackParticleMatching(ctx) : nullptr;
            const MeasurementParticlesMap* measMap =
                m_haveMeasMap ? &m_inputMeasParticlesMap(ctx) : nullptr;

            // Pre-count how many seeds share each middle space point (<= maxSeedsPerSpM).
            std::map<std::size_t, int> midSpCount;
            for (const auto i : std::views::iota(std::size_t{0}, n)) {
                midSpCount[static_cast<std::size_t>(seeds[i].spacePoints()[1].index())]++;
            }

            // Union of truth particles contributing to a space point's measurements.
            auto spParticles = [&](const auto& sp) {
                std::set<SimBarcode> parts;
                if (measMap == nullptr) return parts;
                for (const auto& sl : sp.sourceLinks()) {
                    const auto idx = sl.template get<IndexSourceLink>().index();
                    auto range = measMap->equal_range(idx);
                    for (auto it = range.first; it != range.second; ++it) {
                        parts.insert(it->second);
                    }
                }
                return parts;
            };

            for (const auto i : std::views::iota(std::size_t{0}, n)) {
                const auto& seed = seeds[i];
                const auto& seed_track = seed_tracks.getTrack(i);
                m_seedAnalysisPlotTool.fill(seed, seed_track);

                const auto sps = seed.spacePoints();
                const auto bot = sps[0];
                const auto mid = sps[1];
                const auto top = sps[2];

                m_bEvent = static_cast<int>(ctx.eventNumber);
                m_bSeedQuality = seed.quality();
                m_bVertexZ = seed.vertexZ();
                m_bNSeedsPerSpM = midSpCount[static_cast<std::size_t>(mid.index())];
                m_bBotR = bot.r(); m_bBotZ = bot.z();
                m_bMidR = mid.r(); m_bMidZ = mid.z();
                m_bTopR = top.r(); m_bTopZ = top.z();

                // Recompute the seed-finder's linearized transverse impact
                // estimate 'im' = |(A - B*rM)*rM|, i.e. the exact quantity cut by
                // impactMax, from the 3 space points (conformal U/V transform in
                // the middle-SP frame; cosPhiM = xM/rM, sinPhiM = yM/rM).
                m_bIm = std::numeric_limits<float>::quiet_NaN();
                {
                    const double xM = mid.x(), yM = mid.y(), rM = mid.r();
                    if (rM > 0.0) {
                        const double cosPhiM = xM / rM, sinPhiM = yM / rM;
                        auto uv = [&](const auto& sp, double& u, double& v) -> bool {
                            const double dX = sp.x() - xM, dY = sp.y() - yM;
                            const double dR2 = dX * dX + dY * dY;
                            if (dR2 == 0.0) return false;
                            const double xN = dX * cosPhiM + dY * sinPhiM;
                            const double yN = dY * cosPhiM - dX * sinPhiM;
                            u = xN / dR2; v = yN / dR2;
                            return true;
                        };
                        double uB, vB, uT, vT;
                        if (uv(bot, uB, vB) && uv(top, uT, vT)) {
                            const double dU = uT - uB;
                            if (dU != 0.0) {
                                const double A = (vT - vB) / dU;
                                const double B = vB - A * uB;
                                m_bIm = static_cast<float>(std::abs((A - B * rM) * rM));
                            }
                        }
                    }
                }

                // The seed-track bound params are referenced to the BOTTOM-SP
                // sensor surface, so loc0/loc1 are local hit coords there (not
                // d0/z0). Direction (phi/theta) and momentum are still valid.
                const auto& p = seed_track.parameters();
                const double theta = p[Acts::eBoundTheta];
                m_bBotLoc0 = static_cast<float>(p[Acts::eBoundLoc0]);
                m_bBotLoc1 = static_cast<float>(p[Acts::eBoundLoc1]);
                m_bPhi = static_cast<float>(p[Acts::eBoundPhi]);
                m_bEta = static_cast<float>(-std::log(std::tan(theta / 2.0)));
                m_bPt = static_cast<float>(seed_track.transverseMomentum());

                // Flag A: was the track built from this seed matched to a particle?
                m_bIsFakeSeedTrack = -1;
                if (matching != nullptr) {
                    auto it = matching->find(seed_track.index());
                    const bool truth =
                        it != matching->end() &&
                        (it->second.classification ==
                             TrackMatchClassification::Matched ||
                         it->second.classification ==
                             TrackMatchClassification::Duplicate);
                    m_bIsFakeSeedTrack = truth ? 0 : 1;
                }

                // Flag B: do all three space points share a common truth particle?
                m_bIsFakeShared = -1;
                if (measMap != nullptr) {
                    const auto pb = spParticles(bot);
                    const auto pm = spParticles(mid);
                    const auto pt = spParticles(top);
                    bool shared = false;
                    for (const auto& b : pb) {
                        if (pm.count(b) != 0 && pt.count(b) != 0) { shared = true; break; }
                    }
                    m_bIsFakeShared = shared ? 0 : 1;
                }

                m_tree->Fill();
            }

            return ProcessCode::SUCCESS;
    }

}; // AliceActsTrk
