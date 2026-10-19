#pragma once

#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <memory>
#include <vector>
#include <string>
#include <mutex>


#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/EventData/Seed.hpp"
//#include "Acts/EventData/SeedContainer.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "SeedAnalysisPlotTool.hpp"

#include <mutex>
#include <string>

class TFile;
class TTree;

using namespace ActsExamples;

namespace AliceActsTrk
{

class SeedAnalysisAlg final
  : public WriterT<SeedContainer>
{
 public:
  struct Config {
    /// Input seed collection.
    std::string inputSeeds;
    /// Input Track seed collection
    std::string inputTracks;
    /// Input particles collection.
    std::string inputParticles;
    /// Input track-particle matching.
    std::string inputTrackParticleMatching;
    /// Input measurement -> particles map (for the strict shared-particle flag).
    std::string inputMeasurementParticlesMap;
    /// Output filename.
    std::string filePath = "seed_analysis.root";
    /// Plot tool configurations.
    ///...///
    AliceActsTrk::SeedAnalysisPlotTool::Config seedAnalysisPlotToolConfig;
  };

  /// Construct from configuration and log level.
  /// @param config The configuration
  /// @param level The logger level
  SeedAnalysisAlg(Config config, Acts::Logging::Level level);

  ~SeedAnalysisAlg() override;

  /// Finalize plots.
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const SeedContainer& seeds) override;

  Config m_cfg;

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{
    this, "InputTrackParticleMatching"};
  ReadDataHandle<MeasurementParticlesMap> m_inputMeasParticlesMap{
    this, "InputMeasurementParticlesMap"};
  ReadDataHandle<ConstTrackContainer> m_inputTracks{this,"InputTracks"};

  /// Whether the optional truth inputs were configured (drive the two flags).
  bool m_haveMatching{false};
  bool m_haveMeasMap{false};

  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile* m_outputFile{nullptr};
  /// Per-seed output tree and its branch buffers.
  TTree* m_tree{nullptr};
  int m_bEvent{0};
  // Local coordinates of the bottom space point on its sensor surface -- this is
  // what the seed-track bound parameters are referenced to (NOT d0/z0).
  float m_bBotLoc0{0.f}, m_bBotLoc1{0.f};
  float m_bPt{0.f}, m_bEta{0.f}, m_bPhi{0.f};
  float m_bSeedQuality{0.f}, m_bVertexZ{0.f};
  // Seed-finder linearized transverse impact estimate (the quantity impactMax
  // actually cuts), recomputed from the three space points.
  float m_bIm{0.f};
  int m_bNSeedsPerSpM{0};
  float m_bBotR{0.f}, m_bBotZ{0.f}, m_bMidR{0.f}, m_bMidZ{0.f}, m_bTopR{0.f}, m_bTopZ{0.f};
  int m_bIsFakeSeedTrack{-1};
  int m_bIsFakeShared{-1};
  /// Plot tool for residuals and pulls.
  AliceActsTrk::SeedAnalysisPlotTool m_seedAnalysisPlotTool;
}; //SeedAnalysisAlg

} // namespace AliceActsTrk