#include "SeedAnalysisPlotTool.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <format>
#include <iomanip>
#include <sstream>




using namespace ActsExamples;

namespace AliceActsTrk
{

    SeedAnalysisPlotTool::SeedAnalysisPlotTool(const SeedAnalysisPlotTool::Config& cfg,
        Acts::Logging::Level lvl)
        : m_cfg(cfg), m_logger(Acts::getDefaultLogger("SeedAnalysisPlotTool", lvl))
        {
            const auto& xAxis = m_cfg.varBinning.at("X");
            const auto& yAxis = m_cfg.varBinning.at("Y");
            const auto& zAxis = m_cfg.varBinning.at("Z");
            const auto& rAxis = m_cfg.varBinning.at("R");

            ACTS_DEBUG("Initialize the histograms for seed analysis");

            // middle space space points
            std::vector<std::string> spNames{"botSP","midSP","topSP"};

            for (const auto& spName : spNames) {
            m_RvsZ.emplace(std::format("{}_RvsZ",spName),
                Histogram2(
                        std::format("{}_RvsZ", spName),
                        std::format("{}_RvsZ", spName),
                        std::array{zAxis, rAxis}));
                }
        }
    void SeedAnalysisPlotTool::fill(ConstSeedProxy seed, ConstTrackProxy seed_track) {

        const auto& spacePoints = seed.spacePoints();

        m_RvsZ.at("botSP_RvsZ").fill({spacePoints[0].z(), spacePoints[0].r()});
        m_RvsZ.at("midSP_RvsZ").fill({spacePoints[1].z(), spacePoints[1].r()});
        m_RvsZ.at("topSP_RvsZ").fill({spacePoints[2].z(), spacePoints[2].r()});

    }


} // namespace AliceActsTrk
