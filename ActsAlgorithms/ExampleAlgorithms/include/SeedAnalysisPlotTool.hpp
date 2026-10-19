#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Seed.hpp"
#include "Acts/EventData/SeedContainer.hpp"
#include "ActsExamples/EventData/Track.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace ActsExamples;

namespace AliceActsTrk
{
    class SeedAnalysisPlotTool
{
 public:
  using AxisVariant = Acts::Experimental::AxisVariant;
  using BoostRegularAxis = Acts::Experimental::BoostRegularAxis;
  using Histogram1 = Acts::Experimental::Histogram1;
  using Histogram2 = Acts::Experimental::Histogram2;
  
  
  struct Config {

    /// Binning info for variables
    std::map<std::string, AxisVariant> varBinning = {
      {"X", BoostRegularAxis(100, -500, 500, "X")},
      {"Y", BoostRegularAxis(70, -3.5, 3.5, "Y")},
      {"Z", BoostRegularAxis(600, -3000, 3000, "Z")},
      {"R", BoostRegularAxis(200, 0, 1000, "R")},
    };
  }; // config
    
    
    /// Constructor
    ///
    /// @param cfg Configuration struct
    /// @param level Message level declaration
    SeedAnalysisPlotTool(const Config& cfg, Acts::Logging::Level lvl);

    void fill(ConstSeedProxy seed, ConstTrackProxy seed_track);

    /// @brief Accessors for histograms (const reference)
    const std::map<std::string, Histogram2>& RvsZ() const
  {
    return m_RvsZ;
  }

  private:
  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
  /// R vs Z scatter plot
  std::map<std::string, Histogram2> m_RvsZ;


}; // class SeedAnalysisPlotTool

} // Namespace AliceActsTrk