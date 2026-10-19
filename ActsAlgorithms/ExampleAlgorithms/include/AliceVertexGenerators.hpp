#pragma once

// ----------------------------------------------------------------------------
//  BunchCrossingVertexPositionGenerator
//
//  Drop-in replacement for ActsExamples::GaussianVertexGenerator with a
//  DISCRETE (bunch-crossing) time distribution:
//    - x, y, z are Gaussian (mean +/- stddev), exactly like the Gaussian one;
//    - the time is drawn UNIFORMLY among `nBunchCrossings` equally-probable
//      slots spaced by `bunchSpacing`, starting at mean[3].
//
//  For the defaults (20 slots x 25 ns) the times are
//      t in mean[3] + {0, 25, ..., 475} ns.
// ----------------------------------------------------------------------------

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Utilities/VertexGenerators.hpp"

#include <cstddef>
#include <random>

namespace AliceActsTrk {

struct BunchCrossingVertexPositionGenerator
    : public ActsExamples::PrimaryVertexPositionGenerator {
  /// Spatial spread (x,y,z); the time component (index 3) is ignored.
  Acts::Vector4 stddev = {0.0, 0.0, 0.0, 0.0};
  /// Mean position; mean[3] is an overall time offset of the bunch train.
  Acts::Vector4 mean = {0.0, 0.0, 0.0, 0.0};
  /// Number of bunch crossings (equally-probable time slots).
  std::size_t nBunchCrossings = 20;
  /// Spacing between bunch crossings, in native time units (e.g. 25 * ns).
  double bunchSpacing = 25.0 * Acts::UnitConstants::ns;

  Acts::Vector4 operator()(ActsExamples::RandomEngine& rng,
                           std::size_t /*eventNumber*/) const override {
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_int_distribution<std::size_t> bc(
        0, (nBunchCrossings > 0) ? nBunchCrossings - 1 : 0);
    const std::size_t k = bc(rng);  // pick a bunch slot, uniform probability
    const double t = mean[3] + static_cast<double>(k) * bunchSpacing;
    return {mean[0] + normal(rng) * stddev[0],
            mean[1] + normal(rng) * stddev[1],
            mean[2] + normal(rng) * stddev[2], t};
  }
};

}  // namespace AliceActsTrk
