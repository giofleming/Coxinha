#pragma once

#if __has_include("Acts/EventData/BoundTrackParameters.hpp")
#include "Acts/EventData/BoundTrackParameters.hpp"
#else
#include "Acts/EventData/TrackParameters.hpp"
#endif

#include "Acts/Geometry/GeometryContext.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace standalone_o2
{

enum class PID : std::uint8_t {
  Electron,
  Muon,
  Pion,
  Kaon,
  Proton,
  Deuteron,
  Triton,
  Helium3,
  Alpha
};

inline constexpr std::size_t kNParams = 5;
inline constexpr std::size_t kCovMatSize = 15;
inline constexpr std::size_t kLabCovMatSize = 21;

enum ParameterIndex : std::size_t {
  kY = 0,
  kZ = 1,
  kSnp = 2,
  kTgl = 3,
  kQ2Pt = 4
};

class TrackParametrization
{
 public:
  using value_t = float;
  using dim3_t = std::array<value_t, 3>;
  using params_t = std::array<value_t, kNParams>;

  TrackParametrization() = default;
  TrackParametrization(const dim3_t& xyz, const dim3_t& pxpypz,
                       int charge, bool sectorAlpha = true,
                       const PID pid = PID::Pion);

  value_t getX() const noexcept { return mX; }
  value_t getAlpha() const noexcept { return mAlpha; }
  value_t getY() const noexcept { return mParameters[kY]; }
  value_t getZ() const noexcept { return mParameters[kZ]; }
  value_t getSnp() const noexcept { return mParameters[kSnp]; }
  value_t getTgl() const noexcept { return mParameters[kTgl]; }
  value_t getQ2Pt() const noexcept { return mParameters[kQ2Pt]; }
  const params_t& getParams() const noexcept { return mParameters; }
  int getAbsCharge() const noexcept { return mAbsCharge; }
  int getCharge() const noexcept;
  PID getPID() const noexcept { return mPID; }

  value_t getPt() const noexcept;
  value_t getP() const noexcept;
  void getXYZGlo(dim3_t& xyz) const noexcept;
  bool getPxPyPzGlo(dim3_t& pxpypz) const noexcept;

 protected:
  value_t mX = 0.f;
  value_t mAlpha = 0.f;
  params_t mParameters{};
  int mAbsCharge = 1;
  PID mPID = PID::Pion;
};

class TrackParametrizationWithError : public TrackParametrization
{
 public:
  using covMat_t = std::array<value_t, kCovMatSize>;
  using labCovMat_t = std::array<value_t, kLabCovMatSize>;

  TrackParametrizationWithError() = default;
  TrackParametrizationWithError(const dim3_t& xyz, const dim3_t& pxpypz,
                                const labCovMat_t& globalCovariance, int charge,
                                bool sectorAlpha = true, const PID pid = PID::Pion);

  const covMat_t& getCov() const noexcept { return mCovariance; }
  value_t getCovariance(std::size_t row, std::size_t column) const;

 private:
  covMat_t mCovariance{};
};

TrackParametrization convertActsToO2Track(const Acts::BoundTrackParameters& actsTrack,
                                          const Acts::GeometryContext& geometryContext,
                                          bool sectorAlpha = true,
                                          const PID pid = PID::Pion);

TrackParametrizationWithError convertActsToO2TrackWithCovariance(const Acts::BoundTrackParameters& actsTrack,
                                                                 const Acts::GeometryContext& geometryContext,
                                                                 bool sectorAlpha = true,
                                                                 const PID pid = PID::Pion);

} // namespace standalone_o2
